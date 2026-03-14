// Standalone LoRa receiver - pure C/C++ reimplementation of the gr-lora_sdr
// demodulation chain. No GNU Radio dependency.
//
// Pipeline: RTL-SDR IQ -> Frame Sync -> FFT Demod -> Gray Demap ->
//           Deinterleave -> Hamming Decode -> Header Decode ->
//           Dewhiten -> CRC Verify -> Output
//
// Build:
//   g++ -O2 -o lora_rx lora_rx.cpp kiss_fft.c -lrtlsdr -lpthread -lm
//
// Usage:
//   ./lora_rx [options]
//     -f <freq_hz>      Center frequency (default 869618000)
//     -s <samp_rate>    Sample rate (default 250000)
//     -b <bandwidth>    LoRa bandwidth (default 62500)
//     -S <sf>           Spreading factor 5-12 (default 8)
//     -c <cr>           Coding rate 1-4 (default 4)
//     -w <sync_word>    Sync word hex (default 0x12)
//     -g <gain>         RTL-SDR gain in 0.1 dB (default 490)
//     -p <ppm>          Frequency correction ppm (default -3)
//     -I                Implicit header mode
//     -L <pay_len>      Payload length for implicit header (default 11)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <csignal>
#include <complex>
#include <vector>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <string>
#include <getopt.h>
#include <rtl-sdr.h>

extern "C" {
#include "kiss_fft.h"
}

// ============================================================================
// Types and constants
// ============================================================================

using cx = std::complex<float>;

static constexpr int MIN_SF = 5;
static constexpr int MAX_SF = 12;
static constexpr float LDRO_MAX_DURATION_MS = 16.0f;

static volatile bool g_running = true;
static rtlsdr_dev_t *g_dev = nullptr;

static void signal_handler(int) {
    g_running = false;
    if (g_dev) rtlsdr_cancel_async(g_dev);
}

// ============================================================================
// Whitening sequence (LFSR)
// ============================================================================

static const uint8_t whitening_seq[] = {
    0xFF, 0xFE, 0xFC, 0xF8, 0xF0, 0xE1, 0xC2, 0x85, 0x0B, 0x17, 0x2F, 0x5E, 0xBC, 0x78, 0xF1, 0xE3,
    0xC6, 0x8D, 0x1A, 0x34, 0x68, 0xD0, 0xA0, 0x40, 0x80, 0x01, 0x02, 0x04, 0x08, 0x11, 0x23, 0x47,
    0x8E, 0x1C, 0x38, 0x71, 0xE2, 0xC4, 0x89, 0x12, 0x25, 0x4B, 0x97, 0x2E, 0x5C, 0xB8, 0x70, 0xE0,
    0xC0, 0x81, 0x03, 0x06, 0x0C, 0x19, 0x32, 0x64, 0xC9, 0x92, 0x24, 0x49, 0x93, 0x26, 0x4D, 0x9B,
    0x37, 0x6E, 0xDC, 0xB9, 0x72, 0xE4, 0xC8, 0x90, 0x20, 0x41, 0x82, 0x05, 0x0A, 0x15, 0x2B, 0x56,
    0xAD, 0x5B, 0xB6, 0x6D, 0xDA, 0xB5, 0x6B, 0xD6, 0xAC, 0x59, 0xB2, 0x65, 0xCB, 0x96, 0x2C, 0x58,
    0xB0, 0x61, 0xC3, 0x87, 0x0F, 0x1F, 0x3E, 0x7D, 0xFB, 0xF6, 0xED, 0xDB, 0xB7, 0x6F, 0xDE, 0xBD,
    0x7A, 0xF5, 0xEB, 0xD7, 0xAE, 0x5D, 0xBA, 0x74, 0xE8, 0xD1, 0xA2, 0x44, 0x88, 0x10, 0x21, 0x43,
    0x86, 0x0D, 0x1B, 0x36, 0x6C, 0xD8, 0xB1, 0x63, 0xC7, 0x8F, 0x1E, 0x3C, 0x79, 0xF3, 0xE7, 0xCE,
    0x9C, 0x39, 0x73, 0xE6, 0xCC, 0x98, 0x31, 0x62, 0xC5, 0x8B, 0x16, 0x2D, 0x5A, 0xB4, 0x69, 0xD2,
    0xA4, 0x48, 0x91, 0x22, 0x45, 0x8A, 0x14, 0x29, 0x52, 0xA5, 0x4A, 0x95, 0x2A, 0x54, 0xA9, 0x53,
    0xA7, 0x4E, 0x9D, 0x3B, 0x77, 0xEE, 0xDD, 0xBB, 0x76, 0xEC, 0xD9, 0xB3, 0x67, 0xCF, 0x9E, 0x3D,
    0x7B, 0xF7, 0xEF, 0xDF, 0xBF, 0x7E, 0xFD, 0xFA, 0xF4, 0xE9, 0xD3, 0xA6, 0x4C, 0x99, 0x33, 0x66,
    0xCD, 0x9A, 0x35, 0x6A, 0xD4, 0xA8, 0x51, 0xA3, 0x46, 0x8C, 0x18, 0x30, 0x60, 0xC1, 0x83, 0x07,
    0x0E, 0x1D, 0x3A, 0x75, 0xEA, 0xD5, 0xAA, 0x55, 0xAB, 0x57, 0xAF, 0x5F, 0xBE, 0x7C, 0xF9, 0xF2,
    0xE5, 0xCA, 0x94, 0x28, 0x50, 0xA1, 0x42, 0x84, 0x09, 0x13, 0x27, 0x4F, 0x9F, 0x3F, 0x7F
};

// ============================================================================
// Utility functions
// ============================================================================

static inline long mod(long a, long b) { return ((a % b) + b) % b; }

static inline cx expj(float phase) { return cx(cosf(phase), sinf(phase)); }

static inline int my_roundf(float x) {
    return (int)(x > 0 ? (int)(x + 0.5f) : (int)ceilf(x - 0.5f));
}

static int most_frequent(const int *arr, int n) {
    std::unordered_map<int, int> hash;
    for (int i = 0; i < n; i++) hash[arr[i]]++;
    int max_count = 0, res = -1;
    for (auto &p : hash)
        if (p.second > max_count) { res = p.first; max_count = p.second; }
    return res;
}

// Build modulated upchirp
static void build_upchirp(cx *chirp, uint32_t id, uint8_t sf, uint8_t os_factor = 1) {
    double N = (1 << sf);
    int n_fold = (int)(N * os_factor - id * os_factor);
    for (int n = 0; n < (int)(N * os_factor); n++) {
        double phase;
        if (n < n_fold)
            phase = 2.0 * M_PI * (n * n / (2.0 * N) / (os_factor * os_factor) + ((double)id / N - 0.5) * n / os_factor);
        else
            phase = 2.0 * M_PI * (n * n / (2.0 * N) / (os_factor * os_factor) + ((double)id / N - 1.5) * n / os_factor);
        chirp[n] = cx((float)cos(phase), (float)sin(phase));
    }
}

// Build reference up/down chirps
static void build_ref_chirps(cx *upchirp, cx *downchirp, uint8_t sf, uint8_t os_factor = 1) {
    int N = (1 << sf) * os_factor;
    build_upchirp(upchirp, 0, sf, os_factor);
    for (int i = 0; i < N; i++) downchirp[i] = std::conj(upchirp[i]);
}

// Complex multiply: out = a * b, element-wise
static void cx_multiply(cx *out, const cx *a, const cx *b, int n) {
    for (int i = 0; i < n; i++) out[i] = a[i] * b[i];
}

// FFT magnitude squared, returns argmax
static uint32_t fft_argmax(const cx *samples, const cx *ref_chirp, int N,
                           kiss_fft_cfg cfg, float *energy_out = nullptr) {
    std::vector<cx> dechirped(N);
    cx_multiply(dechirped.data(), samples, ref_chirp, N);

    kiss_fft_cpx *cx_in  = new kiss_fft_cpx[N];
    kiss_fft_cpx *cx_out = new kiss_fft_cpx[N];
    for (int i = 0; i < N; i++) { cx_in[i].r = dechirped[i].real(); cx_in[i].i = dechirped[i].imag(); }
    kiss_fft(cfg, cx_in, cx_out);

    float max_val = -1;
    uint32_t max_idx = 0;
    float tot_en = 0;
    for (int i = 0; i < N; i++) {
        float mag = cx_out[i].r * cx_out[i].r + cx_out[i].i * cx_out[i].i;
        tot_en += mag;
        if (mag > max_val) { max_val = mag; max_idx = i; }
    }
    if (energy_out) *energy_out = tot_en;

    delete[] cx_in;
    delete[] cx_out;
    return (tot_en > 0) ? max_idx : (uint32_t)-1;
}

// ============================================================================
// Hamming Decoder (soft + hard)
// ============================================================================

static const uint8_t cw_LUT[16]     = {0, 23, 45, 58, 78, 89, 99, 116, 139, 156, 166, 177, 197, 210, 232, 255};
static const uint8_t cw_LUT_cr5[16] = {0, 24, 40, 48, 72, 80, 96, 120, 136, 144, 160, 184, 192, 216, 232, 240};

static uint8_t hamming_decode_soft(const double *codeword_LLR, int cr_app) {
    int cw_len = cr_app + 4;
    double cw_proba[16] = {};
    for (int n = 0; n < 16; n++) {
        for (int j = 0; j < cw_len; j++) {
            bool bit = (((cr_app != 1) ? cw_LUT[n] : cw_LUT_cr5[n]) >> (8 - cw_len)) & (1u << (cw_len - 1 - j));
            if ((bit && codeword_LLR[j] > 0) || (!bit && codeword_LLR[j] < 0))
                cw_proba[n] += fabs(codeword_LLR[j]);
            else
                cw_proba[n] -= fabs(codeword_LLR[j]);
        }
    }
    int idx_max = (int)(std::max_element(cw_proba, cw_proba + 16) - cw_proba);
    uint8_t data = cw_LUT[idx_max] >> 4;
    // reverse bit order
    return ((data & 1) << 3) | (((data >> 1) & 1) << 2) | (((data >> 2) & 1) << 1) | ((data >> 3) & 1);
}

static uint8_t hamming_decode_hard(uint8_t cw_byte, int cr_app) {
    int cw_len = cr_app + 4;
    std::vector<bool> codeword(cw_len);
    for (int i = 0; i < cw_len; i++)
        codeword[i] = (cw_byte >> (cw_len - 1 - i)) & 1;

    std::vector<bool> data_nibble = {codeword[3], codeword[2], codeword[1], codeword[0]};
    bool s0, s1, s2;

    switch (cr_app) {
    case 4:
        if (!(std::count(codeword.begin(), codeword.end(), true) % 2))
            break;
        // fallthrough
    case 3:
        s0 = codeword[0] ^ codeword[1] ^ codeword[2] ^ codeword[4];
        s1 = codeword[1] ^ codeword[2] ^ codeword[3] ^ codeword[5];
        s2 = codeword[0] ^ codeword[1] ^ codeword[3] ^ codeword[6];
        {
            int syndrom = s0 + (s1 << 1) + (s2 << 2);
            switch (syndrom) {
            case 5: data_nibble[3].flip(); break;
            case 7: data_nibble[2].flip(); break;
            case 3: data_nibble[1].flip(); break;
            case 6: data_nibble[0].flip(); break;
            default: break;
            }
        }
        break;
    case 2:
        break;
    case 1:
        break;
    }
    uint32_t val = 0;
    for (int i = 0; i < 4; i++) val = (val << 1) | data_nibble[i];
    return (uint8_t)val;
}

// ============================================================================
// CRC-16 CCITT
// ============================================================================

static uint16_t crc16(const uint8_t *data, uint32_t len) {
    uint16_t crc = 0x0000;
    for (uint32_t i = 0; i < len; i++) {
        uint8_t byte = data[i];
        for (int b = 0; b < 8; b++) {
            if (((crc & 0x8000) >> 8) ^ (byte & 0x80))
                crc = (crc << 1) ^ 0x1021;
            else
                crc = (crc << 1);
            byte <<= 1;
        }
    }
    return crc;
}

// ============================================================================
// LoRa Receiver Configuration
// ============================================================================

struct LoRaConfig {
    uint32_t freq           = 869618000;
    uint32_t samp_rate      = 250000;
    uint32_t bw             = 62500;
    uint8_t  sf             = 8;
    uint8_t  cr             = 4;    // coding rate 1-4 (maps to 4/5 .. 4/8)
    uint16_t sync_word      = 0x12;
    bool     impl_head      = false;
    bool     has_crc        = true;
    uint32_t pay_len        = 11;
    int      gain           = 490;  // in 0.1 dB
    int      ppm            = -3;
    uint16_t preamble_len   = 8;
    bool     soft_decoding  = true;

    // derived
    uint32_t n_bins;
    uint32_t samples_per_symbol;
    uint8_t  os_factor;
    uint16_t sync_words[2];

    void compute_derived() {
        os_factor = samp_rate / bw;
        n_bins = 1u << sf;
        samples_per_symbol = n_bins * os_factor;
        sync_words[0] = ((sync_word & 0xF0) >> 4) << 3;
        sync_words[1] = (sync_word & 0x0F) << 3;
    }
};

// ============================================================================
// LoRa Demodulator - complete pipeline
// ============================================================================

class LoRaDemodulator {
public:
    LoRaDemodulator(const LoRaConfig &cfg) : cfg(cfg) {
        N = cfg.n_bins;
        os = cfg.os_factor;
        sps = cfg.samples_per_symbol;

        upchirp.resize(N);
        downchirp.resize(N);
        build_ref_chirps(upchirp.data(), downchirp.data(), cfg.sf);

        fft_cfg = kiss_fft_alloc(N, 0, nullptr, nullptr);

        n_up_req = cfg.preamble_len - 3;
        up_symb_to_use = n_up_req - 1;
        preamb_up_vals.resize(n_up_req, 0);
        preamble_raw.resize(n_up_req * N);
        preamble_raw_up.resize((n_up_req + 3) * sps);
        preamble_upchirps.resize(n_up_req * N);
        net_id_samp.resize((int)(sps * 2.5));
        additional_symbol_samp.resize(2 * sps);
        in_down.resize(N);
        symb_corr.resize(N);
        CFO_frac_correc.resize(N);

        reset_state();
    }

    ~LoRaDemodulator() {
        kiss_fft_free(fft_cfg);
    }

    // Feed oversampled IQ samples. This processes one symbol-worth at a time.
    // Returns true if a complete packet was decoded since last call.
    // Decoded packets are in last_packets.
    struct Packet {
        std::vector<uint8_t> payload;
        bool crc_valid;
        bool has_crc;
        float snr_est;
    };

    std::vector<Packet> last_packets;

    void process_samples(const cx *samples, int n_samples) {
        // Append to ring buffer
        ring.insert(ring.end(), samples, samples + n_samples);

        while ((int)ring.size() >= (int)sps + 2 * os) {
            // Downsample: pick one sample per os period, corrected for STO
            for (uint32_t ii = 0; ii < N; ii++) {
                int idx = (int)(os / 2 + os * ii - my_roundf(m_sto_frac * os));
                if (idx < 0) idx = 0;
                if (idx >= (int)ring.size()) { return; } // need more data
                in_down[ii] = ring[idx];
            }

            int items_to_consume = sps;

            switch (state) {
            case DETECT: {
                uint32_t bin_idx_new = fft_argmax(in_down.data(), downchirp.data(), N, fft_cfg);

                if (std::abs((long)mod(std::abs((long)bin_idx_new - bin_idx) + 1, N) - 1) <= 1 && bin_idx_new != (uint32_t)-1) {
                    if (symbol_cnt == 1 && bin_idx != (uint32_t)-1)
                        preamb_up_vals[0] = bin_idx;
                    preamb_up_vals[symbol_cnt] = bin_idx_new;
                    memcpy(&preamble_raw[N * symbol_cnt], in_down.data(), N * sizeof(cx));
                    memcpy(&preamble_raw_up[sps * symbol_cnt], &ring[os / 2], sps * sizeof(cx));
                    symbol_cnt++;
                } else {
                    memcpy(&preamble_raw[0], in_down.data(), N * sizeof(cx));
                    memcpy(&preamble_raw_up[0], &ring[os / 2], sps * sizeof(cx));
                    symbol_cnt = 1;
                }
                bin_idx = bin_idx_new;

                if (symbol_cnt == n_up_req) {
                    additional_upchirps = 0;
                    state = SYNC;
                    symbol_cnt = 0;
                    cfo_frac_sto_frac_est = false;
                    k_hat = most_frequent((const int *)preamb_up_vals.data(), preamb_up_vals.size());

                    int copy_start = (int)(0.75 * sps - k_hat * os);
                    if (copy_start >= 0 && copy_start + (int)(0.25 * sps) <= (int)ring.size())
                        memcpy(net_id_samp.data(), &ring[copy_start], (int)(0.25 * sps) * sizeof(cx));

                    items_to_consume = os * (N - k_hat);
                }
                break;
            }

            case SYNC: {
                if (!cfo_frac_sto_frac_est) {
                    m_cfo_frac = estimate_CFO_frac_Bernier(&preamble_raw[N - k_hat]);
                    m_sto_frac = estimate_STO_frac();
                    for (uint32_t n = 0; n < N; n++)
                        CFO_frac_correc[n] = expj(-2.0f * (float)M_PI * m_cfo_frac / N * n);
                    cfo_frac_sto_frac_est = true;
                }

                items_to_consume = sps;
                cx_multiply(symb_corr.data(), in_down.data(), CFO_frac_correc.data(), N);
                bin_idx = fft_argmax(symb_corr.data(), downchirp.data(), N, fft_cfg);

                switch (symbol_cnt) {
                case NET_ID1:
                    if (bin_idx == 0 || bin_idx == 1 || bin_idx == N - 1) {
                        // Additional upchirp, not net ID yet
                        if (additional_upchirps < 3) additional_upchirps++;
                    } else {
                        symbol_cnt = NET_ID2;
                        net_ids[0] = bin_idx;
                        if ((int)(0.25 * sps) + (int)sps <= (int)ring.size())
                            memcpy(&net_id_samp[(int)(0.25 * sps)], ring.data(), sps * sizeof(cx));
                    }
                    break;
                case NET_ID2:
                    symbol_cnt = DOWNCHIRP1;
                    net_ids[1] = bin_idx;
                    break;
                case DOWNCHIRP1:
                    symbol_cnt = DOWNCHIRP2;
                    break;
                case DOWNCHIRP2:
                    down_val = fft_argmax(symb_corr.data(), upchirp.data(), N, fft_cfg);
                    memcpy(additional_symbol_samp.data(), ring.data(), sps * sizeof(cx));
                    symbol_cnt = QUARTER_DOWN;
                    break;
                case QUARTER_DOWN: {
                    if ((int)(sps + sps) <= (int)ring.size())
                        memcpy(&additional_symbol_samp[sps], ring.data(), sps * sizeof(cx));

                    if ((uint32_t)down_val < N / 2)
                        m_cfo_int = (int)floor((double)down_val / 2.0);
                    else
                        m_cfo_int = (int)floor((double)(down_val - (int)N) / 2.0);

                    // Validate sync word
                    bool sync_ok = true;
                    if (cfg.sync_words[0] != 0) {
                        // We don't have perfectly corrected net_ids here - use a tolerance check
                        if (std::abs((int)net_ids[0] - (int)cfg.sync_words[0]) > 2) {
                            sync_ok = false;
                        }
                    }

                    if (!sync_ok) {
                        reset_state();
                        break;
                    }

                    // Estimate SFO
                    sfo_hat = (float)(m_cfo_int + m_cfo_frac) * (float)cfg.bw / (float)cfg.freq;
                    m_sto_frac += sfo_hat * cfg.preamble_len;
                    if (fabsf(m_sto_frac) > 0.5f)
                        m_sto_frac += (m_sto_frac > 0 ? -1.0f : 1.0f);
                    m_sto_frac += sfo_hat * 4.25f;
                    sfo_cum = ((m_sto_frac * os) - my_roundf(m_sto_frac * os)) / (float)os;

                    // Build demod downchirp with CFO correction
                    demod_upchirp.resize(N);
                    demod_downchirp.resize(N);
                    build_upchirp(demod_upchirp.data(), mod(m_cfo_int, N), cfg.sf);
                    for (uint32_t i = 0; i < N; i++) demod_downchirp[i] = std::conj(demod_upchirp[i]);
                    for (uint32_t n = 0; n < N; n++)
                        demod_downchirp[n] *= expj(-2.0f * (float)M_PI * m_cfo_frac / N * n);

                    // Estimate SNR from a preamble symbol
                    float snr_e = 0;
                    if (preamble_raw.size() >= N) {
                        float sig_e = 0, tot_e = 0;
                        std::vector<float> fft_mag(N);
                        std::vector<cx> dc(N);
                        cx_multiply(dc.data(), preamble_raw.data(), downchirp.data(), N);
                        kiss_fft_cpx *ci = new kiss_fft_cpx[N];
                        kiss_fft_cpx *co = new kiss_fft_cpx[N];
                        for (uint32_t i = 0; i < N; i++) { ci[i].r = dc[i].real(); ci[i].i = dc[i].imag(); }
                        kiss_fft(fft_cfg, ci, co);
                        for (uint32_t i = 0; i < N; i++) {
                            fft_mag[i] = co[i].r * co[i].r + co[i].i * co[i].i;
                            tot_e += fft_mag[i];
                        }
                        int mi = (int)(std::max_element(fft_mag.begin(), fft_mag.end()) - fft_mag.begin());
                        sig_e = fft_mag[mi];
                        if (tot_e > sig_e) snr_e = 10.0f * log10f(sig_e / (tot_e - sig_e));
                        delete[] ci;
                        delete[] co;
                    }
                    current_snr = snr_e;

                    state = PAYLOAD;
                    payload_symbol_cnt = 0;
                    m_received_head = false;
                    m_ldro = false;
                    frame_cnt++;

                    // Reset payload decode state
                    fft_block.clear();
                    nibbles.clear();
                    is_header = true;
                    m_pay_cr = cfg.cr;
                    m_pay_len = cfg.pay_len;
                    m_pay_has_crc = cfg.has_crc;

                    items_to_consume = sps / 4 + os * m_cfo_int;
                    if (items_to_consume < 0) items_to_consume = 0;

                    fprintf(stderr, "[sync] Frame #%u detected, CFO=%.2f, SNR=%.1f dB\n",
                            frame_cnt, m_cfo_int + m_cfo_frac, snr_e);
                    break;
                }
                default:
                    break;
                }
                break;
            }

            case PAYLOAD: {
                if (payload_symbol_cnt < 8 || ((uint32_t)payload_symbol_cnt < total_payload_symbols && m_received_head)) {
                    // Output downsampled symbol for demodulation
                    demodulate_symbol(in_down.data());
                    items_to_consume = sps;

                    // SFO compensation
                    if (fabsf(sfo_cum) > 1.0f / 2.0f / os) {
                        int sign_val = (sfo_cum >= 0) ? 1 : -1;
                        items_to_consume -= sign_val;
                        sfo_cum -= sign_val * 1.0f / os;
                    }
                    sfo_cum += sfo_hat;
                    payload_symbol_cnt++;
                } else if (!m_received_head) {
                    // Waiting for header decode
                    items_to_consume = 0;
                } else {
                    // Frame done
                    flush_remaining();
                    reset_state();
                }
                break;
            }
            }

            // Consume processed samples
            if (items_to_consume > 0 && items_to_consume <= (int)ring.size())
                ring.erase(ring.begin(), ring.begin() + items_to_consume);
            else if (items_to_consume <= 0)
                break; // waiting state, don't consume
            else
                break; // not enough data
        }
    }

private:
    LoRaConfig cfg;
    uint32_t N;    // number of bins = 2^sf
    uint8_t os;    // oversampling factor
    uint32_t sps;  // samples per symbol

    std::vector<cx> upchirp, downchirp;
    std::vector<cx> demod_upchirp, demod_downchirp;
    kiss_fft_cfg fft_cfg;

    // Ring buffer for incoming samples
    std::vector<cx> ring;

    // Frame sync state
    enum State { DETECT, SYNC, PAYLOAD };
    enum SyncState { NET_ID1 = 0, NET_ID2, DOWNCHIRP1, DOWNCHIRP2, QUARTER_DOWN };

    State state;
    int symbol_cnt;
    uint32_t bin_idx;
    int n_up_req;
    int up_symb_to_use;
    std::vector<int> preamb_up_vals;
    std::vector<cx> preamble_raw;
    std::vector<cx> preamble_raw_up;
    std::vector<cx> preamble_upchirps;
    std::vector<cx> net_id_samp;
    std::vector<cx> additional_symbol_samp;
    std::vector<cx> in_down;
    std::vector<cx> symb_corr;
    std::vector<cx> CFO_frac_correc;
    int net_ids[2];
    int k_hat;
    float m_cfo_frac;
    int m_cfo_int;
    float m_sto_frac;
    float sfo_hat;
    float sfo_cum;
    bool cfo_frac_sto_frac_est;
    uint8_t additional_upchirps;
    int down_val;
    unsigned int frame_cnt;
    float current_snr;

    // Payload demod state
    int payload_symbol_cnt;
    uint32_t total_payload_symbols;
    bool m_received_head;
    bool m_ldro;
    bool is_header;
    uint8_t m_pay_cr;
    uint32_t m_pay_len;
    bool m_pay_has_crc;

    // Demod pipeline buffers
    std::vector<std::vector<double>> fft_block; // LLR blocks for soft decoding
    std::vector<uint16_t> hard_block;           // symbol values for hard decoding
    std::vector<uint8_t> nibbles;               // decoded nibbles

    void reset_state() {
        state = DETECT;
        symbol_cnt = 1;
        bin_idx = 0;
        k_hat = 0;
        m_sto_frac = 0;
        m_cfo_frac = 0;
        m_cfo_int = 0;
        sfo_hat = 0;
        sfo_cum = 0;
        cfo_frac_sto_frac_est = false;
        additional_upchirps = 0;
        down_val = 0;
        payload_symbol_cnt = 0;
        total_payload_symbols = 0;
        m_received_head = false;
        m_ldro = false;
        is_header = true;
        fft_block.clear();
        hard_block.clear();
        nibbles.clear();
    }

    // ---- CFO estimation (Bernier's algorithm) ----
    float estimate_CFO_frac_Bernier(const cx *samples) {
        std::vector<int> k0v(up_symb_to_use);
        std::vector<double> k0_mag(up_symb_to_use);
        std::vector<cx> fft_val(up_symb_to_use * N);
        std::vector<cx> dechirped(N);

        kiss_fft_cpx *ci = new kiss_fft_cpx[N];
        kiss_fft_cpx *co = new kiss_fft_cpx[N];

        for (int i = 0; i < up_symb_to_use; i++) {
            cx_multiply(dechirped.data(), &samples[N * i], downchirp.data(), N);
            for (uint32_t j = 0; j < N; j++) { ci[j].r = dechirped[j].real(); ci[j].i = dechirped[j].imag(); }
            kiss_fft(fft_cfg, ci, co);

            float max_mag = -1;
            int max_idx = 0;
            for (uint32_t j = 0; j < N; j++) {
                float m = co[j].r * co[j].r + co[j].i * co[j].i;
                fft_val[j + i * N] = cx(co[j].r, co[j].i);
                if (m > max_mag) { max_mag = m; max_idx = j; }
            }
            k0v[i] = max_idx;
            k0_mag[i] = max_mag;
        }
        delete[] ci;
        delete[] co;

        int best = (int)(std::max_element(k0_mag.begin(), k0_mag.end()) - k0_mag.begin());
        int idx_max = k0v[best];

        cx four_cum(0, 0);
        for (int i = 0; i < up_symb_to_use - 1; i++)
            four_cum += fft_val[idx_max + N * i] * std::conj(fft_val[idx_max + N * (i + 1)]);

        float cfo_frac = -std::arg(four_cum) / (2.0f * (float)M_PI);

        // Correct CFO in preamble
        std::vector<cx> correc(up_symb_to_use * N);
        for (int n = 0; n < up_symb_to_use * (int)N; n++)
            correc[n] = expj(-2.0f * (float)M_PI * cfo_frac / N * n);
        cx_multiply(preamble_upchirps.data(), samples, correc.data(), up_symb_to_use * N);

        return cfo_frac;
    }

    // ---- STO fractional estimation (RCTSL) ----
    float estimate_STO_frac() {
        std::vector<float> fft_mag_sq(2 * N, 0);
        std::vector<cx> dechirped(N);

        kiss_fft_cpx *ci = new kiss_fft_cpx[2 * N];
        kiss_fft_cpx *co = new kiss_fft_cpx[2 * N];
        kiss_fft_cfg cfg2 = kiss_fft_alloc(2 * N, 0, nullptr, nullptr);

        for (int i = 0; i < up_symb_to_use; i++) {
            cx_multiply(dechirped.data(), &preamble_upchirps[N * i], downchirp.data(), N);
            for (uint32_t j = 0; j < 2 * N; j++) {
                if (j < N) { ci[j].r = dechirped[j].real(); ci[j].i = dechirped[j].imag(); }
                else { ci[j].r = 0; ci[j].i = 0; }
            }
            kiss_fft(cfg2, ci, co);
            for (uint32_t j = 0; j < 2 * N; j++)
                fft_mag_sq[j] += co[j].r * co[j].r + co[j].i * co[j].i;
        }
        free(cfg2);
        delete[] ci;
        delete[] co;

        int k0 = (int)(std::max_element(fft_mag_sq.begin(), fft_mag_sq.end()) - fft_mag_sq.begin());
        double Y_1 = fft_mag_sq[mod(k0 - 1, 2 * N)];
        double Y0 = fft_mag_sq[k0];
        double Y1 = fft_mag_sq[mod(k0 + 1, 2 * N)];

        double u = 64.0 * N / 406.5506497;
        double v = u * 2.4674;
        double wa = (Y1 - Y_1) / (u * (Y1 + Y_1) + v * Y0);
        double ka = wa * N / M_PI;
        double k_residual = fmod((k0 + ka) / 2.0, 1.0);
        float sto_frac = (float)(k_residual - (k_residual > 0.5 ? 1.0 : 0.0));
        return sto_frac;
    }

    // ---- FFT demodulation of one symbol ----
    void demodulate_symbol(const cx *symbol) {
        int block_size = 4 + (is_header ? 4 : m_pay_cr);

        if (cfg.soft_decoding) {
            // Compute LLRs
            std::vector<double> LLRs = compute_LLRs(symbol);
            fft_block.push_back(LLRs);

            if ((int)fft_block.size() == block_size) {
                // Deinterleave + Hamming decode the block
                decode_block_soft();
                fft_block.clear();
            }
        } else {
            // Hard demod
            uint32_t idx = fft_argmax(symbol, demod_downchirp.data(), N, fft_cfg);
            uint16_t val = (uint16_t)mod((long)idx - 1, (1L << cfg.sf));
            if (is_header || m_ldro) val /= 4;
            // Gray demap
            val = val ^ (val >> 1);
            hard_block.push_back(val);

            if ((int)hard_block.size() == block_size) {
                decode_block_hard();
                hard_block.clear();
            }
        }
    }

    // ---- Compute LLRs for one symbol (soft demod) ----
    std::vector<double> compute_LLRs(const cx *samples) {
        std::vector<cx> dechirped(N);
        cx_multiply(dechirped.data(), samples, demod_downchirp.data(), N);

        kiss_fft_cpx *ci = new kiss_fft_cpx[N];
        kiss_fft_cpx *co = new kiss_fft_cpx[N];
        for (uint32_t i = 0; i < N; i++) { ci[i].r = dechirped[i].real(); ci[i].i = dechirped[i].imag(); }
        kiss_fft(fft_cfg, ci, co);

        std::vector<float> fft_mag_sq(N);
        for (uint32_t i = 0; i < N; i++)
            fft_mag_sq[i] = co[i].r * co[i].r + co[i].i * co[i].i;

        delete[] ci;
        delete[] co;

        // SNR estimation
        int symbol_idx = (int)(std::max_element(fft_mag_sq.begin(), fft_mag_sq.end()) - fft_mag_sq.begin());
        double signal_energy = 0, noise_energy = 0;
        for (uint32_t i = 0; i < N; i++) {
            if (mod(std::abs((int)i - symbol_idx), N - 1) < 2)
                signal_energy += fft_mag_sq[i];
            else
                noise_energy += fft_mag_sq[i];
        }
        (void)signal_energy;
        (void)noise_energy;

        // Normalize FFT magnitudes
        for (uint32_t i = 0; i < N; i++)
            fft_mag_sq[i] *= N;

        // Use max-log approximation with |Y[n]|^2 directly (avoids Bessel overflow)
        std::vector<double> LLs(N);
        for (uint32_t n = 0; n < N; n++)
            LLs[n] = (double)fft_mag_sq[n];

        // Compute LLRs (max-log approximation)
        std::vector<double> LLRs(MAX_SF, 0);
        for (uint32_t i = 0; i < cfg.sf; i++) {
            double max_X1 = -1e30, max_X0 = -1e30;
            for (uint32_t n = 0; n < N; n++) {
                uint32_t s = (uint32_t)mod((long)n - 1, (1L << cfg.sf)) / ((is_header || m_ldro) ? 4 : 1);
                s = s ^ (s >> 1); // Gray demap
                if (s & (1u << i)) {
                    if (LLs[n] > max_X1) max_X1 = LLs[n];
                } else {
                    if (LLs[n] > max_X0) max_X0 = LLs[n];
                }
            }
            LLRs[cfg.sf - 1 - i] = max_X1 - max_X0;
        }
        return LLRs;
    }

    // ---- Decode a block (soft) ----
    void decode_block_soft() {
        int sf_app = (is_header || m_ldro) ? cfg.sf - 2 : cfg.sf;
        int cw_len = is_header ? 8 : m_pay_cr + 4;
        int cr_app = is_header ? 4 : m_pay_cr;

        // Deinterleave
        std::vector<std::vector<double>> inter_bin(cw_len, std::vector<double>(sf_app, 0));
        std::vector<std::vector<double>> deinter_bin(sf_app, std::vector<double>(cw_len, 0));

        for (int i = 0; i < cw_len; i++) {
            for (int j = 0; j < sf_app; j++)
                inter_bin[i][j] = fft_block[i][cfg.sf - sf_app + j];
        }

        for (int i = 0; i < cw_len; i++)
            for (int j = 0; j < sf_app; j++)
                deinter_bin[mod(i - j - 1, sf_app)][i] = inter_bin[i][j];

        // Hamming decode each codeword
        for (int i = 0; i < sf_app; i++) {
            uint8_t nibble = hamming_decode_soft(deinter_bin[i].data(), cr_app);
            nibbles.push_back(nibble);
        }

        // Try to decode header or process payload
        process_nibbles();
    }

    // ---- Decode a block (hard) ----
    void decode_block_hard() {
        int sf_app = (is_header || m_ldro) ? cfg.sf - 2 : cfg.sf;
        int cw_len = is_header ? 8 : m_pay_cr + 4;
        int cr_app = is_header ? 4 : m_pay_cr;

        // Deinterleave
        // Convert symbols to binary
        std::vector<std::vector<bool>> inter_bin(cw_len);
        for (int i = 0; i < cw_len; i++) {
            inter_bin[i].resize(sf_app);
            for (int j = 0; j < sf_app; j++)
                inter_bin[i][j] = (hard_block[i] >> (sf_app - 1 - j)) & 1;
        }

        std::vector<std::vector<bool>> deinter_bin(sf_app, std::vector<bool>(cw_len, false));
        for (int i = 0; i < cw_len; i++)
            for (int j = 0; j < sf_app; j++)
                deinter_bin[mod(i - j - 1, sf_app)][i] = inter_bin[i][j];

        // Convert back and hamming decode
        for (int i = 0; i < sf_app; i++) {
            uint8_t cw_byte = 0;
            for (int j = 0; j < cw_len; j++)
                cw_byte = (cw_byte << 1) | deinter_bin[i][j];
            uint8_t nibble = hamming_decode_hard(cw_byte, cr_app);
            nibbles.push_back(nibble);
        }

        process_nibbles();
    }

    // ---- Process accumulated nibbles (header decode / payload assembly) ----
    void process_nibbles() {
        if (is_header && nibbles.size() >= 5 && !cfg.impl_head) {
            // Decode explicit header
            m_pay_len = (nibbles[0] << 4) | nibbles[1];
            m_pay_has_crc = nibbles[2] & 1;
            m_pay_cr = nibbles[2] >> 1;

            if (m_pay_cr < 1) m_pay_cr = 1;
            if (m_pay_cr > 4) m_pay_cr = 4;

            uint8_t header_chk = ((nibbles[3] & 1) << 4) | nibbles[4];

            // Verify header checksum
            bool c4 = ((nibbles[0] >> 3) & 1) ^ ((nibbles[0] >> 2) & 1) ^ ((nibbles[0] >> 1) & 1) ^ (nibbles[0] & 1);
            bool c3 = ((nibbles[0] >> 3) & 1) ^ ((nibbles[1] >> 3) & 1) ^ ((nibbles[1] >> 2) & 1) ^ ((nibbles[1] >> 1) & 1) ^ (nibbles[2] & 1);
            bool c2 = ((nibbles[0] >> 2) & 1) ^ ((nibbles[1] >> 3) & 1) ^ (nibbles[1] & 1) ^ ((nibbles[2] >> 3) & 1) ^ ((nibbles[2] >> 1) & 1);
            bool c1 = ((nibbles[0] >> 1) & 1) ^ ((nibbles[1] >> 2) & 1) ^ (nibbles[1] & 1) ^ ((nibbles[2] >> 2) & 1) ^ ((nibbles[2] >> 1) & 1) ^ (nibbles[2] & 1);
            bool c0 = (nibbles[0] & 1) ^ ((nibbles[1] >> 1) & 1) ^ ((nibbles[2] >> 3) & 1) ^ ((nibbles[2] >> 2) & 1) ^ ((nibbles[2] >> 1) & 1) ^ (nibbles[2] & 1);

            int computed_chk = (c4 << 4) | (c3 << 3) | (c2 << 2) | (c1 << 1) | c0;

            fprintf(stderr, "  Header: pay_len=%u, CRC=%d, CR=%d\n", m_pay_len, m_pay_has_crc, m_pay_cr);

            if (computed_chk != header_chk || m_pay_len == 0) {
                fprintf(stderr, "  Header checksum INVALID!\n");
                reset_state();
                return;
            }

            fprintf(stderr, "  Header checksum valid\n");

            // Determine LDRO
            m_ldro = ((float)(1u << cfg.sf) * 1e3f / cfg.bw) > LDRO_MAX_DURATION_MS;

            // Calculate total payload symbols
            int sf_minus_2ldro = cfg.sf - 2 * m_ldro;
            total_payload_symbols = 8 + (int)ceil((double)(2 * m_pay_len - cfg.sf + 2 + (!cfg.impl_head ? 5 : 0) + (m_pay_has_crc ? 4 : 0)) / sf_minus_2ldro) * (4 + m_pay_cr);

            m_received_head = true;
            is_header = false;

            // Remove header nibbles, keep remaining payload nibbles
            std::vector<uint8_t> remaining(nibbles.begin() + 5, nibbles.end());
            nibbles = remaining;

        } else if (is_header && cfg.impl_head) {
            // Implicit header - parameters from config
            m_ldro = ((float)(1u << cfg.sf) * 1e3f / cfg.bw) > LDRO_MAX_DURATION_MS;
            int sf_minus_2ldro = cfg.sf - 2 * m_ldro;
            total_payload_symbols = 8 + (int)ceil((double)(2 * m_pay_len - cfg.sf + 2 + (m_pay_has_crc ? 4 : 0)) / sf_minus_2ldro) * (4 + m_pay_cr);
            m_received_head = true;
            is_header = false;
        }

        // Check if we have enough nibbles for the full payload + CRC
        uint32_t needed_nibbles = m_pay_len * 2 + (m_pay_has_crc ? 4 : 0);
        if (!is_header && nibbles.size() >= needed_nibbles) {
            assemble_and_output();
        }
    }

    // ---- Dewhiten + CRC check + output ----
    void assemble_and_output() {
        uint32_t total_bytes = m_pay_len + (m_pay_has_crc ? 2 : 0);
        std::vector<uint8_t> bytes(total_bytes);

        for (uint32_t i = 0; i < total_bytes; i++) {
            uint8_t low_nib = nibbles[2 * i];
            uint8_t high_nib = nibbles[2 * i + 1];

            if (i < m_pay_len) {
                // Dewhiten payload
                low_nib ^= (whitening_seq[i] & 0x0F);
                high_nib ^= (whitening_seq[i] >> 4) & 0x0F;
            }
            // CRC bytes are NOT dewhitened

            bytes[i] = (high_nib << 4) | low_nib;
        }

        Packet pkt;
        pkt.payload.assign(bytes.begin(), bytes.begin() + m_pay_len);
        pkt.has_crc = m_pay_has_crc;
        pkt.snr_est = current_snr;
        pkt.crc_valid = false;

        if (m_pay_has_crc && m_pay_len >= 2) {
            uint16_t calc_crc = crc16(pkt.payload.data(), m_pay_len - 2);
            calc_crc = calc_crc ^ pkt.payload[m_pay_len - 1] ^ ((uint16_t)pkt.payload[m_pay_len - 2] << 8);
            uint16_t rx_crc = bytes[m_pay_len] | ((uint16_t)bytes[m_pay_len + 1] << 8);
            pkt.crc_valid = (calc_crc == rx_crc);
        } else if (!m_pay_has_crc) {
            pkt.crc_valid = true; // no CRC to check
        }

        // Print output
        printf("rx msg: ");
        for (uint32_t i = 0; i < m_pay_len; i++) {
            if (i > 0) printf(", ");
            printf("0x%02x", pkt.payload[i]);
        }
        printf("\n");

        // Also print as ASCII if printable
        bool printable = true;
        for (uint32_t i = 0; i < m_pay_len; i++)
            if (pkt.payload[i] < 0x20 || pkt.payload[i] > 0x7e) { printable = false; break; }
        if (printable) {
            printf("rx str: ");
            for (uint32_t i = 0; i < m_pay_len; i++) putchar(pkt.payload[i]);
            printf("\n");
        }

        if (pkt.has_crc) {
            if (pkt.crc_valid)
                printf("CRC valid!\n\n");
            else
                printf("\033[31mCRC invalid\033[0m\n\n");
        } else {
            printf("(no CRC)\n\n");
        }

        fflush(stdout);
        last_packets.push_back(pkt);

        nibbles.clear();
    }

    // Called when frame ends to flush any partial data
    void flush_remaining() {
        uint32_t needed = m_pay_len * 2 + (m_pay_has_crc ? 4 : 0);
        if (!is_header && nibbles.size() >= needed)
            assemble_and_output();
    }
};

// ============================================================================
// RTL-SDR async callback
// ============================================================================

struct RtlContext {
    LoRaDemodulator *demod;
    int n_samples;
};

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *ctx) {
    if (!g_running) return;
    RtlContext *rc = (RtlContext *)ctx;

    int n = len / 2;
    std::vector<cx> samples(n);
    for (int i = 0; i < n; i++) {
        float re = ((float)buf[2 * i]     - 127.5f) / 127.5f;
        float im = ((float)buf[2 * i + 1] - 127.5f) / 127.5f;
        samples[i] = cx(re, im);
    }

    rc->demod->process_samples(samples.data(), n);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char *argv[]) {
    LoRaConfig cfg;
    int opt;

    while ((opt = getopt(argc, argv, "f:s:b:S:c:w:g:p:IL:")) != -1) {
        switch (opt) {
        case 'f': cfg.freq = (uint32_t)atol(optarg); break;
        case 's': cfg.samp_rate = (uint32_t)atol(optarg); break;
        case 'b': cfg.bw = (uint32_t)atol(optarg); break;
        case 'S': cfg.sf = (uint8_t)atoi(optarg); break;
        case 'c': cfg.cr = (uint8_t)atoi(optarg); break;
        case 'w': cfg.sync_word = (uint16_t)strtol(optarg, nullptr, 16); break;
        case 'g': cfg.gain = atoi(optarg); break;
        case 'p': cfg.ppm = atoi(optarg); break;
        case 'I': cfg.impl_head = true; break;
        case 'L': cfg.pay_len = (uint32_t)atoi(optarg); break;
        default:
            fprintf(stderr, "Usage: %s [-f freq] [-s samp_rate] [-b bw] [-S sf] [-c cr] [-w sync_word_hex] [-g gain] [-p ppm] [-I] [-L pay_len]\n", argv[0]);
            return 1;
        }
    }

    if (cfg.sf < MIN_SF || cfg.sf > MAX_SF) {
        fprintf(stderr, "SF must be between %d and %d\n", MIN_SF, MAX_SF);
        return 1;
    }
    if (cfg.samp_rate % cfg.bw != 0) {
        fprintf(stderr, "Warning: samp_rate should be an integer multiple of bandwidth\n");
    }

    cfg.compute_derived();

    fprintf(stderr, "LoRa Standalone RX\n");
    fprintf(stderr, "  Freq:       %u Hz\n", cfg.freq);
    fprintf(stderr, "  Samp rate:  %u Hz\n", cfg.samp_rate);
    fprintf(stderr, "  Bandwidth:  %u Hz\n", cfg.bw);
    fprintf(stderr, "  SF:         %u\n", cfg.sf);
    fprintf(stderr, "  CR:         4/%u\n", cfg.cr + 4);
    fprintf(stderr, "  Sync word:  0x%02X\n", cfg.sync_word);
    fprintf(stderr, "  OS factor:  %u\n", cfg.os_factor);
    fprintf(stderr, "  Bins:       %u\n", cfg.n_bins);
    fprintf(stderr, "  Samples/sym: %u\n", cfg.samples_per_symbol);
    fprintf(stderr, "  Implicit:   %s\n", cfg.impl_head ? "yes" : "no");
    fprintf(stderr, "  Soft decode: %s\n", cfg.soft_decoding ? "yes" : "no");
    fprintf(stderr, "\n");

    // Open RTL-SDR
    rtlsdr_dev_t *dev = nullptr;
    int dev_index = 0;

    int n_devices = rtlsdr_get_device_count();
    if (n_devices == 0) {
        fprintf(stderr, "No RTL-SDR devices found\n");
        return 1;
    }
    fprintf(stderr, "Found %d RTL-SDR device(s)\n", n_devices);

    if (rtlsdr_open(&dev, dev_index) < 0) {
        fprintf(stderr, "Failed to open RTL-SDR device\n");
        return 1;
    }
    g_dev = dev;

    // Configure
    rtlsdr_set_sample_rate(dev, cfg.samp_rate);
    rtlsdr_set_center_freq(dev, cfg.freq);
    rtlsdr_set_freq_correction(dev, cfg.ppm);
    rtlsdr_set_tuner_gain_mode(dev, 1); // manual gain
    rtlsdr_set_tuner_gain(dev, cfg.gain);
    rtlsdr_set_agc_mode(dev, 0);
    rtlsdr_reset_buffer(dev);

    fprintf(stderr, "RTL-SDR configured. Tuner gain: %.1f dB\n", rtlsdr_get_tuner_gain(dev) / 10.0);
    fprintf(stderr, "Listening...\n\n");

    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    LoRaDemodulator demod(cfg);
    RtlContext ctx = { &demod, 0 };

    // Use async read with a buffer size that's a good multiple of samples_per_symbol
    // 16384 bytes = 8192 IQ samples
    int buf_len = 16384;

    rtlsdr_read_async(dev, rtlsdr_callback, &ctx, 0, buf_len);

    rtlsdr_close(dev);

    fprintf(stderr, "\nDone. Received %zu packets.\n", demod.last_packets.size());
    return 0;
}
