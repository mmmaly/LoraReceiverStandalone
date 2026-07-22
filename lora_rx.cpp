// Standalone LoRa receiver - pure C/C++ reimplementation of the gr-lora_sdr
// demodulation chain. No GNU Radio dependency.
//
// Pipeline: RTL-SDR IQ (or IQ file) -> Anti-alias FIR + decimate ->
//           Frame Sync -> FFT Demod -> Gray Demap -> Deinterleave ->
//           Hamming Decode -> Header Decode -> Dewhiten -> CRC Verify -> Output
//
// Multiple demodulators (one per SF/BW combination) can run in parallel on
// the same sample stream (-A auto-scan mode), each on its own worker thread.
//
// Build:
//   make
//
// Usage:
//   ./lora_rx [options]
//     -f <freq_hz>      Center frequency (default 869618000)
//     -s <samp_rate>    Sample rate (default 250000)
//     -b <bandwidth>    LoRa bandwidth (default 62500)
//     -S <sf>           Spreading factor 5-12 (default 8)
//     -c <cr>           Coding rate 1-4 (default 4; explicit header overrides)
//     -w <sync_word>    Sync word hex (default 0x12, 0 = accept any)
//     -g <gain>         RTL-SDR gain in 0.1 dB (default 490)
//     -p <ppm>          Frequency correction ppm (default -3)
//     -I                Implicit header mode
//     -L <pay_len>      Payload length for implicit header (default 11)
//     -A                Auto-scan: run all SF (7-12) and BW combinations in parallel
//     -r <file>         Read IQ from file (rtl_sdr u8 format, "-" = stdin) instead of SDR
//     -D <file>         Dump raw IQ to file while receiving (for later replay with -r)

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <csignal>
#include <complex>
#include <vector>
#include <array>
#include <deque>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <unordered_map>
#include <string>
#include <getopt.h>
#include <rtl-sdr.h>

extern "C" {
#include "kiss_fft.h"
}

#include "lora_common.h"

// ============================================================================
// Types and constants
// ============================================================================

static constexpr int MIN_SF = 5;
static constexpr int MAX_SF = 12;
static constexpr float LDRO_MAX_DURATION_MS = 16.0f;

static volatile bool g_running = true;
static rtlsdr_dev_t *g_dev = nullptr;
static std::mutex g_print_mtx;

static void signal_handler(int) {
    g_running = false;
    if (g_dev) rtlsdr_cancel_async(g_dev);
}

static_assert(sizeof(cx) == sizeof(kiss_fft_cpx),
              "std::complex<float> must be layout-compatible with kiss_fft_cpx");

// ============================================================================
// Utility functions
// ============================================================================

static int most_frequent(const int *arr, int n) {
    std::unordered_map<int, int> hash;
    for (int i = 0; i < n; i++) hash[arr[i]]++;
    int max_count = 0, res = -1;
    for (auto &p : hash)
        if (p.second > max_count) { res = p.first; max_count = p.second; }
    return res;
}

// ============================================================================
// Hamming Decoder (soft + hard)
// ============================================================================

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
// LoRa Receiver Configuration
// ============================================================================

struct LoRaConfig {
    uint32_t freq           = 869618000;
    uint32_t samp_rate      = 250000;
    uint32_t bw             = 62500;
    uint8_t  sf             = 8;
    uint8_t  cr             = 4;    // coding rate 1-4 (maps to 4/5 .. 4/8)
    uint16_t sync_word      = 0x12; // 0 = accept any
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
    struct Packet {
        std::vector<uint8_t> payload;
        bool crc_valid;
        bool has_crc;
        float snr_est;
    };

    std::vector<Packet> last_packets;

    LoRaDemodulator(const LoRaConfig &cfg_, std::string tag_ = "")
        : cfg(cfg_), tag(std::move(tag_)) {
        N = cfg.n_bins;
        os = cfg.os_factor;
        sps = cfg.samples_per_symbol;

        upchirp.resize(N);
        downchirp.resize(N);
        build_ref_chirps(upchirp.data(), downchirp.data(), cfg.sf);

        fft_cfg  = kiss_fft_alloc(N, 0, nullptr, nullptr);
        fft_cfg2 = kiss_fft_alloc(2 * N, 0, nullptr, nullptr);

        // Preallocated FFT scratch (2N covers both the N-point and the
        // zero-padded 2N-point transforms used by the STO estimator)
        dechirp_buf.resize(2 * N);
        fft_out.resize(2 * N);
        fft_mag.resize(2 * N);

        n_up_req = cfg.preamble_len - 3;
        up_symb_to_use = n_up_req - 1;
        preamb_up_vals.resize(n_up_req, 0);
        preamble_raw.resize(n_up_req * N);
        preamble_upchirps.resize(n_up_req * N);
        cfo_fft_val.resize(up_symb_to_use * N);
        in_down.resize(N);
        symb_corr.resize(N);
        CFO_frac_correc.resize(N);

        // Gray demap LUTs: FFT bin -> data symbol, for payload and for
        // header/LDRO (extra /4) modes
        demap_pay.resize(N);
        demap_hdr.resize(N);
        for (uint32_t n = 0; n < N; n++) {
            uint32_t v = (n + N - 1) % N;   // mod(n - 1, N)
            demap_pay[n] = (uint16_t)(v ^ (v >> 1));
            uint32_t h = v / 4;
            demap_hdr[n] = (uint16_t)(h ^ (h >> 1));
        }

        // Anti-alias low-pass FIR applied before decimation (only needed when
        // oversampled). Cutoff slightly above bw/2 so band edges pass intact.
        if (os > 1) {
            int T = 8 * os + 1;
            fir_taps.resize(T);
            fir_delay = T / 2;
            double fc = 0.55 / os;  // as fraction of full sample rate
            double sum = 0;
            for (int k = 0; k < T; k++) {
                double x = k - (T - 1) / 2.0;
                double s = (x == 0) ? 2.0 * M_PI * fc : sin(2.0 * M_PI * fc * x) / x;
                double w = 0.54 - 0.46 * cos(2.0 * M_PI * k / (T - 1));
                fir_taps[k] = (float)(s * w);
                sum += s * w;
            }
            for (auto &t : fir_taps) t = (float)(t / sum);
        } else {
            fir_delay = 0;
        }

        reset_state();
    }

    ~LoRaDemodulator() {
        kiss_fft_free(fft_cfg);
        kiss_fft_free(fft_cfg2);
    }

    // Feed oversampled IQ samples; processes one symbol-worth at a time.
    void process_samples(const cx *samples, size_t n_samples) {
        ring.insert(ring.end(), samples, samples + n_samples);

        const size_t need = sps + 2 * os + fir_delay;
        while (avail() >= need) {
            // Downsample through the anti-alias filter: one output sample per
            // os period, corrected for fractional STO
            for (uint32_t ii = 0; ii < N; ii++) {
                long idx = (long)(os / 2) + (long)os * ii - my_roundf(m_sto_frac * os);
                if (idx < 0) idx = 0;
                if ((size_t)(idx + fir_delay) >= avail()) return; // need more data
                in_down[ii] = sample_at(idx);
            }

            long items_to_consume = sps;

            switch (state) {
            case DETECT: {
                uint32_t bin_idx_new = fft_argmax(in_down.data(), downchirp.data());

                if (std::abs((long)mod(std::abs((long)bin_idx_new - (long)bin_idx) + 1, N) - 1) <= 1 && bin_idx_new != (uint32_t)-1) {
                    if (symbol_cnt == 1 && bin_idx != (uint32_t)-1)
                        preamb_up_vals[0] = bin_idx;
                    preamb_up_vals[symbol_cnt] = bin_idx_new;
                    memcpy(&preamble_raw[N * symbol_cnt], in_down.data(), N * sizeof(cx));
                    symbol_cnt++;
                } else {
                    memcpy(&preamble_raw[0], in_down.data(), N * sizeof(cx));
                    symbol_cnt = 1;
                }
                bin_idx = bin_idx_new;

                if (symbol_cnt == n_up_req) {
                    additional_upchirps = 0;
                    state = SYNC;
                    symbol_cnt = 0;
                    cfo_frac_sto_frac_est = false;
                    k_hat = most_frequent(preamb_up_vals.data(), (int)preamb_up_vals.size());
                    items_to_consume = (long)os * ((long)N - k_hat);
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
                bin_idx = fft_argmax(symb_corr.data(), downchirp.data());

                switch (symbol_cnt) {
                case NET_ID1:
                    if (bin_idx == 0 || bin_idx == 1 || bin_idx == N - 1) {
                        // Additional upchirp, not net ID yet
                        if (additional_upchirps < 3) additional_upchirps++;
                    } else {
                        symbol_cnt = NET_ID2;
                        net_ids[0] = (int)bin_idx;
                    }
                    break;
                case NET_ID2:
                    symbol_cnt = DOWNCHIRP1;
                    net_ids[1] = (int)bin_idx;
                    break;
                case DOWNCHIRP1:
                    symbol_cnt = DOWNCHIRP2;
                    break;
                case DOWNCHIRP2:
                    down_val = (int)fft_argmax(symb_corr.data(), upchirp.data());
                    symbol_cnt = QUARTER_DOWN;
                    break;
                case QUARTER_DOWN: {
                    if ((uint32_t)down_val < N / 2)
                        m_cfo_int = (int)floor((double)down_val / 2.0);
                    else
                        m_cfo_int = (int)floor((double)(down_val - (int)N) / 2.0);

                    // Net IDs corrected for integer CFO -> recovered sync word
                    int nid[2];
                    for (int i = 0; i < 2; i++)
                        nid[i] = (int)mod(net_ids[i] - m_cfo_int, N);
                    uint8_t det_sync = (uint8_t)(((((nid[0] + 4) >> 3) & 0xF) << 4) |
                                                 (((nid[1] + 4) >> 3) & 0xF));

                    // Validate sync word (both symbols) unless accept-any (-w 0)
                    bool sync_ok = true;
                    if (cfg.sync_words[0] != 0 || cfg.sync_words[1] != 0) {
                        for (int i = 0; i < 2; i++) {
                            long d = std::min(mod(nid[i] - cfg.sync_words[i], N),
                                              mod(cfg.sync_words[i] - nid[i], N));
                            if (d > 2) sync_ok = false;
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

                    // Estimate SNR from a preamble symbol (CFO-frac corrected)
                    float snr_e = 0;
                    {
                        cx_multiply(symb_corr.data(), preamble_raw.data(), CFO_frac_correc.data(), N);
                        cx_multiply(dechirp_buf.data(), symb_corr.data(), downchirp.data(), N);
                        fft_N(dechirp_buf.data(), fft_out.data());
                        float tot_e = 0, sig_e = 0;
                        for (uint32_t i = 0; i < N; i++) {
                            float m = std::norm(fft_out[i]);
                            tot_e += m;
                            if (m > sig_e) sig_e = m;
                        }
                        if (tot_e > sig_e && sig_e > 0)
                            snr_e = 10.0f * log10f(sig_e / (tot_e - sig_e));
                    }
                    current_snr = snr_e;

                    state = PAYLOAD;
                    payload_symbol_cnt = 0;
                    m_received_head = false;
                    m_ldro = false;
                    frame_cnt++;

                    // Reset payload decode state
                    fft_block.clear();
                    hard_block.clear();
                    nibbles.clear();
                    is_header = true;
                    m_pay_cr = cfg.cr;
                    m_pay_len = cfg.pay_len;
                    m_pay_has_crc = cfg.has_crc;

                    items_to_consume = sps / 4 + (long)os * m_cfo_int;
                    if (items_to_consume < 0) items_to_consume = 0;

                    {
                        std::lock_guard<std::mutex> lk(g_print_mtx);
                        fprintf(stderr, "%s[sync] Frame #%u detected, CFO=%.2f, sync=0x%02X, SNR=%.1f dB\n",
                                tag.c_str(), frame_cnt, m_cfo_int + m_cfo_frac, det_sync, snr_e);
                    }
                    break;
                }
                default:
                    break;
                }
                break;
            }

            case PAYLOAD: {
                if (payload_symbol_cnt < 8 || ((uint32_t)payload_symbol_cnt < total_payload_symbols && m_received_head)) {
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
            if (items_to_consume > 0 && (size_t)items_to_consume <= avail()) {
                head += items_to_consume;
                compact();
            } else {
                break; // waiting state or not enough data
            }
        }
    }

private:
    LoRaConfig cfg;
    std::string tag;
    uint32_t N;    // number of bins = 2^sf
    uint8_t os;    // oversampling factor
    uint32_t sps;  // samples per symbol

    std::vector<cx> upchirp, downchirp;
    std::vector<cx> demod_upchirp, demod_downchirp;
    kiss_fft_cfg fft_cfg, fft_cfg2;

    // FFT scratch (preallocated, reused every symbol)
    std::vector<cx> dechirp_buf, fft_out;
    std::vector<float> fft_mag;
    std::vector<cx> cfo_fft_val;

    // Gray demap LUTs
    std::vector<uint16_t> demap_pay, demap_hdr;

    // Anti-alias decimation FIR
    std::vector<float> fir_taps;
    int fir_delay = 0;

    // Incoming sample buffer: consumed by advancing `head`, compacted rarely
    std::vector<cx> ring;
    size_t head = 0;

    size_t avail() const { return ring.size() - head; }

    void compact() {
        if (head > (size_t)(1 << 17)) {
            // Keep fir_delay samples of history before head for the FIR
            size_t drop = head - fir_delay;
            ring.erase(ring.begin(), ring.begin() + drop);
            head = fir_delay;
        }
    }

    // Filtered sample at full-rate index idx (relative to head).
    // Applies the anti-alias low-pass; identity when os == 1.
    inline cx sample_at(long idx) const {
        if (fir_taps.empty()) return ring[head + idx];
        long base = (long)head + idx - fir_delay;
        const int T = (int)fir_taps.size();
        int k0 = base < 0 ? (int)(-base) : 0;
        cx acc(0, 0);
        for (int k = k0; k < T; k++)
            acc += fir_taps[k] * ring[base + k];
        return acc;
    }

    inline void fft_N(const cx *in, cx *out) {
        kiss_fft(fft_cfg, (const kiss_fft_cpx *)in, (kiss_fft_cpx *)out);
    }

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
    std::vector<cx> preamble_upchirps;
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
    unsigned int frame_cnt = 0;
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
    std::vector<std::array<double, MAX_SF>> fft_block; // LLR blocks for soft decoding
    std::vector<uint16_t> hard_block;                  // symbol values for hard decoding
    std::vector<uint8_t> nibbles;                      // decoded nibbles

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
        net_ids[0] = net_ids[1] = 0;
        payload_symbol_cnt = 0;
        total_payload_symbols = 0;
        m_received_head = false;
        m_ldro = false;
        is_header = true;
        fft_block.clear();
        hard_block.clear();
        nibbles.clear();
    }

    // Dechirp + FFT, return argmax bin
    uint32_t fft_argmax(const cx *samples, const cx *ref_chirp, float *energy_out = nullptr) {
        cx_multiply(dechirp_buf.data(), samples, ref_chirp, N);
        fft_N(dechirp_buf.data(), fft_out.data());
        float max_val = -1;
        uint32_t max_idx = 0;
        float tot_en = 0;
        for (uint32_t i = 0; i < N; i++) {
            float mag = std::norm(fft_out[i]);
            tot_en += mag;
            if (mag > max_val) { max_val = mag; max_idx = i; }
        }
        if (energy_out) *energy_out = tot_en;
        return (tot_en > 0) ? max_idx : (uint32_t)-1;
    }

    // ---- CFO estimation (Bernier's algorithm) ----
    float estimate_CFO_frac_Bernier(const cx *samples) {
        std::vector<int> k0v(up_symb_to_use);
        std::vector<double> k0_mag(up_symb_to_use);

        for (int i = 0; i < up_symb_to_use; i++) {
            cx_multiply(dechirp_buf.data(), &samples[N * i], downchirp.data(), N);
            fft_N(dechirp_buf.data(), &cfo_fft_val[i * N]);

            float max_mag = -1;
            int max_idx = 0;
            for (uint32_t j = 0; j < N; j++) {
                float m = std::norm(cfo_fft_val[j + i * N]);
                if (m > max_mag) { max_mag = m; max_idx = (int)j; }
            }
            k0v[i] = max_idx;
            k0_mag[i] = max_mag;
        }

        int best = (int)(std::max_element(k0_mag.begin(), k0_mag.end()) - k0_mag.begin());
        int idx_max = k0v[best];

        cx four_cum(0, 0);
        for (int i = 0; i < up_symb_to_use - 1; i++)
            four_cum += cfo_fft_val[idx_max + N * i] * std::conj(cfo_fft_val[idx_max + N * (i + 1)]);

        float cfo_frac = -std::arg(four_cum) / (2.0f * (float)M_PI);

        // Correct CFO in preamble
        for (int n = 0; n < up_symb_to_use * (int)N; n++)
            preamble_upchirps[n] = samples[n] * expj(-2.0f * (float)M_PI * cfo_frac / N * n);

        return cfo_frac;
    }

    // ---- STO fractional estimation (RCTSL) ----
    float estimate_STO_frac() {
        std::fill(fft_mag.begin(), fft_mag.end(), 0.0f);

        for (int i = 0; i < up_symb_to_use; i++) {
            cx_multiply(dechirp_buf.data(), &preamble_upchirps[N * i], downchirp.data(), N);
            std::fill(dechirp_buf.begin() + N, dechirp_buf.begin() + 2 * N, cx(0, 0));
            kiss_fft(fft_cfg2, (const kiss_fft_cpx *)dechirp_buf.data(), (kiss_fft_cpx *)fft_out.data());
            for (uint32_t j = 0; j < 2 * N; j++)
                fft_mag[j] += std::norm(fft_out[j]);
        }

        int k0 = (int)(std::max_element(fft_mag.begin(), fft_mag.end()) - fft_mag.begin());
        double Y_1 = fft_mag[mod(k0 - 1, 2 * N)];
        double Y0 = fft_mag[k0];
        double Y1 = fft_mag[mod(k0 + 1, 2 * N)];

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
            fft_block.push_back(compute_LLRs(symbol));
            if ((int)fft_block.size() == block_size) {
                decode_block_soft();
                fft_block.clear();
            }
        } else {
            uint32_t idx = fft_argmax(symbol, demod_downchirp.data());
            const uint16_t *map = (is_header || m_ldro) ? demap_hdr.data() : demap_pay.data();
            hard_block.push_back(map[idx == (uint32_t)-1 ? 0 : idx]);
            if ((int)hard_block.size() == block_size) {
                decode_block_hard();
                hard_block.clear();
            }
        }
    }

    // ---- Compute LLRs for one symbol (soft demod, max-log approximation) ----
    std::array<double, MAX_SF> compute_LLRs(const cx *samples) {
        cx_multiply(dechirp_buf.data(), samples, demod_downchirp.data(), N);
        fft_N(dechirp_buf.data(), fft_out.data());
        for (uint32_t i = 0; i < N; i++)
            fft_mag[i] = std::norm(fft_out[i]);

        const uint16_t *map = (is_header || m_ldro) ? demap_hdr.data() : demap_pay.data();

        std::array<double, MAX_SF> LLRs{};
        for (uint32_t i = 0; i < cfg.sf; i++) {
            float max_X1 = -1, max_X0 = -1;
            for (uint32_t n = 0; n < N; n++) {
                if (map[n] & (1u << i)) {
                    if (fft_mag[n] > max_X1) max_X1 = fft_mag[n];
                } else {
                    if (fft_mag[n] > max_X0) max_X0 = fft_mag[n];
                }
            }
            LLRs[cfg.sf - 1 - i] = (double)max_X1 - (double)max_X0;
        }
        return LLRs;
    }

    // ---- Decode a block (soft) ----
    void decode_block_soft() {
        int sf_app = (is_header || m_ldro) ? cfg.sf - 2 : cfg.sf;
        int cw_len = is_header ? 8 : m_pay_cr + 4;
        int cr_app = is_header ? 4 : m_pay_cr;

        // Deinterleave
        double deinter[MAX_SF][8];
        for (int i = 0; i < cw_len; i++)
            for (int j = 0; j < sf_app; j++)
                deinter[mod(i - j - 1, sf_app)][i] = fft_block[i][cfg.sf - sf_app + j];

        // Hamming decode each codeword
        for (int i = 0; i < sf_app; i++)
            nibbles.push_back(hamming_decode_soft(deinter[i], cr_app));

        process_nibbles();
    }

    // ---- Decode a block (hard) ----
    void decode_block_hard() {
        int sf_app = (is_header || m_ldro) ? cfg.sf - 2 : cfg.sf;
        int cw_len = is_header ? 8 : m_pay_cr + 4;
        int cr_app = is_header ? 4 : m_pay_cr;

        uint8_t deinter[MAX_SF][8] = {};
        for (int i = 0; i < cw_len; i++)
            for (int j = 0; j < sf_app; j++)
                deinter[mod(i - j - 1, sf_app)][i] = (hard_block[i] >> (sf_app - 1 - j)) & 1;

        for (int i = 0; i < sf_app; i++) {
            uint8_t cw_byte = 0;
            for (int j = 0; j < cw_len; j++)
                cw_byte = (cw_byte << 1) | deinter[i][j];
            nibbles.push_back(hamming_decode_hard(cw_byte, cr_app));
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

            {
                std::lock_guard<std::mutex> lk(g_print_mtx);
                fprintf(stderr, "%s  Header: pay_len=%u, CRC=%d, CR=%d\n", tag.c_str(), m_pay_len, m_pay_has_crc, m_pay_cr);
            }

            if (computed_chk != header_chk || m_pay_len == 0) {
                std::lock_guard<std::mutex> lk(g_print_mtx);
                fprintf(stderr, "%s  Header checksum INVALID!\n", tag.c_str());
                reset_state();
                return;
            }

            {
                std::lock_guard<std::mutex> lk(g_print_mtx);
                fprintf(stderr, "%s  Header checksum valid\n", tag.c_str());
            }

            // Determine LDRO
            m_ldro = ((float)(1u << cfg.sf) * 1e3f / cfg.bw) > LDRO_MAX_DURATION_MS;

            // Calculate total payload symbols
            int sf_minus_2ldro = cfg.sf - 2 * m_ldro;
            total_payload_symbols = 8 + (int)ceil((double)(2 * m_pay_len - cfg.sf + 2 + (!cfg.impl_head ? 5 : 0) + (m_pay_has_crc ? 4 : 0)) / sf_minus_2ldro) * (4 + m_pay_cr);

            m_received_head = true;
            is_header = false;

            // Remove header nibbles, keep remaining payload nibbles
            nibbles.erase(nibbles.begin(), nibbles.begin() + 5);

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

        {
            std::lock_guard<std::mutex> lk(g_print_mtx);

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
        }

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
// Worker: one demodulator on its own thread, fed via a bounded queue
// ============================================================================

class DemodWorker {
public:
    DemodWorker(const LoRaConfig &c, const std::string &tag)
        : demod(c, tag), tag_(tag) {
        th = std::thread([this] { run(); });
    }

    void push(std::shared_ptr<const std::vector<cx>> chunk) {
        std::unique_lock<std::mutex> lk(mtx);
        cv_space.wait(lk, [&] { return q.size() < MAX_Q || stopping; });
        if (stopping) return;
        q.push_back(std::move(chunk));
        cv_data.notify_one();
    }

    // Signals end of stream; drains the queue, then joins the thread.
    void finish() {
        {
            std::lock_guard<std::mutex> lk(mtx);
            stopping = true;
        }
        cv_data.notify_all();
        cv_space.notify_all();
        if (th.joinable()) th.join();
    }

    size_t packet_count() const { return demod.last_packets.size(); }
    const std::string &tag() const { return tag_; }

private:
    void run() {
        for (;;) {
            std::shared_ptr<const std::vector<cx>> chunk;
            {
                std::unique_lock<std::mutex> lk(mtx);
                cv_data.wait(lk, [&] { return !q.empty() || stopping; });
                if (q.empty()) return; // stopping and drained
                chunk = std::move(q.front());
                q.pop_front();
                cv_space.notify_one();
            }
            demod.process_samples(chunk->data(), chunk->size());
        }
    }

    LoRaDemodulator demod;
    std::string tag_;
    std::thread th;
    std::mutex mtx;
    std::condition_variable cv_data, cv_space;
    std::deque<std::shared_ptr<const std::vector<cx>>> q;
    bool stopping = false;
    static constexpr size_t MAX_Q = 64;
};

// ============================================================================
// Sample source -> workers fan-out
// ============================================================================

static float g_u8lut[256];
static std::vector<std::unique_ptr<DemodWorker>> g_workers;
static FILE *g_dump = nullptr;

static void broadcast_u8(const unsigned char *buf, uint32_t len) {
    if (g_dump) fwrite(buf, 1, len, g_dump);

    size_t n = len / 2;
    auto chunk = std::make_shared<std::vector<cx>>(n);
    auto &v = *chunk;
    for (size_t i = 0; i < n; i++)
        v[i] = cx(g_u8lut[buf[2 * i]], g_u8lut[buf[2 * i + 1]]);

    for (auto &w : g_workers)
        w->push(chunk);
}

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *) {
    if (!g_running) return;
    broadcast_u8(buf, len);
}

// ============================================================================
// Main
// ============================================================================

static std::string make_tag(uint8_t sf, uint32_t bw) {
    char buf[32];
    if (bw % 1000 == 0)
        snprintf(buf, sizeof(buf), "[SF%u/%uk] ", sf, bw / 1000);
    else
        snprintf(buf, sizeof(buf), "[SF%u/%.1fk] ", sf, bw / 1000.0);
    return buf;
}

int main(int argc, char *argv[]) {
    LoRaConfig cfg;
    const char *iq_file = nullptr;
    const char *dump_file = nullptr;
    bool auto_mode = false, sf_given = false, bw_given = false;
    int opt;

    while ((opt = getopt(argc, argv, "f:s:b:S:c:w:g:p:IL:r:D:A")) != -1) {
        switch (opt) {
        case 'f': cfg.freq = (uint32_t)atol(optarg); break;
        case 's': cfg.samp_rate = (uint32_t)atol(optarg); break;
        case 'b': cfg.bw = (uint32_t)atol(optarg); bw_given = true; break;
        case 'S': cfg.sf = (uint8_t)atoi(optarg); sf_given = true; break;
        case 'c': cfg.cr = (uint8_t)atoi(optarg); break;
        case 'w': cfg.sync_word = (uint16_t)strtol(optarg, nullptr, 16); break;
        case 'g': cfg.gain = atoi(optarg); break;
        case 'p': cfg.ppm = atoi(optarg); break;
        case 'I': cfg.impl_head = true; break;
        case 'L': cfg.pay_len = (uint32_t)atoi(optarg); break;
        case 'r': iq_file = optarg; break;
        case 'D': dump_file = optarg; break;
        case 'A': auto_mode = true; break;
        default:
            fprintf(stderr, "Usage: %s [-f freq] [-s samp_rate] [-b bw] [-S sf] [-c cr] [-w sync_word_hex]\n"
                            "          [-g gain] [-p ppm] [-I] [-L pay_len] [-A] [-r iq_file] [-D dump_file]\n", argv[0]);
            return 1;
        }
    }

    if (cfg.sf < MIN_SF || cfg.sf > MAX_SF) {
        fprintf(stderr, "SF must be between %d and %d\n", MIN_SF, MAX_SF);
        return 1;
    }
    if (cfg.cr < 1 || cfg.cr > 4) {
        uint8_t clamped = cfg.cr < 1 ? 1 : 4;
        fprintf(stderr, "Warning: CR %u out of range 1-4, using %u (explicit header overrides it anyway)\n",
                cfg.cr, clamped);
        cfg.cr = clamped;
    }
    if (cfg.samp_rate % cfg.bw != 0) {
        fprintf(stderr, "Warning: samp_rate should be an integer multiple of bandwidth\n");
    }

    for (int i = 0; i < 256; i++)
        g_u8lut[i] = ((float)i - 127.5f) / 127.5f;

    // Build the list of (SF, BW) decoder configurations
    std::vector<uint8_t> sfs = {cfg.sf};
    std::vector<uint32_t> bws = {cfg.bw};
    if (auto_mode) {
        if (!sf_given) sfs = {7, 8, 9, 10, 11, 12};
        if (!bw_given) {
            bws.clear();
            for (uint32_t bw : {62500u, 125000u, 250000u, 500000u}) {
                if (bw <= cfg.samp_rate && cfg.samp_rate % bw == 0 && cfg.samp_rate / bw <= 16)
                    bws.push_back(bw);
            }
            if (bws.empty()) bws = {cfg.bw};
        }
    }

    bool multi = sfs.size() * bws.size() > 1;
    for (uint32_t bw : bws) {
        for (uint8_t sf : sfs) {
            LoRaConfig c = cfg;
            c.sf = sf;
            c.bw = bw;
            c.compute_derived();
            g_workers.emplace_back(new DemodWorker(c, multi ? make_tag(sf, bw) : ""));
        }
    }

    cfg.compute_derived();
    fprintf(stderr, "LoRa Standalone RX\n");
    fprintf(stderr, "  Freq:       %u Hz\n", cfg.freq);
    fprintf(stderr, "  Samp rate:  %u Hz\n", cfg.samp_rate);
    if (multi) {
        fprintf(stderr, "  Decoders:   %zu in parallel:", g_workers.size());
        for (auto &w : g_workers) fprintf(stderr, " %s", w->tag().c_str());
        fprintf(stderr, "\n");
    } else {
        fprintf(stderr, "  Bandwidth:  %u Hz\n", cfg.bw);
        fprintf(stderr, "  SF:         %u\n", cfg.sf);
        fprintf(stderr, "  OS factor:  %u\n", cfg.os_factor);
        fprintf(stderr, "  Bins:       %u\n", cfg.n_bins);
        fprintf(stderr, "  Samples/sym: %u\n", cfg.samples_per_symbol);
    }
    fprintf(stderr, "  CR:         4/%u (explicit header overrides)\n", cfg.cr + 4);
    fprintf(stderr, "  Sync word:  0x%02X%s\n", cfg.sync_word, cfg.sync_word == 0 ? " (accept any)" : "");
    fprintf(stderr, "  Implicit:   %s\n", cfg.impl_head ? "yes" : "no");
    fprintf(stderr, "  Soft decode: %s\n", cfg.soft_decoding ? "yes" : "no");
    fprintf(stderr, "\n");

    if (dump_file) {
        g_dump = fopen(dump_file, "wb");
        if (!g_dump) {
            perror("dump file");
            return 1;
        }
        fprintf(stderr, "Dumping raw IQ to %s\n", dump_file);
    }

    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    if (iq_file) {
        // ---- Offline: replay IQ from file (rtl_sdr u8 format) ----
        FILE *f = strcmp(iq_file, "-") == 0 ? stdin : fopen(iq_file, "rb");
        if (!f) {
            perror("iq file");
            return 1;
        }
        fprintf(stderr, "Replaying IQ from %s...\n\n", iq_file);
        std::vector<unsigned char> buf(1 << 18);
        size_t nr;
        while (g_running && (nr = fread(buf.data(), 1, buf.size(), f)) > 0)
            broadcast_u8(buf.data(), (uint32_t)nr);
        if (f != stdin) fclose(f);
    } else {
        // ---- Live: RTL-SDR ----
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

        rtlsdr_set_sample_rate(dev, cfg.samp_rate);
        rtlsdr_set_center_freq(dev, cfg.freq);
        rtlsdr_set_freq_correction(dev, cfg.ppm);
        rtlsdr_set_tuner_gain_mode(dev, 1); // manual gain
        rtlsdr_set_tuner_gain(dev, cfg.gain);
        rtlsdr_set_agc_mode(dev, 0);
        rtlsdr_reset_buffer(dev);

        fprintf(stderr, "RTL-SDR configured. Tuner gain: %.1f dB\n", rtlsdr_get_tuner_gain(dev) / 10.0);
        fprintf(stderr, "Listening...\n\n");

        // 16384 bytes = 8192 IQ samples per callback
        rtlsdr_read_async(dev, rtlsdr_callback, nullptr, 0, 16384);

        rtlsdr_close(dev);
    }

    size_t total = 0;
    for (auto &w : g_workers)
        w->finish();
    fprintf(stderr, "\nDone.\n");
    for (auto &w : g_workers) {
        total += w->packet_count();
        if (multi && w->packet_count() > 0)
            fprintf(stderr, "  %s%zu packet(s)\n", w->tag().c_str(), w->packet_count());
    }
    fprintf(stderr, "Received %zu packet(s) total.\n", total);

    if (g_dump) fclose(g_dump);
    return 0;
}
