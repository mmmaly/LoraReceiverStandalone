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

// Soft-decision Hamming decode. Optionally reports the runner-up data nibble
// and the score margin between best and runner-up (small margin = uncertain
// decision), which feeds the CRC-guided chase recovery on failed payloads.
static uint8_t hamming_decode_soft(const double *codeword_LLR, int cr_app,
                                   uint8_t *alt_nibble = nullptr, double *margin = nullptr) {
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
    auto to_nibble = [](uint8_t data) -> uint8_t {
        // reverse bit order
        return ((data & 1) << 3) | (((data >> 1) & 1) << 2) | (((data >> 2) & 1) << 1) | ((data >> 3) & 1);
    };
    if (alt_nibble || margin) {
        double best2 = -1e300;
        int idx2 = (idx_max + 1) & 15;
        for (int n = 0; n < 16; n++) {
            if (n == idx_max) continue;
            if (cw_proba[n] > best2) { best2 = cw_proba[n]; idx2 = n; }
        }
        if (alt_nibble) *alt_nibble = to_nibble(cw_LUT[idx2] >> 4);
        if (margin) *margin = cw_proba[idx_max] - best2;
    }
    return to_nibble(cw_LUT[idx_max] >> 4);
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
    bool     rtl_agc        = false; // RTL2832 internal digital AGC
    bool     tuner_agc      = false; // R820T analog tuner AGC (instead of manual gain)
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

    LoRaDemodulator(const LoRaConfig &cfg_, std::string tag_ = "", std::string stdout_tag_ = "")
        : cfg(cfg_), tag(std::move(tag_)), stdout_tag(std::move(stdout_tag_)) {
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

                    // Net ID bins already have integer CFO absorbed by the
                    // k_hat timing alignment in DETECT - use them raw.
                    // (Subtracting m_cfo_int here double-corrects and rejects
                    // every frame with |CFO| > 2 bins.)
                    int nid[2];
                    for (int i = 0; i < 2; i++)
                        nid[i] = (int)mod(net_ids[i], N);
                    uint8_t det_sync = (uint8_t)(((((nid[0] + 4) >> 3) & 0xF) << 4) |
                                                 (((nid[1] + 4) >> 3) & 0xF));

                    // Validate sync word unless accept-any (-w 0). Only the
                    // first net-ID symbol is checked: real frames sometimes
                    // mis-latch the second one yet still decode fine, and the
                    // header checksum + CRC reject false syncs anyway.
                    bool sync_ok = true;
                    if (cfg.sync_words[0] != 0 || cfg.sync_words[1] != 0) {
                        long d = std::min(mod(nid[0] - cfg.sync_words[0], N),
                                          mod(cfg.sync_words[0] - nid[0], N));
                        if (d > 2) sync_ok = false;
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
                    // Remember where the frame's symbols start (and the SFO
                    // accumulator state) so a failed frame can be re-run at
                    // shifted timing from the retained ring samples.
                    frame_start_abs = ring_abs0 + head + (size_t)items_to_consume;
                    frame_sfo_cum0 = sfo_cum;

                    {
                        // Signed residual of the first net-ID symbol from its
                        // nominal sync bin: a known-symbol measurement of any
                        // constant bin offset the payload will suffer.
                        long dpos = mod(nid[0] - cfg.sync_words[0], N);
                        long nid_resid = dpos <= (long)N / 2 ? dpos : dpos - (long)N;
                        std::lock_guard<std::mutex> lk(g_print_mtx);
                        fprintf(stderr, "%s[sync] Frame #%u detected, CFO=%.2f, sync=0x%02X, SNR=%.1f dB, nid_resid=%+ld\n",
                                tag.c_str(), frame_cnt, m_cfo_int + m_cfo_frac, det_sync, snr_e, nid_resid);
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
    std::string stdout_tag; // when non-empty, printed as "rx cfg: ..." before each packet
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
    size_t ring_abs0 = 0;        // absolute sample index of ring[0]
    size_t frame_start_abs = 0;  // absolute index of first payload-window sample

    size_t avail() const { return ring.size() - head; }

    void compact() {
        if (head > (size_t)(1 << 17)) {
            // Keep fir_delay samples of history before head for the FIR
            size_t drop = head - fir_delay;
            // While decoding a frame, keep everything from a symbol before its
            // start so a failed frame can be re-demodulated at shifted timing
            if (state == PAYLOAD) {
                long keep_from = (long)(frame_start_abs - ring_abs0) - (long)sps - fir_delay;
                if (keep_from < 0) keep_from = 0;
                if ((size_t)keep_from < drop) drop = (size_t)keep_from;
            }
            if (drop == 0) return;
            ring.erase(ring.begin(), ring.begin() + drop);
            head -= drop;
            ring_abs0 += drop;
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

    // Same, addressed by absolute sample index; false if outside the ring.
    inline bool sample_at_abs(long abs_idx, cx &out) const {
        long rel = abs_idx - (long)ring_abs0;
        if (rel < 0 || (size_t)rel >= ring.size()) return false;
        if (fir_taps.empty()) { out = ring[rel]; return true; }
        long base = rel - fir_delay;
        const int T = (int)fir_taps.size();
        if (base + T > (long)ring.size()) return false;
        int k0 = base < 0 ? (int)(-base) : 0;
        cx acc(0, 0);
        for (int k = k0; k < T; k++)
            acc += fir_taps[k] * ring[base + k];
        out = acc;
        return true;
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
    float frame_sfo_cum0 = 0;

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
    std::vector<uint8_t> nib_alt;                      // runner-up nibble per position
    std::vector<double> nib_margin;                    // decision margin per position

    // Per-block raw FFT magnitudes of the payload symbols, kept for the
    // symbol-level chase: one corrupted symbol (collision burst) damages up to
    // sf_app nibbles at once, a pattern the nibble-level chase cannot reach.
    struct BlockRec {
        int nib_start;                        // first nibble index this block produced
        int cw_len, sf_app, cr_app;
        std::vector<std::vector<float>> mags; // cw_len entries of N bin magnitudes
    };
    std::vector<BlockRec> block_recs;
    std::vector<std::vector<float>> pending_mags;      // symbols of the block in progress

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
        nib_alt.clear();
        nib_margin.clear();
        block_recs.clear();
        pending_mags.clear();
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
            if (!is_header)  // keep raw bin magnitudes for the symbol-level chase
                pending_mags.emplace_back(fft_mag.begin(), fft_mag.begin() + N);
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
        return compute_LLRs_map(samples, is_header || m_ldro);
    }

    // rot: constant bin rotation applied when reading the spectrum - a symbol
    // window shifted by k chips lands every dechirped peak k bins away, so a
    // retimed decode must read the constellation counter-rotated.
    std::array<double, MAX_SF> compute_LLRs_map(const cx *samples, bool hdr_map, int rot = 0) {
        cx_multiply(dechirp_buf.data(), samples, demod_downchirp.data(), N);
        fft_N(dechirp_buf.data(), fft_out.data());
        // Amplitude (not power) bin metric: for noncoherent detection the LLR
        // grows ~linearly with amplitude, and power scaling overweights strong
        // symbols against weak ones inside a deinterleaved codeword.
        for (uint32_t i = 0; i < N; i++)
            fft_mag[i] = sqrtf(std::norm(fft_out[i]));

        const uint16_t *map = hdr_map ? demap_hdr.data() : demap_pay.data();

        std::array<double, MAX_SF> LLRs{};
        for (uint32_t i = 0; i < cfg.sf; i++) {
            float max_X1 = -1, max_X0 = -1;
            for (uint32_t n = 0; n < N; n++) {
                float v = fft_mag[(n + (uint32_t)(rot + (int)N)) % N];
                if (map[n] & (1u << i)) {
                    if (v > max_X1) max_X1 = v;
                } else {
                    if (v > max_X0) max_X0 = v;
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

        // Hamming decode each codeword, keeping the runner-up + margin for
        // possible CRC-guided recovery later
        for (int i = 0; i < sf_app; i++) {
            uint8_t alt;
            double mg;
            nibbles.push_back(hamming_decode_soft(deinter[i], cr_app, &alt, &mg));
            nib_alt.push_back(alt);
            nib_margin.push_back(mg);
        }

        if (!is_header && (int)pending_mags.size() == cw_len) {
            BlockRec rec;
            rec.nib_start = (int)nibbles.size() - sf_app;
            rec.cw_len = cw_len;
            rec.sf_app = sf_app;
            rec.cr_app = cr_app;
            rec.mags = std::move(pending_mags);
            block_recs.push_back(std::move(rec));
        }
        pending_mags.clear();

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
            // no soft info in the hard path: mark as fully confident
            nib_alt.push_back(nibbles.back());
            nib_margin.push_back(1e300);
        }

        process_nibbles();
    }

    // Header checksum over the first 5 nibbles (5-bit result).
    static int header_checksum(const uint8_t *n) {
        bool c4 = ((n[0] >> 3) & 1) ^ ((n[0] >> 2) & 1) ^ ((n[0] >> 1) & 1) ^ (n[0] & 1);
        bool c3 = ((n[0] >> 3) & 1) ^ ((n[1] >> 3) & 1) ^ ((n[1] >> 2) & 1) ^ ((n[1] >> 1) & 1) ^ (n[2] & 1);
        bool c2 = ((n[0] >> 2) & 1) ^ ((n[1] >> 3) & 1) ^ (n[1] & 1) ^ ((n[2] >> 3) & 1) ^ ((n[2] >> 1) & 1);
        bool c1 = ((n[0] >> 1) & 1) ^ ((n[1] >> 2) & 1) ^ (n[1] & 1) ^ ((n[2] >> 2) & 1) ^ ((n[2] >> 1) & 1) ^ (n[2] & 1);
        bool c0 = (n[0] & 1) ^ ((n[1] >> 1) & 1) ^ ((n[2] >> 3) & 1) ^ ((n[2] >> 2) & 1) ^ ((n[2] >> 1) & 1) ^ (n[2] & 1);
        return (c4 << 4) | (c3 << 3) | (c2 << 2) | (c1 << 1) | c0;
    }

    static bool header_consistent(const uint8_t *n) {
        uint8_t chk = ((n[3] & 1) << 4) | n[4];
        uint32_t len = (n[0] << 4) | n[1];
        return header_checksum(n) == chk && len != 0;
    }

    // ---- Process accumulated nibbles (header decode / payload assembly) ----
    void process_nibbles() {
        if (is_header && nibbles.size() >= 5 && !cfg.impl_head) {
            // (Header chase recovery was tried here and removed: with only a
            // 5-bit checksum it accepts bogus headers, and a phantom payload
            // deafens the receiver to following real frames - measured net
            // negative on real captures.)

            // Decode explicit header
            m_pay_len = (nibbles[0] << 4) | nibbles[1];
            m_pay_has_crc = nibbles[2] & 1;
            m_pay_cr = nibbles[2] >> 1;

            if (m_pay_cr < 1) m_pay_cr = 1;
            if (m_pay_cr > 4) m_pay_cr = 4;

            uint8_t header_chk = ((nibbles[3] & 1) << 4) | nibbles[4];
            int computed_chk = header_checksum(nibbles.data());

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
            nib_alt.erase(nib_alt.begin(), nib_alt.begin() + 5);
            nib_margin.erase(nib_margin.begin(), nib_margin.begin() + 5);

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

    // Build dewhitened bytes from a nibble sequence (payload + CRC trailer).
    void nibbles_to_bytes(const uint8_t *nib, uint8_t *out, uint32_t total_bytes) const {
        for (uint32_t i = 0; i < total_bytes; i++) {
            uint8_t low_nib = nib[2 * i];
            uint8_t high_nib = nib[2 * i + 1];
            if (i < m_pay_len) {
                // Dewhiten payload; CRC bytes are NOT dewhitened
                low_nib ^= (whitening_seq[i] & 0x0F);
                high_nib ^= (whitening_seq[i] >> 4) & 0x0F;
            }
            out[i] = (high_nib << 4) | low_nib;
        }
    }

    bool payload_crc_ok(const uint8_t *bytes) const {
        uint16_t calc_crc = crc16(bytes, m_pay_len - 2);
        calc_crc = calc_crc ^ bytes[m_pay_len - 1] ^ ((uint16_t)bytes[m_pay_len - 2] << 8);
        uint16_t rx_crc = bytes[m_pay_len] | ((uint16_t)bytes[m_pay_len + 1] << 8);
        return calc_crc == rx_crc;
    }

    // CRC-guided chase recovery: retry the CRC with the runner-up Hamming
    // codeword substituted at the least-confident nibble positions, testing
    // candidate error patterns in order of increasing implausibility (sum of
    // decision margins). Returns the number of substituted nibbles, 0 if no
    // pattern passed.
    int chase_recover(std::vector<uint8_t> &bytes, uint32_t total_bytes) {
        const uint32_t n_nib = 2 * total_bytes;
        if (nib_margin.size() < n_nib || nib_alt.size() < n_nib) return 0;

        // k least-confident positions
        constexpr int CHASE_K = 16;         // positions considered
        constexpr int CHASE_MAX_SUBST = 4;  // max simultaneous substitutions
        constexpr size_t CHASE_MAX_TRIALS = 2000;
        std::vector<int> pos(n_nib);
        for (uint32_t i = 0; i < n_nib; i++) pos[i] = (int)i;
        std::sort(pos.begin(), pos.end(),
                  [&](int a, int b) { return nib_margin[a] < nib_margin[b]; });
        int k = std::min<int>(CHASE_K, (int)n_nib);

        // Enumerate substitution patterns over the k candidates, cheapest first
        struct Cand { double cost; uint32_t mask; };
        std::vector<Cand> cands;
        for (uint32_t mask = 1; mask < (1u << k); mask++) {
            int nb = __builtin_popcount(mask);
            if (nb > CHASE_MAX_SUBST) continue;
            double cost = 0;
            for (int b = 0; b < k; b++)
                if (mask & (1u << b)) cost += nib_margin[pos[b]];
            cands.push_back({cost, mask});
        }
        std::sort(cands.begin(), cands.end(),
                  [](const Cand &a, const Cand &b) { return a.cost < b.cost; });
        if (cands.size() > CHASE_MAX_TRIALS) cands.resize(CHASE_MAX_TRIALS);

        std::vector<uint8_t> trial_nib(nibbles.begin(), nibbles.begin() + n_nib);
        std::vector<uint8_t> trial_bytes(total_bytes);
        for (const Cand &c : cands) {
            for (int b = 0; b < k; b++)
                if (c.mask & (1u << b)) trial_nib[pos[b]] = nib_alt[pos[b]];
            nibbles_to_bytes(trial_nib.data(), trial_bytes.data(), total_bytes);
            if (payload_crc_ok(trial_bytes.data())) {
                bytes = trial_bytes;
                return __builtin_popcount(c.mask);
            }
            // undo substitutions for the next pattern
            for (int b = 0; b < k; b++)
                if (c.mask & (1u << b)) trial_nib[pos[b]] = nibbles[pos[b]];
        }
        return 0;
    }

    // Re-derive one block's nibbles from stored magnitudes, optionally
    // suppressing the dominant bin of some symbols (collision hypothesis:
    // the top peak was the interferer, the runner-up peak is ours) and/or
    // rotating the bin->symbol mapping (integer CFO mis-split hypothesis).
    void redecode_block(const BlockRec &rec, uint32_t suppress_mask, int rot,
                        std::vector<uint8_t> &out_nib, int companion = 0) const {
        const uint16_t *map = m_ldro ? demap_hdr.data() : demap_pay.data();
        std::array<std::array<double, MAX_SF>, 8> llrs{};
        std::vector<float> mags;
        for (int s = 0; s < rec.cw_len; s++) {
            const float *m = rec.mags[s].data();
            if (suppress_mask & (1u << s)) {
                mags.assign(rec.mags[s].begin(), rec.mags[s].end());
                uint32_t mx = (uint32_t)(std::max_element(mags.begin(), mags.end()) - mags.begin());
                mags[mx] = 0;
                m = mags.data();
            }
            for (uint32_t i = 0; i < cfg.sf; i++) {
                float max_X1 = -1, max_X0 = -1;
                for (uint32_t n = 0; n < N; n++) {
                    float v = m[(n + rot + N) % N];
                    // companion combining: add the energy of a same-data copy
                    // at constant bin offset (simultaneous flood rebroadcast
                    // from a second transmitter, or adjacent-bin smear for
                    // companion=1), turning the collision into diversity
                    if (companion) v += m[(n + rot + N + companion) % N];
                    if (map[n] & (1u << i)) {
                        if (v > max_X1) max_X1 = v;
                    } else {
                        if (v > max_X0) max_X0 = v;
                    }
                }
                llrs[s][cfg.sf - 1 - i] = (double)max_X1 - (double)max_X0;
            }
        }
        double deinter[MAX_SF][8];
        for (int i = 0; i < rec.cw_len; i++)
            for (int j = 0; j < rec.sf_app; j++)
                deinter[mod(i - j - 1, rec.sf_app)][i] = llrs[i][cfg.sf - rec.sf_app + j];
        out_nib.resize(rec.sf_app);
        for (int i = 0; i < rec.sf_app; i++)
            out_nib[i] = hamming_decode_soft(deinter[i], rec.cr_app);
    }

    // Rotation chase: hypothesize the integer/fractional CFO split was off by
    // a whole bin (the frac estimator wraps at +/-0.5, and real dongles sit
    // several bins off), which shifts EVERY payload symbol by a constant
    // offset. The header still decodes because its reduced-rate demap (v/4)
    // tolerates small shifts - so a valid header plus confidently-wrong
    // payload is this failure's signature. Re-decode all blocks from stored
    // magnitudes under a rotated bin mapping; only 4 CRC trials.
    int rotation_chase(std::vector<uint8_t> &bytes, uint32_t total_bytes) {
        const uint32_t n_nib = 2 * total_bytes;
        if (block_recs.empty() || nibbles.size() < n_nib) return 0;
        static const int SHIFTS[] = {1, -1, 2, -2};
        std::vector<uint8_t> trial_nib(nibbles.begin(), nibbles.begin() + n_nib);
        std::vector<uint8_t> trial_bytes(total_bytes), blk_nib;
        for (int shift : SHIFTS) {
            for (const BlockRec &rec : block_recs) {
                if (rec.nib_start >= (int)n_nib) break;
                redecode_block(rec, 0, shift, blk_nib);
                for (int i = 0; i < rec.sf_app; i++) {
                    int p = rec.nib_start + i;
                    if (p >= 0 && p < (int)n_nib) trial_nib[p] = blk_nib[i];
                }
            }
            nibbles_to_bytes(trial_nib.data(), trial_bytes.data(), total_bytes);
            if (payload_crc_ok(trial_bytes.data())) {
                bytes = trial_bytes;
                return shift;
            }
        }
        return 0;
    }

    // Estimate a persistent secondary-peak offset across the frame's payload
    // symbols. Two simultaneous rebroadcasts of the same flood frame (common
    // in MeshCore) appear as one signal plus a same-data copy at a constant
    // bin offset (their CFO difference); this finds that offset. Returns 0 if
    // no offset repeats often enough to be structural.
    int estimate_companion_offset() const {
        std::vector<int> hist(N, 0);
        int nsym = 0;
        for (const auto &rec : block_recs) {
            for (int s = 0; s < rec.cw_len; s++, nsym++) {
                const auto &m = rec.mags[s];
                int b1 = (int)(std::max_element(m.begin(), m.end()) - m.begin());
                float best2 = -1;
                int b2 = -1;
                for (int nn = 0; nn < (int)N; nn++) {
                    int dd = std::abs(((nn - b1 + (int)N / 2) % (int)N + (int)N) % (int)N - (int)N / 2);
                    if (dd <= 2) continue; // skip the main peak's own leakage
                    if (m[nn] > best2) { best2 = m[nn]; b2 = nn; }
                }
                if (b2 >= 0 && best2 > 0.35f * m[b1])
                    hist[((b2 - b1) % (int)N + (int)N) % (int)N]++;
            }
        }
        int G = 0, cnt = 0;
        for (int g = 3; g <= (int)N - 3; g++)
            if (hist[g] > cnt) { cnt = hist[g]; G = g; }
        int thresh = nsym / 5 > 4 ? nsym / 5 : 4;
        return cnt >= thresh ? G : 0;
    }

    // Companion chase: re-decode with each bin's metric summed with its
    // same-data copy. G learned from the frame (flood collision) with G=1
    // (adjacent-bin smear) as fallback; each tried at shifts 0/-1/+1.
    int companion_chase(std::vector<uint8_t> &bytes, uint32_t total_bytes) {
        const uint32_t n_nib = 2 * total_bytes;
        if (block_recs.empty() || nibbles.size() < n_nib) return 0;
        int G_hist = estimate_companion_offset();
        const int GS[2] = {G_hist, 1};
        static const int SHIFTS[] = {0, -1, 1};
        std::vector<uint8_t> trial_nib(nibbles.begin(), nibbles.begin() + n_nib);
        std::vector<uint8_t> trial_bytes(total_bytes), blk_nib;
        for (int G : GS) {
            if (G == 0) continue;
            for (int shift : SHIFTS) {
                for (const BlockRec &rec : block_recs) {
                    if (rec.nib_start >= (int)n_nib) break;
                    redecode_block(rec, 0, shift, blk_nib, G);
                    for (int i = 0; i < rec.sf_app; i++) {
                        int p = rec.nib_start + i;
                        if (p >= 0 && p < (int)n_nib) trial_nib[p] = blk_nib[i];
                    }
                }
                nibbles_to_bytes(trial_nib.data(), trial_bytes.data(), total_bytes);
                if (payload_crc_ok(trial_bytes.data())) {
                    bytes = trial_bytes;
                    return G;
                }
            }
        }
        return 0;
    }

    // Straddle chase: when sync locks between two time-offset copies of the
    // same frame (simultaneous flood rebroadcasts ~ms apart), every payload
    // window straddles two symbols nearly 50/50 - each window's runner-up
    // peak is the NEXT symbol's main peak. A symbol's energy is then split
    // between its own window and the previous one at the SAME bin, so
    // re-deciding on the sum of adjacent windows' magnitudes reconstitutes
    // full symbol energy while the straddle ghost keeps only a partial share.
    // dir=+1 sums with the previous window, dir=-1 with the next.
    int straddle_chase(std::vector<uint8_t> &bytes, uint32_t total_bytes) {
        const uint32_t n_nib = 2 * total_bytes;
        if (block_recs.empty() || nibbles.size() < n_nib) return 0;
        // flat list of all stored payload symbols, in order
        std::vector<const std::vector<float> *> sym;
        for (const auto &rec : block_recs)
            for (int s = 0; s < rec.cw_len; s++) sym.push_back(&rec.mags[s]);
        std::vector<uint8_t> trial_nib(nibbles.begin(), nibbles.begin() + n_nib);
        std::vector<uint8_t> trial_bytes(total_bytes), blk_nib;
        for (int dir : {1, -1}) {
            int base = 0;
            for (const BlockRec &rec : block_recs) {
                if (rec.nib_start >= (int)n_nib) break;
                BlockRec comb;
                comb.nib_start = rec.nib_start;
                comb.cw_len = rec.cw_len;
                comb.sf_app = rec.sf_app;
                comb.cr_app = rec.cr_app;
                comb.mags.resize(rec.cw_len);
                for (int s = 0; s < rec.cw_len; s++) {
                    int gi = base + s, oi = gi - dir;
                    comb.mags[s] = rec.mags[s];
                    if (oi >= 0 && oi < (int)sym.size())
                        for (uint32_t n = 0; n < N; n++) comb.mags[s][n] += (*sym[oi])[n];
                }
                redecode_block(comb, 0, 0, blk_nib);
                for (int i = 0; i < rec.sf_app; i++) {
                    int p = rec.nib_start + i;
                    if (p >= 0 && p < (int)n_nib) trial_nib[p] = blk_nib[i];
                }
                base += rec.cw_len;
            }
            nibbles_to_bytes(trial_nib.data(), trial_bytes.data(), total_bytes);
            if (payload_crc_ok(trial_bytes.data())) {
                bytes = trial_bytes;
                return dir;
            }
        }
        return 0;
    }

    // Symbol-level chase: hypothesize that the strongest FFT peak of a few
    // weak symbols was a colliding transmission, take each one's runner-up
    // peak instead, re-decode the affected block(s) and re-check the CRC.
    // Tries single symbols first, then pairs of the weakest ones.
    int symbol_chase(std::vector<uint8_t> &bytes, uint32_t total_bytes) {
        const uint32_t n_nib = 2 * total_bytes;
        if (block_recs.empty() || nibbles.size() < n_nib) return 0;

        struct SymCand { float dominance; int blk, sym; };
        std::vector<SymCand> sc;
        for (int b = 0; b < (int)block_recs.size(); b++) {
            const BlockRec &rec = block_recs[b];
            if (rec.nib_start >= (int)n_nib) break; // padding blocks
            for (int s = 0; s < rec.cw_len; s++) {
                const auto &m = rec.mags[s];
                float m1 = -1, m2 = -1;
                for (float v : m) {
                    if (v > m1) { m2 = m1; m1 = v; }
                    else if (v > m2) m2 = v;
                }
                if (m1 <= 0) continue;
                sc.push_back({m2 / m1, b, s});
            }
        }
        // weakest decisions first (runner-up peak nearly as strong as winner)
        std::sort(sc.begin(), sc.end(),
                  [](const SymCand &a, const SymCand &b) { return a.dominance > b.dominance; });

        constexpr int SC_SINGLES = 24, SC_PAIRS_OF = 12;
        std::vector<uint8_t> trial_nib(nibbles.begin(), nibbles.begin() + n_nib);
        std::vector<uint8_t> trial_bytes(total_bytes), blk_nib;

        auto apply_block = [&](int b, uint32_t mask) {
            const BlockRec &rec = block_recs[b];
            redecode_block(rec, mask, 0, blk_nib);
            for (int i = 0; i < rec.sf_app; i++) {
                int p = rec.nib_start + i;
                if (p >= 0 && p < (int)n_nib) trial_nib[p] = blk_nib[i];
            }
        };
        auto restore_block = [&](int b) {
            const BlockRec &rec = block_recs[b];
            for (int i = 0; i < rec.sf_app; i++) {
                int p = rec.nib_start + i;
                if (p >= 0 && p < (int)n_nib) trial_nib[p] = nibbles[p];
            }
        };
        auto crc_pass = [&]() {
            nibbles_to_bytes(trial_nib.data(), trial_bytes.data(), total_bytes);
            if (payload_crc_ok(trial_bytes.data())) { bytes = trial_bytes; return true; }
            return false;
        };

        int n_single = std::min<int>(SC_SINGLES, (int)sc.size());
        for (int i = 0; i < n_single; i++) {
            apply_block(sc[i].blk, 1u << sc[i].sym);
            if (crc_pass()) return 1;
            restore_block(sc[i].blk);
        }
        int n_pair = std::min<int>(SC_PAIRS_OF, (int)sc.size());
        for (int i = 0; i < n_pair; i++) {
            for (int j = i + 1; j < n_pair; j++) {
                if (sc[i].blk == sc[j].blk) {
                    apply_block(sc[i].blk, (1u << sc[i].sym) | (1u << sc[j].sym));
                } else {
                    apply_block(sc[i].blk, 1u << sc[i].sym);
                    apply_block(sc[j].blk, 1u << sc[j].sym);
                }
                if (crc_pass()) return 2;
                restore_block(sc[i].blk);
                restore_block(sc[j].blk);
            }
        }
        return 0;
    }

    // Diagnostic (LORA_DEBUG_FAIL=1): per-symbol peak stats of a frame that
    // failed CRC even after chase, to identify the corruption mechanism
    // (fade = dominance collapse, timing slip = fractional offset drift).
    void dump_failed_frame() const {
        std::lock_guard<std::mutex> lk(g_print_mtx);
        fprintf(stderr, "%s  [faildump] pay_len=%u snr=%.1f blocks=%zu companion_G=%d\n",
                tag.c_str(), m_pay_len, current_snr, block_recs.size(),
                estimate_companion_offset());
        int sym_i = 0;
        for (const auto &rec : block_recs) {
            for (int s = 0; s < rec.cw_len; s++, sym_i++) {
                std::vector<float> m = rec.mags[s];
                uint32_t mx = (uint32_t)(std::max_element(m.begin(), m.end()) - m.begin());
                float m1 = m[mx];
                // top-3 peaks with +-2 exclusion zones, like the offline analysis
                int b2 = -1, b3 = -1;
                {
                    std::vector<float> mm = m;
                    for (int k = -2; k <= 2; k++) mm[(mx + N + k) % N] = 0;
                    b2 = (int)(std::max_element(mm.begin(), mm.end()) - mm.begin());
                    float a2 = mm[b2];
                    for (int k = -2; k <= 2; k++) mm[(b2 + N + k) % N] = 0;
                    b3 = (int)(std::max_element(mm.begin(), mm.end()) - mm.begin());
                    fprintf(stderr, "%s  [faildump] sym=%3d bin=%3u dom=%.2f p2=%3d:%.2f p3=%3d:%.2f\n",
                            tag.c_str(), sym_i, mx, a2 / m1, b2, a2 / m1, b3, mm[b3] / m1);
                }
            }
        }
    }

    // Re-demodulate the whole frame (header included) from the retained ring
    // samples with the symbol windows shifted by `delta` samples. Rescues
    // frames where sync locked between two time-offset copies of the same
    // flood frame: every window straddles two symbols, poisoning both the
    // header (wrong pay_len that passes the 5-bit checksum by luck) and the
    // payload. Re-parsing the header under the shifted timing recovers the
    // TRUE frame parameters. Returns true and fills the outputs on CRC pass.
    bool redecode_frame(long delta, std::vector<uint8_t> &out_payload, bool &out_has_crc) {
        std::vector<cx> win(N);
        std::vector<std::array<double, MAX_SF>> blk;
        std::vector<uint8_t> nibs;
        bool hdr = true, ldro = false, pcrc = cfg.has_crc;
        int pcr = cfg.cr;
        uint32_t plen = 0, total_syms = 8;
        float scum = frame_sfo_cum0;
        long start = (long)frame_start_abs + delta;
        long sto_off = my_roundf(m_sto_frac * os);
        for (uint32_t s = 0; s < total_syms && s < 4096; s++) {
            for (uint32_t ii = 0; ii < N; ii++) {
                long idx = start + os / 2 + (long)os * ii - sto_off;
                if (!sample_at_abs(idx, win[ii])) return false;
            }
            blk.push_back(compute_LLRs_map(win.data(), hdr || ldro, (int)(delta / os)));
            int block_size = 4 + (hdr ? 4 : pcr);
            if ((int)blk.size() == block_size) {
                int sf_app = (hdr || ldro) ? cfg.sf - 2 : cfg.sf;
                int cw_len = block_size;
                int cr_app = hdr ? 4 : pcr;
                double deinter[MAX_SF][8];
                for (int i = 0; i < cw_len; i++)
                    for (int j = 0; j < sf_app; j++)
                        deinter[mod(i - j - 1, sf_app)][i] = blk[i][cfg.sf - sf_app + j];
                for (int i = 0; i < sf_app; i++)
                    nibs.push_back(hamming_decode_soft(deinter[i], cr_app));
                blk.clear();
                if (hdr && nibs.size() >= 5) {
                    if (!header_consistent(nibs.data())) return false;
                    plen = (nibs[0] << 4) | nibs[1];
                    pcrc = nibs[2] & 1;
                    pcr = nibs[2] >> 1;
                    if (pcr < 1) pcr = 1;
                    if (pcr > 4) pcr = 4;
                    if (plen == 0 || plen > 255) return false;
                    // Without a payload CRC there is nothing to validate a
                    // retimed hypothesis against - the 5-bit header checksum
                    // alone would accept garbage. Only pursue CRC'd frames.
                    if (!pcrc) return false;
                    ldro = ((float)(1u << cfg.sf) * 1e3f / cfg.bw) > LDRO_MAX_DURATION_MS;
                    int sf2 = cfg.sf - 2 * (ldro ? 1 : 0);
                    total_syms = 8 + (int)ceil((double)(2 * plen - cfg.sf + 2 + 5 + (pcrc ? 4 : 0)) / sf2) * (4 + pcr);
                    hdr = false;
                    nibs.erase(nibs.begin(), nibs.begin() + 5);
                }
            }
            long consume = sps;
            if (fabsf(scum) > 1.0f / 2.0f / os) {
                int sg = (scum >= 0) ? 1 : -1;
                consume -= sg;
                scum -= sg * 1.0f / os;
            }
            scum += sfo_hat;
            start += consume;
        }
        if (hdr) return false;
        uint32_t need = plen * 2 + (pcrc ? 4 : 0);
        if (nibs.size() < need) return false;
        uint32_t tb = plen + (pcrc ? 2 : 0);
        std::vector<uint8_t> bytes(tb);
        for (uint32_t i = 0; i < tb; i++) {
            uint8_t lo = nibs[2 * i], hi = nibs[2 * i + 1];
            if (i < plen) {
                lo ^= whitening_seq[i] & 0x0F;
                hi ^= (whitening_seq[i] >> 4) & 0x0F;
            }
            bytes[i] = (hi << 4) | lo;
        }
        if (pcrc) {
            if (plen < 2) return false;
            uint16_t calc = crc16(bytes.data(), plen - 2);
            calc = calc ^ bytes[plen - 1] ^ ((uint16_t)bytes[plen - 2] << 8);
            uint16_t rx = bytes[plen] | ((uint16_t)bytes[plen + 1] << 8);
            if (calc != rx) return false;
        }
        out_payload.assign(bytes.begin(), bytes.begin() + plen);
        out_has_crc = pcrc;
        return true;
    }

    // Scan retimings, nearest first. Half-symbol reach in steps of N/16 chips.
    long retime_chase(std::vector<uint8_t> &out_payload, bool &out_has_crc) {
        for (uint32_t step = N / 16; step <= N / 2; step += N / 16) {
            for (int sgn : {1, -1}) {
                long delta = (long)sgn * (long)step * os;
                if (redecode_frame(delta, out_payload, out_has_crc))
                    return delta;
            }
        }
        return 0;
    }

    void assemble_and_output() {
        uint32_t total_bytes = m_pay_len + (m_pay_has_crc ? 2 : 0);
        std::vector<uint8_t> bytes(total_bytes);
        nibbles_to_bytes(nibbles.data(), bytes.data(), total_bytes);

        Packet pkt;
        pkt.has_crc = m_pay_has_crc;
        pkt.snr_est = current_snr;
        pkt.crc_valid = false;
        int chase_n = 0;

        if (m_pay_has_crc && m_pay_len >= 2) {
            pkt.crc_valid = payload_crc_ok(bytes.data());
            if (!pkt.crc_valid && cfg.soft_decoding) {
                // Recovery ladder, most-trustworthy hypothesis class first:
                // rotation (4 trials) -> nibble chase (<=2000) -> symbol chase.
                if (int rot = rotation_chase(bytes, total_bytes)) {
                    pkt.crc_valid = true;
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    fprintf(stderr, "%s  Rotation chase: CRC ok with %+d bin shift\n",
                            tag.c_str(), rot);
                } else if ((chase_n = chase_recover(bytes, total_bytes)) > 0) {
                    pkt.crc_valid = true;
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    fprintf(stderr, "%s  Chase recovery: CRC ok after %d nibble substitution(s)\n",
                            tag.c_str(), chase_n);
                } else if (int sym_n = symbol_chase(bytes, total_bytes)) {
                    pkt.crc_valid = true;
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    fprintf(stderr, "%s  Symbol chase: CRC ok after %d symbol substitution(s)\n",
                            tag.c_str(), sym_n);
                } else if (int G = companion_chase(bytes, total_bytes)) {
                    pkt.crc_valid = true;
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    fprintf(stderr, "%s  Companion chase: CRC ok combining bins %+d apart\n",
                            tag.c_str(), G);
                } else if (int sd = straddle_chase(bytes, total_bytes)) {
                    pkt.crc_valid = true;
                    std::lock_guard<std::mutex> lk(g_print_mtx);
                    fprintf(stderr, "%s  Straddle chase: CRC ok summing adjacent windows (dir %+d)\n",
                            tag.c_str(), sd);
                } else {
                    std::vector<uint8_t> rt_payload;
                    bool rt_crc = true;
                    if (long delta = retime_chase(rt_payload, rt_crc)) {
                        // re-demod may have recovered a DIFFERENT (true)
                        // header, so adopt its length/CRC flag wholesale
                        m_pay_len = (uint32_t)rt_payload.size();
                        m_pay_has_crc = rt_crc;
                        bytes = rt_payload;
                        pkt.has_crc = rt_crc;
                        pkt.crc_valid = true;
                        std::lock_guard<std::mutex> lk(g_print_mtx);
                        fprintf(stderr, "%s  Retime chase: CRC ok at %+ld sample shift (pay_len=%u)\n",
                                tag.c_str(), delta, m_pay_len);
                    } else if (getenv("LORA_DEBUG_FAIL")) {
                        dump_failed_frame();
                    }
                }
            }
        } else if (!m_pay_has_crc) {
            pkt.crc_valid = true; // no CRC to check
        }
        pkt.payload.assign(bytes.begin(), bytes.begin() + m_pay_len);

        {
            std::lock_guard<std::mutex> lk(g_print_mtx);

            if (!stdout_tag.empty())
                printf("rx cfg: %s\n", stdout_tag.c_str());
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
                    printf("CRC valid!\n");
                else
                    printf("\033[31mCRC invalid\033[0m\n");
            } else {
                printf("(no CRC)\n");
            }

            // Machine-readable line for CRC-valid packets only: plain hex,
            // ready to pipe into a packet decoder (e.g. meshcore-decoder stream)
            if (pkt.crc_valid) {
                printf("rx ok: ");
                for (uint32_t i = 0; i < m_pay_len; i++)
                    printf("%02x", pkt.payload[i]);
                printf("\n");
            }
            printf("\n");

            fflush(stdout);
        }

        last_packets.push_back(pkt);
        nibbles.clear();
        nib_alt.clear();
        nib_margin.clear();
        block_recs.clear();
        pending_mags.clear();
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
    DemodWorker(const LoRaConfig &c, const std::string &tag, const std::string &stdout_tag)
        : demod(c, tag, stdout_tag), tag_(tag) {
        th = std::thread([this] { run(); });
    }

    // Guard: a worker must never be destroyed with its thread still running,
    // or teardown deadlocks (e.g. early error exits in main).
    ~DemodWorker() { finish(); }

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

// Stop all workers and join their threads. Must run on every exit path once
// workers exist, or the process hangs on the live threads during teardown.
static void finish_workers() {
    for (auto &w : g_workers)
        w->finish();
}
static FILE *g_dump = nullptr;

// One entry per RF channel: frequency-shifts the shared stream to put that
// channel at baseband, then fans out to the workers decoding it. This is how
// one dongle covers several LoRa channels at once (e.g. -C f1,f2 at 1 MS/s).
struct ChannelFanout {
    uint32_t freq = 0;               // channel center frequency (Hz)
    double phase = 0, phase_inc = 0; // oscillator state (radians, radians/sample)
    bool identity = true;            // channel is at tuner center, no shift needed
    std::vector<DemodWorker *> workers;
};
static std::vector<ChannelFanout> g_channels;

static void broadcast_u8(const unsigned char *buf, uint32_t len) {
    if (g_dump && fwrite(buf, 1, len, g_dump) != len) {
        // Disk full or write error: a silently truncated dump is worse than
        // no dump. Warn, stop dumping, keep receiving.
        fprintf(stderr, "Warning: IQ dump write failed (disk full?) - dump disabled, reception continues\n");
        fclose(g_dump);
        g_dump = nullptr;
    }

    size_t n = len / 2;
    auto base = std::make_shared<std::vector<cx>>(n);
    {
        auto &v = *base;
        for (size_t i = 0; i < n; i++)
            v[i] = cx(g_u8lut[buf[2 * i]], g_u8lut[buf[2 * i + 1]]);
    }

    for (auto &ch : g_channels) {
        std::shared_ptr<const std::vector<cx>> chunk;
        if (ch.identity) {
            chunk = base;
        } else {
            auto shifted = std::make_shared<std::vector<cx>>(n);
            const auto &v = *base;
            auto &s = *shifted;
            // Oscillator restarts from the double-precision phase each chunk,
            // so float rounding never accumulates across the stream
            cx osc((float)cos(ch.phase), (float)sin(ch.phase));
            const cx step((float)cos(ch.phase_inc), (float)sin(ch.phase_inc));
            for (size_t i = 0; i < n; i++) {
                s[i] = v[i] * osc;
                osc *= step;
            }
            ch.phase = fmod(ch.phase + ch.phase_inc * (double)n, 2.0 * M_PI);
            chunk = std::move(shifted);
        }
        for (auto *w : ch.workers)
            w->push(chunk);
    }
}

static void rtlsdr_callback(unsigned char *buf, uint32_t len, void *) {
    if (!g_running) return;
    broadcast_u8(buf, len);
}

// ============================================================================
// Main
// ============================================================================

static std::string make_tag(uint32_t chan_freq, uint8_t sf, uint32_t bw) {
    char bwbuf[16], buf[48];
    if (bw % 1000 == 0)
        snprintf(bwbuf, sizeof(bwbuf), "%uk", bw / 1000);
    else
        snprintf(bwbuf, sizeof(bwbuf), "%.1fk", bw / 1000.0);
    if (chan_freq)
        snprintf(buf, sizeof(buf), "[%.4fM SF%u/%s] ", chan_freq / 1e6, sf, bwbuf);
    else
        snprintf(buf, sizeof(buf), "[SF%u/%s] ", sf, bwbuf);
    return buf;
}

// Parse a comma-separated list of integers ("7,8" / "869432000,869618000")
static std::vector<long> parse_list(const char *s) {
    std::vector<long> out;
    while (*s) {
        char *end;
        long v = strtol(s, &end, 10);
        if (end == s) break;
        out.push_back(v);
        s = (*end == ',') ? end + 1 : end;
    }
    return out;
}

int main(int argc, char *argv[]) {
    LoRaConfig cfg;
    const char *iq_file = nullptr;
    const char *dump_file = nullptr;
    bool auto_mode = false, freq_given = false, samp_given = false;
    std::vector<long> sfs_cli, bws_cli, chans_cli;
    int opt;

    while ((opt = getopt(argc, argv, "f:s:b:S:c:w:g:GTp:IL:r:D:AC:")) != -1) {
        switch (opt) {
        case 'f': cfg.freq = (uint32_t)atol(optarg); freq_given = true; break;
        case 's': cfg.samp_rate = (uint32_t)atol(optarg); samp_given = true; break;
        case 'b': bws_cli = parse_list(optarg); break;
        case 'S': sfs_cli = parse_list(optarg); break;
        case 'c': cfg.cr = (uint8_t)atoi(optarg); break;
        case 'w': cfg.sync_word = (uint16_t)strtol(optarg, nullptr, 16); break;
        case 'g': cfg.gain = atoi(optarg); break;
        case 'G': cfg.rtl_agc = true; break;
        case 'T': cfg.tuner_agc = true; break;
        case 'p': cfg.ppm = atoi(optarg); break;
        case 'I': cfg.impl_head = true; break;
        case 'L': cfg.pay_len = (uint32_t)atoi(optarg); break;
        case 'r': iq_file = optarg; break;
        case 'D': dump_file = optarg; break;
        case 'A': auto_mode = true; break;
        case 'C': chans_cli = parse_list(optarg); break;
        default:
            fprintf(stderr, "Usage: %s [-f tuner_freq] [-s samp_rate] [-b bw[,bw..]] [-S sf[,sf..]] [-c cr]\n"
                            "          [-w sync_word_hex] [-g gain] [-G] [-p ppm] [-I] [-L pay_len] [-A]\n"
                            "          [-C chan_freq[,chan_freq..]] [-r iq_file] [-D dump_file]\n"
                            "  -G enables the RTL2832 internal digital AGC\n"
                            "  -T enables the R820T analog tuner AGC (overrides -g)\n", argv[0]);
            return 1;
        }
    }

    for (long sf : sfs_cli) {
        if (sf < MIN_SF || sf > MAX_SF) {
            fprintf(stderr, "SF must be between %d and %d\n", MIN_SF, MAX_SF);
            return 1;
        }
    }
    if (!sfs_cli.empty()) cfg.sf = (uint8_t)sfs_cli[0];
    if (!bws_cli.empty()) cfg.bw = (uint32_t)bws_cli[0];
    if (cfg.cr < 1 || cfg.cr > 4) {
        uint8_t clamped = cfg.cr < 1 ? 1 : 4;
        fprintf(stderr, "Warning: CR %u out of range 1-4, using %u (explicit header overrides it anyway)\n",
                cfg.cr, clamped);
        cfg.cr = clamped;
    }

    // Channel list: default is a single channel at the tuner frequency
    std::vector<uint32_t> channels;
    for (long c : chans_cli) channels.push_back((uint32_t)c);
    if (channels.empty()) channels = {cfg.freq};

    // Multi-channel needs enough span; default to 1 MS/s unless -s was given
    if (channels.size() > 1 && !samp_given)
        cfg.samp_rate = 1000000;

    // Tuner center: explicit -f wins; otherwise midpoint of the channels
    // (which places the RTL-SDR DC spike between them, outside every channel)
    uint32_t tuner_freq;
    if (freq_given || chans_cli.empty()) {
        tuner_freq = cfg.freq;
    } else {
        uint64_t lo = *std::min_element(channels.begin(), channels.end());
        uint64_t hi = *std::max_element(channels.begin(), channels.end());
        tuner_freq = (uint32_t)((lo + hi) / 2);
    }

    if (cfg.samp_rate % cfg.bw != 0) {
        fprintf(stderr, "Warning: samp_rate should be an integer multiple of bandwidth\n");
    }

    for (int i = 0; i < 256; i++)
        g_u8lut[i] = ((float)i - 127.5f) / 127.5f;

    // Build the list of (SF, BW) decoder configurations
    std::vector<uint8_t> sfs;
    for (long s : sfs_cli) sfs.push_back((uint8_t)s);
    std::vector<uint32_t> bws;
    for (long b : bws_cli) bws.push_back((uint32_t)b);
    if (sfs.empty()) sfs = auto_mode ? std::vector<uint8_t>{7, 8, 9, 10, 11, 12} : std::vector<uint8_t>{cfg.sf};
    if (bws.empty()) {
        if (auto_mode) {
            for (uint32_t bw : {62500u, 125000u, 250000u, 500000u}) {
                if (bw <= cfg.samp_rate && cfg.samp_rate % bw == 0 && cfg.samp_rate / bw <= 16)
                    bws.push_back(bw);
            }
        }
        if (bws.empty()) bws = {cfg.bw};
    }

    uint32_t max_bw = *std::max_element(bws.begin(), bws.end());
    bool multi = channels.size() * sfs.size() * bws.size() > 1;
    bool multi_chan = channels.size() > 1;

    for (uint32_t chf : channels) {
        double off = (double)chf - (double)tuner_freq;
        if (fabs(off) + max_bw / 2.0 > 0.48 * cfg.samp_rate)
            fprintf(stderr, "Warning: channel %u Hz (offset %+.0f Hz) is at/beyond the Nyquist edge for samp_rate %u\n",
                    chf, off, cfg.samp_rate);

        g_channels.push_back({});
        ChannelFanout &ch = g_channels.back();
        ch.freq = chf;
        ch.identity = (off == 0.0);
        ch.phase_inc = -2.0 * M_PI * off / (double)cfg.samp_rate;

        for (uint32_t bw : bws) {
            for (uint8_t sf : sfs) {
                LoRaConfig c = cfg;
                c.freq = chf;
                c.sf = sf;
                c.bw = bw;
                c.compute_derived();
                std::string tag, stag;
                if (multi) {
                    tag = make_tag(multi_chan ? chf : 0, sf, bw);
                    char sbuf[64];
                    snprintf(sbuf, sizeof(sbuf), "freq=%u sf=%u bw=%u", chf, sf, bw);
                    stag = sbuf;
                }
                g_workers.emplace_back(new DemodWorker(c, tag, stag));
                ch.workers.push_back(g_workers.back().get());
            }
        }
    }

    cfg.compute_derived();
    fprintf(stderr, "LoRa Standalone RX\n");
    if (multi_chan) {
        fprintf(stderr, "  Tuner freq: %u Hz\n", tuner_freq);
        fprintf(stderr, "  Channels:  ");
        for (uint32_t chf : channels)
            fprintf(stderr, " %u (%+.0f Hz)", chf, (double)chf - (double)tuner_freq);
        fprintf(stderr, "\n");
    } else {
        fprintf(stderr, "  Freq:       %u Hz\n", tuner_freq);
    }
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
            finish_workers();
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
            finish_workers();
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
            finish_workers();
            return 1;
        }
        fprintf(stderr, "Found %d RTL-SDR device(s)\n", n_devices);

        if (rtlsdr_open(&dev, dev_index) < 0) {
            fprintf(stderr, "Failed to open RTL-SDR device\n");
            finish_workers();
            return 1;
        }
        g_dev = dev;

        rtlsdr_set_sample_rate(dev, cfg.samp_rate);
        rtlsdr_set_center_freq(dev, tuner_freq);
        rtlsdr_set_freq_correction(dev, cfg.ppm);
        if (cfg.tuner_agc) {
            rtlsdr_set_tuner_gain_mode(dev, 0); // analog tuner AGC
        } else {
            rtlsdr_set_tuner_gain_mode(dev, 1); // manual gain
            rtlsdr_set_tuner_gain(dev, cfg.gain);
        }
        // RTL2832 internal digital AGC: scales a weak analog signal up to use
        // the 8-bit ADC range. Essential when the tuner output is so low that
        // the baseband RMS is only a couple of LSB (-G to enable).
        rtlsdr_set_agc_mode(dev, cfg.rtl_agc ? 1 : 0);
        rtlsdr_reset_buffer(dev);

        fprintf(stderr, "RTL-SDR configured. Tuner gain: %.1f dB\n", rtlsdr_get_tuner_gain(dev) / 10.0);
        fprintf(stderr, "Listening...\n\n");

        // 16384 bytes = 8192 IQ samples per callback
        rtlsdr_read_async(dev, rtlsdr_callback, nullptr, 0, 16384);

        rtlsdr_close(dev);
    }

    size_t total = 0;
    finish_workers();
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
