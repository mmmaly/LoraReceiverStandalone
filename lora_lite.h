// lora_lite - embedded-friendly LoRa demodulator core.
//
// A no-malloc, no-thread, single-SF/BW reduction of the lora_rx demodulation
// chain, sized entirely with compile-time static buffers so it can run inside
// a Cortex-M4 baseband image (PortaPack/Mayhem) as well as on a host.
//
// It expects IQ already decimated to os = LORA_LITE_OS samples per chip
// (i.e. sample_rate = os * bandwidth), fed as interleaved complex floats.
// Call feed() with each block of samples; decoded packets are delivered via a
// callback. Nothing here allocates or blocks.
//
// Derived from lora_rx.cpp (which is in turn based on gr-lora_sdr). The DSP is
// bit-for-bit the same; only the memory model differs.
#pragma once
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <cmath>

// Defaults chosen to fit the PortaPack M4 baseband SRAM budget (96 KiB data,
// shared with the RTOS). SF up to 8 (N=256) covers MeshCore SF7/SF8; a host
// build can raise LORA_LITE_MAX_SF for wider coverage.
// NOTE: these must be defined before including lora_kiss.h, which sizes its
// twiddle table from LORA_LITE_MAX_SF.
#ifndef LORA_LITE_MAX_SF
#define LORA_LITE_MAX_SF 8
#endif
#ifndef LORA_LITE_OS
#define LORA_LITE_OS 2               // samples per chip after decimation
#endif
#ifndef LORA_LITE_MAX_PAYLOAD
#define LORA_LITE_MAX_PAYLOAD 64     // bytes; MeshCore frames fit easily
#endif

#include "lora_kiss.h"   // tiny self-contained radix-2 FFT

namespace lora_lite {

static constexpr int   MAX_N   = 1 << LORA_LITE_MAX_SF;
static constexpr int   OS      = LORA_LITE_OS;
static constexpr int   MAX_SPS = MAX_N * OS;
static constexpr float PI_F    = 3.14159265358979323846f;
// Preamble capture depth: N_UP_REQ symbols raw, UP_USE for CFO/STO estimation.
static constexpr int   PRE_SYMS = 5;   // == N_UP_REQ
static constexpr int   CFO_SYMS = 4;   // == UP_USE

struct Cx { float r, i; };
static inline Cx cxmul(Cx a, Cx b) { return {a.r * b.r - a.i * b.i, a.r * b.i + a.i * b.r}; }
static inline Cx cxconj(Cx a) { return {a.r, -a.i}; }
static inline float cxnorm(Cx a) { return a.r * a.r + a.i * a.i; }
static inline long lmod(long a, long b) { return ((a % b) + b) % b; }

// Decoded packet (or detection) handed to the caller.
enum PacketStatus : uint8_t {
    PKT_HEADER_ONLY = 0,  // valid header decoded; payload not yet/never assembled
    PKT_CRC_FAIL    = 1,  // full payload assembled but CRC failed
    PKT_DECODED     = 2,  // full payload assembled, CRC valid (or no CRC)
};

struct Packet {
    uint8_t  data[LORA_LITE_MAX_PAYLOAD];
    uint8_t  len;         // payload bytes present in data[] (0 for header-only)
    uint8_t  hdr_len;     // payload length declared by the header
    uint8_t  status;      // PacketStatus
    bool     crc_ok;
    bool     has_crc;
    float    snr;
    float    cfo;         // in bins
    uint8_t  sync;        // detected sync word
};

// ---- Hamming / whitening / CRC tables (shared with lora_rx) ----
static const uint8_t WHITEN[255] = {
    0xFF,0xFE,0xFC,0xF8,0xF0,0xE1,0xC2,0x85,0x0B,0x17,0x2F,0x5E,0xBC,0x78,0xF1,0xE3,
    0xC6,0x8D,0x1A,0x34,0x68,0xD0,0xA0,0x40,0x80,0x01,0x02,0x04,0x08,0x11,0x23,0x47,
    0x8E,0x1C,0x38,0x71,0xE2,0xC4,0x89,0x12,0x25,0x4B,0x97,0x2E,0x5C,0xB8,0x70,0xE0,
    0xC0,0x81,0x03,0x06,0x0C,0x19,0x32,0x64,0xC9,0x92,0x24,0x49,0x93,0x26,0x4D,0x9B,
    0x37,0x6E,0xDC,0xB9,0x72,0xE4,0xC8,0x90,0x20,0x41,0x82,0x05,0x0A,0x15,0x2B,0x56,
    0xAD,0x5B,0xB6,0x6D,0xDA,0xB5,0x6B,0xD6,0xAC,0x59,0xB2,0x65,0xCB,0x96,0x2C,0x58,
    0xB0,0x61,0xC3,0x87,0x0F,0x1F,0x3E,0x7D,0xFB,0xF6,0xED,0xDB,0xB7,0x6F,0xDE,0xBD,
    0x7A,0xF5,0xEB,0xD7,0xAE,0x5D,0xBA,0x74,0xE8,0xD1,0xA2,0x44,0x88,0x10,0x21,0x43,
    0x86,0x0D,0x1B,0x36,0x6C,0xD8,0xB1,0x63,0xC7,0x8F,0x1E,0x3C,0x79,0xF3,0xE7,0xCE,
    0x9C,0x39,0x73,0xE6,0xCC,0x98,0x31,0x62,0xC5,0x8B,0x16,0x2D,0x5A,0xB4,0x69,0xD2,
    0xA4,0x48,0x91,0x22,0x45,0x8A,0x14,0x29,0x52,0xA5,0x4A,0x95,0x2A,0x54,0xA9,0x53,
    0xA7,0x4E,0x9D,0x3B,0x77,0xEE,0xDD,0xBB,0x76,0xEC,0xD9,0xB3,0x67,0xCF,0x9E,0x3D,
    0x7B,0xF7,0xEF,0xDF,0xBF,0x7E,0xFD,0xFA,0xF4,0xE9,0xD3,0xA6,0x4C,0x99,0x33,0x66,
    0xCD,0x9A,0x35,0x6A,0xD4,0xA8,0x51,0xA3,0x46,0x8C,0x18,0x30,0x60,0xC1,0x83,0x07,
    0x0E,0x1D,0x3A,0x75,0xEA,0xD5,0xAA,0x55,0xAB,0x57,0xAF,0x5F,0xBE,0x7C,0xF9,0xF2,
    0xE5,0xCA,0x94,0x28,0x50,0xA1,0x42,0x84,0x09,0x13,0x27,0x4F,0x9F,0x3F,0x7F};

static inline uint16_t crc16_ccitt(const uint8_t* d, int len) {
    uint16_t crc = 0;
    for (int i = 0; i < len; i++) {
        uint8_t b = d[i];
        for (int k = 0; k < 8; k++) {
            if (((crc & 0x8000) >> 8) ^ (b & 0x80)) crc = (crc << 1) ^ 0x1021;
            else crc = (crc << 1);
            b <<= 1;
        }
    }
    return crc;
}

// (A hard-decision Hamming decoder lived here; lora_lite uses only the soft
//  path, so it was removed to keep the M4 baseband image within its 32 KiB.)

// ============================================================================
// The demodulator
// ============================================================================
class Demod {
public:
    typedef void (*PacketCb)(const Packet&, void* user);

    void init(int sf, int cr, bool has_crc, uint16_t sync_word,
              PacketCb cb, void* user, uint32_t bw = 62500, uint32_t freq = 869618000) {
        sf_ = sf; N_ = 1 << sf; sps_ = N_ * OS;
        cr_ = cr; has_crc_ = has_crc;
        bw_ = bw; freq_ = freq;
        sync_words_[0] = ((sync_word & 0xF0) >> 4) << 3;
        sync_words_[1] = (sync_word & 0x0F) << 3;
        cb_ = cb; user_ = user;
        fft_.init(sf + 1);   // one plan covers both N and 2N transforms

        // reference chirps
        build_upchirp(upchirp_, 0);
        for (int i = 0; i < N_; i++) downchirp_[i] = cxconj(upchirp_[i]);

        // Gray demap LUTs
        for (int n = 0; n < N_; n++) {
            uint32_t v = (n + N_ - 1) % N_;
            demap_pay_[n] = v ^ (v >> 1);
            uint32_t h = v / 4;
            demap_hdr_[n] = h ^ (h >> 1);
        }
        reset();
        rlen_ = 0; rhead_ = 0;
    }

    // Feed decimated complex samples (os per chip). Non-blocking.
    // Pumps as the ring fills so no samples are ever silently overwritten,
    // regardless of how large a block the caller passes.
    void feed(const Cx* s, int n) {
        int i = 0;
        while (i < n) {
            while (i < n && rlen_ < RING) {
                ring_[(rhead_ + rlen_) % RING] = s[i++];
                rlen_++;
            }
            pump();
            // If the ring is full and pump() made no room (a rare stall while
            // waiting on a header), drop one symbol to keep moving.
            if (rlen_ >= RING) consume(sps_);
        }
    }

private:
    // ---- ring buffer: two symbols plus a little slack is all the pump needs ----
    static constexpr int RING = MAX_SPS * 2 + 8 * OS;
    Cx ring_[RING];
    int rhead_ = 0, rlen_ = 0;
    int processed_since_ = 0;

    Cx at(int idx) const { return ring_[(rhead_ + idx) % RING]; }
    void consume(int k) {
        if (k > rlen_) k = rlen_;
        rhead_ = (rhead_ + k) % RING; rlen_ -= k;
    }

    // ---- config ----
    int sf_, N_, sps_, cr_;
    uint32_t bw_, freq_;
    bool has_crc_;
    uint16_t sync_words_[2];
    PacketCb cb_; void* user_;

    // ---- tables/scratch ----
    static constexpr int NIB_CAP = 2 * (LORA_LITE_MAX_PAYLOAD + 2) + 16;
    Cx upchirp_[MAX_N], downchirp_[MAX_N];
    Cx demod_down_[MAX_N];
    uint16_t demap_pay_[MAX_N], demap_hdr_[MAX_N];
    Cx dech_[2 * MAX_N], fout_[2 * MAX_N];
    float mag_[2 * MAX_N];
    Cx in_down_[MAX_N], symb_[MAX_N], cfo_corr_[MAX_N];
    // The 18 KiB preamble capture is only live from DETECT through
    // finish_sync(); every recovery buffer below is only live while decoding
    // the payload. They overlay in a union so the recovery feature set adds
    // ZERO bss: the M4's usable local SRAM ends at ~80 KiB (measured on
    // hardware - the linker's 96 KiB region is not all backed), and the
    // pre-recovery Demod size already sat just under that line.
    union Overlay {
        struct {
            Cx raw[PRE_SYMS * MAX_N];
            Cx up[CFO_SYMS * MAX_N];
        } pre;
        struct {
            float llr_block_rot[8][LORA_LITE_MAX_SF];
            uint8_t nibbles_rot[NIB_CAP];
            uint8_t nib_alt[NIB_CAP];
            float nib_margin[NIB_CAP];
            float deinter[LORA_LITE_MAX_SF][8];
            float deinter_rot[LORA_LITE_MAX_SF][8];
            uint8_t trial[NIB_CAP];
            uint8_t found[LORA_LITE_MAX_PAYLOAD + 2];
            uint8_t bytes[LORA_LITE_MAX_PAYLOAD + 2];
            uint8_t rot_bytes[LORA_LITE_MAX_PAYLOAD + 2];
        } pay;
    } u_;
    static_assert(sizeof(Overlay) == PRE_SYMS * MAX_N * sizeof(Cx) + CFO_SYMS * MAX_N * sizeof(Cx),
                  "recovery buffers must fit inside the preamble capture");
    Fft fft_;

    // ---- frame-sync state (same fields as lora_rx) ----
    enum State { DETECT, SYNC, PAYLOAD } state_;
    int symbol_cnt_, preamb_[8], k_hat_, net_ids_[2], down_val_;
    uint32_t bin_idx_;
    float m_cfo_frac_, m_sto_frac_, sfo_hat_, sfo_cum_, cur_snr_;
    int m_cfo_int_, additional_up_;
    bool cfo_est_;
    static constexpr int N_UP_REQ = 5, UP_USE = 4, PREAMBLE_LEN = 8;

    // ---- payload decode state ----
    int payload_sym_cnt_;
    uint32_t total_payload_symbols_;
    bool received_head_, ldro_, is_header_;
    int pay_cr_; uint32_t pay_len_; bool pay_crc_;
    float llr_block_[8][LORA_LITE_MAX_SF];
    int llr_block_n_;
    uint8_t nibbles_[NIB_CAP];
    int nib_n_;
    // CRC-guided recovery state: shadow nibble stream decoded under a bin
    // mapping rotated by the net-ID residual (integer CFO mis-split
    // hypothesis), plus per-nibble runner-up + margin for the nibble chase.
    int nid_resid_;

    void reset() {
        state_ = DETECT; symbol_cnt_ = 1; bin_idx_ = 0; k_hat_ = 0;
        m_sto_frac_ = m_cfo_frac_ = 0; m_cfo_int_ = 0; sfo_hat_ = sfo_cum_ = 0;
        cfo_est_ = false; additional_up_ = 0; down_val_ = 0;
        net_ids_[0] = net_ids_[1] = 0;
        payload_sym_cnt_ = 0; total_payload_symbols_ = 0;
        received_head_ = false; ldro_ = false; is_header_ = true;
        llr_block_n_ = 0; nib_n_ = 0; nid_resid_ = 0;
    }

    void build_upchirp(Cx* out, uint32_t id) {
        float N = N_;
        int nfold = (int)(N - id);
        for (int n = 0; n < N_; n++) {
            float ph = n < nfold
                ? 2.0 * LORA_PI * (n * n / (2.0 * N) + ((float)id / N - 0.5) * n)
                : 2.0 * LORA_PI * (n * n / (2.0 * N) + ((float)id / N - 1.5) * n);
            out[n] = {cosf((float)ph), sinf((float)ph)};
        }
    }

    void fftN(const Cx* in, Cx* out, int n) { fft_.run(in, out, n); }

    uint32_t argmax_dechirp(const Cx* s, const Cx* ref) {
        for (int i = 0; i < N_; i++) dech_[i] = cxmul(s[i], ref[i]);
        fftN(dech_, fout_, N_);
        float mx = -1, tot = 0; uint32_t mi = 0;
        for (int i = 0; i < N_; i++) {
            float m = cxnorm(fout_[i]); tot += m;
            if (m > mx) { mx = m; mi = i; }
        }
        return tot > 0 ? mi : (uint32_t)-1;
    }

    static int most_frequent(const int* a, int n) {
        int best = a[0], bc = 0;
        for (int i = 0; i < n; i++) {
            int c = 0;
            for (int j = 0; j < n; j++) if (a[j] == a[i]) c++;
            if (c > bc) { bc = c; best = a[i]; }
        }
        return best;
    }

    // ---- main pump: process whole symbols while enough samples are buffered ----
    void pump() {
        const int need = sps_ + 2 * OS;
        while (rlen_ >= need) {
            for (int ii = 0; ii < N_; ii++) {
                long idx = OS / 2 + (long)OS * ii - lroundf(m_sto_frac_ * OS);
                if (idx < 0) idx = 0;
                if (idx >= rlen_) return;
                in_down_[ii] = at((int)idx);
            }
            long consume_n = sps_;

            if (state_ == DETECT) consume_n = do_detect();
            else if (state_ == SYNC) consume_n = do_sync();
            else consume_n = do_payload();

            if (consume_n > 0 && consume_n <= rlen_) consume((int)consume_n);
            else break;
        }
    }

    long do_detect() {
        uint32_t nb = argmax_dechirp(in_down_, downchirp_);
        long consume_n = sps_;
        if (labs(lmod(labs((long)nb - (long)bin_idx_) + 1, N_) - 1) <= 1 && nb != (uint32_t)-1) {
            if (symbol_cnt_ == 1 && bin_idx_ != (uint32_t)-1) preamb_[0] = bin_idx_;
            preamb_[symbol_cnt_] = nb;
            memcpy(&u_.pre.raw[N_ * symbol_cnt_], in_down_, N_ * sizeof(Cx));
            symbol_cnt_++;
        } else {
            memcpy(&u_.pre.raw[0], in_down_, N_ * sizeof(Cx));
            symbol_cnt_ = 1;
        }
        bin_idx_ = nb;
        if (symbol_cnt_ == N_UP_REQ) {
            additional_up_ = 0; state_ = SYNC; symbol_cnt_ = 0; cfo_est_ = false;
            k_hat_ = most_frequent(preamb_, N_UP_REQ);
            consume_n = (long)OS * (N_ - k_hat_);
        }
        return consume_n;
    }

    long do_sync() {
        if (!cfo_est_) {
            estimate_cfo_frac(&u_.pre.raw[N_ - k_hat_]);
            estimate_sto_frac();
            for (int n = 0; n < N_; n++) {
                float a = -2.0f * PI_F * m_cfo_frac_ / N_ * n;
                cfo_corr_[n] = {cosf(a), sinf(a)};
            }
            cfo_est_ = true;
        }
        for (int i = 0; i < N_; i++) symb_[i] = cxmul(in_down_[i], cfo_corr_[i]);
        bin_idx_ = argmax_dechirp(symb_, downchirp_);

        switch (symbol_cnt_) {
        case 0: // NET_ID1
            if (bin_idx_ == 0 || bin_idx_ == 1 || (int)bin_idx_ == N_ - 1) {
                if (additional_up_ < 3) additional_up_++;
            } else { symbol_cnt_ = 1; net_ids_[0] = bin_idx_; }
            break;
        case 1: symbol_cnt_ = 2; net_ids_[1] = bin_idx_; break;   // NET_ID2
        case 2: symbol_cnt_ = 3; break;                            // DOWNCHIRP1
        case 3:                                                    // DOWNCHIRP2
            down_val_ = argmax_dechirp(symb_, upchirp_);
            symbol_cnt_ = 4; break;
        case 4: return finish_sync();                              // QUARTER_DOWN
        }
        return sps_;
    }

    long finish_sync() {
        if ((uint32_t)down_val_ < (uint32_t)N_ / 2) m_cfo_int_ = down_val_ / 2;
        else m_cfo_int_ = (down_val_ - N_) / 2;

        int nid0 = (int)lmod(net_ids_[0], N_);
        uint8_t det_sync = (uint8_t)(((((nid0 + 4) >> 3) & 0xF) << 4) |
                                     ((((int)lmod(net_ids_[1], N_) + 4) >> 3) & 0xF));
        if (sync_words_[0] != 0 || sync_words_[1] != 0) {
            long d = lmod(nid0 - sync_words_[0], N_);
            long d2 = lmod(sync_words_[0] - nid0, N_);
            if ((d < d2 ? d : d2) > 2) { reset(); return sps_; }
            // Signed net-ID residual: the first net-ID symbol is a KNOWN
            // value, so its bin offset measures any constant shift the
            // payload will suffer when the CFO int/frac split is off by a
            // whole bin (the frac estimator wraps at +/-0.5). Used to decode
            // a shadow nibble stream tried when the payload CRC fails.
            nid_resid_ = (int)(d <= (long)N_ / 2 ? d : d - (long)N_);
        }
        det_sync_ = det_sync;

        // Per-symbol SFO drift rate; applied incrementally in the payload loop.
        sfo_hat_ = (float)(m_cfo_int_ + m_cfo_frac_) * (float)bw_ / (float)freq_;
        sfo_cum_ = 0;

        // demod downchirp with CFO correction
        build_upchirp(demod_down_, lmod(m_cfo_int_, N_));
        for (int i = 0; i < N_; i++) demod_down_[i] = cxconj(demod_down_[i]);
        for (int n = 0; n < N_; n++) {
            float a = -2.0f * PI_F * m_cfo_frac_ / N_ * n;
            demod_down_[n] = cxmul(demod_down_[n], (Cx){cosf(a), sinf(a)});
        }

        // SNR from CFO-corrected preamble symbol
        for (int i = 0; i < N_; i++) symb_[i] = cxmul(u_.pre.raw[i], cfo_corr_[i]);
        for (int i = 0; i < N_; i++) dech_[i] = cxmul(symb_[i], downchirp_[i]);
        fftN(dech_, fout_, N_);
        float tot = 0, sig = 0;
        for (int i = 0; i < N_; i++) { float m = cxnorm(fout_[i]); tot += m; if (m > sig) sig = m; }
        cur_snr_ = (tot > sig && sig > 0) ? 10.0f * log10f(sig / (tot - sig)) : 0.0f;

        state_ = PAYLOAD; payload_sym_cnt_ = 0; received_head_ = false; ldro_ = false;
        llr_block_n_ = 0; nib_n_ = 0; is_header_ = true;
        pay_cr_ = cr_; pay_len_ = LORA_LITE_MAX_PAYLOAD; pay_crc_ = has_crc_;

        long consume_n = sps_ / 4 + (long)OS * m_cfo_int_;
        if (consume_n < 0) consume_n = 0;
        return consume_n;
    }
    uint8_t det_sync_ = 0;

    void estimate_cfo_frac(const Cx* s) {
        // Pass 1: find each preamble symbol's peak bin (FFT into scratch).
        int k0v[UP_USE]; float k0m[UP_USE];
        for (int i = 0; i < UP_USE; i++) {
            for (int j = 0; j < N_; j++) dech_[j] = cxmul(s[N_ * i + j], downchirp_[j]);
            fftN(dech_, fout_, N_);
            float mx = -1; int mi = 0;
            for (int j = 0; j < N_; j++) { float m = cxnorm(fout_[j]); if (m > mx) { mx = m; mi = j; } }
            k0v[i] = mi; k0m[i] = mx;
        }
        int best = 0; for (int i = 1; i < UP_USE; i++) if (k0m[i] > k0m[best]) best = i;
        int im = k0v[best];
        // Pass 2: recompute the FFTs and keep only the value at bin `im`.
        // (Recomputing is cheaper in RAM than storing every symbol's FFT, and
        //  only runs once per frame at sync.)
        Cx binval[UP_USE];
        for (int i = 0; i < UP_USE; i++) {
            for (int j = 0; j < N_; j++) dech_[j] = cxmul(s[N_ * i + j], downchirp_[j]);
            fftN(dech_, fout_, N_);
            binval[i] = fout_[im];
        }
        Cx acc = {0, 0};
        for (int i = 0; i < UP_USE - 1; i++)
            acc = (Cx){acc.r + (binval[i].r * binval[i + 1].r + binval[i].i * binval[i + 1].i),
                       acc.i + (binval[i].i * binval[i + 1].r - binval[i].r * binval[i + 1].i)};
        m_cfo_frac_ = -atan2f(acc.i, acc.r) / (2.0f * PI_F);
        for (int n = 0; n < UP_USE * N_; n++) {
            float a = -2.0f * PI_F * m_cfo_frac_ / N_ * n;
            u_.pre.up[n] = cxmul(s[n], (Cx){cosf(a), sinf(a)});
        }
    }

    void estimate_sto_frac() {
        int NN = 2 * N_;
        for (int j = 0; j < NN; j++) mag_[j] = 0;
        for (int i = 0; i < UP_USE; i++) {
            for (int j = 0; j < N_; j++) dech_[j] = cxmul(u_.pre.up[N_ * i + j], downchirp_[j]);
            for (int j = N_; j < NN; j++) dech_[j] = (Cx){0, 0};
            fftN(dech_, fout_, NN);
            for (int j = 0; j < NN; j++) mag_[j] += cxnorm(fout_[j]);
        }
        int k0 = 0; for (int j = 1; j < NN; j++) if (mag_[j] > mag_[k0]) k0 = j;
        float Y_1 = mag_[lmod(k0 - 1, NN)], Y0 = mag_[k0], Y1 = mag_[lmod(k0 + 1, NN)];
        float u = 64.0 * N_ / 406.5506497, v = u * 2.4674;
        float wa = (Y1 - Y_1) / (u * (Y1 + Y_1) + v * Y0);
        float ka = wa * N_ / LORA_PI;
        float kr = fmodf((k0 + ka) / 2.0, 1.0);
        m_sto_frac_ = (float)(kr - (kr > 0.5 ? 1.0 : 0.0));
    }

    long do_payload() {
        bool want = payload_sym_cnt_ < 8 ||
                    ((uint32_t)payload_sym_cnt_ < total_payload_symbols_ && received_head_);
        if (want) {
            demod_symbol(in_down_);
            long consume_n = sps_;
            // SFO compensation: nudge the window by a sample when drift exceeds
            // half a chip, keeping long payloads aligned.
            if (fabsf(sfo_cum_) > 1.0f / 2.0f / OS) {
                int sign = (sfo_cum_ >= 0) ? 1 : -1;
                consume_n -= sign;
                sfo_cum_ -= sign * 1.0f / OS;
            }
            sfo_cum_ += sfo_hat_;
            payload_sym_cnt_++;
            return consume_n;
        } else if (!received_head_) {
            return 0; // wait
        } else {
            reset();
            return sps_;
        }
    }

    void demod_symbol(const Cx* sym) {
        int block = 4 + (is_header_ ? 4 : pay_cr_);
        compute_llrs(sym, llr_block_[llr_block_n_],
                     nid_resid_ ? u_.pay.llr_block_rot[llr_block_n_] : nullptr);
        llr_block_n_++;
        if (llr_block_n_ == block) { decode_block(); llr_block_n_ = 0; }
    }

    void compute_llrs(const Cx* sym, float* out, float* out_rot) {
        for (int i = 0; i < N_; i++) dech_[i] = cxmul(sym[i], demod_down_[i]);
        fftN(dech_, fout_, N_);
        // Amplitude (not power) bin metric: closer to the noncoherent LLR and
        // avoids overweighting strong symbols inside a deinterleaved codeword.
        for (int i = 0; i < N_; i++) mag_[i] = sqrtf(cxnorm(fout_[i]));
        const uint16_t* map = (is_header_ || ldro_) ? demap_hdr_ : demap_pay_;
        for (int i = 0; i < sf_; i++) {
            float x1 = -1, x0 = -1, r1 = -1, r0 = -1;
            for (int n = 0; n < N_; n++) {
                bool one = map[n] & (1u << i);
                float v = mag_[n];
                if (one) { if (v > x1) x1 = v; }
                else { if (v > x0) x0 = v; }
                if (out_rot) {
                    float vr = mag_[(n + nid_resid_ + N_) % N_];
                    if (one) { if (vr > r1) r1 = vr; }
                    else { if (vr > r0) r0 = vr; }
                }
            }
            out[sf_ - 1 - i] = x1 - x0;
            if (out_rot) out_rot[sf_ - 1 - i] = r1 - r0;
        }
    }

    // soft Hamming over the LLR codeword; optionally reports the runner-up
    // data nibble and the best-vs-runner-up score margin (small = uncertain)
    static const uint8_t* cwlut();
    uint8_t hamming_soft(const float* llr, int cr_app,
                         uint8_t* alt = nullptr, float* margin = nullptr) {
        static const uint8_t LUT[16]     = {0,23,45,58,78,89,99,116,139,156,166,177,197,210,232,255};
        static const uint8_t LUT5[16]    = {0,24,40,48,72,80,96,120,136,144,160,184,192,216,232,240};
        int cw_len = cr_app + 4;
        float best = -1e30f, best2 = -1e30f; int bi = 0, bi2 = 1;
        for (int n = 0; n < 16; n++) {
            float p = 0;
            for (int j = 0; j < cw_len; j++) {
                bool bit = (((cr_app != 1) ? LUT[n] : LUT5[n]) >> (8 - cw_len)) & (1u << (cw_len - 1 - j));
                float a = llr[j] < 0 ? -llr[j] : llr[j];
                if ((bit && llr[j] > 0) || (!bit && llr[j] < 0)) p += a; else p -= a;
            }
            if (p > best) { best2 = best; bi2 = bi; best = p; bi = n; }
            else if (p > best2) { best2 = p; bi2 = n; }
        }
        auto to_nib = [](uint8_t d) -> uint8_t {
            return ((d & 1) << 3) | (((d >> 1) & 1) << 2) | (((d >> 2) & 1) << 1) | ((d >> 3) & 1);
        };
        if (alt) *alt = to_nib(LUT[bi2] >> 4);
        if (margin) *margin = best - best2;
        return to_nib(LUT[bi] >> 4);
    }

    void decode_block() {
        int sf_app = (is_header_ || ldro_) ? sf_ - 2 : sf_;
        int cw_len = is_header_ ? 8 : pay_cr_ + 4;
        int cr_app = is_header_ ? 4 : pay_cr_;
        for (int i = 0; i < cw_len; i++)
            for (int j = 0; j < sf_app; j++)
                u_.pay.deinter[lmod(i - j - 1, sf_app)][i] = llr_block_[i][sf_ - sf_app + j];
        if (nid_resid_)
            for (int i = 0; i < cw_len; i++)
                for (int j = 0; j < sf_app; j++)
                    u_.pay.deinter_rot[lmod(i - j - 1, sf_app)][i] = u_.pay.llr_block_rot[i][sf_ - sf_app + j];
        for (int i = 0; i < sf_app; i++) {
            if (nib_n_ < NIB_CAP) {
                nibbles_[nib_n_] = hamming_soft(u_.pay.deinter[i], cr_app,
                                                &u_.pay.nib_alt[nib_n_], &u_.pay.nib_margin[nib_n_]);
                u_.pay.nibbles_rot[nib_n_] = nid_resid_
                    ? hamming_soft(u_.pay.deinter_rot[i], cr_app) : nibbles_[nib_n_];
                nib_n_++;
            }
        }
        process_nibbles();
    }

    void process_nibbles() {
        if (is_header_ && nib_n_ >= 5) {
            pay_len_ = (nibbles_[0] << 4) | nibbles_[1];
            pay_crc_ = nibbles_[2] & 1;
            pay_cr_ = nibbles_[2] >> 1;
            if (pay_cr_ < 1) pay_cr_ = 1; if (pay_cr_ > 4) pay_cr_ = 4;
            uint8_t chk = ((nibbles_[3] & 1) << 4) | nibbles_[4];
            uint8_t n0 = nibbles_[0], n1 = nibbles_[1], n2 = nibbles_[2];
            auto b = [](uint8_t v, int k) { return (v >> k) & 1; };
            int c4 = b(n0,3)^b(n0,2)^b(n0,1)^b(n0,0);
            int c3 = b(n0,3)^b(n1,3)^b(n1,2)^b(n1,1)^b(n2,0);
            int c2 = b(n0,2)^b(n1,3)^b(n1,0)^b(n2,3)^b(n2,1);
            int c1 = b(n0,1)^b(n1,2)^b(n1,0)^b(n2,2)^b(n2,1)^b(n2,0);
            int c0 = b(n0,0)^b(n1,1)^b(n2,3)^b(n2,2)^b(n2,1)^b(n2,0);
            int comp = (c4<<4)|(c3<<3)|(c2<<2)|(c1<<1)|c0;
            if (comp != chk || pay_len_ == 0 || pay_len_ > LORA_LITE_MAX_PAYLOAD) { reset(); return; }
            // Low Data Rate Optimize engages when a symbol lasts > 16 ms
            ldro_ = ((float)(1u << sf_) * 1e3f / (float)bw_) > 16.0f;
            int sf2 = sf_ - 2 * (ldro_ ? 1 : 0);
            total_payload_symbols_ = 8 + ((int)ceilf((2.0f * pay_len_ - sf_ + 2 + 5 + (pay_crc_ ? 4 : 0)) / sf2)) * (4 + pay_cr_);
            received_head_ = true; is_header_ = false;
            // drop 5 header nibbles (all parallel streams stay aligned)
            for (int i = 5; i < nib_n_; i++) {
                nibbles_[i - 5] = nibbles_[i];
                u_.pay.nibbles_rot[i - 5] = u_.pay.nibbles_rot[i];
                u_.pay.nib_alt[i - 5] = u_.pay.nib_alt[i];
                u_.pay.nib_margin[i - 5] = u_.pay.nib_margin[i];
            }
            nib_n_ -= 5;
            // A valid header checksum + matching sync is a high-confidence real
            // packet. Emit a detection now so the caller sees it even if the
            // payload never fully assembles (weak signal fades before the end).
            emit_detection();
        }
        uint32_t need = pay_len_ * 2 + (pay_crc_ ? 4 : 0);
        if (!is_header_ && (uint32_t)nib_n_ >= need) emit();
    }

    void emit_detection() {
        Packet pkt{};
        pkt.len = 0;
        pkt.hdr_len = (uint8_t)pay_len_;
        pkt.status = PKT_HEADER_ONLY;
        pkt.has_crc = pay_crc_;
        pkt.crc_ok = false;
        pkt.snr = cur_snr_;
        pkt.cfo = m_cfo_int_ + m_cfo_frac_;
        pkt.sync = det_sync_;
        if (cb_) cb_(pkt, user_);
    }

    void build_bytes(const uint8_t* nib, uint8_t* out, int total) const {
        for (int i = 0; i < total; i++) {
            uint8_t lo = nib[2 * i], hi = nib[2 * i + 1];
            if ((uint32_t)i < pay_len_) { lo ^= WHITEN[i] & 0x0F; hi ^= (WHITEN[i] >> 4) & 0x0F; }
            out[i] = (hi << 4) | lo;
        }
    }

    bool crc_ok(const uint8_t* bytes) const {
        uint16_t calc = crc16_ccitt(bytes, pay_len_ - 2);
        calc ^= bytes[pay_len_ - 1] ^ ((uint16_t)bytes[pay_len_ - 2] << 8);
        uint16_t rx = bytes[pay_len_] | ((uint16_t)bytes[pay_len_ + 1] << 8);
        return calc == rx;
    }

    // Nibble chase: retry the CRC substituting the runner-up nibble at the
    // least-confident positions - singles and pairs of the 12 weakest plus
    // triples of the 8 weakest (~134 CRC trials, a few ms on the M4, only
    // ever runs on an already-lost frame). ALL patterns are tested; the
    // frame is adopted only when EXACTLY ONE pattern passes. With multiple
    // passers the 16-bit CRC cannot tell which is real, and a handheld UI
    // showing a green "decoded" packet must not be reporting a coin flip.
    bool nibble_chase(uint8_t* bytes, int total) {
        const int n_nib = 2 * total;
        constexpr int K = 12, K3 = 8;
        int k = n_nib < K ? n_nib : K;
        int k3 = k < K3 ? k : K3;
        int pos[K];
        for (int i = 0; i < k; i++) pos[i] = -1;
        for (int p = 0; p < n_nib; p++) {          // partial selection sort: k smallest margins
            for (int i = 0; i < k; i++) {
                if (pos[i] < 0 || u_.pay.nib_margin[p] < u_.pay.nib_margin[pos[i]]) {
                    for (int j = k - 1; j > i; j--) pos[j] = pos[j - 1];
                    pos[i] = p;
                    break;
                }
            }
        }
        uint8_t* trial = u_.pay.trial;
        memcpy(trial, nibbles_, n_nib);
        uint8_t* found = u_.pay.found;
        int n_passed = 0;
        uint8_t* tb = u_.pay.bytes;
        auto test = [&]() {
            build_bytes(trial, tb, total);
            if (crc_ok(tb)) {
                if (n_passed == 0) memcpy(found, tb, total);
                n_passed++;
            }
        };
        for (int i = 0; i < k && n_passed < 2; i++) {
            trial[pos[i]] = u_.pay.nib_alt[pos[i]];
            test();
            for (int j = i + 1; j < k && n_passed < 2; j++) {
                trial[pos[j]] = u_.pay.nib_alt[pos[j]];
                test();
                if (i < k3 && j < k3) {
                    for (int c = j + 1; c < k3 && n_passed < 2; c++) {
                        trial[pos[c]] = u_.pay.nib_alt[pos[c]];
                        test();
                        trial[pos[c]] = nibbles_[pos[c]];
                    }
                }
                trial[pos[j]] = nibbles_[pos[j]];
            }
            trial[pos[i]] = nibbles_[pos[i]];
        }
        if (n_passed == 1) {
            memcpy(bytes, found, total);
            return true;
        }
        return false;
    }

    void emit() {
        Packet pkt;
        int total = pay_len_ + (pay_crc_ ? 2 : 0);
        uint8_t bytes[LORA_LITE_MAX_PAYLOAD + 2];
        build_bytes(nibbles_, bytes, total);
        pkt.crc_ok = true;
        if (pay_crc_ && pay_len_ >= 2) {
            pkt.crc_ok = crc_ok(bytes);
            if (!pkt.crc_ok && nid_resid_) {
                // rotation hypothesis: shadow stream decoded under the
                // net-ID-residual bin shift
                build_bytes(u_.pay.nibbles_rot, u_.pay.rot_bytes, total);
                if (crc_ok(u_.pay.rot_bytes)) {
                    memcpy(bytes, u_.pay.rot_bytes, total);
                    pkt.crc_ok = true;
                }
            }
            if (!pkt.crc_ok && nibble_chase(bytes, total))
                pkt.crc_ok = true;
        }
        for (uint32_t i = 0; i < pay_len_; i++) pkt.data[i] = bytes[i];
        pkt.len = pay_len_; pkt.has_crc = pay_crc_; pkt.snr = cur_snr_;
        pkt.cfo = m_cfo_int_ + m_cfo_frac_; pkt.sync = det_sync_;
        pkt.hdr_len = (uint8_t)pay_len_;
        pkt.status = pkt.crc_ok ? PKT_DECODED : PKT_CRC_FAIL;
        if (cb_) cb_(pkt, user_);
        nib_n_ = 0;
        reset();
    }
};

} // namespace lora_lite
