// LoRa frame synthesis: payload bytes -> complex baseband samples.
//
// The exact inverse of the lora_rx decode chain (explicit header with
// checksum, whitening, CRC-16, Hamming, diagonal interleaver, Gray), shared
// by lora_tx_gen (writes IQ files) and lora_tx (transmits over a HackRF) so
// there is exactly one encoder to keep in step with the receiver.
#pragma once

#include <cstdint>
#include <cstddef>
#include <complex>
#include <vector>

#include "lora_common.h"

struct FrameParams {
    int sf = 8;
    int cr = 1;                 // 1..4 => 4/5..4/8
    uint16_t sync_word = 0x12;
    int os = 4;                 // samp_rate / bw
    int id_offset = 1;          // chirp id = data value + this (rx convention)
    size_t pre_silence = 0;     // samples of silence before the preamble
    size_t tail_silence = 0;    // samples of silence after the last symbol
};

struct FrameInfo {
    size_t payload_symbols = 0;
    bool ldro = false;
};

inline uint8_t lora_reverse4(uint8_t d) {
    return (uint8_t)(((d & 1) << 3) | (((d >> 1) & 1) << 2) | (((d >> 2) & 1) << 1) | ((d >> 3) & 1));
}

// Encode payload into chirp ids (header + whitened payload + CRC).
// Payload must be 2..255 bytes: the CRC convention below reads the last two
// bytes directly, which is what the receiver checks.
inline std::vector<uint32_t> lora_encode_ids(const uint8_t *payload, size_t paylen,
                                             const FrameParams &p, bool ldro) {
    uint32_t N = 1u << p.sf;

    std::vector<uint8_t> nib;
    uint8_t n0 = (uint8_t)((paylen >> 4) & 0xF), n1 = (uint8_t)(paylen & 0xF);
    uint8_t n2 = (uint8_t)((p.cr << 1) | 1); // CRC present
    auto b = [](uint8_t v, int k) { return (v >> k) & 1; };
    int c4 = b(n0, 3) ^ b(n0, 2) ^ b(n0, 1) ^ b(n0, 0);
    int c3 = b(n0, 3) ^ b(n1, 3) ^ b(n1, 2) ^ b(n1, 1) ^ b(n2, 0);
    int c2 = b(n0, 2) ^ b(n1, 3) ^ b(n1, 0) ^ b(n2, 3) ^ b(n2, 1);
    int c1 = b(n0, 1) ^ b(n1, 2) ^ b(n1, 0) ^ b(n2, 2) ^ b(n2, 1) ^ b(n2, 0);
    int c0 = b(n0, 0) ^ b(n1, 1) ^ b(n2, 3) ^ b(n2, 2) ^ b(n2, 1) ^ b(n2, 0);
    nib.push_back(n0);
    nib.push_back(n1);
    nib.push_back(n2);
    nib.push_back((uint8_t)c4);
    nib.push_back((uint8_t)((c3 << 3) | (c2 << 2) | (c1 << 1) | c0));

    for (size_t i = 0; i < paylen; i++) {
        uint8_t w = (uint8_t)(payload[i] ^ whitening_seq[i]);
        nib.push_back(w & 0xF);
        nib.push_back(w >> 4);
    }

    uint16_t crc = crc16(payload, (uint32_t)(paylen - 2));
    crc = (uint16_t)(crc ^ payload[paylen - 1] ^ ((uint16_t)payload[paylen - 2] << 8));
    nib.push_back(crc & 0xF);
    nib.push_back((crc >> 4) & 0xF);
    nib.push_back((crc >> 8) & 0xF);
    nib.push_back((crc >> 12) & 0xF);

    std::vector<uint32_t> ids;
    size_t pos = 0;
    bool first = true;
    while (pos < nib.size()) {
        int sf_app = (first || ldro) ? p.sf - 2 : p.sf;
        int cw_len = first ? 8 : p.cr + 4;
        int cr_app = first ? 4 : p.cr;

        std::vector<uint8_t> rows(sf_app);
        for (int i = 0; i < sf_app; i++) {
            uint8_t d = pos < nib.size() ? nib[pos++] : 0;
            uint8_t rev = lora_reverse4(d);
            rows[i] = (uint8_t)((cr_app != 1 ? cw_LUT[rev] : cw_LUT_cr5[rev]) >> (8 - cw_len));
        }

        for (int i = 0; i < cw_len; i++) {
            uint32_t s = 0;
            for (int j = 0; j < sf_app; j++) {
                int r = (int)mod(i - j - 1, sf_app);
                int bit = (rows[r] >> (cw_len - 1 - i)) & 1;
                s |= (uint32_t)bit << (sf_app - 1 - j);
            }
            uint32_t v = s;
            v ^= v >> 1; v ^= v >> 2; v ^= v >> 4; v ^= v >> 8;
            uint32_t val = (first || ldro) ? v * 4 : v;
            ids.push_back((uint32_t)mod((long)val + p.id_offset, N));
        }
        first = false;
    }
    return ids;
}

// Full frame as unit-amplitude complex baseband: silence, 8 preamble upchirps,
// 2 sync-word symbols, 2.25 downchirps, payload symbols, silence.
inline std::vector<cx> build_lora_frame(const uint8_t *payload, size_t paylen,
                                        const FrameParams &p, uint32_t bw,
                                        FrameInfo *info = nullptr) {
    uint32_t N = 1u << p.sf;
    uint32_t sps = N * (uint32_t)p.os;
    bool ldro = (double)N * 1000.0 / bw > LORA_LDRO_MAX_MS;

    std::vector<uint32_t> ids = lora_encode_ids(payload, paylen, p, ldro);
    if (info) { info->payload_symbols = ids.size(); info->ldro = ldro; }

    std::vector<cx> sym(sps), sig;
    sig.reserve(p.pre_silence + (13 + ids.size()) * sps + p.tail_silence);
    sig.resize(p.pre_silence, cx(0, 0));

    auto up = [&](uint32_t id) {
        build_upchirp(sym.data(), id, (uint8_t)p.sf, (uint8_t)p.os);
        sig.insert(sig.end(), sym.begin(), sym.end());
    };

    for (int i = 0; i < 8; i++) up(0);                              // preamble
    up((uint32_t)(((p.sync_word & 0xF0) >> 4) << 3));               // net id 1
    up((uint32_t)((p.sync_word & 0x0F) << 3));                      // net id 2
    build_upchirp(sym.data(), 0, (uint8_t)p.sf, (uint8_t)p.os);     // 2.25 downchirps
    for (auto &z : sym) z = std::conj(z);
    sig.insert(sig.end(), sym.begin(), sym.end());
    sig.insert(sig.end(), sym.begin(), sym.end());
    sig.insert(sig.end(), sym.begin(), sym.begin() + sps / 4);
    for (uint32_t id : ids) up(id);                                 // payload
    sig.resize(sig.size() + p.tail_silence, cx(0, 0));
    return sig;
}
