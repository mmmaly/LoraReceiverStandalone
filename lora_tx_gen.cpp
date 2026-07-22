// Offline LoRa frame synthesizer. Produces rtl_sdr-format u8 IQ files so
// lora_rx can be tested end-to-end without any hardware:
//
//   ./lora_tx_gen -o test.iq -S 8 -b 62500 -s 250000 -c 1 -m "HELLO"
//   ./lora_rx -r test.iq -S 8 -b 62500 -s 250000
//
// Encodes exactly the inverse of the lora_rx decode chain (explicit header
// with checksum, whitening, CRC-16, Hamming, diagonal interleaver, Gray).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <complex>
#include <getopt.h>

#include "lora_common.h"

static uint8_t reverse4(uint8_t d) {
    return (uint8_t)(((d & 1) << 3) | (((d >> 1) & 1) << 2) | (((d >> 2) & 1) << 1) | ((d >> 3) & 1));
}

int main(int argc, char *argv[]) {
    uint32_t samp_rate = 250000, bw = 62500;
    int sf = 8, cr = 1;
    uint16_t sync_word = 0x12;
    std::string payload = "HELLO LORA";
    const char *out_path = "test.iq";
    float noise = 0.02f, amp = 0.5f, cfo_hz = 0.0f;
    int id_offset = 1;      // chirp id = data value + this (matches rx convention)
    int pre_silence = 12345;

    int opt;
    while ((opt = getopt(argc, argv, "s:b:S:c:w:m:o:N:O:a:i:")) != -1) {
        switch (opt) {
        case 's': samp_rate = (uint32_t)atol(optarg); break;
        case 'b': bw = (uint32_t)atol(optarg); break;
        case 'S': sf = atoi(optarg); break;
        case 'c': cr = atoi(optarg); break;
        case 'w': sync_word = (uint16_t)strtol(optarg, nullptr, 16); break;
        case 'm': payload = optarg; break;
        case 'o': out_path = optarg; break;
        case 'N': noise = (float)atof(optarg); break;
        case 'O': cfo_hz = (float)atof(optarg); break;
        case 'a': amp = (float)atof(optarg); break;
        case 'i': id_offset = atoi(optarg); break;
        default:
            fprintf(stderr, "Usage: %s [-s samp_rate] [-b bw] [-S sf] [-c cr] [-w sync_hex]\n"
                            "          [-m message] [-o out.iq] [-N noise] [-O cfo_hz] [-a amp]\n", argv[0]);
            return 1;
        }
    }

    if (samp_rate % bw != 0) { fprintf(stderr, "samp_rate must be an integer multiple of bw\n"); return 1; }
    if (sf < 5 || sf > 12) { fprintf(stderr, "SF must be 5-12\n"); return 1; }
    if (cr < 1 || cr > 4) { fprintf(stderr, "CR must be 1-4\n"); return 1; }
    if (payload.size() < 2 || payload.size() > 255) { fprintf(stderr, "payload must be 2-255 bytes\n"); return 1; }

    int os = samp_rate / bw;
    uint32_t N = 1u << sf;
    uint32_t sps = N * os;
    bool ldro = (double)N * 1000.0 / bw > 16.0;
    uint32_t paylen = (uint32_t)payload.size();

    // ---- Nibble stream: explicit header (5) + whitened payload + CRC ----
    std::vector<uint8_t> nib;
    uint8_t n0 = (paylen >> 4) & 0xF, n1 = paylen & 0xF;
    uint8_t n2 = (uint8_t)((cr << 1) | 1); // CRC present
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

    for (uint32_t i = 0; i < paylen; i++) {
        uint8_t w = (uint8_t)payload[i] ^ whitening_seq[i];
        nib.push_back(w & 0xF);
        nib.push_back(w >> 4);
    }

    uint16_t crc = crc16((const uint8_t *)payload.data(), paylen - 2);
    crc = crc ^ (uint8_t)payload[paylen - 1] ^ ((uint16_t)(uint8_t)payload[paylen - 2] << 8);
    nib.push_back(crc & 0xF);
    nib.push_back((crc >> 4) & 0xF);
    nib.push_back((crc >> 8) & 0xF);
    nib.push_back((crc >> 12) & 0xF);

    // ---- Encode blocks -> chirp ids ----
    std::vector<uint32_t> ids;
    size_t pos = 0;
    bool first = true;
    while (pos < nib.size()) {
        int sf_app = (first || ldro) ? sf - 2 : sf;
        int cw_len = first ? 8 : cr + 4;
        int cr_app = first ? 4 : cr;

        // Hamming-encode sf_app nibbles into codeword rows (pad with 0)
        std::vector<uint8_t> rows(sf_app);
        for (int i = 0; i < sf_app; i++) {
            uint8_t d = pos < nib.size() ? nib[pos++] : 0;
            uint8_t rev = reverse4(d);
            rows[i] = (uint8_t)((cr_app != 1 ? cw_LUT[rev] : cw_LUT_cr5[rev]) >> (8 - cw_len));
        }

        // Diagonal interleave + inverse Gray -> chirp id per symbol
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
            ids.push_back((uint32_t)mod((long)val + id_offset, N));
        }
        first = false;
    }

    // ---- Synthesize samples ----
    std::vector<cx> sym(sps), sig;
    sig.reserve(pre_silence + (13 + ids.size()) * sps);
    sig.resize(pre_silence, cx(0, 0));

    auto up = [&](uint32_t id) {
        build_upchirp(sym.data(), id, (uint8_t)sf, (uint8_t)os);
        sig.insert(sig.end(), sym.begin(), sym.end());
    };

    for (int i = 0; i < 8; i++) up(0);                          // preamble
    up((uint32_t)(((sync_word & 0xF0) >> 4) << 3));             // net id 1
    up((uint32_t)((sync_word & 0x0F) << 3));                    // net id 2
    build_upchirp(sym.data(), 0, (uint8_t)sf, (uint8_t)os);     // 2.25 downchirps
    for (auto &z : sym) z = std::conj(z);
    sig.insert(sig.end(), sym.begin(), sym.end());
    sig.insert(sig.end(), sym.begin(), sym.end());
    sig.insert(sig.end(), sym.begin(), sym.begin() + sps / 4);
    for (uint32_t id : ids) up(id);                             // payload
    sig.resize(sig.size() + 4 * sps, cx(0, 0));                 // tail silence

    // ---- CFO, noise, quantize to u8, write ----
    FILE *f = fopen(out_path, "wb");
    if (!f) { perror("fopen"); return 1; }

    uint32_t rng = 0xC0FFEE42;
    auto frand = [&]() {
        rng = rng * 1664525u + 1013904223u;
        return (float)(rng >> 8) * (1.0f / 16777216.0f) - 0.5f;
    };

    std::vector<uint8_t> out(sig.size() * 2);
    for (size_t n = 0; n < sig.size(); n++) {
        cx z = sig[n] * amp;
        if (cfo_hz != 0)
            z *= expj(2.0f * (float)M_PI * cfo_hz / (float)samp_rate * (float)n);
        float re = z.real() + noise * (frand() + frand() + frand()) * 2.0f;
        float im = z.imag() + noise * (frand() + frand() + frand()) * 2.0f;
        int r8 = (int)lrintf(127.5f + 127.5f * re);
        int i8 = (int)lrintf(127.5f + 127.5f * im);
        out[2 * n]     = (uint8_t)(r8 < 0 ? 0 : (r8 > 255 ? 255 : r8));
        out[2 * n + 1] = (uint8_t)(i8 < 0 ? 0 : (i8 > 255 ? 255 : i8));
    }
    fwrite(out.data(), 1, out.size(), f);
    fclose(f);

    fprintf(stderr, "wrote %s: %zu IQ samples, %zu payload symbols (SF%d, BW %u, os %d%s, CR 4/%d)\n",
            out_path, sig.size(), ids.size(), sf, bw, os, ldro ? ", LDRO" : "", cr + 4);
    return 0;
}
