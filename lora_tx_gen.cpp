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
#include "lora_frame.h"

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

    FrameParams fp;
    fp.sf = sf;
    fp.cr = cr;
    fp.sync_word = sync_word;
    fp.os = samp_rate / bw;
    fp.id_offset = id_offset;
    fp.pre_silence = pre_silence;
    fp.tail_silence = 4 * (1u << sf) * fp.os;

    FrameInfo fi;
    std::vector<cx> sig = build_lora_frame((const uint8_t *)payload.data(),
                                           payload.size(), fp, bw, &fi);


    // ---- CFO, noise, quantize to u8, write ----
    FILE *f = fopen(out_path, "wb");
    if (!f) { perror("fopen"); return 1; }

    uint32_t rng = 0xC0FFEE42;
    auto frand = [&]() {
        rng = rng * 1664525u + 1013904223u;
        return (float)(rng >> 8) * (1.0f / 16777216.0f) - 0.5f;
    };

    std::vector<uint8_t> out(sig.size() * 2);
    // Incremental oscillator with double phase: float phase accumulation would
    // lose precision after ~1e6 samples at large offsets
    double ph = 0.0;
    const double ph_inc = 2.0 * M_PI * (double)cfo_hz / (double)samp_rate;
    for (size_t n = 0; n < sig.size(); n++) {
        cx z = sig[n] * amp;
        if (cfo_hz != 0) {
            z *= cx((float)cos(ph), (float)sin(ph));
            ph += ph_inc;
            if (ph > M_PI) ph -= 2.0 * M_PI;
            else if (ph < -M_PI) ph += 2.0 * M_PI;
        }
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
            out_path, sig.size(), fi.payload_symbols, sf, bw, fp.os, fi.ldro ? ", LDRO" : "", cr + 4);
    return 0;
}
