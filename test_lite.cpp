// Host test: run lora_lite over an rtl_sdr u8 IQ file and print decoded packets.
// Verifies the embedded core against real captures before it goes to the M4.
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "lora_lite.h"
using namespace lora_lite;

static int g_count = 0, g_crcok = 0;
static void on_pkt(const Packet& p, void*) {
    g_count++;
    if (p.crc_ok) g_crcok++;
    printf("%s len=%2d snr=%+.1f cfo=%+.1f sync=0x%02x  ",
           p.crc_ok ? "OK " : "BAD", p.len, p.snr, p.cfo, p.sync);
    for (int i = 0; i < p.len; i++) printf("%02x", p.data[i]);
    printf("\n");
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s file.iq sf os_in_file\n", argv[0]); return 1; }
    int sf = atoi(argv[2]);
    int os_in = atoi(argv[3]);          // os of the file (250k/62.5k = 4)
    static Demod demod;
    demod.init(sf, 1, true, 0x12, on_pkt, nullptr);

    FILE* f = fopen(argv[1], "rb");
    if (!f) { perror("open"); return 1; }
    std::vector<unsigned char> buf(1 << 16);
    std::vector<Cx> samp;
    // The M4 will hand lora_lite os=LORA_LITE_OS. If the file already has that
    // os, feed directly; else this simple harness just requires os_in==OS.
    if (os_in != OS) { fprintf(stderr, "note: file os=%d, core os=%d; feeding as-is\n", os_in, OS); }
    size_t nr;
    while ((nr = fread(buf.data(), 1, buf.size(), f)) > 0) {
        int n = nr / 2;
        samp.resize(n);
        for (int i = 0; i < n; i++)
            samp[i] = Cx{((float)buf[2*i] - 127.5f) / 127.5f,
                         ((float)buf[2*i+1] - 127.5f) / 127.5f};
        demod.feed(samp.data(), n);
    }
    fclose(f);
    fprintf(stderr, "\nlora_lite: %d packets, %d CRC-valid\n", g_count, g_crcok);
    return 0;
}
