// LoRa transmitter for the HackRF. The counterpart to lora_rx: same PHY
// parameters, same option letters, and it eats the very lines lora_rx prints,
// so a captured log can be replayed on air.
//
//   ./lora_tx -f 869618000 -S 7 -b 62500 -m "HELLO MESH"
//   ./lora_tx -f 869618000 -S 7 -x 2e0092293a8d
//   ./lora_tx -f 869618000 -S 7 -r captured.txt -n 3
//   grep "rx ok" obe-pakety.txt | ./lora_tx -f 869618000 -S 7 -r -
//
// Packet file format (one packet per line, "-" reads stdin):
//   2e0092293a8d...        hex payload (spaces and colons ignored)
//   rx ok: 2e009229...     a line straight out of lora_rx (prefix stripped)
//   @2.5                   wait 2.5 s before the next packet
//   # comment              ignored, as are blank lines
// With -t, each line is instead sent as literal text.
//
// Transmitting is regulated. The EU 868 MHz ISM band caps how long a device
// may occupy a channel, so this tool enforces a duty cycle (-y, default 1%,
// the conservative sub-band limit) by idling after each frame; it never
// leaves the PA keyed between frames.

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <csignal>
#include <string>
#include <vector>
#include <complex>
#include <getopt.h>
#include <unistd.h>

#include "lora_common.h"
#include "lora_frame.h"

#ifdef HAVE_HACKRF
#include <libhackrf/hackrf.h>
#endif

static volatile bool g_running = true;
static void on_signal(int) { g_running = false; }

// A packet to send, or a pause when `delay` is set (payload then empty)
struct TxItem {
    std::vector<uint8_t> payload;
    double delay = 0.0;
};

static bool parse_hex(const std::string &s, std::vector<uint8_t> &out) {
    int hi = -1;
    for (char c : s) {
        if (c == ' ' || c == ':' || c == ',' || c == '\t' || c == '_') continue;
        int v;
        if (c >= '0' && c <= '9') v = c - '0';
        else if (c >= 'a' && c <= 'f') v = c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v = c - 'A' + 10;
        else return false;
        if (hi < 0) hi = v;
        else { out.push_back((uint8_t)((hi << 4) | v)); hi = -1; }
    }
    return hi < 0;   // an odd number of nibbles is a malformed line
}

// Read the packet list. Understands lora_rx's own "rx ok:" output so a
// capture log can be piped straight back on air.
static bool load_packets(const char *path, bool as_text, std::vector<TxItem> &items) {
    FILE *f = strcmp(path, "-") == 0 ? stdin : fopen(path, "r");
    if (!f) { perror(path); return false; }
    char line[8192];
    int lineno = 0;
    bool ok = true;
    while (fgets(line, sizeof(line), f)) {
        lineno++;
        std::string s(line);
        while (!s.empty() && (s.back() == '\n' || s.back() == '\r' || s.back() == ' ')) s.pop_back();
        size_t b = s.find_first_not_of(" \t");
        if (b == std::string::npos) continue;
        s = s.substr(b);
        if (s[0] == '#') continue;
        if (s[0] == '@') {
            TxItem it;
            it.delay = atof(s.c_str() + 1);
            if (it.delay > 0) items.push_back(it);
            continue;
        }
        if (s.rfind("rx ok:", 0) == 0) s = s.substr(6);
        b = s.find_first_not_of(" \t");
        if (b == std::string::npos) continue;
        s = s.substr(b);

        TxItem it;
        if (as_text) {
            it.payload.assign(s.begin(), s.end());
        } else if (!parse_hex(s, it.payload)) {
            fprintf(stderr, "%s:%d: not valid hex (use -t for text input): %s\n", path, lineno, s.c_str());
            ok = false;
            continue;
        }
        if (it.payload.size() < 2 || it.payload.size() > 255) {
            fprintf(stderr, "%s:%d: payload is %zu bytes, must be 2-255\n", path, lineno, it.payload.size());
            ok = false;
            continue;
        }
        items.push_back(it);
    }
    if (f != stdin) fclose(f);
    return ok;
}

#ifdef HAVE_HACKRF
// libhackrf keeps TRANSFER_COUNT (4) buffers of TRANSFER_BUFFER_SIZE (256 KB)
// in flight; at 2 MS/s that is a quarter second of samples queued ahead of
// whatever the callback has just handed over.
static constexpr size_t HACKRF_QUEUE_BYTES = 4 * 262144;

struct TxState {
    const int8_t *buf = nullptr;
    size_t len = 0;
    size_t pos = 0;
};

static int tx_callback(hackrf_transfer *t) {
    TxState *st = (TxState *)t->tx_ctx;
    size_t want = (size_t)t->valid_length;
    size_t avail = st->len > st->pos ? st->len - st->pos : 0;
    size_t n = want < avail ? want : avail;
    if (n) memcpy(t->buffer, st->buf + st->pos, n);
    // Past the end of the frame the PA must fall silent, not repeat samples
    if (n < want) memset(t->buffer + n, 0, want - n);
    st->pos += n;
    return 0;
}
#endif

int main(int argc, char *argv[]) {
    uint32_t freq = 869618000, samp_rate = 2000000, bw = 62500;
    int sf = 8, cr = 1, ppm = 0;
    uint16_t sync_word = 0x12;
    long vga = 20;
    bool amp = false;
    double level = 0.7;
    double gap = 1.0, duty = 1.0;
    double offset_hz = 0.0;
    // -P raises the whole chain to maximum, but only for the knobs the user
    // did not set by hand, so an explicit -g/-L still wins whatever the order
    bool vga_set = false, level_set = false, offset_set = false, max_power = false;
    long repeats = 1;
    bool as_text = false, dry_run = false;
    const char *pkt_file = nullptr, *out_path = nullptr;
    std::vector<TxItem> items;

    int opt;
    while ((opt = getopt(argc, argv, "f:s:b:S:c:w:p:g:aL:m:x:r:tn:G:y:O:o:NP")) != -1) {
        switch (opt) {
        case 'f': freq = (uint32_t)atol(optarg); break;
        case 's': samp_rate = (uint32_t)atol(optarg); break;
        case 'b': bw = (uint32_t)atol(optarg); break;
        case 'S': sf = atoi(optarg); break;
        case 'c': cr = atoi(optarg); break;
        case 'w': sync_word = (uint16_t)strtol(optarg, nullptr, 16); break;
        case 'p': ppm = atoi(optarg); break;
        case 'g': vga = atol(optarg); vga_set = true; break;
        case 'a': amp = true; break;
        case 'L': level = atof(optarg); level_set = true; break;
        case 'P': max_power = true; break;
        case 'm': {
            TxItem it;
            it.payload.assign(optarg, optarg + strlen(optarg));
            items.push_back(it);
            break;
        }
        case 'x': {
            TxItem it;
            if (!parse_hex(optarg, it.payload)) { fprintf(stderr, "-x: not valid hex\n"); return 1; }
            items.push_back(it);
            break;
        }
        case 'r': pkt_file = optarg; break;
        case 't': as_text = true; break;
        case 'n': repeats = atol(optarg); break;
        case 'G': gap = atof(optarg); break;
        case 'y': duty = atof(optarg); break;
        case 'O': offset_hz = atof(optarg); offset_set = true; break;
        case 'o': out_path = optarg; break;
        case 'N': dry_run = true; break;
        default:
            fprintf(stderr,
                "Usage: %s [-f freq] [-s samp_rate] [-b bw] [-S sf] [-c cr] [-w sync_hex] [-p ppm]\n"
                "          [-m text | -x hex | -r pktfile] [-t] [-n repeats] [-G gap_s] [-y duty_pct]\n"
                "          [-g tx_vga_db] [-a] [-P] [-L level] [-O offset_hz] [-o file.iq] [-N]\n"
                "  PHY options match lora_rx: -f -s -b -S -c -w -p\n"
                "  -m <text>   send this text as one packet (repeatable)\n"
                "  -x <hex>    send this hex payload as one packet (repeatable)\n"
                "  -r <file>   packet list, \"-\" = stdin; one packet per line as hex or as a\n"
                "              \"rx ok: <hex>\" line copied from lora_rx, \"@<sec>\" to pause,\n"
                "              \"#\" comments. With -t each line is literal text instead.\n"
                "  -n <count>  repeat the whole list count times, 0 = until interrupted (default 1)\n"
                "  -G <sec>    gap between packets (default 1.0)\n"
                "  -y <pct>    duty cycle cap: idle after each frame so on-air time stays under\n"
                "              this share of the channel (default 1.0 = the conservative EU 868\n"
                "              sub-band limit; 869.4-869.65 MHz allows 10; 0 disables the cap)\n"
                "  -g <0-47>   HackRF TX VGA gain in dB (default 20)\n"
                "  -a          enable the HackRF RF amp (+11 dB) - off by default\n"
                "  -P          maximum power: -g 47 -a -L 1.0, for whichever of those you did\n"
                "              not set yourself. LoRa is constant-envelope, so driving the\n"
                "              chain into compression costs no signal quality - but it does\n"
                "              raise harmonics and spurs, which are separately regulated.\n"
                "  -L <0-1>    digital amplitude (default 0.7). A LoRa chirp has 0 dB PAPR, so\n"
                "              1.0 is the full DAC scale without clipping and is worth +3.1 dB\n"
                "  -O <hz>     transmit offset from the LO so its leakage lands outside the\n"
                "              channel (default 2x bw, min 100 kHz, capped at samp_rate/4;\n"
                "              0 tunes the LO onto the channel). Larger offsets cost signal\n"
                "              to the DAC sinc roll-off, so this stays no wider than it must\n"
                "  -o <file>   write baseband IQ (rtl_sdr u8) instead of transmitting; decode it\n"
                "              back with lora_rx -r <file>. No LO offset is applied.\n"
                "  -N          dry run: encode and report airtime, transmit nothing\n", argv[0]);
            return 1;
        }
    }

    if (pkt_file && !load_packets(pkt_file, as_text, items)) return 1;
    if (items.empty()) {
        fprintf(stderr, "Nothing to send: give -m, -x or -r (see -h)\n");
        return 1;
    }
    if (sf < 5 || sf > 12) { fprintf(stderr, "SF must be 5-12\n"); return 1; }
    if (cr < 1 || cr > 4) { fprintf(stderr, "CR must be 1-4\n"); return 1; }
    if (samp_rate % bw != 0) { fprintf(stderr, "samp_rate must be an integer multiple of bw\n"); return 1; }
    if (vga < 0 || vga > 47) { fprintf(stderr, "-g: HackRF TX VGA gain is 0-47 dB\n"); return 1; }
    if (level <= 0 || level > 1.0) { fprintf(stderr, "-L: level must be in (0,1]\n"); return 1; }
    if (duty < 0 || duty > 100) { fprintf(stderr, "-y: duty cycle is 0-100 percent\n"); return 1; }
    for (const TxItem &it : items) {
        if (it.delay == 0.0 && (it.payload.size() < 2 || it.payload.size() > 255)) {
            fprintf(stderr, "payload is %zu bytes, must be 2-255\n", it.payload.size());
            return 1;
        }
    }
    if (!out_path && !dry_run && samp_rate < 2000000) {
        fprintf(stderr, "HackRF minimum sample rate is 2000000 (use -o for a lower-rate file)\n");
        return 1;
    }
    if (max_power) {
        if (!vga_set) vga = 47;
        if (!level_set) level = 1.0;
        amp = true;
    }
    // The offset only has to clear the receiver's channel filter, and every Hz
    // of it costs signal: the DAC's zero-order hold rolls off as sinc(f/fs) and
    // the baseband filter adds its own skirt. The old samp_rate/4 default threw
    // away 0.9 dB at 2 MS/s to move the LO spur 8x further out than it needed
    // to be. Two channel widths is ample, capped at the old value so a wide bw
    // at a low sample rate is no worse off than before.
    if (!offset_set) {
        offset_hz = 2.0 * bw;
        if (offset_hz < 100000.0) offset_hz = 100000.0;
        if (offset_hz > samp_rate / 4.0) offset_hz = samp_rate / 4.0;
    }
    if (fabs(offset_hz) > 0.45 * samp_rate) {
        fprintf(stderr, "-O %.0f Hz is outside the usable span at samp_rate %u\n", offset_hz, samp_rate);
        return 1;
    }

    FrameParams fp;
    fp.sf = sf;
    fp.cr = cr;
    fp.sync_word = sync_word;
    fp.os = (int)(samp_rate / bw);
    // A little lead-in and tail keeps the PA ramp away from the preamble and
    // lets the final symbol clear the filter before the carrier drops
    fp.pre_silence = samp_rate / 100;
    fp.tail_silence = samp_rate / 100;

    // Encode every packet once; the same frame may be sent many times
    std::vector<std::vector<cx>> frames(items.size());
    double total_air = 0.0;
    size_t n_pkts = 0;
    for (size_t i = 0; i < items.size(); i++) {
        if (items[i].delay > 0.0) continue;
        FrameInfo fi;
        frames[i] = build_lora_frame(items[i].payload.data(), items[i].payload.size(), fp, bw, &fi);
        total_air += (double)frames[i].size() / samp_rate;
        n_pkts++;
        if (i == 0)
            fprintf(stderr, "SF%d BW %u os %d%s CR 4/%d sync 0x%02X\n",
                    sf, bw, fp.os, fi.ldro ? " LDRO" : "", cr + 4, sync_word);
    }

    double mean_air = total_air / n_pkts;
    fprintf(stderr, "%zu packet(s), %.3f s on air per pass (%.3f s per frame)\n",
            n_pkts, total_air, mean_air);
    if (duty > 0)
        fprintf(stderr, "duty cycle cap %.2f%% -> idling >= %.2f s after each frame\n",
                duty, mean_air * (100.0 / duty - 1.0));
    else
        fprintf(stderr, "duty cycle cap disabled (-y 0): only the %.2f s gap paces transmission\n", gap);

    if (!out_path) {
        double sinc_db = 0.0;
        if (offset_hz != 0.0) {
            double x = M_PI * offset_hz / samp_rate;
            sinc_db = 20.0 * log10(fabs(sin(x) / x));
        }
        fprintf(stderr, "VGA %ld dB, amp %s, level %.2f (%+.1f dBFS), LO offset %+.0f Hz (%.2f dB sinc loss)\n",
                vga, amp ? "ON" : "off", level, 20.0 * log10(level), offset_hz, sinc_db);
        if (max_power)
            fprintf(stderr, "-P: chain at maximum. LoRa is constant-envelope so compression costs no\n"
                            "    signal quality, but harmonics and spurious emissions are limited\n"
                            "    separately from output power - filter the output if you mean it.\n");
    }

    // ---- File output: plain baseband, decodable with lora_rx -r ----
    if (out_path) {
        FILE *f = fopen(out_path, "wb");
        if (!f) { perror(out_path); return 1; }
        size_t gap_samples = (size_t)(gap * samp_rate);
        std::vector<uint8_t> silence(2 * gap_samples, 127);
        for (size_t i = 0; i < items.size(); i++) {
            if (items[i].delay > 0.0) continue;
            std::vector<uint8_t> out(frames[i].size() * 2);
            for (size_t n = 0; n < frames[i].size(); n++) {
                cx z = frames[i][n] * (float)level;
                int r8 = (int)lrintf(127.5f + 127.5f * z.real());
                int i8 = (int)lrintf(127.5f + 127.5f * z.imag());
                out[2 * n]     = (uint8_t)(r8 < 0 ? 0 : (r8 > 255 ? 255 : r8));
                out[2 * n + 1] = (uint8_t)(i8 < 0 ? 0 : (i8 > 255 ? 255 : i8));
            }
            fwrite(out.data(), 1, out.size(), f);
            if (!silence.empty()) fwrite(silence.data(), 1, silence.size(), f);
        }
        fclose(f);
        fprintf(stderr, "wrote %s (decode: ./lora_rx -r %s -s %u -b %u -S %d)\n",
                out_path, out_path, samp_rate, bw, sf);
        return 0;
    }

    if (dry_run) {
        fprintf(stderr, "dry run: nothing transmitted\n");
        return 0;
    }

#ifndef HAVE_HACKRF
    fprintf(stderr, "This build has no HackRF support (libhackrf was not found at build time).\n"
                    "Use -o to write an IQ file instead.\n");
    return 1;
#else
    // ---- Quantize to the HackRF's signed 8-bit IQ, at the LO offset ----
    // The frame steps from silence to full amplitude in one sample. That edge
    // rings in the HackRF's analog reconstruction filter -- tolerable at -L 0.7,
    // but at full scale the overshoot clips, and either way the step splatters
    // energy into the neighbouring channels. A 20 us raised-cosine ramp on each
    // end removes both, at the cost of 1% of one of the eight preamble symbols.
    // The -o file output keeps its hard edges: the ringing is a property of the
    // radio, not of the waveform, and the tests decode that file byte-exact.
    const size_t ramp = samp_rate / 50000;

    std::vector<std::vector<int8_t>> tx(items.size());
    for (size_t i = 0; i < items.size(); i++) {
        if (items[i].delay > 0.0) continue;
        tx[i].resize(frames[i].size() * 2);
        const size_t r0 = fp.pre_silence;                        // first live sample
        const size_t r1 = frames[i].size() - fp.tail_silence;    // one past the last
        double ph = 0.0, ph_inc = 2.0 * M_PI * offset_hz / (double)samp_rate;
        for (size_t n = 0; n < frames[i].size(); n++) {
            cx z = frames[i][n] * (float)level;
            if (ramp && n >= r0 && n < r1) {
                size_t into = n - r0, from_end = r1 - 1 - n;
                size_t k = into < from_end ? into : from_end;
                if (k < ramp) z *= (float)(0.5 - 0.5 * cos(M_PI * (k + 0.5) / ramp));
            }
            if (offset_hz != 0.0) {
                z *= cx((float)cos(ph), (float)sin(ph));
                ph += ph_inc;
                if (ph > M_PI) ph -= 2.0 * M_PI;
            }
            int re = (int)lrintf(127.0f * z.real());
            int im = (int)lrintf(127.0f * z.imag());
            tx[i][2 * n]     = (int8_t)(re < -127 ? -127 : (re > 127 ? 127 : re));
            tx[i][2 * n + 1] = (int8_t)(im < -127 ? -127 : (im > 127 ? 127 : im));
        }
        frames[i].clear();
        frames[i].shrink_to_fit();
    }

    signal(SIGINT, on_signal);
    signal(SIGTERM, on_signal);

    hackrf_device *dev = nullptr;
    if (hackrf_init() != HACKRF_SUCCESS || hackrf_open(&dev) != HACKRF_SUCCESS) {
        fprintf(stderr, "Failed to open HackRF\n");
        return 1;
    }
    // libhackrf has no ppm knob, so -p (positive = crystal fast, as in rtlsdr)
    // is folded into the tune frequency, as it is on the receive side
    double lo = (double)freq - offset_hz;
    uint64_t tune = (uint64_t)llround(lo * (1.0 - ppm * 1e-6));
    hackrf_set_sample_rate(dev, samp_rate);
    hackrf_set_baseband_filter_bandwidth(dev, hackrf_compute_baseband_filter_bw_round_down_lt(samp_rate));
    hackrf_set_freq(dev, tune);
    hackrf_set_txvga_gain(dev, (uint32_t)vga);
    hackrf_set_amp_enable(dev, amp ? 1 : 0);

    fprintf(stderr, "HackRF TX on %.6f MHz (LO %.6f MHz)\n", freq / 1e6, tune / 1e6);

    long sent = 0;
    for (long pass = 0; g_running && (repeats == 0 || pass < repeats); pass++) {
        for (size_t i = 0; g_running && i < items.size(); i++) {
            if (items[i].delay > 0.0) {
                usleep((useconds_t)(items[i].delay * 1e6));
                continue;
            }
            TxState st;
            st.buf = tx[i].data();
            st.len = tx[i].size();
            st.pos = 0;
            if (hackrf_start_tx(dev, tx_callback, &st) != HACKRF_SUCCESS) {
                fprintf(stderr, "hackrf_start_tx failed\n");
                break;
            }
            // start_tx flips the transceiver mode, which reconfigures the RF
            // path. The gain and amp settings are meant to survive that; a
            // control transfer per frame is cheap enough not to depend on it.
            hackrf_set_txvga_gain(dev, (uint32_t)vga);
            hackrf_set_amp_enable(dev, amp ? 1 : 0);
            while (g_running && st.pos < st.len) usleep(2000);
            // The callback handing over the last sample only means libhackrf
            // has it, not that it has been radiated: a full USB transfer queue
            // still sits ahead of it. The callback feeds silence from here on,
            // so wait out that backlog before dropping the carrier. Stopping
            // early truncates the tail of the frame -- which is the payload,
            // leaving a receiver that locks on and decodes the header cleanly
            // and then sees nothing but noise where the payload should be.
            double drain = (double)(HACKRF_QUEUE_BYTES / 2) / samp_rate + 0.1;
            usleep((useconds_t)(drain * 1e6));
            hackrf_stop_tx(dev);
            sent++;

            double air = (double)(tx[i].size() / 2) / samp_rate;
            fprintf(stderr, "tx %ld: %zu bytes, %.3f s air\n", sent, items[i].payload.size(), air);

            double wait = gap;
            if (duty > 0) {
                double need = air * (100.0 / duty - 1.0);
                if (need > wait) wait = need;
            }
            bool last = (i + 1 == items.size()) && (repeats != 0 && pass + 1 == repeats);
            if (!last)
                for (double slept = 0; g_running && slept < wait; slept += 0.05) usleep(50000);
        }
    }

    hackrf_stop_tx(dev);
    hackrf_close(dev);
    hackrf_exit();
    fprintf(stderr, "\nDone: %ld frame(s) transmitted.\n", sent);
    return 0;
#endif
}
