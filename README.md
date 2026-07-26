# LoraReceiverStandalone

Standalone LoRa packet receiver in pure C/C++. No GNU Radio dependency.

Connects directly to an RTL-SDR dongle (or replays a recorded IQ file) and implements the full LoRa PHY demodulation chain:

- Anti-alias low-pass FIR before decimation (~5 dB sensitivity gain at 4x oversampling)
- Preamble detection with CFO/STO estimation (Bernier + RCTSL algorithms)
- FFT-based chirp de-spreading with soft-decision LLR computation
- Gray demapping, deinterleaving, Hamming FEC decoding
- Explicit header parsing with checksum verification
- LFSR dewhitening and CRC-16 CCITT verification
- Parallel auto-scan mode: all SF/BW combinations decoded simultaneously
- Multi-channel: one dongle receives several LoRa channels at once (`-C f1,f2,...`)

Based on the DSP algorithms from [gr-lora_sdr](https://github.com/tapparelj/gr-lora_sdr) (EPFL/TCL), reimplemented as a single self-contained application.

## Dependencies

- `librtlsdr` — RTL-SDR driver library
- A C++17 compiler (g++ or clang++)
- [KissFFT](https://github.com/mborgerding/kissfft) (included in source)

### Linux (Debian/Ubuntu)

```bash
sudo apt install librtlsdr-dev build-essential
```

### macOS

```bash
brew install librtlsdr
```

Xcode command line tools provide the compiler (`xcode-select --install` if needed). The Makefile auto-detects macOS and adds the correct Homebrew paths for both Apple Silicon and Intel Macs.

## Build

```bash
git clone https://github.com/mmmaly/LoraReceiverStandalone.git
cd LoraReceiverStandalone
make            # builds lora_rx
make check      # builds the offline frame generator and runs end-to-end tests (no SDR needed)
```

## Usage

```bash
./lora_rx [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-f <hz>` | Tuner center frequency | 869618000 |
| `-C <hz,hz,..>` | Channel frequencies (multi-channel mode; tuner auto-centers between them) | `-f` value |
| `-s <hz>` | Sample rate | 250000 (1000000 with multi `-C`) |
| `-b <hz,hz,..>` | LoRa bandwidth(s) | 62500 |
| `-S <5-12,..>` | Spreading factor(s) | 8 |
| `-c <1-4>` | Coding rate (1=4/5 .. 4=4/8); explicit header overrides | 4 |
| `-w <hex>` | Sync word; `0` = accept any | 0x12 |
| `-g <0.1dB>` | RTL-SDR tuner gain | 490 |
| `-p <ppm>` | Frequency correction | -3 |
| `-I` | Implicit header mode | off |
| `-L <bytes>` | Payload length (implicit mode) | 11 |
| `-A` | Auto-scan: run all SF (7-12) x BW decoders in parallel | off |
| `-W <lo>,<hi>[,<dwell>]` | Survey sweep: retune through the band in overlapping chunks (central 75% of the sample rate), running the full channel x SF x BW demod fan-out on each for `dwell` seconds (default 30), cycling forever. Two stages per chunk: a **coarse** pass on a raster of the smallest BW detects frames, then each detection's CFO gives the transmitter's true centre and a **refine** pass decodes there (real networks are rarely on a raster - MeshCore SK at 869.618 MHz is 29 bins off the 62.5 kHz grid). The tuner parks off-grid so no channel lands on the DC spike. Defaults: BW 125k, SF 7-12, 2 MS/s (auto-lowered so every BW ratio stays <= 16). `DETECT`/`HIT`/`REFINED HIT` summaries on stderr; normal `rx cfg`/`rx ok` stdout. | — |
| `-M <n>` | Max decoders running at once in `-W` mode (default: CPU cores). Bigger sets are split into sequential dwells at the same tuner position - oversubscribing the CPU makes the reader too slow and the dongle silently drops samples, which shreds frames mid-payload. A sweep warns when it keeps <95% of samples. | cores |
| `-r <file>` | Replay IQ from file (rtl_sdr u8 format, `-` = stdin) instead of SDR; a `.C16` file is read as PortaPack/Mayhem int16 IQ, with tuner freq and sample rate defaulted from its `.TXT` sidecar | — |
| `-D <file>` | In `-W` mode this becomes a prefix: one file per tuner position (`<file>-<tuner_hz>.iq`), each with a printed ready-to-run replay command - a single file spanning retunes has no well-defined centre frequency and cannot be replayed. Otherwise: Dump raw IQ to file while receiving (replay later with `-r`) | — |

### Examples

```bash
# EU 868 MHz, SF8, BW 62.5 kHz
./lora_rx -f 869618000 -b 62500 -S 8

# US 915 MHz, SF7, BW 125 kHz
./lora_rx -f 915000000 -b 125000 -S 7 -s 500000

# Unknown network: scan all SF/BW combos at once, accept any sync word.
# The [SFx/BWk] tag on the sync line tells you the network's parameters.
./lora_rx -f 869618000 -A -w 0

# Record while receiving, then replay later (no SDR needed for replay)
./lora_rx -D capture.iq
./lora_rx -r capture.iq -A

# One dongle, two meshes: Czech (869.432 MHz, SF7) and Slovak (869.618 MHz, SF8)
# received simultaneously. The tuner auto-centers between the channels at 1 MS/s,
# which also puts the RTL-SDR DC spike harmlessly between them.
./lora_rx -C 869432000,869618000 -S 7,8
```

`-S` and `-b` accept comma-separated lists; each channel gets a decoder for every SF x BW combination. In auto-scan mode they act as filters: `-A -S 9` scans only bandwidths at SF9; `-A -b 125000` scans only SFs at 125 kHz. `-A` combines with `-C` (all SF/BW on every channel).

When more than one decoder runs, each packet on stdout is preceded by a `rx cfg: freq=... sf=... bw=...` line identifying which channel/decoder produced it (the `rx msg` / `rx str` / CRC line format itself is unchanged).

CRC-valid packets additionally emit a machine-readable `rx ok: <hex>` line (plain hex, one line per packet), designed for piping into a packet decoder. For MeshCore, [meshcore-cpp-decoder](https://github.com/mmmaly/meshcore-cpp-decoder) consumes it directly with its `stream` command:

```bash
# Receive both meshes and decode/decrypt MeshCore packets in one pipeline
./lora_rx -C 869432000,869618000 -S 7,8 | meshcore-decoder stream -K ~/.config/meshcore/keys.txt
```

Non-packet lines pass through the decoder untouched, so channel attribution and CRC status stay visible in the combined output.

### Output

```
[SF8/62.5k] [sync] Frame #1 detected, CFO=0.67, sync=0x12, SNR=-14.8 dB
[SF8/62.5k]   Header: pay_len=111, CRC=1, CR=1
[SF8/62.5k]   Header checksum valid
rx msg: 0x11, 0x04, 0x67, ...
CRC valid!
```

Diagnostics (with the decoder tag) go to stderr; `rx msg` / `rx str` / CRC lines go to stdout, so `lora_rx | grep "rx msg"` pipelines keep working in every mode.

The `sync=0x..` field reports the sync word actually detected in the frame — with `-w 0` this identifies an unknown network's sync word.

## Offline testing

`lora_tx_gen` synthesizes a complete LoRa frame (header, whitening, Hamming, interleaving, CRC) as an rtl_sdr-format IQ file:

```bash
./lora_tx_gen -o test.iq -S 9 -b 125000 -s 250000 -c 2 -m "HELLO" -N 0.3 -O 300
./lora_rx -r test.iq -A          # decodes it, reports [SF9/125k]
```

Options: `-N <sigma>` additive noise, `-O <hz>` carrier frequency offset, `-a <amp>` amplitude, `-w <hex>` sync word. `./run_tests.sh` (or `make check`) runs an end-to-end suite covering SF7-SF12, LDRO, CFO, noise, auto-scan, and sync-word detection.

## License

KissFFT: BSD-3-Clause (see kiss_fft.h).
DSP algorithms based on gr-lora_sdr (GPL-3.0).
