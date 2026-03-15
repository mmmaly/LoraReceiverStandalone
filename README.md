# LoraReceiverStandalone

Standalone LoRa packet receiver in pure C/C++. No GNU Radio dependency.

Connects directly to an RTL-SDR dongle and implements the full LoRa PHY demodulation chain:

- Preamble detection with CFO/STO estimation (Bernier + RCTSL algorithms)
- FFT-based chirp de-spreading with soft-decision LLR computation
- Gray demapping, deinterleaving, Hamming FEC decoding
- Explicit header parsing with checksum verification
- LFSR dewhitening and CRC-16 CCITT verification

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
make
```

## Usage

```bash
./lora_rx [options]
```

| Option | Description | Default |
|--------|-------------|---------|
| `-f <hz>` | Center frequency | 869618000 |
| `-s <hz>` | Sample rate | 250000 |
| `-b <hz>` | LoRa bandwidth | 62500 |
| `-S <5-12>` | Spreading factor | 8 |
| `-c <1-4>` | Coding rate (1=4/5 .. 4=4/8) | 4 |
| `-w <hex>` | Sync word | 0x12 |
| `-g <0.1dB>` | RTL-SDR tuner gain | 490 |
| `-p <ppm>` | Frequency correction | -3 |
| `-I` | Implicit header mode | off |
| `-L <bytes>` | Payload length (implicit mode) | 11 |

### Example

```bash
# EU 868 MHz, SF8, BW 62.5 kHz
./lora_rx -f 869618000 -b 62500 -S 8

# US 915 MHz, SF7, BW 125 kHz
./lora_rx -f 915000000 -b 125000 -S 7 -s 500000
```

### Output

```
[sync] Frame #1 detected, CFO=0.67, SNR=-14.8 dB
  Header: pay_len=111, CRC=1, CR=1
  Header checksum valid
rx msg: 0x11, 0x04, 0x67, ...
CRC valid!
```

## License

KissFFT: BSD-3-Clause (see kiss_fft.h).
DSP algorithms based on gr-lora_sdr (GPL-3.0).
