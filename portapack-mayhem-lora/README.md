# LoRa RX app for PortaPack / Mayhem (HackRF)

A native LoRa receiver app for the HackRF + PortaPack running [Mayhem firmware](https://github.com/portapack-mayhem/mayhem-firmware) — the long-standing community request ([#400](https://github.com/portapack-mayhem/mayhem-firmware/issues/400), [#2400](https://github.com/portapack-mayhem/mayhem-firmware/issues/2400)).

It runs the demodulator entirely on the device: the DSP lives in an **M4 baseband image** (`proc_lorarx`) built around `lora_lite` — a no-malloc, no-thread, static-buffer reduction of this repo's `lora_rx` chain — and an **M0 UI app** (`lorarx`) shows decoded packets with SNR on the display. No host PC needed.

Not a sensitivity contender (the HackRF front end and 8-bit ADC are why), but a genuinely useful field tool: sniff LoRa / MeshCore traffic, confirm a channel is active, read packet hex and SNR, all handheld.

## What's here

| File | Role |
|------|------|
| `0001-lora-rx-app.patch` | Complete patch against Mayhem `master` — adds every file and wiring below |
| `proc_lorarx.hpp/.cpp` | M4 baseband processor (÷8 decimation → `lora_lite::Demod`) |
| `app_ui/` | M0 external app: `main.cpp`, `ui_lorarx.hpp/.cpp` |

The `lora_lite.h` / `lora_kiss.h` core lives in the parent repo (and is copied into `firmware/baseband/` by the patch). It is verified on real captures by `run_tests.sh` (`lite` cases) here.

## Capabilities & limits

- **SF7, SF8** (compile-time `LORA_LITE_MAX_SF`, N≤256 — fits the 96 KiB M4 baseband SRAM at ~70 KB).
- **BW 62.5 kHz and 125 kHz** (÷8 decimation from 2.0 / 4.0 MHz, os = 4).
- Explicit-header mode: coding rate, length, CRC read from each packet's header.
- Sync word 0x12 (public), 0x34, or accept-any.
- Payloads up to 64 bytes (`LORA_LITE_MAX_PAYLOAD`).

Wider SF/BW is a matter of raising `LORA_LITE_MAX_SF` and adding decimation stages, bounded by baseband SRAM.

## Verification status

- ✅ **Verified on hardware** (2026-07-27, HackRF One + PortaPack, Mayhem v2.2.0 1 MB build): receives live MeshCore traffic off-air — yellow detections and full CRC-valid decodes at SF7/62.5 kHz with an external antenna. Follow `BUILD_AND_FLASH.md` for the exact working build recipe (pinned gcc 9.2.1, no Docker needed).
- ✅ **DSP correctness**: `lora_lite` decodes real off-air MeshCore captures (see the parent repo's regression suite; `make check` runs the `lite` cases) and clean synthetic frames SF7–SF10.
- ⚠️ **Sensitivity is modest — set expectations accordingly.** The HackRF's front end and 8-bit ADC give up real dB against both a dedicated LoRa chip and an RTL-SDR running the host `lora_rx` (the parent repo measured the RTL-SDR itself ~4–5 dB behind an SX1262; the HackRF sits below that). It reliably hears strong-to-moderate local traffic; weak frames show as yellow detections (valid header, failed payload CRC) or not at all. A good antenna matters more here than anywhere else in this project.

Practical settings that worked on air: correct SF first (networks differ — SF7 vs SF8 is the difference between traffic and silence), amp ON with LNA 32–40 / VGA ~40 on an external antenna. The RF path tolerates hot gain (no observed clipping in captures), but back off if the spectrum view saturates red.

## Build (official Mayhem dev container)

```bash
git clone --recurse-submodules https://github.com/portapack-mayhem/mayhem-firmware.git
cd mayhem-firmware
git apply /path/to/0001-lora-rx-app.patch     # or 'git am' to keep the commit
docker run --rm -v "$PWD":/havoc -w /havoc portapackmayhem/mayhem-firmware-dev \
    bash -c "mkdir -p build && cd build && cmake .. && make firmware -j$(nproc)"
```

The result `build/firmware/portapack-mayhem*.bin` (or the `.ppfw.tar` / `.tar` bundle) flashes the whole firmware. The LoRa app appears under the **RX** menu as "LoRa RX".

To integrate without the patch, copy `proc_lorarx.*` into `firmware/baseband/`, `app_ui/*` into `firmware/application/external/lorarx/`, `lora_lite.h`+`lora_kiss.h` into `firmware/baseband/`, then apply the small edits the patch makes to `message.hpp`, `spi_image.hpp`, `baseband_api.*`, `external.cmake`, `external.ld`, and `baseband/CMakeLists.txt`.

## On-device tuning to expect

If packets don't decode on first flash, the usual suspects (in order):
1. **Decimation level / scaling** — `proc_lorarx.cpp` normalizes the ÷8 decimator output by `1/32768`; the actual c16 magnitude may want a different scale for the FFT to peak cleanly.
2. **Sample-rate/bandwidth mapping** — `configure_baseband()` sets `fs = bw * 4 * 8`; confirm the radio accepts and the front-end filter matches.
3. **os** — the core is compiled `LORA_LITE_OS = 4`; the decimation must deliver exactly 4 samples per chip.

All three are observable on-device by adding a bin-peak/`feed_channel_stats` readout to the console.
