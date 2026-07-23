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

## Verification status (important, read before flashing)

- ✅ **DSP correctness**: `lora_lite` decodes real off-air MeshCore captures (63 CRC-valid packets from `capture.iq`, matching the full `lora_rx` where alignment allows) and clean synthetic frames SF7–SF10. Run `make check` in the parent repo.
- ✅ **Compiles for the target**: both `lora_lite` and the baseband processor `proc_lorarx.cpp` compile cleanly with `arm-none-eabi-g++ 14.2` for `cortex-m4 -mfloat-abi=hard -mfpu=fpv4-sp-d16` against the real Mayhem headers.
- ⚠️ **Not yet run on hardware, and the full firmware image was not linked here.** Mayhem's build is pinned to its own Docker toolchain (arm-none-eabi 9.2.1); a from-scratch image build outside that environment hit *pre-existing* CMake issues that also affect an unmodified `master` checkout (empty-sources on `application.elf` / `sd_over_usb`), i.e. unrelated to this app. Build the image in the official Mayhem dev container (below), where those targets resolve normally.

So: treat this as a compile-verified, DSP-proven app that needs an on-device shakedown. The UI app follows the `afsk_rx` template closely; the risk is in on-hardware integration (decimation levels, RSSI/AGC), not the demodulator.

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
