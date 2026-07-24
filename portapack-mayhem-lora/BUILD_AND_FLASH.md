# Building & flashing the LoRa RX firmware (what actually worked)

This records the exact, reproducible recipe used to build a flashable Mayhem
firmware image with the LoRa RX app, on macOS (Apple Silicon), **without Docker**.

## IMPORTANT: match your device's flash size

A **stock HackRF One has 1 MB of SPI flash**. Mayhem `master`/v2.4.0 is a **2 MB**
build and will NOT flash on it — the on-device Flash Utility rejects any file
larger than the flash ("BAD FIRMWARE FILE OR W/R ERR"). The target device here
runs **v2.2.0 on 1 MB flash**, so the deliverable is a **1 MB v2.2.0 image**:

- `portapack-mayhem_v2.2.0_LORA.bin` (exactly 1 MB, checksum-valid) — flash this.
- Patch: `0001-lora-rx-app-v2.2.0.patch` (against the `v2.2.0` tag).

If your device has 2 MB+ flash, the master build (`0001-lora-rx-app.patch`,
`-DFLASH_MB_SIZE=2 -DFLASH_MB_LIMIT_SIZE=2`) also works and keeps all stock apps.

### The gotcha that made the 1 MB build fit

Do **NOT** pass `-DCMAKE_BUILD_TYPE=Release`. Mayhem's app is compiled `-Os`
(`USE_OPT`), but `Release` appends `-O3 -DNDEBUG` which *overrides* `-Os` and
bloats the M0 app by ~66 KB — enough to overflow 1 MB. Configure with a bare
`cmake ..` (no build type).

## Toolchain (all no-sudo, downloaded to a scratch dir)

| Tool | Version | Why |
|------|---------|-----|
| arm-none-eabi-gcc | **9.2.1** (9-2019-q4-major) | Mayhem's *pinned* compiler. Newer GCC (14.2) bloats code ~10-15% and overflows flash. Use exactly this. |
| CMake | 3.28.3 | 4.x is fine too; 3.x used here. |
| Python + PyYAML | any 3.x | libopencm3 header generator. |
| dfu-util (`dfu-suffix`) | any | packages the HackRF USB DFU (part of `make firmware`). |
| `readelf` | from the arm toolchain | the external-app packager calls plain `readelf`; symlink `arm-none-eabi-readelf` → `readelf` on PATH. |

## Recipe

```bash
git clone --recurse-submodules https://github.com/portapack-mayhem/mayhem-firmware.git
cd mayhem-firmware
git submodule update --init --recursive     # CRUCIAL: nested hackrf/libopencm3, or you get
                                            # confusing "No SOURCES given to target" errors
git apply /path/to/0001-lora-rx-app.patch

# toolchain + shims on PATH (9.2.1 gcc, python with yaml, a readelf shim)
export PATH="/path/to/gcc-arm-none-eabi-9-2019-q4-major/bin:/path/to/python-with-yaml:$PATH"
ln -sf "$(command -v arm-none-eabi-readelf)" /tmp/shims/readelf && export PATH="/tmp/shims:$PATH"

mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DFLASH_MB_SIZE=2 -DFLASH_MB_LIMIT_SIZE=2   # 2 MB = official release config
make firmware -j$(sysctl -n hw.ncpu 2>/dev/null || nproc)
# -> build/firmware/portapack-mayhem-firmware.bin
# -> build/firmware/application/lorarx.ppma
```

### Gotchas that cost real time (so you don't repeat them)

1. **Nested submodules**: `hackrf` has its own `libopencm3` submodule. Without a
   *recursive* init, CMake generate fails with "No SOURCES given to target:
   application.elf / sd_over_usb" — misleading; it's really missing source files.
2. **Compiler version**: with GCC 14.2 the *core* M0 app overflowed flash by ~100 KB.
   The pinned 9.2.1 is not optional for a flashable image.
3. **Flash limit**: default is `FLASH_MB_LIMIT_SIZE=1`; the app is ~1.15 MB, so you
   must pass `=2` (this is exactly what the official nightly/stable CI does).
4. **Oversized stock apps on this master commit**: 14 stock external apps
   (jammer, sstvrx, random_password, remote, fmradio, scanner, waterfall_designer,
   battleship, epirb_rx, flex_rx, subcarrx, morse_radio, flex_tx, two_tone_rx) plus
   epirb_tx exceed the 32 KB external-app slot on this HEAD (e.g. battleship is only
   192 B over) and were **disabled** to let the image link. They are unrelated to
   LoRa. If you build on a tagged stable release instead of bleeding master, they
   should fit and you can re-enable them.

## LoRa-specific size work

The M4 baseband image must fit a 32 KB slot **combined with the M0 UI** (the .ppma
bundles both). Two changes got it there:
- `-Os` on `proc_lorarx.cpp` (baseband is `-O3` globally).
- **All `double` → `float`** in `lora_lite`: the Cortex-M4F has only a
  single-precision FPU, so every `double` op was software-emulated (kilobytes of
  `__aeabi_d*`). Final combined app: ~32.3 KB. Host decode is unchanged (still 63
  CRC-valid on `capture.iq`, 14/14 tests).

## Flashing — two options

### Option A: SD-card app only (no firmware flash, lowest risk) — try this first

`lorarx.ppma` is self-contained (bundles its baseband). On a PortaPack already
running Mayhem **v2.4.0** (or a close nightly):

1. Copy `lorarx.ppma` to the SD card under `/APPS/`.
2. On the device: main menu → **App Loader** (the "..." / external-apps menu) →
   **LoRa RX**.

If it refuses to load (firmware ABI mismatch), use Option B.

### Option B: flash the full firmware (guaranteed) 

Flash `portapack-mayhem-firmware.bin` by any normal Mayhem method:

- **On-device (easiest)**: copy the `.bin` to the SD card, then on the PortaPack go
  **Utilities → Flash Utility**, select the file, confirm. The device reflashes and
  reboots.
- **Web flasher**: https://flasher.mayhem.mnmlabs.dev/ (or the official one linked
  from the Mayhem wiki) — connect over USB, pick the `.bin`.
- **USB DFU / CLI**:
  ```bash
  hackrf_spiflash -w portapack-mayhem-firmware.bin      # HackRF tools
  ```
  (hold the HackRF in DFU if required; see the Mayhem "Flashing" wiki page.)

⚠️ Back up your current firmware first if you can. This image is built from bleeding
master with 15 stock apps omitted (see gotcha 4); it is meant for trying the LoRa
app, not as a daily driver. Once the LoRa app is validated on hardware, a clean
build on a stable tag (with all apps) is the right production step.

## Using the app

Main menu → RX → **LoRa RX** (or via App Loader if using the .ppma). Set frequency,
SF (7/8), bandwidth (62.5/125 kHz), sync word (0x12 public). Decoded packets scroll
with hex payload + SNR. First on-air test: point it at your Slovak MeshCore channel
(869.618 MHz, SF8, BW 62.5 kHz, sync 0x12) — expect the same packets your SDR sees,
just fewer (HackRF sensitivity).
