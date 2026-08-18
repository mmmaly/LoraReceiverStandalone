#!/bin/bash
# End-to-end offline tests: synthesize LoRa frames with lora_tx_gen, decode
# them with lora_rx via file replay, check payload and CRC. No SDR needed.
set -e
cd "$(dirname "$0")"

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
FAIL=0

run_case() {
    local desc="$1"; shift
    local genargs="$1"; shift
    local rxargs="$1"; shift
    local expect="$1"; shift

    eval "./lora_tx_gen -o '$TMP/t.iq' $genargs" 2>/dev/null
    local out
    out=$(eval "./lora_rx -r '$TMP/t.iq' $rxargs" 2>/dev/null || true)
    local expect_hex
    expect_hex=$(printf '%s' "$expect" | od -An -tx1 | tr -d ' \n')
    if echo "$out" | grep -q "CRC valid" && echo "$out" | grep -qF "rx str: $expect" \
       && echo "$out" | grep -qF "rx ok: $expect_hex"; then
        echo "PASS: $desc"
    else
        echo "FAIL: $desc"
        echo "$out" | sed 's/^/    /'
        FAIL=1
    fi
}

# Single-config decodes across SF / BW / CR (SF12@62.5k exercises LDRO)
run_case "SF7  BW62.5k CR4/5"  "-S 7  -b 62500  -s 250000 -c 1 -m 'AHOJ SF7'"   "-S 7  -b 62500  -s 250000" "AHOJ SF7"
run_case "SF8  BW62.5k CR4/8"  "-S 8  -b 62500  -s 250000 -c 4 -m 'AHOJ SF8'"   "-S 8  -b 62500  -s 250000" "AHOJ SF8"
run_case "SF9  BW125k  CR4/6"  "-S 9  -b 125000 -s 250000 -c 2 -m 'AHOJ SF9'"   "-S 9  -b 125000 -s 250000" "AHOJ SF9"
run_case "SF12 BW62.5k LDRO"   "-S 12 -b 62500  -s 250000 -c 1 -m 'AHOJ SF12'"  "-S 12 -b 62500  -s 250000" "AHOJ SF12"

# CFO robustness (300 Hz offset ~ a few bins at SF8/62.5k)
run_case "SF8 with 300Hz CFO"  "-S 8 -b 62500 -s 250000 -c 1 -m 'CFO TEST' -O 300" "-S 8 -b 62500 -s 250000" "CFO TEST"

# Large CFO: ~3.3 bins, integer part 3. Regression test for the sync check --
# real RTL-SDR dongles routinely sit several bins off even after ppm correction
run_case "SF8 with 800Hz CFO (3+ bins)" "-S 8 -b 62500 -s 250000 -c 1 -m 'BIG CFO TEST' -O 800" "-S 8 -b 62500 -s 250000" "BIG CFO TEST"

# Noise robustness
run_case "SF8 noisy (sigma 0.15)" "-S 8 -b 62500 -s 250000 -c 1 -m 'NOISY TEST' -N 0.15" "-S 8 -b 62500 -s 250000" "NOISY TEST"

# Auto-scan: RX not told SF or BW, must find the right decoder
run_case "auto-scan finds SF9/125k" "-S 9 -b 125000 -s 250000 -c 2 -m 'AUTO TEST'" "-A -s 250000" "AUTO TEST"
run_case "auto-scan finds SF10/62.5k" "-S 10 -b 62500 -s 250000 -c 1 -m 'AUTO TEST 2'" "-A -s 250000" "AUTO TEST 2"

# Unknown sync word: -w 0 accepts and reports it
run_case "accept-any sync (-w 0)" "-S 8 -b 62500 -s 250000 -c 1 -w 34 -m 'SYNC TEST'" "-S 8 -b 62500 -s 250000 -w 0" "SYNC TEST"

# Dual-channel: one 1 MS/s stream carrying two meshes 186 kHz apart
# (Czech 869.432 MHz SF7 + Slovak 869.618 MHz SF8), transmitted
# simultaneously. Tuner auto-centers at the midpoint; offsets are +/-93 kHz.
./lora_tx_gen -o "$TMP/cz.iq" -S 7 -b 62500 -s 1000000 -c 1 -m 'CESKA SIT' -O -93000 2>/dev/null
./lora_tx_gen -o "$TMP/sk.iq" -S 8 -b 62500 -s 1000000 -c 1 -m 'SLOVENSKA SIET' -O 93000 2>/dev/null
python3 - "$TMP/cz.iq" "$TMP/sk.iq" "$TMP/mix.iq" <<'PY'
import sys
a = open(sys.argv[1], 'rb').read()
b = open(sys.argv[2], 'rb').read()
n = max(len(a), len(b))
a += bytes([128]) * (n - len(a))
b += bytes([128]) * (n - len(b))
open(sys.argv[3], 'wb').write(bytes(min(255, max(0, x + y - 128)) for x, y in zip(a, b)))
PY
out=$(./lora_rx -r "$TMP/mix.iq" -s 1000000 -C 869432000,869618000 -S 7,8 -b 62500 2>/dev/null || true)
if echo "$out" | grep -qF "rx str: CESKA SIT" && echo "$out" | grep -qF "rx str: SLOVENSKA SIET"; then
    echo "PASS: dual-channel (one dongle, two meshes, simultaneous frames)"
else
    echo "FAIL: dual-channel (one dongle, two meshes, simultaneous frames)"
    echo "$out" | sed 's/^/    /'
    FAIL=1
fi

# ---- lora_tx (transmitter) ----
# Same encoder as lora_tx_gen, driven through the transmitter's own input
# parsing: a hex line, a "rx ok:" line copied from lora_rx, and a text line.
# -o writes the baseband it would have sent, so the loopback needs no radio.
tx_case() {
    local desc="$1" txargs="$2" rxargs="$3" expect_hex="$4"
    eval "./lora_tx -o '$TMP/tx.iq' -G 0.02 $txargs" >/dev/null 2>&1
    local out
    out=$(eval "./lora_rx -r '$TMP/tx.iq' $rxargs" 2>/dev/null || true)
    if echo "$out" | grep -qF "rx ok: $expect_hex"; then
        echo "PASS: lora_tx $desc"
    else
        echo "FAIL: lora_tx $desc"
        echo "$out" | sed 's/^/    /'
        FAIL=1
    fi
}

tx_case "text payload (-m)" "-S 7 -b 62500 -s 250000 -m 'TX HELLO'" \
        "-S 7 -b 62500 -s 250000" "$(printf 'TX HELLO' | od -An -tx1 | tr -d ' \n')"
tx_case "hex payload (-x)" "-S 8 -b 62500 -s 250000 -c 2 -x 2e0092293a8dc83d" \
        "-S 8 -b 62500 -s 250000" "2e0092293a8dc83d"

# The transmitter must accept lora_rx's own log lines, so a capture replays
printf '# comment\nrx ok: aabbccddeeff\n@0.1\n0102030405\n' > "$TMP/pkts.txt"
out=$(./lora_tx -o "$TMP/tx2.iq" -G 0.02 -S 7 -b 62500 -s 250000 -r "$TMP/pkts.txt" >/dev/null 2>&1 && \
      ./lora_rx -r "$TMP/tx2.iq" -S 7 -b 62500 -s 250000 2>/dev/null || true)
if echo "$out" | grep -qF "rx ok: aabbccddeeff" && echo "$out" | grep -qF "rx ok: 0102030405"; then
    echo "PASS: lora_tx packet file (rx-ok replay + hex + comments)"
else
    echo "FAIL: lora_tx packet file (rx-ok replay + hex + comments)"
    echo "$out" | sed 's/^/    /'
    FAIL=1
fi

# ---- lora_lite (embedded core, used by the PortaPack/Mayhem app) ----
# Same synthesized frames, decoded by the no-malloc/no-thread static-buffer core.
if command -v c++ >/dev/null 2>&1; then
    lite_case() {
        local desc="$1" sf="$2" os="$3" bin="$TMP/lite_$2_$3"
        c++ -O2 -std=c++17 -DLORA_LITE_OS=$os -DLORA_LITE_MAX_SF=8 -I. -o "$bin" test_lite.cpp 2>/dev/null
        ./lora_tx_gen -o "$TMP/l.iq" -S $sf -b 62500 -s $((62500*os)) -c 1 -m "LITE SF$sf" -N 0.02 2>/dev/null
        if "$bin" "$TMP/l.iq" $sf $os 2>/dev/null | grep -q "^OK.*$(printf 'LITE SF%s' $sf | od -An -tx1 | tr -d ' \n')"; then
            echo "PASS: lora_lite $desc"
        else
            echo "FAIL: lora_lite $desc"; FAIL=1
        fi
    }
    lite_case "SF7 os2 (M4 default)" 7 2
    lite_case "SF8 os2 (M4 default)" 8 2
    lite_case "SF8 os4 (firmware config)" 8 4
else
    echo "SKIP: lora_lite (no host c++ compiler)"
fi

exit $FAIL
