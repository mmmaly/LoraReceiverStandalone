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
    if echo "$out" | grep -q "CRC valid" && echo "$out" | grep -qF "rx str: $expect"; then
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

# Noise robustness
run_case "SF8 noisy (sigma 0.15)" "-S 8 -b 62500 -s 250000 -c 1 -m 'NOISY TEST' -N 0.15" "-S 8 -b 62500 -s 250000" "NOISY TEST"

# Auto-scan: RX not told SF or BW, must find the right decoder
run_case "auto-scan finds SF9/125k" "-S 9 -b 125000 -s 250000 -c 2 -m 'AUTO TEST'" "-A -s 250000" "AUTO TEST"
run_case "auto-scan finds SF10/62.5k" "-S 10 -b 62500 -s 250000 -c 1 -m 'AUTO TEST 2'" "-A -s 250000" "AUTO TEST 2"

# Unknown sync word: -w 0 accepts and reports it
run_case "accept-any sync (-w 0)" "-S 8 -b 62500 -s 250000 -c 1 -w 34 -m 'SYNC TEST'" "-S 8 -b 62500 -s 250000 -w 0" "SYNC TEST"

exit $FAIL
