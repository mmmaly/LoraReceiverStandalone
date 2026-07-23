// Tiny self-contained radix-2 iterative FFT for lora_lite.
// No malloc, no recursion; twiddles precomputed once into a static table sized
// for the largest transform (2*MAX_N). Forward transform only, matching the
// sign convention of kiss_fft (which lora_rx used), so host and target agree.
#pragma once
#include <cstdint>
#include <cmath>

namespace lora_lite {
struct Cx;  // fwd (defined in lora_lite.h)

// LORA_PI is not guaranteed under strict ISO C++ (e.g. newlib on the M4 target).
static constexpr double LORA_PI = 3.14159265358979323846;

class Fft {
public:
    // max_log2 = log2 of the largest transform size that will be requested.
    void init(int max_log2) {
        maxlog_ = max_log2;
        int maxN = 1 << max_log2;
        for (int i = 0; i < maxN; i++) {
            double a = -2.0 * LORA_PI * i / maxN;   // full-size twiddles; subsampled per stage
            tw_[i][0] = (float)cos(a);
            tw_[i][1] = (float)sin(a);
        }
        maxN_ = maxN;
    }

    // out[k] = sum_n in[n] exp(-j 2pi n k / n)   (kiss_fft forward convention)
    void run(const Cx* in, Cx* out, int n) {
        const float* pin = reinterpret_cast<const float*>(in);
        float* pout = reinterpret_cast<float*>(out);
        int log2n = 0; while ((1 << log2n) < n) log2n++;

        // bit-reversal copy
        for (int i = 0; i < n; i++) {
            int j = bitrev(i, log2n);
            pout[2 * j]     = pin[2 * i];
            pout[2 * j + 1] = pin[2 * i + 1];
        }
        // Danielson-Lanczos stages
        for (int s = 1; s <= log2n; s++) {
            int m = 1 << s;          // butterfly span
            int half = m >> 1;
            int twstep = maxN_ / m;  // subsample the full twiddle table
            for (int k = 0; k < n; k += m) {
                for (int j = 0; j < half; j++) {
                    float wr = tw_[j * twstep][0];
                    float wi = tw_[j * twstep][1];
                    float* a = &pout[2 * (k + j)];
                    float* b = &pout[2 * (k + j + half)];
                    float tr = wr * b[0] - wi * b[1];
                    float ti = wr * b[1] + wi * b[0];
                    b[0] = a[0] - tr; b[1] = a[1] - ti;
                    a[0] = a[0] + tr; a[1] = a[1] + ti;
                }
            }
        }
    }

private:
    static int bitrev(int x, int bits) {
        int r = 0;
        for (int i = 0; i < bits; i++) { r = (r << 1) | (x & 1); x >>= 1; }
        return r;
    }
#ifndef LORA_LITE_MAX_SF
#define LORA_LITE_MAX_SF 10
#endif
    static constexpr int TWMAX = 1 << (LORA_LITE_MAX_SF + 1);
    float tw_[TWMAX][2];
    int maxlog_ = 0, maxN_ = 0;
};

} // namespace lora_lite
