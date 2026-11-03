# Core FFT Algorithms

Ten 1-D kernels, all compiled into `libfft_lib_c.a`. Each source file also keeps
a demonstration `main()` behind `#ifndef LIB_BUILD`; those demos are not
self-contained (they call each other), so the CMake build does not compile them.
Use `test_all` or `benchmark_all` instead.

| File | Lengths | Scratch | Notes |
|------|---------|---------|-------|
| `radix2_dit.c` | power of 2 | none | In place; bit-reverses the input before the butterflies. The baseline for every ratio on the performance page, and the planner's choice up to 256 points. |
| `radix2_dif.c` | power of 2 | none | Bit-reversal at the end instead of the start. Same cost and natural-order output as DIT; choose by data flow, not speed. |
| `radix4.c` | any power of 2 | O(n) digit reversal | A leading radix-2 stage covers an odd `log2(n)`. Saves 25% of the multiplications; the planner picks it from 257 to 4095 points, or when AVX2/OpenMP are unavailable. `FFT_CONSERVE_MEMORY` selects DIT instead. |
| `radix4_avx2.c` | any power of 2 | O(n) reorder + per-stage twiddle table | AVX2 version of the radix-4 stage, two complex values per register. Runtime dispatch; selected from 4096 points when the CPU and OS support AVX2. Falls back to scalar below 4096. |
| `radix4_omp.c` | any power of 2 | O(n) reorder + per-stage twiddle table | OpenMP version, built only with `-DFFT_ENABLE_OPENMP=ON`. Selected from 524288 points when more than one thread is available; scalar fallback below that or with one thread. |
| `four_step.c` | any power of 2 | n + about 2*sqrt(n) | Cache-blocked N1*N2 decomposition with radix-4 sub-transforms and a tiled transpose. Correctness-tested but not selected by the size heuristic: it measured slower than AVX2/OpenMP at 1M-4M. `FFT_MEASURE` can still choose it. |
| `split_radix.c` | power of 2 | 4n pool, peak about 3n | Lowest multiplication count. The recursion reads strided input instead of gathering at each level; still slower than radix-4 because of the memory traffic per stage. |
| `bluestein.c` | any length up to 2^29 | about 9n | Chirp z-transform; pads the convolution to `next_power_of_two(2n-1)`. The default for non-power-of-two lengths, and the only kernel that handles primes at O(n log n). |
| `mixed_radix.c` | any composite | O(n) | Recursive Cooley-Tukey factorisation with natural-order output. Prime lengths fall back to a direct O(p^2) DFT, which is why the planner sends primes to Bluestein. `FFT_CONSERVE_MEMORY` picks it for highly composite lengths. |
| `recursive_fft.c` | power of 2 | O(n) peak, about 3n | Divide-and-conquer reference implementation, kept for clarity. Not in the benchmark table. |

The relative speeds are the 3.4.0 interleaved figures from
[the performance guide](../../docs/performance.md#power-of-two-sizes); other
machines and compilers will move them. The planner rules are summarised under
[Choosing a kernel](../../docs/algorithms.md#choosing-a-kernel).

## Usage

```c
#include "fft_common.h"
#include "fft_algorithms.h"

int main(void) {
    enum { n = 1024 };
    complex_t *signal = allocate_complex_array(n);
    if (!signal)
        return 1;

    for (int i = 0; i < n; i++)
        signal[i] = sin(2.0 * PI * 10.0 * i / n);

    if (radix2_dit_fft(signal, n, FFT_FORWARD) != FFT_OK) {
        free_complex_array(signal);
        return 1;
    }

    /* spectrum in signal[0..n/2] */

    if (radix2_dit_fft(signal, n, FFT_INVERSE) != FFT_OK) {
        free_complex_array(signal);
        return 1;
    }

    /* the inverse divides by n, so this is the original signal again */
    free_complex_array(signal);
    return 0;
}
```

## Building and testing

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build
./build/bin/test_all        # exits non-zero if any check fails
```

Add `-DFFT_ENABLE_OPENMP=ON` to build `radix4_omp.c` and the `test_threads`
concurrency test.
