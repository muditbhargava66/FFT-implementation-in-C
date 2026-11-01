# Core FFT Algorithms

Seven 1-D kernels, all compiled into `libfft_lib_c.a`. Each source file also
keeps a demonstration `main()` behind `#ifndef LIB_BUILD`; those demos are not
self-contained (they call each other), so the CMake build does not compile them.
Use `test_all` or `benchmark_all` instead.

| File | Lengths | Scratch | Notes |
|------|---------|---------|-------|
| `radix2_dit.c` | power of 2 | none | In place; bit-reverses the input before the butterflies. The baseline for every ratio on the performance page. |
| `radix2_dif.c` | power of 2 | none | Bit-reversal at the end instead of the start. Same cost and natural-order output as DIT; choose by data flow, not speed. |
| `radix4.c` | power of 2, not only power of 4 | O(n) digit reversal | A leading radix-2 stage covers an odd `log2(n)`. Saves 25% of the multiplications; measured faster than DIT from n=512, slower at 16-256. `FFT_CONSERVE_MEMORY` selects DIT instead. |
| `split_radix.c` | power of 2 | O(n) | Lowest multiplication count, but 1.6-2.9x slower than DIT in the measured runs: its recursion moves about 2.5n elements per level. |
| `bluestein.c` | any length up to 2^29 | O(n) | Chirp z-transform; pads the convolution to `next_power_of_two(2n-1)`. About 10x DIT at a power of two, and the choice for prime lengths. |
| `mixed_radix.c` | any composite | O(n) | Recursive Cooley-Tukey factorisation with natural-order output. Prime lengths fall back to a direct O(p^2) DFT, which is why the planner sends primes to Bluestein. |
| `recursive_fft.c` | power of 2 | O(n) peak | Divide-and-conquer reference implementation, kept for clarity. Not in the benchmark table. |

The relative speeds are the Windows 11 figures from
[the performance guide](../../docs/performance.md#power-of-two-sizes); other
machines and compilers will move them.

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
