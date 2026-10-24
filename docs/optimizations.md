# Optimization Notes

This page describes what the library does to go fast, the measurements behind
the choices, and what is deliberately absent.

## Algorithm work

**Radix-4 twiddles.** The radix-4 kernel used to call `cexp()` three times per
butterfly. Twiddles now advance multiplicatively per stage, as in radix-2 DIT.
That moved radix-4 from 3.7-4.0x slower than radix-2 DIT to 0.83-0.97x of it at
every power-of-two length from 512 up.

**Split-radix.** The previous implementation was radix-2 DIT under a
split-radix name. It now performs the one `n/2` plus two `n/4` decomposition,
using a bump-allocated scratch pool instead of a `malloc` per recursion node.
The operation count is lower, but wall time is 1.6-2.9x slower than radix-2 DIT
because the recursion moves roughly 2.5n elements per level. The measurements
are in the [Performance Guide](performance.md).

**Mixed-radix.** Non-power-of-two lengths previously fell through to a direct
O(n²) DFT. The kernel is now a recursive Cooley-Tukey implementation for
composite lengths. That took N=4095 from 1513 ms to 5.90 ms; the 1513 ms figure
predates a timing-harness correction, so the ratio is the meaningful part.

**Bit reversal.** `bit_reverse()` uses a four-stage swap network for
`log2n <= 8` and a shift loop above that. All four stages run unconditionally.
An earlier guard skipped the byte swap for `log2n <= 4`, which produced a
partial permutation and wrong spectra at N = 4, 8 and 16.

```c
static inline unsigned int bit_reverse(unsigned int x, int log2n) {
    if (log2n <= 0) return 0;
    if (log2n <= 8) {
        x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
        x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
        x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
        x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);
        return x >> (16 - log2n);
    }
    unsigned int r = 0;
    for (int i = 0; i < log2n; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}
```

**Plans allocate less.** `fft_plan_dft_1d()` used to build an n-entry twiddle
table and an n-entry bit-reversal table that no kernel read. They cost two
allocations and n `cexp()` calls per plan and have been removed. The kernels
compute twiddles as they run.

**Selection threshold.** The radix-2/radix-4 crossover moved from 64 to 256
based on the measured ratios. With `FFT_MEASURE`, planning times the candidates
on a private buffer, interleaves the samples, and requires a repeatable win
before replacing the heuristic's choice. The batch target is 20 ms of wall
time, expressed as `0.020 * CLOCKS_PER_SEC`; a fixed tick count was not portable
because `clock()` counts milliseconds on Windows and microseconds on POSIX.

## Measured result

Cost relative to radix-2 DIT at the same N; below 1.00 is faster.

| Algorithm | N=1024 | N=4096 | N=16384 |
| --- | --- | --- | --- |
| Radix-2 DIF | 0.97 | 1.02 | 1.00 |
| Radix-4 | 0.97 | 0.89 | 0.93 |
| Split-Radix | 2.30 | 1.99 | 2.86 |
| Bluestein | 11.41 | 9.54 | 9.68 |

The full table, the method, and the noise floor are in the
[Performance Guide](performance.md).

## Memory

- `fft_alloc_complex()` and `fft_alloc_real()` return 64-byte aligned buffers.
  The alignment is for cache-line behaviour; there is no SIMD path that needs
  it.
- `FFT_CONSERVE_MEMORY` selects in-place radix-2 DIT for power-of-two lengths
  and mixed-radix for highly composite lengths, avoiding radix-4's O(n)
  digit-reversal buffer and Bluestein's roughly 9n complex numbers of scratch.
- The 2-D plan allocates one `rows * cols` transpose buffer per execution. It is
  not cached in the plan.

## Compiler flags

The CMake build uses:

```cmake
add_compile_options(-Wall -Wextra -O3)
add_compile_definitions(_GNU_SOURCE _POSIX_C_SOURCE=200809L)
```

`-march=native` is behind `-DFFT_NATIVE_ARCH=ON` and is off by default. It makes
the binary non-portable, and with MinGW-w64 GCC 13.2.0 it triggered a
misaligned 32-byte AVX store that made `realtime_analyzer` crash on a fraction
of launches. The option prints a warning when enabled.

`-ffast-math` is not used. It permits reassociation of butterfly sums and
assumes no NaN or Inf, both of which are load-bearing for an FFT.

Feature-test macros are defined in `fft_common.h` as well as on the command
line, so strict `-std=c11` and `-std=c99` builds see `posix_memalign` and
`sysconf`. GCC 14 and later treat the resulting implicit declarations as hard
errors.

For a manual build:

```bash
gcc -O3 -std=c11 -D_POSIX_C_SOURCE=200809L -DLIB_BUILD -Iinclude \
    your_program.c \
    algorithms/core/*.c algorithms/dft/*.c algorithms/auto/*.c \
    utils/*.c -lm -o your_program
```

`-DLIB_BUILD` is required when compiling the algorithm sources yourself: each
one keeps a demonstration `main()` unless that macro is defined, and the link
otherwise fails with multiple definitions of `main`.

## What is not implemented

- **SIMD.** The former `optimizations/` directory held SSE2/AVX butterflies,
  a pthreads/OpenMP stage splitter and Q15 kernels. No header declared them and
  nothing called them, so they were compiled into the archive and never
  executed. They have been deleted. A replacement has to be reachable from
  `select_algorithm()` and verified against the reference DFT.
- **Threading.** The planner and wisdom cache use unlocked globals. The library
  is not thread-safe, and no threaded kernel exists.
- **Fixed point.** There is no Q15 or integer path.
- **GPU execution.** The CUDA and Metal sources are outside the default build;
  `fft_gpu_*` resolves to CPU stubs, and no planner path reaches a GPU kernel.

## Remaining work

In rough order of expected value:

1. Make algorithm selection and the wisdom cache thread-safe, then add a
   threaded kernel that the planner can choose.
2. Write a SIMD kernel with the same test and selection integration as the
   scalar kernels.
3. Reduce split-radix's memory traffic; the lower flop count does not pay off
   with the current layout.
4. Add cache blocking for transforms that exceed L2.
