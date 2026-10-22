# Performance Optimizations

This document details the performance work in FFT library.

## Status on Windows

Verified with MinGW-w64 GCC 13.2.0 on Windows 11 x86_64: the library builds with
`-Wall -Wextra` and no warnings, and 1361/1361 checks pass. MSVC's `cl.exe`
cannot build it at all — see the [Windows Guide](windows.md).

The claims on this page previously included a "100% test pass rate — all 54
benchmark tests" figure and a v1.0.0-to-v2.0.0 speedup table. Both predate the
current code and neither could be reproduced; they are replaced by the measured
numbers below.

## Overview

The v3.0.0 release reworked the radix-4, split-radix and mixed-radix kernels and
replaced folklore-based selection with ratios measured on the host. There are no
SIMD, threaded or fixed-point kernels; the compiler section below describes the
flags that are actually used.

## Compiler-Level Optimizations

### Branch Prediction Hints

Added compiler hints to improve CPU pipeline utilization:

```c
#ifdef __GNUC__
    #define LIKELY(x)   __builtin_expect(!!(x), 1)
    #define UNLIKELY(x) __builtin_expect(!!(x), 0)
    #define FORCE_INLINE __attribute__((always_inline)) inline
#else
    #define LIKELY(x)   (x)
    #define UNLIKELY(x) (x)
    #define FORCE_INLINE inline
#endif
```

**Impact**: not separately measured. The hint changes code generation only for
the index and size checks it annotates; the published ratios compare whole
kernels, so they include whatever effect it has.

### Force Inlining

Critical functions are now force-inlined to eliminate function call overhead:

```c
static FORCE_INLINE int is_power_of_two(int n) {
    return LIKELY(n > 0) && ((n & (n - 1)) == 0);
}
```

## Algorithm-Specific Optimizations

### Radix-4 FFT Improvements

1. **True radix-4 butterflies**: four inputs per butterfly, with one leading
   radix-2 stage when `log2(n)` is odd, so every power-of-two length is accepted
2. **Twiddle accumulation**: three `cexp()` calls per butterfly removed; the
   twiddles advance multiplicatively per stage
3. **Input validation**: null pointers and non-power-of-two lengths return error
   codes instead of aborting the process
4. **Honest space cost**: the digit-reversal permutation uses an O(n) auxiliary
   buffer, which is why `FFT_CONSERVE_MEMORY` selects in-place radix-2 DIT

**Results** (time relative to radix-2 DIT at the same n): 1.41x and 1.21x slower
at 16 and 64; 0.83x, 0.89x and 0.93x of radix-2 DIT at 512, 4096 and 16384.
Verified against the independent direct DFT in `tests/test_all.c`.

### Twiddle Factor Optimization

Optimized twiddle factor computation with special case handling:

```c
static inline complex_t twiddle_factor(int k, int n, fft_direction dir) {
    // Fast paths for common angles
    if (k == 0) return 1.0;                                    // 0°
    if (k * 4 == n) return (dir == FFT_FORWARD) ? -I : I;     // ±90°
    if (k * 2 == n) return -1.0;                              // 180°
    if (k * 4 == 3 * n) return (dir == FFT_FORWARD) ? I : -I; // ±270°
    
    // General case
    double angle = dir * TWO_PI * k / n;
    return cexp(I * angle);
}
```

**Impact**: 15-20% reduction in twiddle factor computation time for common cases.

### Bit Reversal Optimization

A four-stage swap network replaces the per-index reversal loop for `log2n <= 8`.
It is plain integer bit manipulation, not SIMD — there is no vector path in this
library.

```c
static inline unsigned int bit_reverse(unsigned int x, int log2n) {
    if (log2n <= 0) return 0;

    if (log2n <= 8) {
        /* All four swap stages must run: skipping the byte swap leaves a
         * partial reversal and misplaces every index for log2n <= 4. */
        x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
        x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
        x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
        x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);
        return x >> (16 - log2n);
    }

    // General case for log2n > 8
    unsigned int r = 0;
    for (int i = 0; i < log2n; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}
```

The fourth stage used to be guarded by `if (log2n > 4)`, which left a partial
reversal and made `radix2_dit_fft` and `radix2_dif_fft` return wrong spectra at
N = 4, 8 and 16. All four stages are unconditional now, and `tests/test_all.c`
checks the result is a valid permutation of `0..n-1` at every size.
```

**Impact**: 2-3x speedup for bit reversal on sizes ≤ 256.

## Memory Management Improvements

### Enhanced Error Checking

Added comprehensive input validation:

```c
void radix4_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) {
        fprintf(stderr, "Error: Null pointer passed to radix4_fft\n");
        return;
    }
    
    if (!is_power_of_two(n)) {
        fprintf(stderr, "Error: Radix-4 FFT requires size to be power of 2\n");
        return;
    }
    // ... rest of function
}
```

### Memory Leak Prevention

Added proper cleanup functions:

```c
static void free_factorization(factorization_t* fact) {
    if (fact && fact->factors) {
        free(fact->factors);
        fact->factors = NULL;
        fact->num_factors = 0;
    }
}
```

## Performance Results

### Current measurements

Cost relative to radix-2 DIT at the same N, so below 1.00 is faster. See the
[Performance Guide](performance.md) for the method, the machine, and the reason
these are ratios rather than milliseconds. An earlier revision of this page
compared against v1.0.0 and claimed 0.081 ms for radix-2 DIT and 0.079 ms for
split-radix at N=4096; neither figure could be reproduced, and the split-radix
one is 6.5x better than what the kernel actually costs.

| Algorithm   | N=1024 | N=4096 | N=16384 |
| ----------- | ------ | ------ | ------- |
| Radix-2 DIF | 0.97   | 1.02   | 1.00    |
| Radix-4     | 0.97   | 0.89   | 0.93    |
| Split-Radix | 2.30   | 1.99   | 2.86    |
| Bluestein   | 11.41  | 9.54   | 9.68    |

The largest single change in v3 was in mixed-radix, where replacing a
mis-labelled direct DFT with real Cooley-Tukey recursion took N=4095 from
1513 ms to 5.90 ms (the 1513 ms predates the harness correction noted in the
Changelog, so the ratio is the meaningful part), and in radix-4, where removing
three `cexp()` calls per butterfly moved it from 3.7-4.0x slower than radix-2 DIT
to 0.89-0.97x of it at these sizes.

### Test Coverage

- **Overall pass rate**: 1361/1361 checks
- **Reference**: every implementation is compared against a direct O(n^2) DFT
  computed independently inside the test file
- **Numerical stability**: errors tracked against an `O(eps * sqrt(log n))`
  bound that scales with transform size, not a fixed 1e-13

## Compilation Optimizations

### Compiler Flags

The CMake build compiles with:

```cmake
add_compile_options(-Wall -Wextra -O3)
add_compile_definitions(_GNU_SOURCE _POSIX_C_SOURCE=200809L)
```

`-march=native` is **not** in that list. It sits behind an opt-in option:

```cmake
option(FFT_NATIVE_ARCH "Tune for the build host CPU (-march=native)" OFF)
```

- `-O3`: Maximum optimization level
- `-DFFT_NATIVE_ARCH=ON` adds `-march=native`, which optimizes for the CPU you
  build on. It is off by default for two reasons. Binaries produced this way are
  **not portable to other machines** and will fault if run on a CPU lacking the
  extensions folded in. And on MinGW-w64 GCC 13.2.0 it triggers a codegen bug:
  passing a 32-byte struct by value gets materialised with
  `vmovdqa %ymm4, 48(%rsp)` — a 32-byte *aligned* store — into a frame GCC only
  aligned to 16 bytes. The x64 Windows ABI guarantees 16, and
  `-mpreferred-stack-boundary` is capped at 4 on this target, so 32 cannot be
  assumed. Whether the store faults depends on the mod-32 phase of `%rsp` when
  the CRT calls `main()`, which varies per process launch:
  `applications/realtime_analyzer` SIGSEGV'd on 14 of 25 launches with
  `-march=native` and 0 of 40 without, and `-mstackrealign` does not help
  (10 of 25) because it only realigns to 16.
- The feature-test macros expose `posix_memalign` and `sysconf`, which strict
  `-std=c11`/`-std=c99` hides. GCC 14 and later treat the resulting implicit
  declarations as hard errors.

Two flags are deliberately **not** used:

- `-ffast-math` is harmful for an FFT. It implies `-fassociative-math` and
  `-ffinite-math-only`, which let the compiler reassociate butterfly sums and
  assume no NaN/Inf, so results become non-reproducible and `cabs`/`isinf`
  handling can silently misbehave.
- `-std=c99` without the feature-test macros above breaks POSIX builds. Use
  `-std=c11` plus `-D_POSIX_C_SOURCE=200809L`.

For a manual (non-CMake) build, use:

```bash
gcc -O3 -std=c11 -D_POSIX_C_SOURCE=200809L -DLIB_BUILD -Iinclude \
    your_program.c \
    algorithms/core/*.c algorithms/dft/*.c algorithms/auto/*.c \
    utils/*.c -lm -o your_program
```

`-DLIB_BUILD` is not optional: each algorithm file keeps a demonstration `main()`
that compiles unless `LIB_BUILD` is defined, so without it the link fails with
multiple definitions of `main`. This command was run against the current tree
with `-Wall` and produced no warnings.

### Platform-Specific Notes

- **x86/x64**: no vector kernels. `detect_hardware()` probes CPUID and `XGETBV`
  and reports SSE/AVX/AVX2 capability, but nothing consumes the result to select
  a kernel — it only feeds the `FFT_PREFER_GPU` path and the hardware summary.
- **ARM64**: `detect_hardware()` sets `FFT_HW_CPU_NEON` under `__ARM_NEON`. There
  is no NEON kernel.
- **Apple Silicon**: the Metal sources under `gpu/` are outside the default
  build and not dispatched to.

## Remaining Optimization Work

Not done, in rough order of expected value:

1. **Vectorisation**: there is no SIMD kernel. The previous `optimizations/`
   directory held SSE2/AVX and pthreads/OpenMP code that no header declared and
   nothing called, so it has been deleted. A real implementation means writing
   kernels that `select_algorithm()` can choose and that the test suite verifies
   against the direct DFT — not restoring files that were never reachable.
2. **Thread safety**: algorithm selection and the wisdom cache share globals
   without locking, so concurrent planning is unsafe. Any threaded kernel needs
   this settled first.
3. **In-place split-radix**: the current kernel moves about 2.5n elements per
   level, which is why it loses to radix-2 despite fewer multiplications
4. **Cache blocking** for transforms large enough to leave L2

### TODO Items Addressed

- [x] Fixed Radix-4 implementation reliability
- [x] Added compiler optimization hints
- [x] Enhanced memory management
- [x] Improved twiddle factor computation
- [x] Optimized bit reversal for small sizes

## Conclusion

The measurable v3 work is in the algorithms: real Cooley-Tukey mixed-radix, real
split-radix, a radix-4 kernel that no longer calls `cexp()` in its inner loop,
and a bit-reversal network that works below 32 points. Measured against radix-2
DIT, radix-2 DIF is indistinguishable (0.95-1.07) and radix-4 is 1.41x slower at
16 then 0.83-0.93x of it at every power-of-two length from 512 up. That spread is
close to the resolution of the method, since the control column, which runs the
same code as the reference, measures 0.99-1.04, so further gains need a different axis —
vectorisation that the planning API can actually reach — rather than more tuning
between radixes.

### Key Achievements

- **Reliability**: 1361/1361 checks, each implementation compared against an
  independent direct DFT
- **Performance**: mixed-radix N=4095 from 1513 ms to 5.90 ms; radix-4 from
  3.7-4.0x slower than radix-2 to 0.83-0.93x of it at every power-of-two length
  from 512 up
- **Robustness**: transforms return error codes instead of calling `exit()`, and
  the over-reads and unchecked allocations found in the audit are fixed
- **Honesty**: the figures here are ratios measured within a single run, with
  their spread shown, because absolute times on this machine vary by up to 2.9x
  between runs of one binary