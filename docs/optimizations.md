# Optimization Notes

This page describes what the library does to go fast, the measurements behind
the choices, and what is deliberately absent.

## Algorithm work

The radix-4 kernel used to call `cexp()` three times per butterfly. Twiddles now
advance multiplicatively per stage, as in radix-2 DIT. That moved radix-4 from
3.7-4.0x slower than radix-2 DIT to 0.76x at 4096 points and 0.32-0.51x from
16384 up on the 3.4.0 timing host.

Split-radix used to be radix-2 DIT under a split-radix name. It now performs the
one `n/2` plus two `n/4` decomposition, using a bump-allocated scratch pool
instead of a `malloc` per recursion node. In 3.4.0 the recursion also reads its
input with a stride rather than gathering even and odd samples into scratch at
every level, and the two combine passes share one loop. Peak scratch is about 3n
(pool sized 4n). That cut wall time from 2.0-2.9x radix-2 DIT to 1.1-1.4x, but it
still loses to radix-4. The measurements are in the
[Performance Guide](performance.md).

Mixed-radix previously fell through to a direct O(n²) DFT for non-power-of-two
lengths. The kernel is now a recursive Cooley-Tukey implementation for composite
lengths. That took N=4095 from 1513 ms to 5.90 ms; the 1513 ms figure predates a
timing-harness correction, so the ratio is the meaningful part.

For bit reversal, `bit_reverse()` uses a four-stage swap network for
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

Plans used to build an n-entry twiddle table and an n-entry bit-reversal table
that no kernel read. They cost two allocations and n `cexp()` calls per plan and
have been removed. The kernels compute twiddles as they run.

The radix-2/radix-4 selection threshold moved from 64 to 256 based on the
measured ratios. With `FFT_MEASURE`, planning times the candidates on a private
buffer, interleaves the samples, and requires a repeatable win before replacing
the heuristic's choice. The batch target is 20 ms of wall time, expressed as
`0.020 * CLOCKS_PER_SEC`; a fixed tick count was not portable because `clock()`
counts milliseconds on Windows and microseconds on POSIX.

## AVX2 radix-4

`algorithms/core/radix4_avx2.c` vectorizes the radix-4 stage over two complex
values per 256-bit register. The dispatch and portability details:

- The file is compiled with GCC/Clang `__attribute__((target("avx2")))`; the
  rest of the library keeps the baseline ISA. `-march=native` is not required.
- `radix4_avx2_available()` calls `__builtin_cpu_init()` and
  `__builtin_cpu_supports("avx2")` once per dispatch decision. Non-x86 targets
  and compilers without the builtin get a stub that returns 0.
- The scalar kernel stays in charge below 4096 points, where the per-stage
  twiddle table and vector setup cost more than the two-wide arithmetic saves.
- `complex_t` is two consecutive doubles on every target this file is built
  for. `cmul()` uses `addsub` to get the real and imaginary products from one
  multiply pair; `mul_minus_i` and `mul_plus_i` are shuffle-and-swap operations.

## OpenMP radix-4

`algorithms/core/radix4_omp.c` is built only with `-DFFT_ENABLE_OPENMP=ON`. It
reorders in place, runs the leading radix-2 stage (when `log2(n)` is odd) in a
parallel loop, and then parallelizes the butterfly loop of each radix-4 stage
with `#pragma omp parallel for collapse(2) schedule(static)`. One twiddle table
is built per stage and shared; each butterfly derives `w2 = w1*w1` and
`w3 = w2*w1` instead of reading two more tables.

Below 524288 points, or with one thread, the kernel calls scalar radix-4. That
keeps small transforms out of a parallel region whose launch costs more than the
work. On the 6-core/12-thread verification host it beat AVX2 by about 24% at
524288 and tied it from 1M up. Scaling past one socket was not measured.

## Four-step

`algorithms/core/four_step.c` implements the cache-blocked decomposition for
large power-of-two lengths. With `N = N1*N2` and factors near `sqrt(N)`, the
input is treated as an N1-by-N2 row-major matrix: columns are gathered into
cache-sized scratch and transformed, the result is multiplied by
`W_N^(k1*c)`, rows are transformed, and a 16x16 tiled transpose returns natural
order. Sub-transforms use `radix4_fft_scratch()` and share one scratch buffer, so
blocking does not add an allocation per row or column.

The kernel is a correctness-tested candidate for `FFT_MEASURE`, not a heuristic
choice. On the 3.4.0 hosts it measured 17-30% slower than AVX2 at 1M-4M points:
the column gather, twiddle pass and transpose add roughly 4n element transfers
that the in-cache radix kernels do not pay. It stays available for hosts where
that balance differs.

## Planning state

The hardware probe and the wisdom cache are process-global, so 3.4.0 guards them
with spin locks built on C11 `atomic_flag`. The probe uses a double-checked
acquire/release pair; the cache locks each lookup and store. Export copies a
snapshot of the cache under the lock and formats outside it; import parses into a
local array and replaces the cache under the lock. That keeps a concurrent
export from observing a half-written entry.

A plan is not a synchronization object. The locks cover planning state, not a
plan's buffers, so separate plans can be created and executed concurrently while
a single plan stays thread-confined.

## Measured result

Cost relative to radix-2 DIT at the same N; below 1.00 is faster. One
interleaved run on the WSL host described in the
[Performance Guide](performance.md).

| Algorithm | N=4096 | N=16384 | N=1048576 |
| --- | --- | --- | --- |
| Radix-4 | 0.76 | 0.51 | 0.41 |
| Radix-4 AVX2 | 0.26 | 0.24 | 0.37 |
| Radix-4 OpenMP | — | — | 0.35 |
| Four-Step | 1.12 | 0.63 | 0.45 |
| Split-Radix | 1.35 | 1.32 | 1.31 |

The OpenMP kernel is timed only from 524288 points; below that it calls scalar
radix-4.

## Memory

- `fft_alloc_complex()` and `fft_alloc_real()` return 64-byte aligned buffers.
  The alignment is a cache-line optimization; the AVX2 kernel uses unaligned
  loads because the stage stride does not preserve 32-byte boundaries.
- `FFT_CONSERVE_MEMORY` selects in-place radix-2 DIT for power-of-two lengths
  and mixed-radix for highly composite lengths, avoiding radix-4's O(n)
  digit-reversal buffer and Bluestein's roughly 9n complex numbers of scratch.
- Split-radix allocates a 4n pool and peaks near 3n live.
- Four-step allocates n elements for the reordered output plus about `2*sqrt(n)`
  for the column and sub-transform buffers.
- The 2-D plan allocates one `rows * cols` transpose buffer per execution. It is
  not cached in the plan.

## Compiler flags

The CMake build uses:

```cmake
add_compile_options(-Wall -Wextra -O3)
add_compile_definitions(_GNU_SOURCE _POSIX_C_SOURCE=200809L)
```

The AVX2 kernel is compiled with a function-level target attribute, so AVX2
instructions do not leak into the baseline build and `-march=native` is not
needed for it. The OpenMP kernel is opt-in:

```bash
cmake -S . -B build-omp -DFFT_ENABLE_OPENMP=ON
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

Add `-fopenmp` and define `FFT_ENABLE_OPENMP` through your build system if you
want the threaded kernel. `-DLIB_BUILD` is required when compiling the algorithm
sources yourself: each one keeps a demonstration `main()` unless that macro is
defined, and the link otherwise fails with multiple definitions of `main`.

## What is not implemented

- **AVX-512 / NEON.** The AVX2 kernel is the only SIMD path. A NEON port would
  need its own target attribute and availability check, and NEON has no
  `addsub` equivalent, so the complex multiply needs a different shuffle
  sequence.
- **Fixed point.** There is no Q15 or integer path.
- **Automatic threading.** The OpenMP kernel is opt-in and the planner does not
  change its choice based on the thread count beyond the availability check. A
  program that sets `OMP_NUM_THREADS=1` gets scalar radix-4.
- **Native Windows CUDA.** The CUDA backend reaches the planner and is tested
  on Linux/WSL. It cannot be built natively on Windows: nvcc requires MSVC as
  its host compiler, and MSVC cannot compile the library's C11 `_Complex`
  kernels. Metal is macOS-only and float32; CUDA is double precision.

## Remaining work

1. Measure OpenMP scaling past one socket and decide whether the fixed 524288
   threshold should become a thread-count-aware rule.
2. Fuse the four-step twiddle pass into the row FFT's first stage, or drop the
   kernel if that does not close the 17-30% gap to AVX2.
3. Add an AVX-512 path behind the same runtime-dispatch pattern as AVX2.
4. Add real-input kernels that exploit Hermitian symmetry instead of promoting
   to complex and transforming an all-zero imaginary part.
