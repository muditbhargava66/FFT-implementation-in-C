# Performance Guide

The numbers on this page were measured on the Windows machine described below.
They describe that machine and that build; they are not a cross-platform
comparison. Correctness was re-checked on Linux and macOS (see the README), but
those runs did not produce the table below.

## How the numbers were measured

Each figure is a cost relative to radix-2 DIT at the same `n`. Ratios are taken
inside one run, so both kernels see the same CPU state. Below 1.00 means faster
than radix-2 DIT.

A batch of transforms is timed; if `clock()` cannot resolve the interval, the
batch grows and is retried until it spans at least 20 ms. The minimum of three
batches is kept, because the minimum discards scheduler interference that a mean
would include. Candidates are visited in rotation so a clock-frequency change
moves them together. Before each repetition the work buffer is restored from a
seed: a forward transform scales its input by `n`, so repeatedly transforming
one buffer in place overflows it to `Inf` and then `NaN`.

Absolute milliseconds are not tabulated because three runs of one probe, minutes
apart, put radix-2 DIT at n=16384 at 0.501, 1.372 and 1.433 ms — a 2.9x spread
for identical code. Ratios stay steady across those runs.

Environment: Intel Core i7-9750H (6 cores, 2.6 GHz base), Windows 11,
MinGW-w64 GCC 13.2.0, default build options (`-O3`, no `-march=native`).

Two limits worth remembering:

- `clock()` advances about once a millisecond on Windows. A figure measured from
  a single small transform is quantisation, not timing.
- Repeated runs of one kernel move by a few percent. Gaps smaller than the
  displayed range are not measurements of the algorithms.

## Power-of-two sizes

Range in parentheses is the observed min-max across the sampled runs.

| n | Radix-2 DIF | Radix-4 | Split-Radix | Mixed-Radix | Bluestein |
| --- | --- | --- | --- | --- | --- |
| 16 | 1.06 (0.93-1.15) | 1.41 (1.38-1.43) | 1.64 (1.57-1.71) | 0.99 (0.93-1.02) | 11.11 (10.93-11.22) |
| 64 | 1.06 (1.05-1.08) | 1.21 (1.18-1.26) | 2.44 (2.10-2.67) | 1.02 (0.96-1.10) | 11.30 (10.99-11.85) |
| 256 | 1.00 (0.96-1.03) | 1.10 (0.95-1.26) | 2.38 (2.31-2.50) | 0.99 (0.94-1.02) | 11.90 (10.84-12.44) |
| 512 | 0.95 (0.92-0.99) | **0.83 (0.81-0.84)** | 2.00 (1.95-2.08) | 1.01 (0.96-1.09) | 10.34 (10.01-10.85) |
| 1024 | 0.97 (0.94-1.02) | 0.97 (0.88-1.13) | 2.30 (2.25-2.37) | 1.04 (1.01-1.06) | 11.41 (10.44-12.22) |
| 4096 | 1.02 (0.97-1.05) | **0.89 (0.83-0.96)** | 1.99 (1.91-2.09) | 1.03 (0.97-1.08) | 9.54 (9.42-9.73) |
| 16384 | 1.00 (0.95-1.05) | **0.93 (0.85-0.99)** | 2.86 (2.58-3.15) | 1.01 (0.96-1.09) | 9.68 (9.62-9.78) |

What the table shows:

- **Mixed-radix is the control.** `mixed_radix_fft()` forwards power-of-two
  lengths to `radix2_dit_fft()`, so its column must read 1.00. It measures
  0.99-1.04, which is the resolution of the method.
- **Radix-2 DIT and DIF are indistinguishable** at these sizes, 0.95-1.07.
- **Radix-4 crosses over between 256 and 512.** It is slower at 16 (1.41x) and
  64 (1.21x), a wash at 256, and faster from 512 up (0.83x, 0.89x, 0.93x). That
  is where `select_algorithm()` puts the threshold. Below it, the O(n)
  digit-reversal buffer costs more than the saved multiplications.
- **Split-radix is slower in wall time**, 1.6-2.9x, despite the lowest
  multiplication count. Its recursion moves roughly 2.5n elements per level, and
  memory traffic is the constraint here.
- **Bluestein costs about 10x at a power of two** because it pads the
  convolution to `next_power_of_two(2n-1)` and runs three transforms of that
  length.

## Non-power-of-two sizes

Ratio of Bluestein to mixed-radix at the same length; below 1.00 means Bluestein
is faster.

| n | Factorisation | Bluestein / Mixed-Radix | Mixed-Radix slower by |
| --- | --- | --- | --- |
| 97 | prime | 0.042x | 24x |
| 251 | prime | 0.016x | 64x |
| 360 | 2^3 * 3^2 * 5 | 0.479x | 2.1x |
| 1000 | 2^3 * 5^3 | 0.313x | 3.2x |
| 1021 | prime | 0.004x | 241x |
| 4095 | 3 * 5 * 7 * 13 | 0.168x | 6.0x |
| 6561 | 3^8 | 0.280x | 3.6x |

Bluestein wins at every length measured. Mixed-radix has no factor to split on
at a prime and falls back to a direct O(n²) DFT, which accounts for the largest
gaps. It remains useful under `FFT_CONSERVE_MEMORY`: Bluestein needs about 9n
complex numbers of scratch, mixed-radix about n.

## What is not measured

There is no SIMD, threaded or fixed-point kernel to measure. The GPU sources
under `gpu/` are not in the default build, and a default build resolves every
`fft_gpu_*` symbol to a CPU stub, so no GPU timing exists here.

`fft_recursive`, `naive_dft` and `optimized_dft` are excluded from the ratio
table. They are reference implementations, and the DFTs are O(n²): a naive
16384-point transform performs about 268 million complex exponentials. The
benchmark does time them up to n=4096 and skips them above that.

## Operation counts

| Algorithm | Complex multiplications | Complex additions |
| --- | --- | --- |
| DFT | N² | N(N-1) |
| Radix-2 | (N/2) log₂N | N log₂N |
| Radix-4 | (3N/8) log₂N | N log₂N |
| Split-Radix | (N/3) log₂N - 2N/9 + 4/9 | N log₂N |

Operation counts predict relative speed only until memory traffic dominates.
Split-radix is the clearest example: it has the fewest multiplications and is
the slowest of the power-of-two kernels here.

## Build flags

The CMake build uses `-Wall -Wextra -O3`. Two flags are deliberately not used:

- `-march=native` is opt-in (`-DFFT_NATIVE_ARCH=ON`). It produces binaries that
  only run on the build CPU, and with MinGW-w64 GCC 13.2.0 it triggered
  misaligned 32-byte AVX stores that crashed `realtime_analyzer`
  intermittently. The published numbers were measured without it.
- `-ffast-math` is unsuitable for an FFT. It allows reassociation of butterfly
  sums and assumes no NaN or Inf, which makes results non-reproducible.

On POSIX builds the project defines `_POSIX_C_SOURCE=200809L` (through
`fft_common.h`), so `posix_memalign` and `sysconf` are declared under strict
C11.

## Practical advice

- Use a power of two when the application allows it.
- Use `fft_plan_r2c_1d()` for real input: it returns `n/2+1` bins and avoids the
  work of transforming an all-zero imaginary part.
- Reuse a plan for repeated transforms of the same size; planning with
  `FFT_MEASURE` is deliberately not free.
- Use `FFT_CONSERVE_MEMORY` when the working set matters more than speed.
- Avoid timing a single small transform. Time a batch, or use
  `QueryPerformanceCounter` / `clock_gettime` directly.
- For 2-D work, both dimensions must be powers of two. The transform transposes
  between row and column passes.

## Profiling

On Linux:

```bash
perf record ./build/bin/benchmark_all
perf report
perf stat -e cache-misses,cache-references ./build/bin/test_all
```

On macOS:

```bash
instruments -t "Time Profiler" ./build/bin/benchmark_all
instruments -t "Allocations" ./build/bin/benchmark_all
```

The built-in timer (`fft_timer_t`, `timer_start`/`timer_stop`) reports
milliseconds using `clock()`. It is sufficient for batches, not for one small
transform.
