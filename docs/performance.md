# Performance Guide

The 3.4.0 power-of-two ratios were measured on the WSL host described below with
the new kernels. The non-power-of-two and GPU tables come from earlier runs on
the same laptop; those kernels did not change in 3.4.0, so their ratios still
describe the code. None of this is a cross-platform comparison. Correctness was
re-checked on Linux, macOS and Windows (see the README).

## How the numbers were measured

For the 3.4.0 table, a small probe times one batch of each kernel per round, in
rotation, so a clock-frequency change moves all candidates together. Before each
batch the work buffer is reset to an impulse: a forward transform scales its
input by `n`, and repeatedly transforming arbitrary data in place would overflow
it to `Inf` and then `NaN`. Each batch is calibrated to span about 150 ms, and
the best of five batches is kept. The clock is `omp_get_wtime()`.

The earlier tables used the same batching idea with `clock()`, keeping the
minimum of three batches. On Windows `clock()` advances about once a millisecond,
so a figure measured from a single small transform is quantisation, not timing.

Absolute milliseconds vary with thermal state and scheduler load. Ratios taken
inside one run are steadier. As a scale reference, one run of the probe put a
1048576-point forward transform at 172 ms for radix-2 DIT, 70 ms for radix-4,
63 ms for AVX2 radix-4, 60 ms for OpenMP radix-4, 77 ms for four-step and 225 ms
for split-radix.

Environment for the 3.4.0 ratios: Intel Core i7-9750H (6 cores, 12 threads,
2.6 GHz base), Ubuntu 22.04 under WSL, GCC 11.4.0, `-O3`, no `-march=native`,
`FFT_ENABLE_OPENMP=ON` so every kernel is compiled in. The OpenMP kernel used 12
threads.

Two limits worth remembering:

- A gap smaller than about 5% is inside the run-to-run spread of this machine.
- The threshold table is a heuristic for one class of machine. `FFT_MEASURE`
  times the candidates on the actual host.

## Power-of-two sizes

Cost relative to radix-2 DIT at the same `n`; below 1.00 is faster. One
representative interleaved run. The OpenMP kernel is only timed from 524288
points up; below that it calls scalar radix-4, so the Radix-4 column applies.

| n | Radix-4 | AVX2 | OpenMP | Four-Step | Split-Radix |
| --- | --- | --- | --- | --- | --- |
| 4096 | 0.76 | **0.26** | — | 1.12 | 1.35 |
| 16384 | 0.51 | **0.24** | — | 0.63 | 1.32 |
| 65536 | 0.32 | **0.20** | — | 0.49 | 1.30 |
| 262144 | 0.33 | **0.24** | — | 0.45 | 1.35 |
| 524288 | 0.38 | 0.29 | **0.22** | 0.40 | 1.13 |
| 1048576 | 0.41 | 0.37 | **0.35** | 0.45 | 1.31 |
| 2097152 | 0.47 | 0.42 | 0.42 | 0.49 | 1.35 |
| 4194304 | 0.49 | 0.42 | 0.42 | 0.54 | 1.41 |

What the table says:

- AVX2 radix-4 is the fastest single-threaded kernel from 4096 to 262144
  points, at 3.8x to 5x the scalar radix-2 speed. Its advantage over scalar
  radix-4 shrinks from about 3x at 4096 to 10-15% at 1M, where memory traffic
  binds.
- OpenMP radix-4 wins at 524288 by about 24% over AVX2 on this 12-thread host,
  then ties it from 1M up. Scaling beyond one socket was not measured.
- Four-step was 17-30% slower than AVX2 on this host at 1M-4M. It lost to AVX2
  on the Windows host as well, and to OpenMP where that kernel was compiled in,
  so the planner leaves it to `FFT_MEASURE` rather than selecting it by size.
- Split-radix now measures 1.1-1.4x radix-2 DIT, down from 2.0-2.9x before the
  strided-input rewrite. It still loses to radix-4 because the recursion carries
  more memory traffic per stage.
- Radix-4 beats radix-2 DIT at every size in the table, by 1.3x at 4096, about
  3x at 65536-262144, and 2.0-2.6x above that. The planner still switches at
  256 points because the crossover moves with the compiler and the
  digit-reversal buffer's cost.

## Non-power-of-two sizes

Ratio of Bluestein to mixed-radix at the same length; below 1.00 means Bluestein
is faster. Measured on the Windows 3.0.0 machine; both kernels are unchanged
since then.

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

The OpenMP kernel's scaling past 12 threads, AVX-512, fixed-point arithmetic and
a sustained batched GPU comparison are not measured. The opt-in Metal and CUDA
backends are measured separately below and are not part of the CPU ratio tables.

`fft_recursive`, `naive_dft` and `optimized_dft` are excluded from the ratio
table. They are reference implementations, and the DFTs are O(n²): a naive
16384-point transform performs about 268 million complex exponentials. The
benchmark does time them up to n=4096 and skips them above that.

## Metal backend

`gpu_benchmark` in a `-DFFT_ENABLE_METAL=ON` build times one forward transform
per size, including host upload, graph execution and host read-back. On the
M1 Pro used for the macOS verification, three runs of that program gave:

| n | CPU (radix-2 DIT) | GPU (MPSGraph) | Ratio |
| --- | --- | --- | --- |
| 1024 | 20-22 us | 7.0-7.9 ms | <0.01x |
| 4096 | 96-101 us | 6.9-8.3 ms | 0.01-0.02x |
| 65536 | 2.01-2.07 ms | 7.3-11.0 ms | 0.19-0.28x |
| 1048576 | 48.0-60.1 ms | 13.1-16.0 ms | 3.3-3.8x |

The fixed cost of an MPSGraph call is several milliseconds, so a single
small transform is slower on the GPU. The GPU only wins on this machine at
1048576 points, and the ratio is for one transform, not a sustained batch.

The backend computes in float32 while the CPU kernels use double. Maximum
relative error against the CPU result was 4.5e-6 at 1024 and 5.3e-4 at 1048576.
These are single-machine one-shot measurements; they are not comparable to the
WSL ratios above.

## CUDA backend

`gpu_benchmark` in a `-DFFT_ENABLE_CUDA=ON` build times one forward transform
per size, including host upload, cuFFT execution and host read-back. On the
GTX 1050 used for verification, under WSL with CUDA 11.5.119:

| n | CPU (radix-2 DIT) | GPU (cuFFT) | Ratio |
| --- | --- | --- | --- |
| 1024 | 0.090 ms | 0.406 ms | 0.22x |
| 4096 | 0.375 ms | 0.558 ms | 0.67x |
| 65536 | 7.308 ms | 2.051 ms | 3.56x |
| 1048576 | 192.278 ms | 14.340 ms | 13.41x |

The GPU is slower for small transforms because of the transfer and launch
overhead, and because the GTX 1050 has little FP64 throughput. It wins from
65536 points on this machine and by 13.4x at 1048576. The CUDA path stays in
double precision; maximum relative error against the CPU kernel was 1.4e-13 at
1024, 1.2e-12 at 4096, 3.6e-11 at 65536 and 3.8e-9 at 1048576. These are
one-shot measurements, not a sustained batched comparison, and they use a
consumer Pascal card rather than a datacenter GPU.

## Operation counts

| Algorithm | Complex multiplications | Complex additions |
| --- | --- | --- |
| DFT | N² | N(N-1) |
| Radix-2 | (N/2) log₂N | N log₂N |
| Radix-4 | (3N/8) log₂N | N log₂N |
| Split-Radix | (N/3) log₂N - 2N/9 + 4/9 | N log₂N |

AVX2, OpenMP and four-step use the radix-4 count; the first two change how the
same arithmetic is scheduled, not how much of it there is. Operation counts
predict relative speed only until memory traffic dominates. Split-radix is the
clearest example: it has the fewest multiplications and is still slower than
radix-4 here.

## Build flags

The CMake build uses `-Wall -Wextra -O3`. Three flags need attention:

- `-march=native` is opt-in (`-DFFT_NATIVE_ARCH=ON`). It produces binaries that
  only run on the build CPU, and with MinGW-w64 GCC 13.2.0 it triggered
  misaligned 32-byte AVX stores that crashed `realtime_analyzer`
  intermittently. The AVX2 kernel does not need it: it is compiled with
  `__attribute__((target("avx2")))` and dispatched at runtime.
- `-ffast-math` is unsuitable for an FFT. It allows reassociation of butterfly
  sums and assumes no NaN or Inf, which makes results non-reproducible.
- `-DFFT_ENABLE_OPENMP=ON` builds the threaded radix-4 kernel. It is off by
  default; the CMake option finds OpenMP and links it into the library.

On POSIX builds the project defines `_POSIX_C_SOURCE=200809L` (through
`fft_common.h`), so `posix_memalign` and `sysconf` are declared under strict
C11.

## Practical advice

- Use a power of two when the application allows it.
- Build with OpenMP and plan once for large transforms; below 524288 points the
  threaded kernel falls back, so the option costs nothing there.
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
