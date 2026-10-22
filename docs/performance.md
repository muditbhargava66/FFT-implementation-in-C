# Performance Guide

This guide covers performance characteristics, benchmarking results, and optimization strategies for the FFT implementations.

**Release 3.0.0 — 23 October 2026.** Numbers below were measured on the
reference Windows machine described in *How these were measured*.

## Benchmark Results

### How these were measured

Each figure is a cost **relative to radix-2 DIT at the same n**, so smaller than
1.00 means faster. Ratios are taken inside a single run, where both kernels see
the same machine state.

Publishing ratios rather than milliseconds is a deliberate choice. Three runs of
one probe, minutes apart, put radix-2 DIT at n=16384 at **0.501, 1.372 and
1.433 ms** — a 2.9x spread for byte-identical code, which is this laptop's
power and thermal behaviour rather than anything about the library. An absolute
table would have been a description of the machine on the day it was measured.
The ratios above are steady across runs, and each cell shows the min-max range
observed so the size of the residual noise is visible.

Within a run: a batch is timed, and if `clock()` did not resolve it the batch is
grown and retried until it spans at least 100 ms. The minimum of three such
batches is kept, because the minimum discards scheduler interference while a mean
would keep a share of it, and the candidates are visited in rotation so drift
moves them together. Before every repetition the work buffer is restored from a
seed: a forward transform scales its input by n, so re-transforming one buffer in
place overflows it to `Inf` and then `NaN` within a few dozen reps, which measured
about 2x slower than finite data and was enough to invert the radix-2 DIT versus
radix-4 ranking in earlier revisions of this page.

Environment: Intel Core i7-9750H (6 cores, 2.6 GHz base), Windows 11,
MinGW-w64 GCC 13.2.0, built with the project's default options — `-O3` and no
`-march=native`. Native tuning is opt-in (`-DFFT_NATIVE_ARCH=ON`) and is off by
default; see the build section of `docs/optimizations.md` for why, which includes
a MinGW-w64 GCC codegen bug that makes AVX builds crash intermittently.

Two caveats that belong with the numbers rather than in a footnote:

- `clock()` advances about once a millisecond on Windows. Any figure that was
  measured from a single transform, or a batch of a handful, is quantisation
  rather than timing. Earlier revisions of this page contained such figures.
- Repeated runs of *one* kernel move by a few percent even inside one run. Gaps
  smaller than the displayed range are not measurements of the algorithms; they
  are measurements of the machine.

### Power-of-two sizes: cost relative to radix-2 DIT

| n | Radix-2 DIF | Radix-4 | Split-Radix | Mixed-Radix | Bluestein |
| --- | --- | --- | --- | --- | --- |
| 16 | 1.06 (0.93-1.15) | 1.41 (1.38-1.43) | 1.64 (1.57-1.71) | 0.99 (0.93-1.02) | 11.11 (10.93-11.22) |
| 64 | 1.06 (1.05-1.08) | 1.21 (1.18-1.26) | 2.44 (2.10-2.67) | 1.02 (0.96-1.10) | 11.30 (10.99-11.85) |
| 256 | 1.00 (0.96-1.03) | 1.10 (0.95-1.26) | 2.38 (2.31-2.50) | 0.99 (0.94-1.02) | 11.90 (10.84-12.44) |
| 512 | 0.95 (0.92-0.99) | **0.83 (0.81-0.84)** | 2.00 (1.95-2.08) | 1.01 (0.96-1.09) | 10.34 (10.01-10.85) |
| 1024 | 0.97 (0.94-1.02) | 0.97 (0.88-1.13) | 2.30 (2.25-2.37) | 1.04 (1.01-1.06) | 11.41 (10.44-12.22) |
| 4096 | 1.02 (0.97-1.05) | **0.89 (0.83-0.96)** | 1.99 (1.91-2.09) | 1.03 (0.97-1.08) | 9.54 (9.42-9.73) |
| 16384 | 1.00 (0.95-1.05) | **0.93 (0.85-0.99)** | 2.86 (2.58-3.15) | 1.01 (0.96-1.09) | 9.68 (9.62-9.78) |

Reading the table:

- **Mixed-radix is the control.** `mixed_radix_fft()` forwards power-of-two
  lengths straight to `radix2_dit_fft()`, so that column *must* read 1.00. It
  comes out at 0.99-1.04 across seven sizes, which is the size of the residual
  error in the whole method. Any difference smaller than that is not measurable
  here.
- **Radix-2 DIT and DIF are the same speed.** Their column sits at 0.95-1.07 for
  every size, inside the control's own error band.
- **Radix-4 crosses over between 256 and 512.** It is clearly slower at 16
  (1.41x) and 64 (1.21x), a wash at 256 (1.10x, but a 0.95-1.26 range), and
  clearly faster from 512 up (0.83x, 0.89x, 0.93x). `select_algorithm()` uses
  exactly that threshold. Below it the radix-4 kernel's O(n) digit-reversal
  buffer costs more than the flops it saves.
- **Split-radix is reproducibly slower**, 1.6-2.9x, in exchange for the lowest
  multiplication count. Its recursion moves roughly 2.5n elements per level, and
  memory traffic, not flops, is the constraint here.
- **Bluestein costs about 10x at a power of two** because it pads the convolution
  to `next_power_of_two(2n-1)` — about 4n — and runs three transforms of that
  length. It is the most stable column in the table.

### Non-power-of-two: Bluestein against mixed-radix

Ratio of Bluestein's time to mixed-radix's at the same length, so below 1.00
means Bluestein is faster. Absolute milliseconds are deliberately not tabulated
here, for the reason given above.

| n | factorisation | Bluestein / Mixed-Radix | Mixed-Radix slower by |
| --- | --- | --- | --- |
| 97 | prime | 0.042x | 24x |
| 251 | prime | 0.016x | 64x |
| 360 | 2^3 * 3^2 * 5 | 0.479x | 2.1x |
| 1000 | 2^3 * 5^3 | 0.313x | 3.2x |
| 1021 | prime | 0.004x | 241x |
| 4095 | 3 * 5 * 7 * 13 | 0.168x | 6.0x |
| 6561 | 3^8 | 0.280x | 3.6x |

Bluestein is faster at every length measured. Mixed-radix has no factor to
divide by once the length is prime and falls back to a direct O(n^2) DFT, which
is where the two-order-of-magnitude gaps come from — note that the three primes
are the three largest ratios, and the composite lengths the three smallest. It
remains the choice when `FFT_CONSERVE_MEMORY` is set and the length factors into
small primes: Bluestein keeps about 9n complex numbers of scratch, mixed-radix
about n. At power-of-two lengths the same flag keeps radix-2 DIT, which is in
place, instead of radix-4, which digit-reverses through an O(n) buffer.

### What is not measured

There is no SIMD, threaded or fixed-point kernel to measure. An `optimizations/`
directory holding SSE2/AVX butterflies, a pthreads and OpenMP radix-2 stage
splitter, a four-step cache-blocked transform and a Q15 fixed-point path used to
be compiled into `libfft_lib_c.a`, but nothing declared those functions in a
header and nothing called them, so no plan could reach them. They have been
deleted rather than left as unreachable objects. An earlier revision of this page
tabulated "SIMD (AVX2)", "Parallel (4 cores)" and GPU speedups; none of those
configurations could be reproduced here, and a benchmark that cannot be run is
not a result.

The GPU sources under `gpu/` are still on disk. The dispatch in `fft_auto.c` is
real, but a default build resolves every `fft_gpu_*` symbol to
`utils/fft_gpu_stubs.c`, where `fft_gpu_available()` returns 0, so nothing
executes on a GPU and no GPU figure here is measured.

`fft_recursive`, `naive_dft` and `optimized_dft` are excluded from the ratio
table for the opposite reason: they are reference implementations, and the two
DFTs are O(n^2) — a naive 16384-point transform performs 268 million complex
exponentials, which is a demonstration of complexity class, not a benchmark.
`benchmark_all` does time them up to n=4096 and skips them above it.

## Operation Count Analysis

### Theoretical Complexity

| Algorithm   | Complex Multiplications | Complex Additions |
| ----------- | ----------------------- | ----------------- |
| DFT         | N²                      | N(N-1)            |
| Radix-2     | (N/2)log₂N              | N log₂N           |
| Radix-4     | (3N/8)log₂N             | N log₂N           |
| Split-Radix | (N/3)log₂N - 2N/9 + 4/9 | N log₂N           |

### Memory Access Patterns

```
Radix-2 DIT Memory Access Pattern (N=16):

Stage 1: Stride = 8
[0]←→[8], [1]←→[9], [2]←→[10], [3]←→[11], ...

Stage 2: Stride = 4  
[0]←→[4], [1]←→[5], [8]←→[12], [9]←→[13], ...

Stage 3: Stride = 2
[0]←→[2], [1]←→[3], [4]←→[6], [5]←→[7], ...

Stage 4: Stride = 1
[0]←→[1], [2]←→[3], [4]←→[5], [6]←→[7], ...
```

## v3.0.0 Performance Optimizations (23 October 2026)

### Compiler-Level Optimizations

The library now includes several low-level optimizations:

```c
// Branch prediction hints for better CPU pipeline utilization
if (LIKELY(i < j)) {
    // Common case - swap elements
    complex_t temp = x[i];
    x[i] = x[j];
    x[j] = temp;
}

// Force inlining for critical functions
static FORCE_INLINE complex_t twiddle_factor(int k, int n, fft_direction dir) {
    // Optimized special cases
    if (k == 0) return 1.0;
    if (k * 4 == n) return (dir == FFT_FORWARD) ? -I : I;
    // ... general case
}
```

### Bit Reversal Optimization

`bit_reverse()` in `include/fft_common.h` uses a four-stage swap network for
`log2n <= 8` and a shift loop above that. It is plain integer bit manipulation;
there is no SIMD anywhere in this library.

```c
static inline unsigned int bit_reverse(unsigned int x, int log2n) {
    if (log2n <= 0) return 0;
    if (log2n <= 8) {
        x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
        x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
        x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
        x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);   /* unconditional */
        return x >> (16 - log2n);
    }
    /* general case: log2n shift-and-accumulate steps */
}
```

The fourth stage used to be guarded by `if (log2n > 4)`. That guard left a
partial reversal, so `radix2_dit_fft` and `radix2_dif_fft` returned wrong spectra
at N = 4, 8 and 16, and Bluestein and mixed-radix inherited the fault at those
sizes. All four stages run unconditionally now.

### Memory Management Improvements

- Enhanced error checking and validation
- Proper cleanup functions to prevent memory leaks
- Input validation for robustness

## Optimization Strategies

### 1. Algorithm Selection

Choose the right algorithm for your use case:

```c
// Decision tree for algorithm selection
fft_algorithm_t select_fft_algorithm(int size, bool real_time, bool accuracy_critical) {
    if (!is_power_of_two(size)) {
        if (is_prime(size)) {
            return BLUESTEIN_FFT;  // Only option for prime sizes
        } else {
            return MIXED_RADIX_FFT;  // Better for composite sizes
        }
    }
    
    if (real_time && size <= 1024) {
        return RADIX_4_FFT;  // Good balance
    } else if (size >= 16384) {
        return RADIX_2_DIT_FFT;  // Better cache behavior
    } else {
        return SPLIT_RADIX_FFT;  // Minimum operations
    }
}
```

### 2. Memory Optimization

#### Cache-Friendly Access
```c
// Block processing for large FFTs
void cache_optimized_fft(complex_t* x, int n) {
    const int CACHE_LINE = 64;  // bytes
    const int COMPLEX_SIZE = sizeof(complex_t);
    const int BLOCK_SIZE = CACHE_LINE / COMPLEX_SIZE;
    
    // Process in cache-friendly blocks
    for (int block = 0; block < n; block += BLOCK_SIZE) {
        // Process block...
    }
}
```

#### Memory Reuse
```c
// Reuse twiddle factors
typedef struct {
    int size;
    complex_t* twiddles;
} twiddle_cache_t;

twiddle_cache_t* create_twiddle_cache(int n) {
    twiddle_cache_t* cache = malloc(sizeof(twiddle_cache_t));
    cache->size = n;
    cache->twiddles = malloc(n/2 * sizeof(complex_t));
    
    // Precompute all twiddle factors
    for (int k = 0; k < n/2; k++) {
        cache->twiddles[k] = cexp(-2.0 * PI * I * k / n);
    }
    
    return cache;
}
```

### 3. Vectorisation and parallelism — not implemented

There is no SIMD kernel and no threaded kernel in this library. Both are
sketches below, kept because the shape of the work is useful to whoever picks it
up, and labelled so nobody mistakes them for something they can call.

An `optimizations/` directory used to hold SSE2 and AVX butterflies, a pthreads
and OpenMP radix-2 stage splitter, a four-step cache-blocked transform and a Q15
fixed-point path. All of it was compiled into `libfft_lib_c.a`, none of it was
declared in a header, and nothing called it, so no plan could ever reach it. It
has been deleted. Restoring it means more than restoring the files: a kernel has
to be selectable from `select_algorithm()`, verified against the direct DFT by
`tests/test_all.c`, and measured well enough to justify being chosen — and the
library has to become thread-safe first, because selection and the wisdom cache
share unlocked globals.

The AVX2 butterfly would look roughly like this, in struct-of-arrays form so the
loads are contiguous:

```c
/* SKETCH — not in the library, not compiled, not measured. */
void butterfly_avx2(double *ar, double *ai, double *br, double *bi,
                    double wr, double wi) {
    __m256d a_r = _mm256_load_pd(ar), a_i = _mm256_load_pd(ai);
    __m256d b_r = _mm256_load_pd(br), b_i = _mm256_load_pd(bi);
    __m256d w_r = _mm256_set1_pd(wr), w_i = _mm256_set1_pd(wi);

    __m256d t_r = _mm256_sub_pd(_mm256_mul_pd(b_r, w_r),
                                _mm256_mul_pd(b_i, w_i));
    __m256d t_i = _mm256_add_pd(_mm256_mul_pd(b_r, w_i),
                                _mm256_mul_pd(b_i, w_r));

    _mm256_store_pd(ar, _mm256_add_pd(a_r, t_r));
    _mm256_store_pd(ai, _mm256_add_pd(a_i, t_i));
    _mm256_store_pd(br, _mm256_sub_pd(a_r, t_r));
    _mm256_store_pd(bi, _mm256_sub_pd(a_i, t_i));
}
```

Note that this needs the four arrays 32-byte aligned, and that the current
kernels store complex numbers interleaved rather than split, so adopting it
means changing the data layout — which is the real cost, not the intrinsics.

The threaded stage would be the radix-2 butterfly loop over blocks, which is
embarrassingly parallel within a stage because block `k` touches only
`x[k .. k+m-1]`:

```c
/* SKETCH — not in the library, not compiled, not measured. */
#pragma omp parallel for
for (int k = 0; k < n; k += m) {
    complex_t w = 1.0;
    for (int j = 0; j < half_m; j++) {
        int t = k + j, u = t + half_m;
        complex_t temp = x[u] * w;
        x[u] = x[t] - temp;
        x[t] = x[t] + temp;
        w *= w_m;
    }
}
```

The build does not enable OpenMP, so `-fopenmp` would be needed as well.


## Profiling and Analysis

### Using perf (Linux)

```bash
# Profile FFT execution
perf record ./build/bin/benchmark_all
perf report

# Measure cache misses
perf stat -e cache-misses,cache-references ./build/bin/test_all
```

### Using Instruments (macOS)

```bash
# Time Profiler
instruments -t "Time Profiler" ./bin/benchmark_all

# Memory analysis
instruments -t "Allocations" ./bin/benchmark_all
```

### Built-in Profiling

```c
// Use the built-in timer
void profile_fft_algorithms() {
    int sizes[] = {256, 1024, 4096, 16384};
    fft_timer_t timer;
    
    for (int i = 0; i < 4; i++) {
        int n = sizes[i];
        complex_t* data = allocate_complex_array(n);
        
        // Generate test data
        generate_random_signal(data, n);
        
        // Profile each algorithm
        timer_start(&timer);
        radix2_dit_fft(data, n, FFT_FORWARD);
        timer_stop(&timer);
        
        printf("Size %d: %.3f ms\n", n, timer.elapsed_ms);
        
        free_complex_array(data);
    }
}
```

## Performance Tips

### 1. General Guidelines

- **Use power-of-2 sizes** whenever possible
- **Align data to cache boundaries** (64-byte alignment)
- **Minimize memory allocations** in hot paths
- **Precompute twiddle factors** for repeated transforms
- **Use appropriate data types** (float vs double)

### 2. Platform-Specific

#### x86/x64
- No SIMD kernel ships. Enabling AVX2/AVX-512 here would require writing and
  testing a vector kernel that `select_algorithm()` can reach; see the
  vectorisation note in [docs/optimizations.md](optimizations.md).
- `-march=native` is available here as `-DFFT_NATIVE_ARCH=ON`, but it is off by
  default and the numbers above were measured without it. On MinGW-w64 GCC it
  currently produces intermittently crashing binaries; see
  [docs/optimizations.md](optimizations.md#compiler-flags) before turning it on.
- Consider Intel MKL for production

#### ARM
- No NEON kernel ships; a NEON implementation would need the same planner and
  test integration as any other new kernel
- Consider ARM Performance Libraries for production

#### Embedded Systems
- Fixed-point arithmetic would be the usual answer here, but this library has no
  fixed-point path; a Q15 kernel that used to sit in `optimizations/` was
  unreachable from the planning API and has been deleted
- Minimize memory usage
- Consider lookup tables for twiddle factors

### 3. Real-Time Considerations

For real-time applications:

```c
// Pre-allocate all memory
typedef struct {
    complex_t* workspace;
    complex_t* twiddles;
    int size;
} realtime_fft_t;

realtime_fft_t* init_realtime_fft(int size) {
    realtime_fft_t* rt = malloc(sizeof(realtime_fft_t));
    rt->size = size;
    rt->workspace = allocate_complex_array(size);
    rt->twiddles = precompute_twiddles(size);
    return rt;
}

// Lock memory pages (Linux)
#include <sys/mman.h>
void lock_memory(void* addr, size_t size) {
    mlockall(MCL_CURRENT | MCL_FUTURE);
    mlock(addr, size);
}
```

## Benchmarking Your Code

### Simple Benchmark Template

```c
void benchmark_my_application() {
    const int NUM_RUNS = 1000;
    const int WARMUP_RUNS = 100;
    
    // Warmup
    for (int i = 0; i < WARMUP_RUNS; i++) {
        run_my_fft();
    }
    
    // Actual benchmark
    fft_timer_t timer;
    timer_start(&timer);
    
    for (int i = 0; i < NUM_RUNS; i++) {
        run_my_fft();
    }
    
    timer_stop(&timer);
    
    double avg_time = timer.elapsed_ms / NUM_RUNS;
    printf("Average time: %.3f ms\n", avg_time);
    printf("Throughput: %.2f transforms/second\n", 1000.0 / avg_time);
}
```

## Conclusion

Performance optimization is a balance between:
- Algorithm complexity
- Memory access patterns  
- Hardware capabilities
- Implementation complexity

Start with the right algorithm, then optimize based on profiling results. Remember that premature optimization is the root of all evil - profile first, optimize second!
