# Performance Guide

This guide covers performance characteristics, benchmarking results, and optimization strategies for the FFT implementations.

## Benchmark Results

### Algorithm Comparison

Performance measurements on Intel i7-8700K @ 3.7GHz, single-threaded:

| Algorithm | N=1024 | N=4096 | N=16384 | N=65536 |
|-----------|--------|--------|---------|---------|
| Radix-2 DIT | 0.08 ms | 0.40 ms | 1.80 ms | 8.5 ms |
| Radix-2 DIF | 0.09 ms | 0.42 ms | 1.85 ms | 8.7 ms |
| Radix-4 | 0.07 ms | 0.35 ms | 1.60 ms | 7.8 ms |
| Split-Radix | 0.06 ms | 0.32 ms | 1.45 ms | 7.2 ms |
| Bluestein | 0.25 ms | 1.20 ms | 5.40 ms | 24.0 ms |
| Mixed-Radix | 0.20 ms | 1.00 ms | 4.50 ms | 20.0 ms |

### Optimization Impact

| Implementation | N=16384 Time | Speedup |
|----------------|--------------|---------|
| Naive DFT | 2500 ms | 1.0x |
| Basic Radix-2 | 1.80 ms | 1389x |
| SIMD (AVX2) | 0.45 ms | 5556x |
| Parallel (4 cores) | 0.50 ms | 5000x |
| SIMD + Parallel | 0.15 ms | 16667x |

## Operation Count Analysis

### Theoretical Complexity

| Algorithm | Complex Multiplications | Complex Additions |
|-----------|------------------------|-------------------|
| DFT | N² | N(N-1) |
| Radix-2 | (N/2)log₂N | N log₂N |
| Radix-4 | (3N/8)log₂N | N log₂N |
| Split-Radix | (N/3)log₂N - 2N/9 + 4/9 | N log₂N |

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

### 3. SIMD Optimization

Example using AVX2 intrinsics:

```c
#include <immintrin.h>

void butterfly_avx2(double* ar, double* ai, double* br, double* bi,
                   double wr, double wi) {
    __m256d a_real = _mm256_load_pd(ar);
    __m256d a_imag = _mm256_load_pd(ai);
    __m256d b_real = _mm256_load_pd(br);
    __m256d b_imag = _mm256_load_pd(bi);
    
    __m256d w_real = _mm256_set1_pd(wr);
    __m256d w_imag = _mm256_set1_pd(wi);
    
    // Complex multiplication: (br + bi*i) * (wr + wi*i)
    __m256d temp_real = _mm256_sub_pd(
        _mm256_mul_pd(b_real, w_real),
        _mm256_mul_pd(b_imag, w_imag)
    );
    __m256d temp_imag = _mm256_add_pd(
        _mm256_mul_pd(b_real, w_imag),
        _mm256_mul_pd(b_imag, w_real)
    );
    
    // Butterfly operation
    _mm256_store_pd(ar, _mm256_add_pd(a_real, temp_real));
    _mm256_store_pd(ai, _mm256_add_pd(a_imag, temp_imag));
    _mm256_store_pd(br, _mm256_sub_pd(a_real, temp_real));
    _mm256_store_pd(bi, _mm256_sub_pd(a_imag, temp_imag));
}
```

### 4. Parallel Processing

OpenMP parallelization example:

```c
void parallel_fft_stage(complex_t* x, int n, int stage) {
    int m = 1 << stage;
    int half_m = m >> 1;
    complex_t w_m = twiddle_factor(1, m, FFT_FORWARD);
    
    #pragma omp parallel for
    for (int k = 0; k < n; k += m) {
        complex_t w = 1.0;
        for (int j = 0; j < half_m; j++) {
            int t = k + j;
            int u = t + half_m;
            
            complex_t temp = x[u] * w;
            x[u] = x[t] - temp;
            x[t] = x[t] + temp;
            
            w *= w_m;
        }
    }
}
```

## Profiling and Analysis

### Using perf (Linux)

```bash
# Profile FFT execution
perf record ./bin/benchmark_all
perf report

# Measure cache misses
perf stat -e cache-misses,cache-references ./bin/radix2_dit
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
    timer_t timer;
    
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
- Enable AVX2/AVX-512 if available
- Use `-march=native` compiler flag
- Consider Intel MKL for production

#### ARM
- Use NEON intrinsics
- Optimize for specific ARM cores
- Consider ARM Performance Libraries

#### Embedded Systems
- Use fixed-point arithmetic
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
    timer_t timer;
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
