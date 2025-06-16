#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"
#include <immintrin.h>  // Intel intrinsics

/**
 * SIMD Optimized FFT Implementation
 * 
 * Uses SSE/AVX intrinsics for vectorized computation of FFT.
 * Processes multiple complex numbers simultaneously.
 * 
 * Benefits:
 * - 2-4x speedup on modern processors
 * - Better utilization of CPU resources
 * - Reduced memory bandwidth requirements
 * 
 * Requirements:
 * - CPU with SSE2/AVX support
 * - Proper memory alignment
 */

// Check CPU features at runtime
#ifdef __x86_64__
#include <cpuid.h>

typedef struct {
    int sse2;
    int sse3;
    int ssse3;
    int sse41;
    int sse42;
    int avx;
    int avx2;
    int fma;
} cpu_features_t;

cpu_features_t detect_cpu_features() {
    cpu_features_t features = {0};
    
    unsigned int eax, ebx, ecx, edx;
    
    // Check for SSE/AVX support
    __cpuid(1, eax, ebx, ecx, edx);
    
    features.sse2 = (edx >> 26) & 1;
    features.sse3 = ecx & 1;
    features.ssse3 = (ecx >> 9) & 1;
    features.sse41 = (ecx >> 19) & 1;
    features.sse42 = (ecx >> 20) & 1;
    features.avx = (ecx >> 28) & 1;
    features.fma = (ecx >> 12) & 1;
    
    // Check for AVX2
    __cpuid_count(7, 0, eax, ebx, ecx, edx);
    features.avx2 = (ebx >> 5) & 1;
    
    return features;
}
#endif

// Aligned memory allocation
void* aligned_alloc_wrapper(size_t alignment, size_t size) {
    #ifdef _WIN32
        return _aligned_malloc(size, alignment);
    #else
        void* ptr;
        if (posix_memalign(&ptr, alignment, size) != 0) {
            return NULL;
        }
        return ptr;
    #endif
}

void aligned_free_wrapper(void* ptr) {
    #ifdef _WIN32
        _aligned_free(ptr);
    #else
        free(ptr);
    #endif
}

// Complex number structure for SIMD
typedef struct {
    float* real;
    float* imag;
} complex_float_split_t;

// Allocate aligned complex arrays for SIMD
complex_float_split_t* allocate_simd_complex(int n) {
    complex_float_split_t* arr = (complex_float_split_t*)malloc(sizeof(complex_float_split_t));
    arr->real = (float*)aligned_alloc_wrapper(32, n * sizeof(float));
    arr->imag = (float*)aligned_alloc_wrapper(32, n * sizeof(float));
    return arr;
}

void free_simd_complex(complex_float_split_t* arr) {
    aligned_free_wrapper(arr->real);
    aligned_free_wrapper(arr->imag);
    free(arr);
}

#ifdef __SSE2__
// SSE2 optimized butterfly operation
void butterfly_sse2(float* ar, float* ai, float* br, float* bi, 
                   float wr, float wi) {
    __m128 a_real = _mm_load_ps(ar);
    __m128 a_imag = _mm_load_ps(ai);
    __m128 b_real = _mm_load_ps(br);
    __m128 b_imag = _mm_load_ps(bi);
    
    __m128 w_real = _mm_set1_ps(wr);
    __m128 w_imag = _mm_set1_ps(wi);
    
    // Complex multiplication: (br + i*bi) * (wr + i*wi)
    __m128 temp_real = _mm_sub_ps(_mm_mul_ps(b_real, w_real), 
                                  _mm_mul_ps(b_imag, w_imag));
    __m128 temp_imag = _mm_add_ps(_mm_mul_ps(b_real, w_imag), 
                                  _mm_mul_ps(b_imag, w_real));
    
    // Butterfly computation
    __m128 new_a_real = _mm_add_ps(a_real, temp_real);
    __m128 new_a_imag = _mm_add_ps(a_imag, temp_imag);
    __m128 new_b_real = _mm_sub_ps(a_real, temp_real);
    __m128 new_b_imag = _mm_sub_ps(a_imag, temp_imag);
    
    _mm_store_ps(ar, new_a_real);
    _mm_store_ps(ai, new_a_imag);
    _mm_store_ps(br, new_b_real);
    _mm_store_ps(bi, new_b_imag);
}

// SSE2 optimized radix-2 FFT
void fft_radix2_sse2(complex_float_split_t* x, int n, fft_direction dir) {
    if (!is_power_of_two(n)) {
        fprintf(stderr, "SSE2 FFT requires power-of-2 size\n");
        return;
    }
    
    int log2n = log2_int(n);
    
    // Bit reversal (scalar for now)
    for (int i = 0; i < n; i++) {
        int j = bit_reverse(i, log2n);
        if (i < j) {
            float temp = x->real[i];
            x->real[i] = x->real[j];
            x->real[j] = temp;
            
            temp = x->imag[i];
            x->imag[i] = x->imag[j];
            x->imag[j] = temp;
        }
    }
    
    // FFT computation with SSE2
    for (int stage = 1; stage <= log2n; stage++) {
        int m = 1 << stage;
        int half_m = m >> 1;
        
        float angle = dir * TWO_PI / m;
        float w_real = cos(angle);
        float w_imag = sin(angle);
        
        for (int k = 0; k < n; k += m) {
            float wr = 1.0f, wi = 0.0f;
            
            // Process 4 butterflies at a time when possible
            int j = 0;
            for (; j < half_m - 3; j += 4) {
                butterfly_sse2(&x->real[k + j], &x->imag[k + j],
                             &x->real[k + j + half_m], &x->imag[k + j + half_m],
                             wr, wi);
                
                // Update twiddle factors
                float new_wr = wr * w_real - wi * w_imag;
                float new_wi = wr * w_imag + wi * w_real;
                wr = new_wr;
                wi = new_wi;
            }
            
            // Handle remaining butterflies
            for (; j < half_m; j++) {
                int idx1 = k + j;
                int idx2 = idx1 + half_m;
                
                float tr = x->real[idx2] * wr - x->imag[idx2] * wi;
                float ti = x->real[idx2] * wi + x->imag[idx2] * wr;
                
                x->real[idx2] = x->real[idx1] - tr;
                x->imag[idx2] = x->imag[idx1] - ti;
                x->real[idx1] = x->real[idx1] + tr;
                x->imag[idx1] = x->imag[idx1] + ti;
                
                float new_wr = wr * w_real - wi * w_imag;
                float new_wi = wr * w_imag + wi * w_real;
                wr = new_wr;
                wi = new_wi;
            }
        }
    }
    
    // Scale for inverse
    if (dir == FFT_INVERSE) {
        float scale = 1.0f / n;
        __m128 scale_vec = _mm_set1_ps(scale);
        
        int i = 0;
        for (; i < n - 3; i += 4) {
            __m128 real = _mm_load_ps(&x->real[i]);
            __m128 imag = _mm_load_ps(&x->imag[i]);
            _mm_store_ps(&x->real[i], _mm_mul_ps(real, scale_vec));
            _mm_store_ps(&x->imag[i], _mm_mul_ps(imag, scale_vec));
        }
        
        for (; i < n; i++) {
            x->real[i] *= scale;
            x->imag[i] *= scale;
        }
    }
}
#endif

#ifdef __AVX__
// AVX optimized butterfly (processes 8 complex numbers)
void butterfly_avx(float* ar, float* ai, float* br, float* bi, 
                  float wr, float wi) {
    __m256 a_real = _mm256_load_ps(ar);
    __m256 a_imag = _mm256_load_ps(ai);
    __m256 b_real = _mm256_load_ps(br);
    __m256 b_imag = _mm256_load_ps(bi);
    
    __m256 w_real = _mm256_set1_ps(wr);
    __m256 w_imag = _mm256_set1_ps(wi);
    
    // Complex multiplication
    __m256 temp_real = _mm256_sub_ps(_mm256_mul_ps(b_real, w_real), 
                                     _mm256_mul_ps(b_imag, w_imag));
    __m256 temp_imag = _mm256_add_ps(_mm256_mul_ps(b_real, w_imag), 
                                     _mm256_mul_ps(b_imag, w_real));
    
    // Butterfly
    _mm256_store_ps(ar, _mm256_add_ps(a_real, temp_real));
    _mm256_store_ps(ai, _mm256_add_ps(a_imag, temp_imag));
    _mm256_store_ps(br, _mm256_sub_ps(a_real, temp_real));
    _mm256_store_ps(bi, _mm256_sub_ps(a_imag, temp_imag));
}
#endif

// Benchmark SIMD vs scalar
void benchmark_simd_fft() {
    printf("\nSIMD FFT Performance Comparison:\n");
    printf("================================\n");
    
    #ifdef __x86_64__
    cpu_features_t cpu = detect_cpu_features();
    printf("CPU Features: SSE2=%d, AVX=%d, AVX2=%d, FMA=%d\n",
           cpu.sse2, cpu.avx, cpu.avx2, cpu.fma);
    #endif
    
    int sizes[] = {1024, 4096, 16384};
    
    for (int s = 0; s < 3; s++) {
        int n = sizes[s];
        
        // Scalar version
        complex_t* scalar_data = allocate_complex_array(n);
        for (int i = 0; i < n; i++) {
            scalar_data[i] = ((double)rand() / RAND_MAX) + 
                            I * ((double)rand() / RAND_MAX);
        }
        
        timer_t timer;
        timer_start(&timer);
        radix2_dit_fft(scalar_data, n, FFT_FORWARD);
        timer_stop(&timer);
        double scalar_time = timer.elapsed_ms;
        
        #ifdef __SSE2__
        // SIMD version
        complex_float_split_t* simd_data = allocate_simd_complex(n);
        for (int i = 0; i < n; i++) {
            simd_data->real[i] = creal(scalar_data[i]);
            simd_data->imag[i] = cimag(scalar_data[i]);
        }
        
        timer_start(&timer);
        fft_radix2_sse2(simd_data, n, FFT_FORWARD);
        timer_stop(&timer);
        double simd_time = timer.elapsed_ms;
        
        printf("\nSize %d:\n", n);
        printf("  Scalar: %.3f ms\n", scalar_time);
        printf("  SIMD:   %.3f ms (%.2fx speedup)\n", 
               simd_time, scalar_time / simd_time);
        
        free_simd_complex(simd_data);
        #else
        printf("\nSize %d: SIMD not available on this platform\n", n);
        #endif
        
        free_complex_array(scalar_data);
    }
}

// Main demonstration
int main() {
    printf("SIMD Optimized FFT\n");
    printf("==================\n");
    
    // Test basic functionality
    #ifdef __SSE2__
    printf("\nTesting SSE2 FFT correctness:\n");
    
    int n = 16;
    complex_float_split_t* simd_data = allocate_simd_complex(n);
    complex_t* reference = allocate_complex_array(n);
    
    // Generate test data
    for (int i = 0; i < n; i++) {
        float real = sin(2 * PI * 3 * i / n);
        float imag = 0;
        simd_data->real[i] = real;
        simd_data->imag[i] = imag;
        reference[i] = real + I * imag;
    }
    
    // Compute FFT
    fft_radix2_sse2(simd_data, n, FFT_FORWARD);
    radix2_dit_fft(reference, n, FFT_FORWARD);
    
    // Compare results
    double max_error = 0;
    for (int i = 0; i < n; i++) {
        complex_t simd_result = simd_data->real[i] + I * simd_data->imag[i];
        double error = cabs(simd_result - reference[i]);
        if (error > max_error) max_error = error;
    }
    
    printf("Maximum error: %.2e %s\n", max_error,
           max_error < 1e-5 ? "(PASS)" : "(FAIL)");
    
    free_simd_complex(simd_data);
    free_complex_array(reference);
    #endif
    
    // Performance benchmarks
    benchmark_simd_fft();
    
    // SIMD optimization tips
    printf("\n\nSIMD Optimization Tips:\n");
    printf("=======================\n");
    printf("1. Ensure data alignment (16-byte for SSE, 32-byte for AVX)\n");
    printf("2. Use structure-of-arrays layout for complex numbers\n");
    printf("3. Process multiple elements per instruction\n");
    printf("4. Minimize data movement between registers\n");
    printf("5. Use FMA instructions when available\n");
    printf("6. Consider cache blocking for large FFTs\n");
    printf("7. Profile to identify bottlenecks\n");
    
    // Architecture-specific notes
    printf("\n\nArchitecture Notes:\n");
    printf("===================\n");
    printf("- SSE2: 4 single-precision floats per register\n");
    printf("- AVX: 8 single-precision floats per register\n");
    printf("- AVX-512: 16 single-precision floats per register\n");
    printf("- ARM NEON: 4 single-precision floats per register\n");
    
    return 0;
}
