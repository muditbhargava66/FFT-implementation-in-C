#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
    #include <immintrin.h>  // For SIMD on x86/x64 architectures
#elif defined(__arm__) || defined(__aarch64__)
    #include <arm_neon.h>   // For SIMD on ARM architectures
#endif

/**
 * Optimized DFT Implementation
 * 
 * Applies various optimization techniques to the naive DFT:
 * - Precomputed twiddle factors
 * - Loop unrolling
 * - Cache-friendly memory access
 * - Strength reduction
 * - Real input optimization
 * 
 * Still O(n²) but with better constants
 */

// Precompute and cache twiddle factors
typedef struct {
    complex_t* factors;
    int size;
} twiddle_cache_t;

static twiddle_cache_t* create_twiddle_cache(int n, fft_direction dir) {
    twiddle_cache_t* cache = (twiddle_cache_t*)malloc(sizeof(twiddle_cache_t));
    cache->size = n;
    cache->factors = allocate_complex_array(n * n);
    
    // Precompute all twiddle factors
    for (int k = 0; k < n; k++) {
        for (int j = 0; j < n; j++) {
            cache->factors[k * n + j] = cexp(dir * I * TWO_PI * j * k / n);
        }
    }
    
    return cache;
}

static void free_twiddle_cache(twiddle_cache_t* cache) {
    free_complex_array(cache->factors);
    free(cache);
}

// Optimized DFT with precomputed twiddles
void optimized_dft_cached(complex_t* x, int n, twiddle_cache_t* cache) {
    complex_t* temp = allocate_complex_array(n);
    memcpy(temp, x, n * sizeof(complex_t));
    
    // Use precomputed twiddle factors
    for (int k = 0; k < n; k++) {
        complex_t sum = 0;
        complex_t* twiddle_row = &cache->factors[k * n];
        
        // Unroll by 4 for better performance
        int j = 0;
        for (; j < n - 3; j += 4) {
            sum += temp[j] * twiddle_row[j];
            sum += temp[j+1] * twiddle_row[j+1];
            sum += temp[j+2] * twiddle_row[j+2];
            sum += temp[j+3] * twiddle_row[j+3];
        }
        
        // Handle remaining elements
        for (; j < n; j++) {
            sum += temp[j] * twiddle_row[j];
        }
        
        x[k] = sum;
    }
    
    free_complex_array(temp);
}

// Optimized DFT for real input (exploits symmetry)
void optimized_dft_real(double* x_real, complex_t* X, int n) {
    // For real input, X[k] = conj(X[N-k])
    // So we only need to compute half the outputs
    
    for (int k = 0; k <= n/2; k++) {
        complex_t sum = 0;
        
        // First element (no multiplication needed)
        sum = x_real[0];
        
        // Middle elements (use symmetry)
        for (int j = 1; j < n; j++) {
            double angle = -TWO_PI * j * k / n;
            sum += x_real[j] * (cos(angle) + I * sin(angle));
        }
        
        X[k] = sum;
        
        // Use symmetry for the second half
        if (k > 0 && k < n/2) {
            X[n-k] = conj(sum);
        }
    }
}

// Goertzel algorithm for computing single DFT bin
complex_t goertzel_single_bin(complex_t* x, int n, int k) {
    double omega = TWO_PI * k / n;
    double cos_omega = cos(omega);
    double sin_omega = sin(omega);
    double coeff = 2 * cos_omega;
    
    double s0 = 0, s1 = 0, s2 = 0;
    
    // Process input
    for (int i = 0; i < n; i++) {
        s0 = creal(x[i]) + coeff * s1 - s2;
        s2 = s1;
        s1 = s0;
    }
    
    // Final calculation
    double real = s1 - s2 * cos_omega;
    double imag = s2 * sin_omega;
    
    return real - I * imag;
}

// Optimized DFT using various techniques
void optimized_dft(complex_t* x, int n, fft_direction dir) {
    if (n <= 32) {
        // For small sizes, use cached twiddle factors
        twiddle_cache_t* cache = create_twiddle_cache(n, dir);
        optimized_dft_cached(x, n, cache);
        free_twiddle_cache(cache);
    } else {
        // For larger sizes, compute on the fly with optimization
        complex_t* temp = allocate_complex_array(n);
        memcpy(temp, x, n * sizeof(complex_t));
        
        // Use symmetry and reduce redundant calculations
        complex_t* cos_table = allocate_complex_array(n);
        complex_t* sin_table = allocate_complex_array(n);
        
        for (int j = 0; j < n; j++) {
            double angle = dir * TWO_PI * j / n;
            cos_table[j] = cos(angle);
            sin_table[j] = sin(angle);
        }
        
        for (int k = 0; k < n; k++) {
            complex_t sum = temp[0];  // j=0 case
            
            // Use recurrence relation for twiddle factors
            complex_t w = cos_table[k] + I * sin_table[k];
            complex_t w_j = w;
            
            for (int j = 1; j < n; j++) {
                sum += temp[j] * w_j;
                w_j *= w;  // Incremental computation
            }
            
            x[k] = sum;
        }
        
        free_complex_array(temp);
        free_complex_array(cos_table);
        free_complex_array(sin_table);
    }
    
    // Scale for inverse
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

// Benchmark different optimizations
void benchmark_optimizations() {
    printf("\nOptimization Techniques Comparison:\n");
    printf("==================================\n");
    
    int sizes[] = {16, 32, 64, 128};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    for (int s = 0; s < num_sizes; s++) {
        int n = sizes[s];
        printf("\nSize N = %d:\n", n);
        
        complex_t* x_naive = allocate_complex_array(n);
        complex_t* x_opt = allocate_complex_array(n);
        complex_t* x_cached = allocate_complex_array(n);
        double* x_real = (double*)malloc(n * sizeof(double));
        complex_t* X_real = allocate_complex_array(n);
        
        // Generate test data
        for (int i = 0; i < n; i++) {
            double val = sin(2 * PI * 3 * i / n) + 0.5 * cos(2 * PI * 7 * i / n);
            x_naive[i] = x_opt[i] = x_cached[i] = val;
            x_real[i] = val;
        }
        
        timer_t timer;
        
        // Naive DFT
        timer_start(&timer);
        naive_dft(x_naive, n, FFT_FORWARD);
        timer_stop(&timer);
        double time_naive = timer.elapsed_ms;
        
        // Optimized DFT
        timer_start(&timer);
        optimized_dft(x_opt, n, FFT_FORWARD);
        timer_stop(&timer);
        double time_opt = timer.elapsed_ms;
        
        // Cached twiddle DFT
        twiddle_cache_t* cache = create_twiddle_cache(n, FFT_FORWARD);
        timer_start(&timer);
        optimized_dft_cached(x_cached, n, cache);
        timer_stop(&timer);
        double time_cached = timer.elapsed_ms;
        free_twiddle_cache(cache);
        
        // Real-input DFT
        timer_start(&timer);
        optimized_dft_real(x_real, X_real, n);
        timer_stop(&timer);
        double time_real = timer.elapsed_ms;
        
        printf("  Naive DFT:        %.3f ms (baseline)\n", time_naive);
        printf("  Optimized DFT:    %.3f ms (%.1fx speedup)\n", 
               time_opt, time_naive / time_opt);
        printf("  Cached Twiddle:   %.3f ms (%.1fx speedup)\n", 
               time_cached, time_naive / time_cached);
        printf("  Real-input DFT:   %.3f ms (%.1fx speedup)\n", 
               time_real, time_naive / time_real);
        
        // Verify correctness
        double error = 0;
        for (int i = 0; i < n; i++) {
            error += cabs(x_naive[i] - x_opt[i]);
        }
        printf("  Error: %.2e\n", error);
        
        free_complex_array(x_naive);
        free_complex_array(x_opt);
        free_complex_array(x_cached);
        free(x_real);
        free_complex_array(X_real);
    }
}

// Main demonstration
int main() {
    printf("Optimized DFT Implementation\n");
    printf("============================\n");
    
    // Test Goertzel algorithm for single bin
    printf("\nGoertzel Algorithm Demo (single bin computation):\n");
    printf("------------------------------------------------\n");
    
    int n = 32;
    complex_t* signal = allocate_complex_array(n);
    
    // Generate signal with known frequency content
    int freq_bin = 5;
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * freq_bin * i / n);
    }
    
    printf("Signal: Sinusoid at bin %d\n", freq_bin);
    printf("Computing DFT bins using Goertzel:\n");
    
    for (int k = 0; k < 10; k++) {
        complex_t bin_value = goertzel_single_bin(signal, n, k);
        double magnitude = cabs(bin_value);
        printf("  Bin %d: magnitude = %.3f %s\n", 
               k, magnitude, k == freq_bin ? "<-- Peak" : "");
    }
    
    // Benchmark optimizations
    benchmark_optimizations();
    
    // Memory access pattern analysis
    printf("\n\nMemory Access Pattern Analysis:\n");
    printf("===============================\n");
    printf("Naive DFT: Strided access pattern (poor cache usage)\n");
    printf("Optimized: Sequential access with unrolling (better cache usage)\n");
    printf("Cached:    Pre-computed table lookup (trades memory for speed)\n");
    
    // Show optimization benefits for different scenarios
    printf("\n\nOptimization Recommendations:\n");
    printf("=============================\n");
    printf("- Small N (< 32):     Use cached twiddle factors\n");
    printf("- Medium N (32-512):  Use optimized with unrolling\n");
    printf("- Real input only:    Use real-input optimization (2x speedup)\n");
    printf("- Single bins:        Use Goertzel algorithm\n");
    printf("- Large N:            Consider FFT algorithms instead\n");
    
    free_complex_array(signal);
    
    return 0;
}
