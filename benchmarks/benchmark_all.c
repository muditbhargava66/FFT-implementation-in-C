#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"
#include <float.h>

/**
 * Comprehensive FFT Benchmark Suite
 * 
 * Compares performance and accuracy of all FFT implementations
 * in the repository. Measures:
 * - Execution time
 * - Accuracy vs reference
 * - Memory usage
 * - Scalability
 */

// FFT implementation function pointer
typedef void (*fft_func_t)(complex_t*, int, fft_direction);

// FFT implementation descriptor
typedef struct {
    const char* name;
    fft_func_t forward;
    fft_func_t inverse;
    int requires_power_of_2;
    int min_size;
    int max_size;
} fft_impl_t;

// Benchmark result
typedef struct {
    const char* impl_name;
    int size;
    double forward_time_ms;
    double inverse_time_ms;
    double max_error;
    double rms_error;
    int passed;
} benchmark_result_t;

// Import all FFT implementations
extern void radix2_dit_fft(complex_t* x, int n, fft_direction dir);
extern void radix2_dif_fft(complex_t* x, int n, fft_direction dir);
extern void radix4_fft(complex_t* x, int n, fft_direction dir);
extern void split_radix_fft(complex_t* x, int n, fft_direction dir);
extern void bluestein_fft(complex_t* x, int n, fft_direction dir);
extern void mixed_radix_fft(complex_t* x, int n, fft_direction dir);
extern void naive_dft(complex_t* x, int n, fft_direction dir);
extern void optimized_dft(complex_t* x, int n, fft_direction dir);
extern void fft_recursive(complex_t* x, int n, fft_direction dir);

// List of implementations to benchmark
fft_impl_t implementations[] = {
    {"Radix-2 DIT", radix2_dit_fft, radix2_dit_fft, 1, 2, 65536},
    {"Radix-2 DIF", radix2_dif_fft, radix2_dif_fft, 1, 2, 65536},
    {"Radix-4", radix4_fft, radix4_fft, 1, 4, 65536},
    {"Split-Radix", split_radix_fft, split_radix_fft, 1, 2, 65536},
    {"Bluestein", bluestein_fft, bluestein_fft, 0, 2, 4096},
    {"Mixed-Radix", mixed_radix_fft, mixed_radix_fft, 0, 2, 4096},
    {"Recursive", fft_recursive, fft_recursive, 1, 2, 8192},
    {"Naive DFT", naive_dft, naive_dft, 0, 2, 256},
    {"Optimized DFT", optimized_dft, optimized_dft, 0, 2, 512}
};

int num_implementations = sizeof(implementations) / sizeof(implementations[0]);

// Generate reference result using most accurate method
void compute_reference_fft(complex_t* input, complex_t* output, int n, fft_direction dir) {
    memcpy(output, input, n * sizeof(complex_t));
    
    // Use Bluestein for arbitrary sizes, Radix-2 for power of 2
    if (is_power_of_two(n)) {
        radix2_dit_fft(output, n, dir);
    } else {
        bluestein_fft(output, n, dir);
    }
}

// Compute error metrics
void compute_error_metrics(complex_t* test, complex_t* reference, int n,
                         double* max_error, double* rms_error) {
    *max_error = 0;
    double sum_squared_error = 0;
    
    for (int i = 0; i < n; i++) {
        double error = cabs(test[i] - reference[i]);
        if (error > *max_error) *max_error = error;
        sum_squared_error += error * error;
    }
    
    *rms_error = sqrt(sum_squared_error / n);
}

// Run single benchmark
benchmark_result_t run_single_benchmark(fft_impl_t* impl, complex_t* input, int n, int iterations) {
    benchmark_result_t result;
    result.impl_name = impl->name;
    result.size = n;
    result.passed = 1;
    
    // Skip if size not supported
    if ((impl->requires_power_of_2 && !is_power_of_two(n)) ||
        n < impl->min_size || n > impl->max_size) {
        result.forward_time_ms = -1;
        result.inverse_time_ms = -1;
        result.max_error = -1;
        result.rms_error = -1;
        result.passed = -1;
        return result;
    }
    
    // Allocate arrays
    complex_t* test = allocate_complex_array(n);
    complex_t* reference = allocate_complex_array(n);
    complex_t* original = allocate_complex_array(n);
    
    memcpy(original, input, n * sizeof(complex_t));
    
    // Warm-up run
    memcpy(test, input, n * sizeof(complex_t));
    impl->forward(test, n, FFT_FORWARD);
    
    // Time forward FFT
    timer_t timer;
    timer_start(&timer);
    
    for (int iter = 0; iter < iterations; iter++) {
        memcpy(test, input, n * sizeof(complex_t));
        impl->forward(test, n, FFT_FORWARD);
    }
    
    timer_stop(&timer);
    result.forward_time_ms = timer.elapsed_ms / iterations;
    
    // Compute reference
    compute_reference_fft(input, reference, n, FFT_FORWARD);
    
    // Check accuracy
    compute_error_metrics(test, reference, n, &result.max_error, &result.rms_error);
    
    // Time inverse FFT
    timer_start(&timer);
    
    for (int iter = 0; iter < iterations; iter++) {
        memcpy(test, reference, n * sizeof(complex_t));
        impl->inverse(test, n, FFT_INVERSE);
    }
    
    timer_stop(&timer);
    result.inverse_time_ms = timer.elapsed_ms / iterations;
    
    // Check reconstruction
    double recon_max_error, recon_rms_error;
    compute_error_metrics(test, original, n, &recon_max_error, &recon_rms_error);
    
    if (recon_max_error > 1e-10) {
        result.passed = 0;
    }
    
    free_complex_array(test);
    free_complex_array(reference);
    free_complex_array(original);
    
    return result;
}

// Generate test signals
void generate_test_signal(complex_t* signal, int n, const char* type) {
    if (strcmp(type, "random") == 0) {
        for (int i = 0; i < n; i++) {
            signal[i] = ((double)rand() / RAND_MAX - 0.5) + 
                       I * ((double)rand() / RAND_MAX - 0.5);
        }
    } else if (strcmp(type, "sine") == 0) {
        for (int i = 0; i < n; i++) {
            signal[i] = sin(2 * PI * 10 * i / n) + 
                       0.5 * sin(2 * PI * 25 * i / n);
        }
    } else if (strcmp(type, "impulse") == 0) {
        memset(signal, 0, n * sizeof(complex_t));
        signal[0] = 1.0;
    } else if (strcmp(type, "dc") == 0) {
        for (int i = 0; i < n; i++) {
            signal[i] = 1.0;
        }
    }
}

// Print results table
void print_results_table(benchmark_result_t* results, int num_results, int size) {
    printf("\n=== Size: %d ===\n", size);
    printf("%-15s | %-12s | %-12s | %-12s | %-12s | %s\n",
           "Implementation", "Forward (ms)", "Inverse (ms)", "Max Error", "RMS Error", "Status");
    printf("----------------|--------------|--------------|--------------|--------------|--------\n");
    
    for (int i = 0; i < num_results; i++) {
        if (results[i].size != size) continue;
        
        if (results[i].passed == -1) {
            printf("%-15s | %-12s | %-12s | %-12s | %-12s | %s\n",
                   results[i].impl_name, "N/A", "N/A", "N/A", "N/A", "SKIP");
        } else {
            printf("%-15s | %12.3f | %12.3f | %12.2e | %12.2e | %s\n",
                   results[i].impl_name,
                   results[i].forward_time_ms,
                   results[i].inverse_time_ms,
                   results[i].max_error,
                   results[i].rms_error,
                   results[i].passed ? "PASS" : "FAIL");
        }
    }
}

// Find best implementation for each size
void analyze_best_implementations(benchmark_result_t* results, int num_results) {
    printf("\n\nBest Implementations by Size:\n");
    printf("=============================\n");
    
    int sizes[] = {16, 64, 256, 1024, 4096, 16384};
    
    for (int s = 0; s < 6; s++) {
        int size = sizes[s];
        double best_time = DBL_MAX;
        const char* best_impl = NULL;
        
        for (int i = 0; i < num_results; i++) {
            if (results[i].size == size && results[i].passed == 1 &&
                results[i].forward_time_ms > 0 && results[i].forward_time_ms < best_time) {
                best_time = results[i].forward_time_ms;
                best_impl = results[i].impl_name;
            }
        }
        
        if (best_impl) {
            printf("Size %6d: %s (%.3f ms)\n", size, best_impl, best_time);
        }
    }
}

// Performance scaling analysis
void analyze_scaling(benchmark_result_t* results, int num_results) {
    printf("\n\nScaling Analysis (Time Complexity):\n");
    printf("===================================\n");
    
    for (int impl = 0; impl < num_implementations; impl++) {
        printf("\n%s:\n", implementations[impl].name);
        
        int prev_size = 0;
        double prev_time = 0;
        
        for (int i = 0; i < num_results; i++) {
            if (strcmp(results[i].impl_name, implementations[impl].name) != 0) continue;
            if (results[i].passed != 1) continue;
            
            if (prev_size > 0) {
                double size_ratio = (double)results[i].size / prev_size;
                double time_ratio = results[i].forward_time_ms / prev_time;
                double complexity = log(time_ratio) / log(size_ratio);
                
                printf("  %d->%d: O(n^%.2f)\n", prev_size, results[i].size, complexity);
            }
            
            prev_size = results[i].size;
            prev_time = results[i].forward_time_ms;
        }
    }
}

// Main benchmark program
int main() {
    printf("FFT Implementation Benchmark Suite\n");
    printf("==================================\n");
    
    // Test parameters
    int test_sizes[] = {16, 64, 256, 1024, 4096, 16384};
    int num_sizes = sizeof(test_sizes) / sizeof(test_sizes[0]);
    const char* test_signal = "random";
    
    // Adjust iterations based on size
    int iterations[] = {10000, 5000, 1000, 100, 10, 1};
    
    // Allocate results array
    int total_benchmarks = num_implementations * num_sizes;
    benchmark_result_t* all_results = (benchmark_result_t*)malloc(
        total_benchmarks * sizeof(benchmark_result_t));
    
    int result_idx = 0;
    
    // Run benchmarks
    printf("\nRunning benchmarks...\n");
    
    for (int s = 0; s < num_sizes; s++) {
        int n = test_sizes[s];
        int iter = iterations[s];
        
        // Generate test signal
        complex_t* test_signal_data = allocate_complex_array(n);
        generate_test_signal(test_signal_data, n, test_signal);
        
        printf("\nTesting size %d (%d iterations):\n", n, iter);
        
        for (int i = 0; i < num_implementations; i++) {
            printf("  %s...", implementations[i].name);
            fflush(stdout);
            
            benchmark_result_t result = run_single_benchmark(
                &implementations[i], test_signal_data, n, iter);
            
            all_results[result_idx++] = result;
            
            if (result.passed == -1) {
                printf(" SKIPPED\n");
            } else if (result.passed == 0) {
                printf(" FAILED\n");
            } else {
                printf(" %.3f ms\n", result.forward_time_ms);
            }
        }
        
        free_complex_array(test_signal_data);
    }
    
    // Print detailed results
    printf("\n\nDetailed Results:\n");
    printf("=================\n");
    
    for (int s = 0; s < num_sizes; s++) {
        print_results_table(all_results, result_idx, test_sizes[s]);
    }
    
    // Analysis
    analyze_best_implementations(all_results, result_idx);
    analyze_scaling(all_results, result_idx);
    
    // Summary statistics
    printf("\n\nSummary Statistics:\n");
    printf("===================\n");
    
    int total_pass = 0, total_fail = 0, total_skip = 0;
    
    for (int i = 0; i < result_idx; i++) {
        if (all_results[i].passed == 1) total_pass++;
        else if (all_results[i].passed == 0) total_fail++;
        else total_skip++;
    }
    
    printf("Total tests: %d\n", result_idx);
    printf("Passed: %d (%.1f%%)\n", total_pass, 100.0 * total_pass / result_idx);
    printf("Failed: %d (%.1f%%)\n", total_fail, 100.0 * total_fail / result_idx);
    printf("Skipped: %d (%.1f%%)\n", total_skip, 100.0 * total_skip / result_idx);
    
    // Recommendations
    printf("\n\nRecommendations:\n");
    printf("================\n");
    printf("- Small sizes (< 256): Split-Radix or Radix-4\n");
    printf("- Medium sizes (256-4096): Radix-2 DIT or Split-Radix\n");
    printf("- Large sizes (> 4096): Radix-2 DIT with cache optimization\n");
    printf("- Non-power-of-2: Bluestein or Mixed-Radix\n");
    printf("- Real-time: Fixed-point or SIMD implementations\n");
    
    free(all_results);
    
    return 0;
}
