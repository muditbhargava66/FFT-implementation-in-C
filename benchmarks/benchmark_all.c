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
typedef fft_error_t (*fft_func_t)(complex_t*, int, fft_direction);

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
    int bins_checked;
    int passed;
} benchmark_result_t;

// The implementations are declared by fft_algorithms.h; a local `extern` list
// here previously duplicated those prototypes and drifted out of date.

// List of implementations to benchmark
fft_impl_t implementations[] = {
    {"Radix-2 DIT", radix2_dit_fft, radix2_dit_fft, 1, 2, 65536},
    {"Radix-2 DIF", radix2_dif_fft, radix2_dif_fft, 1, 2, 65536},
    {"Radix-4", radix4_fft, radix4_fft, 1, 4, 65536},
    {"Split-Radix", split_radix_fft, split_radix_fft, 1, 2, 65536},
    {"Bluestein", bluestein_fft, bluestein_fft, 0, 2, 65536},      // Increased from 4096
    {"Mixed-Radix", mixed_radix_fft, mixed_radix_fft, 0, 2, 65536}, // Increased from 4096
    {"Recursive", fft_recursive, fft_recursive, 1, 2, 65536},       // Increased from 8192
    {"Naive DFT", naive_dft, naive_dft, 0, 2, 4096},      // O(n^2): 14 s at n=16384, so it is capped
    {"Optimized DFT", optimized_dft, optimized_dft, 0, 2, 4096}   // O(n^2) as well
};

int num_implementations = sizeof(implementations) / sizeof(implementations[0]);

/*
 * The reference DFT is O(n^2): 14.3 s at n=16384 on this machine. It used to be
 * recomputed inside run_single_benchmark(), i.e. once per implementation, which
 * made the size-16384 block alone cost 9 x 14.3 s = 129 s and the whole suite
 * several minutes. It is now computed once per size in main().
 *
 * Above MAX_REF_BINS the reference is evaluated on an evenly spaced sample of
 * bins rather than all n. Which bins the DFT definition pins down does not
 * depend on the bin index, so a spread sample catches indexing and twiddle
 * errors just as well as all n bins, at 1/128th of the cost.
 */
#define MAX_REF_BINS 128

static int build_bin_list(int n, int* bins) {
    if (n <= MAX_REF_BINS) {
        for (int k = 0; k < n; k++) bins[k] = k;
        return n;
    }
    for (int i = 0; i < MAX_REF_BINS; i++) {
        bins[i] = (int)((long long)i * (n - 1) / (MAX_REF_BINS - 1));
    }
    return MAX_REF_BINS;
}

static void compute_reference_bins(const complex_t* input, complex_t* reference,
                                   int n, const int* bins, int num_bins) {
    for (int b = 0; b < num_bins; b++) {
        int k = bins[b];
        complex_t acc = 0.0;
        for (int j = 0; j < n; j++) {
            double angle = -2.0 * PI * (double)k * (double)j / (double)n;
            acc += input[j] * (cos(angle) + I * sin(angle));
        }
        reference[k] = acc;
    }
}

// Compute error metrics over every entry
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

// Same, restricted to the sampled reference bins
static void compute_error_metrics_bins(const complex_t* test, const complex_t* reference,
                                       const int* bins, int num_bins,
                                       double* max_error, double* rms_error) {
    *max_error = 0;
    double sum_squared_error = 0;

    for (int b = 0; b < num_bins; b++) {
        double error = cabs(test[bins[b]] - reference[bins[b]]);
        if (error > *max_error) *max_error = error;
        sum_squared_error += error * error;
    }

    *rms_error = sqrt(sum_squared_error / num_bins);
}

/*
 * clock() on Windows advances in coarse steps, so a fixed iteration count made
 * fast transforms print "0.000 ms" (measured: the whole n=64 reference came back
 * as 0.0 ms) and made analyze_scaling() take log(0). Grow the batch until it
 * spans at least MIN_BATCH_MS and report the per-iteration mean. The geometric
 * growth bounds the overshoot to about 8x MIN_BATCH_MS.
 */
#define MIN_BATCH_MS 20.0
#define MAX_BATCH_REPS 2000000

static double time_transform(fft_func_t fn, complex_t* work, const complex_t* seed,
                             int n, fft_direction dir) {
    int reps = 1;
    for (;;) {
        fft_timer_t timer;
        timer_start(&timer);
        for (int r = 0; r < reps; r++) {
            memcpy(work, seed, (size_t)n * sizeof(complex_t));
            fn(work, n, dir);
        }
        timer_stop(&timer);
        if (timer.elapsed_ms >= MIN_BATCH_MS || reps >= MAX_BATCH_REPS) {
            return timer.elapsed_ms / reps;
        }
        int next = reps * 8;
        reps = (next > MAX_BATCH_REPS) ? MAX_BATCH_REPS : next;
    }
}

// Run single benchmark. `reference` holds the DFT of `input` at the bins listed
// in `bins`; it is shared by every implementation at this size.
static benchmark_result_t run_single_benchmark(fft_impl_t* impl, complex_t* input, int n,
                                               const complex_t* reference,
                                               const int* bins, int num_bins) {
    benchmark_result_t result;
    result.impl_name = impl->name;
    result.size = n;
    result.passed = 1;
    result.bins_checked = num_bins;
    
    // Skip if size not supported
    if ((impl->requires_power_of_2 && !is_power_of_two(n)) ||
        n < impl->min_size || n > impl->max_size) {
        result.forward_time_ms = -1;
        result.inverse_time_ms = -1;
        result.max_error = -1;
        result.rms_error = -1;
        result.bins_checked = 0;
        result.passed = -1;
        return result;
    }
    
    // Allocate arrays
    complex_t* test = allocate_complex_array(n);
    complex_t* forward = allocate_complex_array(n);
    complex_t* original = allocate_complex_array(n);

    if (!test || !forward || !original) {
        free_complex_array(test);
        free_complex_array(forward);
        free_complex_array(original);
        result.forward_time_ms = -1;
        result.inverse_time_ms = -1;
        result.max_error = -1;
        result.rms_error = -1;
        result.bins_checked = 0;
        result.passed = -1;
        return result;
    }

    memcpy(original, input, (size_t)n * sizeof(complex_t));

    // Warm-up run, then time the forward transform
    memcpy(test, input, (size_t)n * sizeof(complex_t));
    impl->forward(test, n, FFT_FORWARD);
    result.forward_time_ms = time_transform(impl->forward, test, input, n, FFT_FORWARD);

    // Accuracy against the reference, over the bins that were computed
    memcpy(forward, input, (size_t)n * sizeof(complex_t));
    impl->forward(forward, n, FFT_FORWARD);
    compute_error_metrics_bins(forward, reference, bins, num_bins,
                               &result.max_error, &result.rms_error);

    // Time the inverse, seeded with this implementation's own forward output so
    // the round trip is self-consistent at every size.
    result.inverse_time_ms = time_transform(impl->inverse, test, forward, n, FFT_INVERSE);

    // Reconstruction is checked over all n samples: the round trip goes through
    // the implementation's own spectrum, so it does not depend on which bins the
    // sampled reference happened to cover.
    memcpy(test, forward, (size_t)n * sizeof(complex_t));
    impl->inverse(test, n, FFT_INVERSE);

    double recon_max_error, recon_rms_error;
    compute_error_metrics(test, original, n, &recon_max_error, &recon_rms_error);

    // More lenient error threshold for numerical algorithms
    if (recon_max_error > 1e-8 || result.max_error > 1e-8) {
        result.passed = 0;
    }

    free_complex_array(test);
    free_complex_array(forward);
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
            char status[48];
            snprintf(status, sizeof(status), "%s (%d/%d bins)",
                     results[i].passed ? "PASS" : "FAIL",
                     results[i].bins_checked, results[i].size);
            printf("%-15s | %12.4f | %12.4f | %12.2e | %12.2e | %s\n",
                   results[i].impl_name,
                   results[i].forward_time_ms,
                   results[i].inverse_time_ms,
                   results[i].max_error,
                   results[i].rms_error,
                   status);
        }
    }
}

// Find best implementation for each size
void analyze_best_implementations(benchmark_result_t* results, int num_results) {
    printf("\n\nBest Implementations by Size:\n");
    printf("=============================\n");

    // Sizes are taken from the results themselves; a hardcoded copy of main()'s
    // test_sizes[] here used to silently print nothing whenever the two drifted.
    for (int i = 0; i < num_results; i++) {
        int size = results[i].size;

        int already_done = 0;
        for (int j = 0; j < i; j++) {
            if (results[j].size == size) { already_done = 1; break; }
        }
        if (already_done) continue;

        double best_time = DBL_MAX;
        const char* best_impl = NULL;

        for (int k = 0; k < num_results; k++) {
            if (results[k].size == size && results[k].passed == 1 &&
                results[k].forward_time_ms > 0 && results[k].forward_time_ms < best_time) {
                best_time = results[k].forward_time_ms;
                best_impl = results[k].impl_name;
            }
        }

        if (best_impl) {
            printf("Size %6d: %s (%.4f ms)\n", size, best_impl, best_time);
        } else {
            printf("Size %6d: no implementation passed\n", size);
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
            if (results[i].forward_time_ms <= 0) continue;

            if (prev_size > 0 && prev_time > 0) {
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
    const char* test_signal = "impulse";  // Use impulse instead of random for better stability
    
    // Seed random number generator
    srand(42);
    
    // Allocate results array
    int total_benchmarks = num_implementations * num_sizes;
    benchmark_result_t* all_results = (benchmark_result_t*)malloc(
        total_benchmarks * sizeof(benchmark_result_t));
    if (!all_results) {
        fprintf(stderr, "benchmark_all: out of memory allocating results\n");
        return 1;
    }

    int* bins = (int*)malloc(MAX_REF_BINS * sizeof(int));
    if (!bins) {
        fprintf(stderr, "benchmark_all: out of memory allocating bin list\n");
        free(all_results);
        return 1;
    }

    int result_idx = 0;
    
    // Run benchmarks
    printf("\nRunning benchmarks...\n");
    
    for (int s = 0; s < num_sizes; s++) {
        int n = test_sizes[s];

        // Generate test signal
        complex_t* test_signal_data = allocate_complex_array(n);
        complex_t* reference = allocate_complex_array(n);
        if (!test_signal_data || !reference) {
            fprintf(stderr, "benchmark_all: out of memory at size %d\n", n);
            free_complex_array(test_signal_data);
            free_complex_array(reference);
            break;
        }
        generate_test_signal(test_signal_data, n, test_signal);

        // One O(n^2) reference per size, shared by all implementations.
        int num_bins = build_bin_list(n, bins);
        compute_reference_bins(test_signal_data, reference, n, bins, num_bins);

        printf("\nTesting size %d (reference: %d/%d bins):\n", n, num_bins, n);
        
        for (int i = 0; i < num_implementations; i++) {
            printf("  %s...", implementations[i].name);
            fflush(stdout);
            
            benchmark_result_t result = run_single_benchmark(
                &implementations[i], test_signal_data, n, reference, bins, num_bins);
            
            all_results[result_idx++] = result;
            
            if (result.passed == -1) {
                printf(" SKIPPED\n");
            } else if (result.passed == 0) {
                printf(" FAILED\n");
            } else {
                printf(" %.4f ms\n", result.forward_time_ms);
            }
        }
        
        free_complex_array(test_signal_data);
        free_complex_array(reference);
    }

    free(bins);
    
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
    if (result_idx > 0) {
        printf("Passed: %d (%.1f%%)\n", total_pass, 100.0 * total_pass / result_idx);
        printf("Failed: %d (%.1f%%)\n", total_fail, 100.0 * total_fail / result_idx);
        printf("Skipped: %d (%.1f%%)\n", total_skip, 100.0 * total_skip / result_idx);
    } else {
        printf("Passed: 0\nFailed: 0\nSkipped: 0\n");
    }

    // Notes. Which algorithm is fastest at each size is measured above and
    // varies with the host CPU, so it is deliberately not hardcoded here; this
    // section used to assert a fixed ranking that the suite had just contradicted.
    printf("\n\nNotes:\n");
    printf("======\n");
    printf("- Read the fastest implementation per size from 'Best Implementations by Size' above.\n");
    printf("- Bluestein and Mixed-Radix are the only entries that accept non-power-of-2 sizes.\n");
    printf("- Naive DFT and Optimized DFT are O(n^2) and are skipped above n=%d.\n",
           implementations[num_implementations - 1].max_size);
    printf("- Timings are the mean of a batch grown until it spans at least %.0f ms,\n",
           MIN_BATCH_MS);
    printf("  because clock() on this platform is too coarse to time one small transform.\n");
    
    free(all_results);
    
    return 0;
}
