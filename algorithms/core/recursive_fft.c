#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Recursive FFT Implementation
 * 
 * Educational implementation showing the divide-and-conquer nature
 * of the FFT algorithm. Not optimized for performance but excellent
 * for understanding the algorithm.
 * 
 * Time Complexity: O(n log n)
 * Space Complexity: O(n log n) due to recursion
 */

// Recursive FFT function
void recursive_fft(complex_t* x, complex_t* result, int n, int stride, fft_direction dir) {
    if (n == 1) {
        result[0] = x[0];
        return;
    }
    
    int half_n = n / 2;
    
    // Allocate temporary arrays for even and odd results
    complex_t* even_result = allocate_complex_array(half_n);
    complex_t* odd_result = allocate_complex_array(half_n);
    
    // Recursively compute FFT of even-indexed elements
    recursive_fft(x, even_result, half_n, stride * 2, dir);
    
    // Recursively compute FFT of odd-indexed elements
    recursive_fft(x + stride, odd_result, half_n, stride * 2, dir);
    
    // Combine results
    for (int k = 0; k < half_n; k++) {
        complex_t t = twiddle_factor(k, n, dir) * odd_result[k];
        result[k] = even_result[k] + t;
        result[k + half_n] = even_result[k] - t;
    }
    
    free_complex_array(even_result);
    free_complex_array(odd_result);
}

// Main recursive FFT wrapper
void fft_recursive(complex_t* x, int n, fft_direction dir) {
    CHECK_POWER_OF_TWO(n);
    
    complex_t* result = allocate_complex_array(n);
    recursive_fft(x, result, n, 1, dir);
    
    // Copy result back to input array
    memcpy(x, result, n * sizeof(complex_t));
    free_complex_array(result);
    
    // Scale for inverse FFT
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

// Wrapper functions
void fft_recursive_forward(complex_t* x, int n) {
    fft_recursive(x, n, FFT_FORWARD);
}

void ifft_recursive(complex_t* x, int n) {
    fft_recursive(x, n, FFT_INVERSE);
}

// Visualize recursion tree
void print_recursion_tree(int n, int level, const char* prefix) {
    if (n == 1) {
        printf("%s└─ Base case (n=1)\n", prefix);
        return;
    }
    
    printf("%s├─ FFT(n=%d)\n", prefix, n);
    
    char new_prefix[256];
    snprintf(new_prefix, sizeof(new_prefix), "%s│  ", prefix);
    
    printf("%s├─ Even subset:\n", prefix);
    print_recursion_tree(n/2, level + 1, new_prefix);
    
    printf("%s└─ Odd subset:\n", prefix);
    snprintf(new_prefix, sizeof(new_prefix), "%s   ", prefix);
    print_recursion_tree(n/2, level + 1, new_prefix);
}

// Count recursion depth and calls
typedef struct {
    int max_depth;
    int total_calls;
    int* calls_per_level;
} recursion_stats_t;

void count_recursion_stats(int n, int depth, recursion_stats_t* stats) {
    stats->total_calls++;
    stats->calls_per_level[depth]++;
    
    if (depth > stats->max_depth) {
        stats->max_depth = depth;
    }
    
    if (n > 1) {
        count_recursion_stats(n/2, depth + 1, stats);
        count_recursion_stats(n/2, depth + 1, stats);
    }
}

// Demonstration and analysis
int main() {
    printf("Recursive FFT Implementation\n");
    printf("============================\n\n");
    
    // Visualize recursion tree
    printf("Recursion Tree Visualization:\n");
    printf("-----------------------------\n");
    printf("FFT with n=8:\n");
    print_recursion_tree(8, 0, "");
    
    printf("\nFFT with n=16:\n");
    print_recursion_tree(16, 0, "");
    
    // Analyze recursion statistics
    printf("\n\nRecursion Statistics:\n");
    printf("--------------------\n");
    printf("Size\tDepth\tTotal Calls\tCalls per Level\n");
    printf("----\t-----\t-----------\t---------------\n");
    
    for (int n = 4; n <= 128; n *= 2) {
        recursion_stats_t stats;
        stats.max_depth = 0;
        stats.total_calls = 0;
        stats.calls_per_level = (int*)calloc(20, sizeof(int));
        
        count_recursion_stats(n, 0, &stats);
        
        printf("%d\t%d\t%d\t\t", n, stats.max_depth, stats.total_calls);
        for (int i = 0; i <= stats.max_depth; i++) {
            printf("%d ", stats.calls_per_level[i]);
        }
        printf("\n");
        
        free(stats.calls_per_level);
    }
    
    // Performance comparison
    printf("\n\nPerformance Comparison:\n");
    printf("----------------------\n");
    printf("Size\tRecursive\tIterative\tRatio\n");
    printf("----\t---------\t---------\t-----\n");
    
    for (int n = 64; n <= 1024; n *= 2) {
        complex_t* x_rec = allocate_complex_array(n);
        complex_t* x_iter = allocate_complex_array(n);
        
        // Generate identical test data
        for (int i = 0; i < n; i++) {
            double val = sin(2 * PI * 5 * i / n);
            x_rec[i] = x_iter[i] = val;
        }
        
        timer_t timer;
        
        // Time recursive FFT
        timer_start(&timer);
        fft_recursive_forward(x_rec, n);
        timer_stop(&timer);
        double time_rec = timer.elapsed_ms;
        
        // Time iterative FFT
        timer_start(&timer);
        radix2_dit_fft(x_iter, n, FFT_FORWARD);
        timer_stop(&timer);
        double time_iter = timer.elapsed_ms;
        
        printf("%d\t%.3f ms\t%.3f ms\t%.2fx\n",
               n, time_rec, time_iter, time_rec / time_iter);
        
        free_complex_array(x_rec);
        free_complex_array(x_iter);
    }
    
    // Demonstrate algorithm correctness
    printf("\n\nCorrectness Test:\n");
    printf("-----------------\n");
    
    int n = 16;
    complex_t* signal = allocate_complex_array(n);
    complex_t* reference = allocate_complex_array(n);
    
    // Generate test signal
    printf("Input signal: ");
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 2 * i / n) + 0.5 * cos(2 * PI * 5 * i / n);
        reference[i] = signal[i];
        if (i < 4) printf("%.3f ", creal(signal[i]));
    }
    printf("...\n");
    
    // Apply recursive FFT
    fft_recursive_forward(signal, n);
    
    // Apply reference FFT
    radix2_dit_fft(reference, n, FFT_FORWARD);
    
    // Compare results
    double max_error = 0;
    for (int i = 0; i < n; i++) {
        double error = cabs(signal[i] - reference[i]);
        if (error > max_error) max_error = error;
    }
    
    printf("Maximum error vs reference: %.2e %s\n",
           max_error, max_error < 1e-10 ? "(PASS)" : "(FAIL)");
    
    // Show magnitude spectrum
    printf("\nMagnitude spectrum: ");
    for (int i = 0; i < n/2; i++) {
        printf("%.1f ", cabs(signal[i]));
    }
    printf("\n");
    
    // Memory usage analysis
    printf("\n\nMemory Usage Analysis:\n");
    printf("---------------------\n");
    printf("Recursive FFT allocates O(n log n) temporary arrays\n");
    printf("For n=%d: ~%d complex numbers allocated\n", 
           n, n * (log2_int(n) + 1));
    printf("Memory overhead: %.1fx compared to in-place iterative\n",
           (double)(log2_int(n) + 1));
    
    free_complex_array(signal);
    free_complex_array(reference);
    
    return 0;
}
