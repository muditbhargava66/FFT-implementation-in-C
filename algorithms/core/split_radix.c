#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Split-Radix FFT Implementation
 * 
 * Combines radix-2 and radix-4 decompositions to minimize the number
 * of arithmetic operations. Uses radix-2 for even indices and radix-4
 * for odd indices.
 * 
 * Benefits:
 * - Lowest operation count among power-of-2 FFT algorithms
 * - About 33% fewer multiplications than radix-2
 * - More complex implementation but better performance
 * 
 * Time Complexity: O(n log n) with the best constants
 */

// Forward declaration for recursion
void split_radix_fft_recursive(complex_t* x, int n, int stride, fft_direction dir);

// Split-radix butterfly operations
static void split_radix_butterflies(complex_t* x, int n, int stride, fft_direction dir) {
    if (n <= 2) {
        if (n == 2) {
            complex_t t = x[0];
            x[0] = t + x[stride];
            x[stride] = t - x[stride];
        }
        return;
    }
    
    int n2 = n / 2;
    int n4 = n / 4;
    
    // Recursive calls
    split_radix_fft_recursive(x, n2, stride * 2, dir);              // Even indices
    split_radix_fft_recursive(x + stride, n4, stride * 4, dir);     // Indices ≡ 1 (mod 4)
    split_radix_fft_recursive(x + 3 * stride, n4, stride * 4, dir); // Indices ≡ 3 (mod 4)
    
    // Combine results with twiddle factors
    for (int k = 0; k < n4; k++) {
        int k1 = 2 * k * stride;
        int k2 = k1 + stride;
        int k3 = k2 + 2 * stride;
        int k4 = k3 + stride;
        
        complex_t w1 = twiddle_factor(k, n, dir);
        complex_t w3 = twiddle_factor(3 * k, n, dir);
        
        complex_t t1 = x[k2] * w1;
        complex_t t2 = x[k4] * w3;
        
        complex_t u1 = t1 + t2;
        complex_t u2 = t1 - t2;
        
        if (dir == FFT_FORWARD) {
            u2 *= -I;
        } else {
            u2 *= I;
        }
        
        t1 = x[k1];
        x[k1] = t1 + u1;
        x[k2] = t1 - u1;
        
        t2 = x[k3];
        x[k3] = t2 + u2;
        x[k4] = t2 - u2;
    }
}

void split_radix_fft_recursive(complex_t* x, int n, int stride, fft_direction dir) {
    if (n <= 1) return;
    split_radix_butterflies(x, n, stride, dir);
}

// Main split-radix FFT function
void split_radix_fft(complex_t* x, int n, fft_direction dir) {
    CHECK_POWER_OF_TWO(n);
    
    // Create working array for bit-reversed data
    complex_t* work = allocate_complex_array(n);
    
    // Bit-reversal permutation
    int log2n = log2_int(n);
    for (int i = 0; i < n; i++) {
        int j = bit_reverse(i, log2n);
        work[j] = x[i];
    }
    
    // Copy back to original array
    memcpy(x, work, n * sizeof(complex_t));
    free_complex_array(work);
    
    // Perform split-radix FFT
    split_radix_fft_recursive(x, n, 1, dir);
    
    // Scale for inverse FFT
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

// Wrapper functions
void fft_split_radix(complex_t* x, int n) {
    split_radix_fft(x, n, FFT_FORWARD);
}

void ifft_split_radix(complex_t* x, int n) {
    split_radix_fft(x, n, FFT_INVERSE);
}

// Operation counting for comparison
typedef struct {
    int adds;
    int muls;
} op_count_t;

op_count_t count_split_radix_ops(int n) {
    op_count_t count = {0, 0};
    if (n <= 1) return count;
    
    // Based on theoretical analysis
    int log2n = log2_int(n);
    count.muls = (n * log2n - 3 * n + 4) / 3;
    count.adds = n * log2n;
    
    return count;
}

op_count_t count_radix2_ops(int n) {
    op_count_t count = {0, 0};
    int log2n = log2_int(n);
    count.muls = n * log2n / 2;
    count.adds = n * log2n;
    return count;
}

// Performance demonstration
int main() {
    printf("Split-Radix FFT Implementation\n");
    printf("==============================\n\n");
    
    // Operation count comparison
    printf("Theoretical Operation Counts:\n");
    printf("Size\tRadix-2 Muls\tSplit-Radix Muls\tReduction\n");
    printf("----\t------------\t----------------\t---------\n");
    
    for (int n = 8; n <= 1024; n *= 2) {
        op_count_t r2_ops = count_radix2_ops(n);
        op_count_t sr_ops = count_split_radix_ops(n);
        double reduction = (1.0 - (double)sr_ops.muls / r2_ops.muls) * 100;
        
        printf("%d\t%d\t\t%d\t\t\t%.1f%%\n", 
               n, r2_ops.muls, sr_ops.muls, reduction);
    }
    
    // Performance comparison
    printf("\n\nPerformance Comparison:\n");
    printf("Size\tRadix-2\t\tRadix-4\t\tSplit-Radix\n");
    printf("----\t-------\t\t-------\t\t-----------\n");
    
    int sizes[] = {128, 512, 2048, 8192};
    int num_sizes = sizeof(sizes) / sizeof(sizes[0]);
    
    for (int s = 0; s < num_sizes; s++) {
        int n = sizes[s];
        complex_t* x_r2 = allocate_complex_array(n);
        complex_t* x_r4 = allocate_complex_array(n);
        complex_t* x_sr = allocate_complex_array(n);
        
        // Generate identical random data
        for (int i = 0; i < n; i++) {
            double real = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double imag = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            x_r2[i] = x_r4[i] = x_sr[i] = real + I * imag;
        }
        
        timer_t timer;
        double time_r2, time_r4, time_sr;
        
        // Time radix-2
        timer_start(&timer);
        radix2_dit_fft(x_r2, n, FFT_FORWARD);
        timer_stop(&timer);
        time_r2 = timer.elapsed_ms;
        
        // Time radix-4
        timer_start(&timer);
        fft_radix4(x_r4, n);
        timer_stop(&timer);
        time_r4 = timer.elapsed_ms;
        
        // Time split-radix
        timer_start(&timer);
        fft_split_radix(x_sr, n);
        timer_stop(&timer);
        time_sr = timer.elapsed_ms;
        
        printf("%d\t%.3f ms\t%.3f ms\t%.3f ms\n",
               n, time_r2, time_r4, time_sr);
        
        // Verify correctness
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = cabs(x_r2[i] - x_sr[i]);
            if (diff > max_diff) max_diff = diff;
        }
        
        if (max_diff > 1e-10) {
            printf("  WARNING: Results differ by %.2e\n", max_diff);
        }
        
        free_complex_array(x_r2);
        free_complex_array(x_r4);
        free_complex_array(x_sr);
    }
    
    // Accuracy test
    printf("\n\nAccuracy Test:\n");
    int n = 32;
    complex_t* signal = allocate_complex_array(n);
    complex_t* original = allocate_complex_array(n);
    
    // Generate test signal
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 3 * i / n) + 0.5 * cos(2 * PI * 7 * i / n);
        original[i] = signal[i];
    }
    
    // Forward and inverse transform
    fft_split_radix(signal, n);
    ifft_split_radix(signal, n);
    
    // Check reconstruction
    double error = 0.0;
    for (int i = 0; i < n; i++) {
        error += cabs(signal[i] - original[i]);
    }
    
    printf("Reconstruction error: %.2e %s\n", error / n,
           error / n < 1e-14 ? "(PASS)" : "(FAIL)");
    
    free_complex_array(signal);
    free_complex_array(original);
    
    return 0;
}
