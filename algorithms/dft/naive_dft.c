#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * @file naive_dft.c
 * @brief Direct DFT Implementation (O(n²) Reference)
 * 
 * @details
 * This file implements the Discrete Fourier Transform using the direct
 * mathematical definition. While computationally expensive (O(n²)), it serves
 * as a reference implementation for verifying the correctness of faster algorithms.
 * 
 * Mathematical Definition:
 * Forward DFT: X[k] = Σ(n=0 to N-1) x[n] * e^(-j2πkn/N)
 * Inverse DFT: x[n] = (1/N) * Σ(k=0 to N-1) X[k] * e^(j2πkn/N)
 * 
 * @author FFT Study Repository
 * @date 2024
 * 
 * Time Complexity: O(n²)
 * Space Complexity: O(n)
 * 
 * Applications:
 * - Reference implementation for testing
 * - Small data sets where n < 100
 * - Educational purposes
 * - When simplicity is more important than speed
 */

/**
 * @brief Compute DFT using direct definition
 * 
 * @details
 * This function implements the textbook definition of the DFT.
 * For each output frequency bin k, it computes the sum of all
 * input samples multiplied by the corresponding complex exponentials.
 * 
 * Data Structure Used:
 * - Temporary array to hold input while computing output in-place
 * 
 * Algorithm:
 * 1. Copy input to temporary array
 * 2. For each output bin k:
 *    a. Initialize sum to 0
 *    b. For each input sample j:
 *       - Compute complex exponential e^(dir*2πi*jk/n)
 *       - Multiply with input sample and add to sum
 *    c. Store sum in output[k]
 * 3. If inverse transform, scale by 1/n
 * 
 * @param x Input/output array of complex numbers
 * @param n Length of the array
 * @param dir Transform direction (FFT_FORWARD or FFT_INVERSE)
 */
void naive_dft(complex_t* x, int n, fft_direction dir) {
    /* Input validation */
    if (!x || n <= 0) {
        fprintf(stderr, "Error: Invalid input to naive_dft\n");
        return;
    }
    
    /* Allocate temporary storage */
    complex_t* temp = allocate_complex_array(n);
    if (!temp) {
        fprintf(stderr, "Error: Memory allocation failed in naive_dft\n");
        return;
    }
    
    /* Copy input to temporary array */
    memcpy(temp, x, n * sizeof(complex_t));
    
    /* 
     * Direct DFT computation
     * For each frequency bin k, compute the weighted sum
     */
    for (int k = 0; k < n; k++) {
        x[k] = 0;  /* Initialize accumulator */
        
        /* Sum over all time samples */
        for (int j = 0; j < n; j++) {
            /* Compute phase: e^(dir * 2πi * jk / n) */
            double angle = dir * TWO_PI * j * k / n;
            complex_t twiddle = cexp(I * angle);
            
            /* Accumulate: X[k] += x[j] * e^(phase) */
            x[k] += temp[j] * twiddle;
        }
        
        /* Scale for inverse transform */
        if (dir == FFT_INVERSE) {
            x[k] /= n;
        }
    }
    
    /* Clean up */
    free_complex_array(temp);
}

/**
 * @brief Compute forward DFT
 * @param x Input/output array
 * @param n Array length
 */
void dft_forward(complex_t* x, int n) {
    naive_dft(x, n, FFT_FORWARD);
}

/**
 * @brief Compute inverse DFT
 * @param x Input/output array
 * @param n Array length
 */
void dft_inverse(complex_t* x, int n) {
    naive_dft(x, n, FFT_INVERSE);
}

/**
 * @brief Operation counting structure for complexity analysis
 */
typedef struct {
    long long adds;    /* Number of complex additions */
    long long muls;    /* Number of complex multiplications */
    long long exps;    /* Number of exponential calculations */
} dft_op_count_t;

/**
 * @brief Count operations required for naive DFT
 * 
 * @param n Transform size
 * @return Structure containing operation counts
 */
dft_op_count_t count_dft_operations(int n) {
    dft_op_count_t count;
    count.adds = (long long)n * (n - 1);  /* n-1 additions per output, n outputs */
    count.muls = (long long)n * n;        /* n multiplications per output */
    count.exps = (long long)n * n;        /* n exponentials per output */
    return count;
}

#ifndef LIB_BUILD

/**
 * @brief Analyze DFT computational complexity
 */
void analyze_dft_performance() {
    printf("\nDFT Computational Complexity Analysis:\n");
    printf("=====================================\n");
    printf("N\tAdditions\tMultiplies\tExponentials\tTotal Ops\n");
    printf("---\t---------\t----------\t------------\t---------\n");
    
    for (int n = 2; n <= 1024; n *= 2) {
        dft_op_count_t ops = count_dft_operations(n);
        long long total = ops.adds + ops.muls + ops.exps;
        printf("%d\t%lld\t\t%lld\t\t%lld\t\t%lld\n", 
               n, ops.adds, ops.muls, ops.exps, total);
    }
    
    printf("\nComparison with FFT (Radix-2):\n");
    printf("N\tDFT Ops\t\tFFT Ops\t\tRatio\n");
    printf("---\t-------\t\t-------\t\t-----\n");
    
    for (int n = 16; n <= 1024; n *= 2) {
        long long dft_ops = 3LL * n * n;  /* Simplified total */
        long long fft_ops = 5LL * n * log2_int(n);  /* Approximate */
        printf("%d\t%lld\t\t%lld\t\t%.1f\n", 
               n, dft_ops, fft_ops, (double)dft_ops / fft_ops);
    }
}

/**
 * @brief Test fundamental DFT properties
 * 
 * Tests:
 * 1. Linearity: DFT(ax + by) = a*DFT(x) + b*DFT(y)
 * 2. Parseval's theorem: Energy conservation
 * 3. Circular shift property
 */
void test_dft_properties() {
    printf("\n\nDFT Properties Verification:\n");
    printf("============================\n");
    
    int n = 8;
    complex_t* x = allocate_complex_array(n);
    complex_t* X = allocate_complex_array(n);
    
    /* Test 1: Linearity */
    printf("\n1. Linearity Test:\n");
    complex_t* a = allocate_complex_array(n);
    complex_t* b = allocate_complex_array(n);
    complex_t* sum = allocate_complex_array(n);
    
    /* Generate test signals */
    for (int i = 0; i < n; i++) {
        a[i] = sin(2 * PI * i / n);
        b[i] = cos(2 * PI * i / n);
        sum[i] = 2 * a[i] + 3 * b[i];
    }
    
    /* DFT of linear combination */
    memcpy(x, sum, n * sizeof(complex_t));
    dft_forward(x, n);
    
    /* Linear combination of DFTs */
    dft_forward(a, n);
    dft_forward(b, n);
    
    /* Verify linearity */
    double error = 0;
    for (int i = 0; i < n; i++) {
        complex_t expected = 2 * a[i] + 3 * b[i];
        error += cabs(x[i] - expected);
    }
    printf("Linearity error: %.2e %s\n", error, 
           error < 1e-10 ? "(PASS)" : "(FAIL)");
    
    /* Test 2: Parseval's theorem */
    printf("\n2. Parseval's Theorem Test:\n");
    
    /* Generate random signal */
    for (int i = 0; i < n; i++) {
        x[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 + 
               I * ((double)rand() / RAND_MAX) * 2.0 - 1.0;
    }
    memcpy(X, x, n * sizeof(complex_t));
    
    /* Time domain energy */
    double time_energy = 0;
    for (int i = 0; i < n; i++) {
        time_energy += cabs(x[i]) * cabs(x[i]);
    }
    
    /* Frequency domain energy */
    dft_forward(X, n);
    double freq_energy = 0;
    for (int i = 0; i < n; i++) {
        freq_energy += cabs(X[i]) * cabs(X[i]);
    }
    freq_energy /= n;  /* Account for DFT scaling */
    
    printf("Time domain energy: %.6f\n", time_energy);
    printf("Freq domain energy: %.6f\n", freq_energy);
    printf("Relative error: %.2e %s\n", 
           fabs(time_energy - freq_energy) / time_energy,
           fabs(time_energy - freq_energy) / time_energy < 1e-10 ? "(PASS)" : "(FAIL)");
    
    /* Test 3: Circular shift property */
    printf("\n3. Circular Shift Property Test:\n");
    
    /* Generate square wave */
    generate_square_wave(x, n, 2, n);
    memcpy(X, x, n * sizeof(complex_t));
    dft_forward(X, n);
    
    /* Create circularly shifted version */
    complex_t* x_shifted = allocate_complex_array(n);
    int shift = 2;
    for (int i = 0; i < n; i++) {
        x_shifted[i] = x[(i - shift + n) % n];
    }
    
    /* DFT of shifted signal */
    complex_t* X_shifted = allocate_complex_array(n);
    memcpy(X_shifted, x_shifted, n * sizeof(complex_t));
    dft_forward(X_shifted, n);
    
    /* Verify phase shift property */
    error = 0;
    for (int k = 0; k < n; k++) {
        complex_t expected = X[k] * cexp(-I * TWO_PI * shift * k / n);
        error += cabs(X_shifted[k] - expected);
    }
    printf("Shift property error: %.2e %s\n", error,
           error < 1e-10 ? "(PASS)" : "(FAIL)");
    
    /* Clean up */
    free_complex_array(a);
    free_complex_array(b);
    free_complex_array(sum);
    free_complex_array(x);
    free_complex_array(X);
    free_complex_array(x_shifted);
    free_complex_array(X_shifted);
}

/**
 * @brief Main demonstration program
 */
int main() {
    printf("Naive DFT Implementation\n");
    printf("========================\n");
    
    /* Basic functionality test */
    int n = 16;
    complex_t* signal = allocate_complex_array(n);
    complex_t* original = allocate_complex_array(n);
    
    /* Generate test signal with known frequency content */
    printf("\nTest signal: sum of 3 sinusoids\n");
    printf("Frequencies: 2, 5, and 7 cycles per window\n");
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 2 * i / n) + 
                   0.5 * sin(2 * PI * 5 * i / n) + 
                   0.25 * sin(2 * PI * 7 * i / n);
        original[i] = signal[i];
    }
    
    /* Forward DFT */
    fft_timer_t timer;
    timer_start(&timer);
    dft_forward(signal, n);
    timer_stop(&timer);
    
    printf("\nForward DFT completed in %.3f ms\n", timer.elapsed_ms);
    
    /* Display magnitude spectrum */
    printf("\nMagnitude spectrum:\n");
    printf("Bin\tFreq\tMagnitude\tDetection\n");
    printf("---\t----\t---------\t---------\n");
    for (int k = 0; k < n/2; k++) {
        double mag = cabs(signal[k]);
        printf("%2d\t%d\t%.3f\t\t", k, k, mag);
        if (mag > 0.1) {
            printf("* Peak detected");
        }
        printf("\n");
    }
    
    /* Inverse DFT */
    timer_start(&timer);
    dft_inverse(signal, n);
    timer_stop(&timer);
    
    printf("\nInverse DFT completed in %.3f ms\n", timer.elapsed_ms);
    
    /* Verify perfect reconstruction */
    double error = 0;
    for (int i = 0; i < n; i++) {
        error += cabs(signal[i] - original[i]);
    }
    printf("Reconstruction error: %.2e %s\n", error / n,
           error / n < 1e-14 ? "(Perfect)" : "(Good)");
    
    /* Performance analysis */
    analyze_dft_performance();
    
    /* Property verification */
    test_dft_properties();
    
    /* Execution time scaling */
    printf("\n\nExecution Time Scaling:\n");
    printf("======================\n");
    printf("N\tTime (ms)\tTime/N²\t\tScaling\n");
    printf("---\t---------\t--------\t-------\n");
    
    double prev_time = 0;
    for (int size = 8; size <= 256; size *= 2) {
        complex_t* test = allocate_complex_array(size);
        
        /* Random data */
        for (int i = 0; i < size; i++) {
            test[i] = ((double)rand() / RAND_MAX) * 2.0 - 1.0 +
                     I * ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
        
        timer_start(&timer);
        dft_forward(test, size);
        timer_stop(&timer);
        
        double time_per_n2 = timer.elapsed_ms / (size * size);
        printf("%d\t%.3f\t\t%.6f", size, timer.elapsed_ms, time_per_n2);
        
        if (prev_time > 0) {
            printf("\t%.1fx", timer.elapsed_ms / prev_time);
        }
        printf("\n");
        
        prev_time = timer.elapsed_ms;
        free_complex_array(test);
    }
    
    /* Summary */
    printf("\n\nAlgorithm Characteristics:\n");
    printf("=========================\n");
    printf("✓ Works for any input size\n");
    printf("✓ Numerically stable\n");
    printf("✓ Simple to understand and implement\n");
    printf("✗ O(n²) complexity makes it impractical for large n\n");
    printf("✓ Excellent for verification and small datasets\n");
    
    /* Clean up */
    free_complex_array(signal);
    free_complex_array(original);
    
    return 0;
}

#endif /* LIB_BUILD */
