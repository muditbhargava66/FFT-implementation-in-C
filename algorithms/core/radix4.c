#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * @file radix4.c
 * @brief Radix-4 FFT Implementation
 * 
 * @details
 * The Radix-4 FFT algorithm processes 4 points at a time instead of 2,
 * reducing the number of stages. log2(n) may be odd, in which case one
 * leading radix-2 stage handles the leftover factor of 2, so any power-of-two
 * length is accepted (not only powers of four).
 * 
 * Mathematical Foundation:
 * The DFT is decomposed into 4 interleaved subsequences:
 *   X[k] = Σ(r=0 to 3) W_N^(rk) * DFT_{N/4}[x[4n+r]]
 * 
 * The radix-4 butterfly combines 4 (pre-twiddled) inputs using:
 *   [A]   [1  1  1  1] [a]
 *   [B] = [1 -j -1  j] [b]
 *   [C]   [1 -1  1 -1] [c]
 *   [D]   [1  j -1 -j] [d]
 * 
 * @author FFT Study Repository
 * @date 2024
 * 
 * Time Complexity: O(n log n). The radix-4 butterfly itself needs 25% fewer
 * twiddle multiplications than radix-2, but this implementation still recomputes
 * the digit-reversal permutation with an O(n) auxiliary buffer per call, so the
 * realised advantage over radix-2 is smaller than the butterfly arithmetic alone
 * implies.
 * Space Complexity: O(n) auxiliary (digit-reversal scratch), not in-place.
 */

/**
 * @brief Radix-4 butterfly operation (reference illustration only)
 * 
 * @details
 * This helper exists to show the 4-point DFT matrix and its output twiddles.
 * It is NOT used by radix4_fft(), which pre-multiplies its inputs rather than
 * post-multiplying its outputs.
 * 
 * @param a,b,c,d Pointers to the four complex values
 * @param w1,w2,w3 Twiddle factors W^k, W^2k, W^3k
 */
static inline void radix4_butterfly_demo(complex_t* a, complex_t* b, complex_t* c, complex_t* d,
                                         complex_t w1, complex_t w2, complex_t w3) {
    /* Store original values */
    complex_t a0 = *a, b0 = *b, c0 = *c, d0 = *d;
    
    /* Standard radix-4 butterfly computation */
    complex_t t0 = a0 + c0;
    complex_t t1 = a0 - c0;
    complex_t t2 = b0 + d0;
    complex_t t3 = (b0 - d0) * (-I);  /* Multiply by -j */
    
    /* Row B of the matrix above is [1 -j -1 j] = t1 + t3, row D is t1 - t3. */
    *a = t0 + t2;                    /* No twiddle factor */
    *b = (t1 + t3) * w1;            /* Apply W^k */
    *c = (t0 - t2) * w2;            /* Apply W^2k */
    *d = (t1 - t3) * w3;            /* Apply W^3k */
}

/**
 * @brief Main Radix-4 FFT implementation
 * 
 * @details
 * Implements the complete radix-4 FFT algorithm with the following steps:
 * 1. Validate input size is power of 2
 * 2. Perform bit reversal permutation
 * 3. Execute stages of radix-4 butterflies where possible
 * 4. Handle remaining radix-2 stage if needed
 * 5. Scale for inverse transform
 * 
 * @param x Input/output array of complex numbers
 * @param n Array length (must be power of 2)
 * @param dir Transform direction (FFT_FORWARD or FFT_INVERSE)
 */
fft_error_t radix4_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    if (!is_power_of_two(n)) return FFT_ERR_SIZE;

    if (n <= 1) return FFT_OK;

    int log2n = log2_int(n);
    
    // Custom permutation for Radix-4 (and Mixed Radix 2x4x4...)
    // We need to reverse the digits.
    // If N = 4^k, digits are base 4.
    // If N = 2*4^k, first digit is base 2, others base 4.
    
    int is_power_of_4 = (log2n % 2 == 0);
    
    // We'll use a temporary array for the permutation to avoid complex in-place logic
    // (In-place digit reversal is hard, O(N) auxiliary space is acceptable here)
    complex_t* temp = (complex_t*)malloc((size_t)n * sizeof(complex_t));
    if (!temp) return FFT_ERR_ALLOC;
    
    for (int i = 0; i < n; i++) {
        int weight = n;
        int temp_i = i;
        int accumulated_rev = 0;
        
        // If odd power of 2, first factor is 2
        if (!is_power_of_4) {
             int radix = 2;
             weight /= radix;
             int digit = temp_i % radix;
             temp_i /= radix;
             accumulated_rev += digit * weight;
        }
        
        // Remaining factors are 4
        while (weight > 1) {
            int radix = 4;
            weight /= radix;
            int digit = temp_i % radix;
            temp_i /= radix;
            accumulated_rev += digit * weight;
        }
        
        temp[i] = x[accumulated_rev];
    }
    
    // Copy back
    for (int i = 0; i < n; i++) {
        x[i] = temp[i];
    }
    free(temp);
    
    // Determine start length for Radix-4 stages
    int start_length = 4;
    
    // If not power of 4 (odd power of 2), perform one radix-2 stage first
    if (!is_power_of_4) {
        int m = 2;
        // Twiddle for size 2 is -1 (for forward), but w=1 always for first element
        
        for (int k = 0; k < n; k += m) {
            complex_t t = x[k + 1];
            complex_t u = x[k];
            x[k] = u + t;
            x[k + 1] = u - t;
        }
        // After Radix-2 stage, we have DFTs of size 2.
        // Next stage should combine 4 DFTs of size 2 -> size 8.
        start_length = 8;
    }
    
    // Remaining stages are Radix-4
    for (int m = start_length; m <= n; m *= 4) {
        int m4 = m / 4;          /* Quarter of current DFT size */
        int m2 = 2 * m4;         /* Half of current DFT size */
        int m3 = 3 * m4;         /* Three-quarters */

        /* Advance the three twiddles multiplicatively within a block instead of
         * calling cexp() per butterfly. */
        complex_t s1 = twiddle_factor(1, m, dir);
        complex_t s2 = twiddle_factor(2, m, dir);
        complex_t s3 = twiddle_factor(3, m, dir);

        for (int k = 0; k < n; k += m) {
            complex_t w1 = 1.0, w2 = 1.0, w3 = 1.0;
            for (int j = 0; j < m4; j++) {
                // Indices
                int i0 = k + j;
                int i1 = i0 + m4;
                int i2 = i0 + m2;
                int i3 = i0 + m3;
                
                // Read inputs, applying twiddles to the upper three branches
                complex_t x0 = x[i0];
                complex_t x1 = x[i1] * w1;
                complex_t x2 = x[i2] * w2;
                complex_t x3 = x[i3] * w3;
                
                complex_t t0 = x0 + x2;
                complex_t t1 = x0 - x2;
                complex_t t2 = x1 + x3;
                complex_t t3;
                
                if (dir == FFT_FORWARD) {
                    t3 = (x1 - x3) * (-I);
                } else {
                    t3 = (x1 - x3) * (I);
                }
                
                x[i0] = t0 + t2;
                x[i1] = t1 + t3;
                x[i2] = t0 - t2;
                x[i3] = t1 - t3;

                w1 *= s1;
                w2 *= s2;
                w3 *= s3;
            }
        }
    }
    
    // Scale for inverse transform
    if (dir == FFT_INVERSE) {
        double scale = 1.0 / n;
        for (int i = 0; i < n; i++) {
            x[i] *= scale;
        }
    }
    return FFT_OK;
}

/**
 * @brief Compute forward FFT using Radix-4
 * @param x Input/output array
 * @param n Array length (should be power of 4 for best performance)
 */
fft_error_t fft_radix4(complex_t* x, int n) {
    return radix4_fft(x, n, FFT_FORWARD);
}

/**
 * @brief Compute inverse FFT using Radix-4
 * @param x Input/output array
 * @param n Array length
 */
fft_error_t ifft_radix4(complex_t* x, int n) {
    return radix4_fft(x, n, FFT_INVERSE);
}

#ifndef LIB_BUILD

/**
 * @brief Demonstrate radix-4 butterfly operation
 */
static void demonstrate_butterfly() {
    printf("\nRadix-4 Butterfly Demonstration:\n");
    printf("================================\n");
    
    /* Create test values */
    complex_t a = 1.0 + 0.0 * I;
    complex_t b = 0.0 + 1.0 * I;
    complex_t c = -1.0 + 0.0 * I;
    complex_t d = 0.0 - 1.0 * I;
    
    printf("Input values:\n");
    printf("  a = "); print_complex(a); printf("\n");
    printf("  b = "); print_complex(b); printf(" = j\n");
    printf("  c = "); print_complex(c); printf("\n");
    printf("  d = "); print_complex(d); printf(" = -j\n");
    
    /* Apply butterfly with unity twiddle factors */
    radix4_butterfly_demo(&a, &b, &c, &d, 1.0, 1.0, 1.0);
    
    printf("\nOutput values (with W=1):\n");
    printf("  A = "); print_complex(a); printf("\n");
    printf("  B = "); print_complex(b); printf("\n");
    printf("  C = "); print_complex(c); printf("\n");
    printf("  D = "); print_complex(d); printf("\n");
    
    /* Mathematical verification */
    printf("\nThis implements the 4-point DFT matrix multiplication\n");
}

/**
 * @brief Analyze operation count savings
 */
static void analyze_operation_count() {
    printf("\n\nOperation Count Analysis:\n");
    printf("========================\n");
    printf("N\tRadix-2 Muls\tRadix-4 Muls\tSavings\n");
    printf("----\t------------\t------------\t-------\n");
    
    for (int n = 16; n <= 4096; n *= 4) {
        int log2n = log2_int(n);
        
        /* Radix-2: (N/2) × log₂(N) complex multiplications */
        int radix2_muls = (n / 2) * log2n;
        
        /* Radix-4: (3N/8) × log₂(N) complex multiplications */
        int radix4_muls = (3 * n / 8) * log2n;
        
        double savings = (1.0 - (double)radix4_muls / radix2_muls) * 100;
        
        printf("%d\t%d\t\t%d\t\t%.1f%%\n", n, radix2_muls, radix4_muls, savings);
    }
    
    printf("\nRadix-4 uses 25%% fewer multiplications than Radix-2\n");
}

/**
 * @brief Main test and demonstration program
 */
int main() {
    printf("Radix-4 FFT Implementation\n");
    printf("==========================\n");
    
    /* Test with various sizes */
    int test_sizes[] = {16, 64, 256, 1024, 4096};
    int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);
    
    printf("\nPerformance Comparison:\n");
    printf("======================\n");
    
    for (int t = 0; t < num_tests; t++) {
        int n = test_sizes[t];
        printf("\nTesting size n = %d:\n", n);
        
        /* Allocate test arrays */
        complex_t* signal_r2 = allocate_complex_array(n);
        complex_t* signal_r4 = allocate_complex_array(n);
        complex_t* original = allocate_complex_array(n);
        
        /* Generate random signal */
        for (int i = 0; i < n; i++) {
            double real = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double imag = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            signal_r2[i] = signal_r4[i] = original[i] = real + I * imag;
        }
        
        /* Time Radix-2 */
        fft_timer_t timer;
        timer_start(&timer);
        radix2_dit_fft(signal_r2, n, FFT_FORWARD);
        timer_stop(&timer);
        double time_r2 = timer.elapsed_ms;
        
        /* Time Radix-4 */
        timer_start(&timer);
        fft_radix4(signal_r4, n);
        timer_stop(&timer);
        double time_r4 = timer.elapsed_ms;
        
        printf("  Radix-2: %.3f ms\n", time_r2);
        printf("  Radix-4: %.3f ms\n", time_r4);
        printf("  Speedup: %.2fx\n", time_r2 / time_r4);
        
        /* Verify correctness */
        double max_error = 0;
        for (int i = 0; i < n; i++) {
            double error = cabs(signal_r2[i] - signal_r4[i]);
            if (error > max_error) max_error = error;
        }
        printf("  Max Error: %.2e\n", max_error);
        
        free_complex_array(signal_r2);
        free_complex_array(signal_r4);
        free_complex_array(original);
    }
    
    demonstrate_butterfly();
    analyze_operation_count();
    
    return 0;
}
#endif /* LIB_BUILD */
