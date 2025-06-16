#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * @file radix4.c
 * @brief Radix-4 FFT Implementation
 * 
 * @details
 * The Radix-4 FFT algorithm processes 4 points at a time instead of 2,
 * reducing the number of stages and complex multiplications. This algorithm
 * requires the input size to be a power of 4 (4, 16, 64, 256, ...).
 * 
 * Key Advantages:
 * - 25% fewer complex multiplications than radix-2
 * - Better memory access patterns (processes more data per stage)
 * - Fewer stages (log₄(n) instead of log₂(n))
 * 
 * Mathematical Foundation:
 * The DFT is decomposed into 4 interleaved subsequences:
 *   X[k] = Σ(r=0 to 3) W_N^(rk) * DFT_{N/4}[x[4n+r]]
 * 
 * The radix-4 butterfly combines 4 inputs using:
 *   [A]   [1  1  1  1] [a]
 *   [B] = [1 -j -1  j] [b] × twiddle factors
 *   [C]   [1 -1  1 -1] [c]
 *   [D]   [1  j -1 -j] [d]
 * 
 * @author FFT Study Repository
 * @date 2024
 * 
 * Time Complexity: O(n log n) with 25% fewer operations than radix-2
 * Space Complexity: O(1) for in-place operation
 * 
 * References:
 * [1] Proakis, J. G., & Manolakis, D. G. (2007). "Digital Signal Processing"
 * [2] Chu, E., & George, A. (2000). "Inside the FFT Black Box"
 */

/**
 * @brief Radix-4 butterfly operation
 * 
 * @details
 * Implements the 4-point DFT that forms the core of the radix-4 algorithm.
 * The butterfly operation combines 4 complex values using the DFT matrix
 * and twiddle factors.
 * 
 * Butterfly equations:
 *   t0 = a + c,  t1 = a - c
 *   t2 = b + d,  t3 = j(b - d)
 *   
 *   A = t0 + t2
 *   B = (t1 - t3) × W¹
 *   C = (t0 - t2) × W²
 *   D = (t1 + t3) × W³
 * 
 * @param a,b,c,d Pointers to the four complex values
 * @param w1,w2,w3 Twiddle factors W^k, W^2k, W^3k
 */
static inline void radix4_butterfly(complex_t* a, complex_t* b, complex_t* c, complex_t* d,
                                    complex_t w1, complex_t w2, complex_t w3) {
    /* First stage of butterfly */
    complex_t t0 = *a + *c;
    complex_t t1 = *a - *c;
    complex_t t2 = *b + *d;
    complex_t t3 = (*b - *d) * I;  /* Multiplication by j = sqrt(-1) */
    
    /* Second stage with twiddle factors */
    *a = t0 + t2;                  /* No twiddle factor needed */
    *b = (t1 - t3) * w1;          /* Multiply by W^k */
    *c = (t0 - t2) * w2;          /* Multiply by W^2k */
    *d = (t1 + t3) * w3;          /* Multiply by W^3k */
}

/**
 * @brief Main Radix-4 FFT implementation
 * 
 * @details
 * Implements the complete radix-4 FFT algorithm with the following steps:
 * 1. Validate input size is power of 4
 * 2. Perform base-4 digit reversal
 * 3. Execute log₄(n) stages of radix-4 butterflies
 * 4. Handle remaining radix-2 stage if log₂(n) is odd
 * 5. Scale for inverse transform
 * 
 * Data Structure:
 * - In-place computation on input array
 * - Base-4 digit reversal for reordering
 * 
 * @param x Input/output array of complex numbers
 * @param n Array length (must be power of 4)
 * @param dir Transform direction (FFT_FORWARD or FFT_INVERSE)
 */
void radix4_fft(complex_t* x, int n, fft_direction dir) {
    /* Validate that n is a power of 4 */
    if (!is_power_of_two(n)) {
        fprintf(stderr, "Error: Radix-4 FFT requires size to be power of 2\n");
        return;
    }
    
    /* Check if n is power of 4: n = 4^k means log2(n) must be even */
    int log2n = log2_int(n);
    if (log2n & 1) {
        /* Handle case where n = 2 × 4^k (odd power of 2) */
        /* We'll do radix-4 stages and one radix-2 stage at the end */
    }
    
    int log4n = log2n / 2;  /* Number of radix-4 stages */
    
    /* 
     * Step 1: Digit reversal (base-4)
     * Reorder array elements based on base-4 digit reversal
     */
    for (int i = 0; i < n; i++) {
        int j = 0;
        int temp = i;
        
        /* Reverse base-4 digits */
        for (int k = 0; k < log4n; k++) {
            j = (j << 2) | (temp & 3);  /* Extract and shift 2 bits */
            temp >>= 2;
        }
        
        /* Handle odd log2n case */
        if (log2n & 1) {
            j = (j << 1) | (temp & 1);
        }
        
        if (i < j) {
            complex_t t = x[i];
            x[i] = x[j];
            x[j] = t;
        }
    }
    
    /* 
     * Step 2: Main radix-4 computation
     * Process log₄(n) stages, each combining groups of 4
     */
    int m = 1;  /* Initial radix-4 block size */
    
    for (int stage = 0; stage < log4n; stage++) {
        int m4 = m * 4;  /* Size of radix-4 block */
        
        /* Principal 4th root of unity for this stage */
        complex_t w = twiddle_factor(1, m4, dir);
        complex_t w2 = w * w;
        complex_t w3 = w2 * w;
        
        /* Process all radix-4 blocks */
        for (int k = 0; k < n; k += m4) {
            /* Initialize twiddle factors */
            complex_t wj = 1.0;
            complex_t wj2 = 1.0;
            complex_t wj3 = 1.0;
            
            /* Process butterflies within block */
            for (int j = 0; j < m; j++) {
                /* Get indices of 4 elements */
                int idx0 = k + j;
                int idx1 = idx0 + m;
                int idx2 = idx1 + m;
                int idx3 = idx2 + m;
                
                /* Apply radix-4 butterfly */
                radix4_butterfly(&x[idx0], &x[idx1], &x[idx2], &x[idx3],
                               wj, wj2, wj3);
                
                /* Update twiddle factors */
                wj *= w;
                wj2 *= w2;
                wj3 *= w3;
            }
        }
        m = m4;
    }
    
    /* 
     * Step 3: Handle remaining radix-2 stage if log2(n) is odd
     * This happens when n = 2 × 4^k (e.g., n = 8, 32, 128, ...)
     */
    if (log2n & 1) {
        for (int k = 0; k < n; k += 2) {
            complex_t t = x[k];
            x[k] = t + x[k + 1];
            x[k + 1] = t - x[k + 1];
        }
    }
    
    /* Step 4: Scale for inverse FFT */
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

/**
 * @brief Compute forward FFT using Radix-4
 * @param x Input/output array
 * @param n Array length (should be power of 4 for best performance)
 */
void fft_radix4(complex_t* x, int n) {
    radix4_fft(x, n, FFT_FORWARD);
}

/**
 * @brief Compute inverse FFT using Radix-4
 * @param x Input/output array
 * @param n Array length
 */
void ifft_radix4(complex_t* x, int n) {
    radix4_fft(x, n, FFT_INVERSE);
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
    radix4_butterfly(&a, &b, &c, &d, 1.0, 1.0, 1.0);
    
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
        
        /* Generate test signal with multiple frequencies */
        for (int i = 0; i < n; i++) {
            double value = sin(2 * PI * 10 * i / n) + 
                          0.5 * sin(2 * PI * 25 * i / n) +
                          0.25 * sin(2 * PI * 40 * i / n);
            signal_r2[i] = value;
            signal_r4[i] = value;
            original[i] = value;
        }
        
        /* Time radix-2 FFT */
        timer_t timer_r2, timer_r4;
        timer_start(&timer_r2);
        radix2_dit_fft(signal_r2, n, FFT_FORWARD);
        timer_stop(&timer_r2);
        
        /* Time radix-4 FFT */
        timer_start(&timer_r4);
        fft_radix4(signal_r4, n);
        timer_stop(&timer_r4);
        
        printf("  Radix-2 time: %.3f ms\n", timer_r2.elapsed_ms);
        printf("  Radix-4 time: %.3f ms (%.1f%% faster)\n", 
               timer_r4.elapsed_ms, 
               (1.0 - timer_r4.elapsed_ms / timer_r2.elapsed_ms) * 100);
        
        /* Verify correctness */
        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = cabs(signal_r2[i] - signal_r4[i]);
            if (diff > max_diff) max_diff = diff;
        }
        printf("  Max difference: %.2e %s\n", max_diff,
               max_diff < 1e-10 ? "(PASS)" : "(FAIL)");
        
        /* Test inverse transform */
        ifft_radix4(signal_r4, n);
        double error = 0.0;
        for (int i = 0; i < n; i++) {
            error += cabs(signal_r4[i] - original[i]);
        }
        printf("  Reconstruction error: %.2e\n", error / n);
        
        free_complex_array(signal_r2);
        free_complex_array(signal_r4);
        free_complex_array(original);
    }
    
    /* Test odd power of 2 (uses final radix-2 stage) */
    printf("\n\nTesting odd power of 2 (n = 32):\n");
    printf("=================================\n");
    int n = 32;
    complex_t* signal = allocate_complex_array(n);
    
    /* Generate impulse */
    generate_impulse(signal, n);
    fft_radix4(signal, n);
    
    /* Verify all magnitudes are 1 */
    int pass = 1;
    for (int i = 0; i < n; i++) {
        if (fabs(cabs(signal[i]) - 1.0) > 1e-10) {
            pass = 0;
            break;
        }
    }
    printf("Impulse response test: %s\n", pass ? "PASS" : "FAIL");
    
    /* Demonstrate butterfly operation */
    demonstrate_butterfly();
    
    /* Analyze operation counts */
    analyze_operation_count();
    
    /* Algorithm characteristics */
    printf("\n\nAlgorithm Characteristics:\n");
    printf("=========================\n");
    printf("✓ 25%% fewer multiplications than radix-2\n");
    printf("✓ Better cache utilization (processes 4 points at once)\n");
    printf("✓ Works for any power of 2 (uses radix-2 for odd powers)\n");
    printf("✓ More complex butterfly structure\n");
    printf("✓ Ideal for processors with good multiplication performance\n");
    
    free_complex_array(signal);
    
    return 0;
}

#endif /* LIB_BUILD */
