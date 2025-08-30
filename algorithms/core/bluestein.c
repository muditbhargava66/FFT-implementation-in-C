#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * @file bluestein.c
 * @brief Bluestein's Algorithm (Chirp Z-Transform) Implementation
 * 
 * @details
 * Bluestein's algorithm computes the DFT for arbitrary lengths (not just powers of 2)
 * by converting the DFT into a convolution using a chirp sequence. This allows us to
 * use fast convolution via FFT to compute DFTs of any size.
 * 
 * Mathematical Foundation:
 * The DFT can be rewritten as:
 *   X[k] = sum(n=0 to N-1) x[n] * exp(-2πi*kn/N)
 *        = sum(n=0 to N-1) x[n] * exp(-πi*n²/N) * exp(-πi*k²/N) * exp(πi*(n-k)²/N)
 * 
 * This is a convolution of x[n]*exp(-πi*n²/N) with exp(-πi*n²/N).
 * 
 * Algorithm Steps:
 * 1. Compute chirp sequence: exp(-πi*k²/N)
 * 2. Multiply input by conjugate chirp: a[k] = x[k] * conj(chirp[k])
 * 3. Prepare convolution kernel from chirp sequence
 * 4. Zero-pad both sequences to next power of 2 ≥ 2N-1
 * 5. Compute FFT-based convolution
 * 6. Multiply result by conjugate chirp to get final DFT
 * 
 * @author FFT Study Repository
 * @date 2024
 * 
 * Time Complexity: O(n log n) for any n
 * Space Complexity: O(n)
 * 
 * References:
 * [1] Bluestein, L. (1970). "A linear filtering approach to the computation of 
 *     discrete Fourier transform"
 * [2] Rabiner, L. R., et al. (1969). "The chirp z-transform algorithm"
 */

/**
 * @brief Compute chirp sequence for Bluestein's algorithm
 * 
 * @details
 * The chirp sequence is defined as: exp(-i * π * k² / n) for forward transform
 * or exp(+i * π * k² / n) for inverse transform.
 * 
 * @param n Length of the sequence
 * @param dir Transform direction (FFT_FORWARD or FFT_INVERSE)
 * @return Pointer to allocated chirp sequence array
 */
static complex_t* compute_chirp_sequence(int n, fft_direction dir) {
    complex_t* chirp = allocate_complex_array(n);
    if (!chirp) {
        fprintf(stderr, "Error: Failed to allocate chirp sequence\n");
        return NULL;
    }
    
    /* Compute exp(-dir * i * π * k² / n) for each k */
    for (int k = 0; k < n; k++) {
        double phase = -dir * PI * k * k / n;
        chirp[k] = cexp(I * phase);
    }
    
    return chirp;
}

/**
 * @brief Main Bluestein FFT algorithm implementation
 * 
 * @details
 * This function implements the Chirp Z-Transform algorithm to compute the DFT
 * of arbitrary length sequences. It converts the DFT into a convolution problem
 * which can be solved efficiently using FFT.
 * 
 * @param x Input/output array of complex numbers
 * @param n Length of the array (can be any positive integer)
 * @param dir Transform direction (FFT_FORWARD or FFT_INVERSE)
 */
void bluestein_fft(complex_t* x, int n, fft_direction dir) {
    /* Input validation */
    if (!x || n <= 0) {
        fprintf(stderr, "Error: Invalid input to bluestein_fft\n");
        return;
    }
    
    /* Find next power of 2 >= 2n - 1 for circular convolution */
    int m = next_power_of_two(2 * n - 1);
    
    /* Allocate working arrays */
    complex_t* a = allocate_complex_array(m);
    complex_t* b = allocate_complex_array(m);
    complex_t* chirp = compute_chirp_sequence(n, dir);
    
    if (!a || !b || !chirp) {
        fprintf(stderr, "Error: Memory allocation failed in bluestein_fft\n");
        goto cleanup;
    }
    
    /* Initialize arrays to zero */
    memset(a, 0, m * sizeof(complex_t));
    memset(b, 0, m * sizeof(complex_t));
    
    /* 
     * Step 1: Prepare sequence a[k] = x[k] * conj(chirp[k])
     * This modulates the input with the conjugate chirp
     */
    for (int k = 0; k < n; k++) {
        a[k] = x[k] * conj(chirp[k]);
    }
    
    /* 
     * Step 2: Prepare convolution kernel b[k]
     * b[k] = chirp[k] for k in [0, n-1]
     * b[m-k] = chirp[k] for k in [1, n-1] (wrap around for circular convolution)
     */
    for (int k = 0; k < n; k++) {
        b[k] = chirp[k];
        if (k > 0) {
            b[m - k] = chirp[k];
        }
    }
    
    /* Step 3: Compute FFT of both sequences */
    radix2_dit_fft(a, m, FFT_FORWARD);
    radix2_dit_fft(b, m, FFT_FORWARD);
    
    /* Step 4: Pointwise multiplication in frequency domain */
    for (int k = 0; k < m; k++) {
        a[k] *= b[k];
    }
    
    /* Step 5: Inverse FFT to get convolution result */
    radix2_dit_fft(a, m, FFT_INVERSE);
    
    /* 
     * Step 6: Extract result and demodulate
     * y[k] = a[k] * conj(chirp[k])
     */
    for (int k = 0; k < n; k++) {
        x[k] = a[k] * conj(chirp[k]);
    }
    
    /* Scale for inverse transform */
    if (dir == FFT_INVERSE) {
        for (int k = 0; k < n; k++) {
            x[k] /= n;
        }
    }
    
cleanup:
    /* Free allocated memory */
    if (a) free_complex_array(a);
    if (b) free_complex_array(b);
    if (chirp) free_complex_array(chirp);
}

/**
 * @brief Compute forward FFT using Bluestein's algorithm
 * @param x Input/output array
 * @param n Array length
 */
void fft_bluestein(complex_t* x, int n) {
    bluestein_fft(x, n, FFT_FORWARD);
}

/**
 * @brief Compute inverse FFT using Bluestein's algorithm
 * @param x Input/output array
 * @param n Array length
 */
void ifft_bluestein(complex_t* x, int n) {
    bluestein_fft(x, n, FFT_INVERSE);
}

#ifndef LIB_BUILD
/**
 * @brief Performance and correctness testing
 * 
 * This main function demonstrates:
 * 1. Arbitrary-length FFT capability
 * 2. Performance comparison with naive DFT
 * 3. Accuracy verification
 * 4. Prime-length FFT handling
 */
int main() {
    printf("Bluestein's Algorithm (Chirp Z-Transform)\n");
    printf("=========================================\n\n");
    
    /* Test with various non-power-of-2 sizes */
    int test_sizes[] = {13, 17, 23, 30, 50, 100, 127, 255, 500, 1000};
    int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);
    
    printf("Testing arbitrary-length FFTs:\n");
    printf("Size\tBluestein\tNaive DFT\tSpeedup\tError\n");
    printf("----\t---------\t---------\t-------\t-----\n");
    
    for (int t = 0; t < num_tests; t++) {
        int n = test_sizes[t];
        
        /* Allocate test arrays */
        complex_t* x_blue = allocate_complex_array(n);
        complex_t* x_naive = allocate_complex_array(n);
        complex_t* original = allocate_complex_array(n);
        
        /* Generate test signal with multiple frequency components */
        for (int i = 0; i < n; i++) {
            double value = sin(2 * PI * 3 * i / n) + 
                          0.5 * cos(2 * PI * 7 * i / n);
            x_blue[i] = x_naive[i] = original[i] = value;
        }
        
        /* Time Bluestein's algorithm */
        fft_timer_t timer;
        timer_start(&timer);
        fft_bluestein(x_blue, n);
        timer_stop(&timer);
        double time_blue = timer.elapsed_ms;
        
        /* Time naive DFT for small sizes only */
        double time_naive = -1;
        if (n <= 100) {
            timer_start(&timer);
            naive_dft(x_naive, n, FFT_FORWARD);
            timer_stop(&timer);
            time_naive = timer.elapsed_ms;
        }
        
        /* Verify correctness by comparing with naive DFT */
        double max_error = 0;
        if (time_naive > 0) {
            for (int i = 0; i < n; i++) {
                double error = cabs(x_blue[i] - x_naive[i]);
                if (error > max_error) max_error = error;
            }
        }
        
        /* Test inverse transform */
        ifft_bluestein(x_blue, n);
        double recon_error = 0;
        for (int i = 0; i < n; i++) {
            recon_error += cabs(x_blue[i] - original[i]);
        }
        recon_error /= n;
        
        /* Print results */
        printf("%d\t%.3f ms\t", n, time_blue);
        if (time_naive > 0) {
            printf("%.3f ms\t%.1fx\t%.2e\n", 
                   time_naive, time_naive / time_blue, max_error);
        } else {
            printf("--\t\t--\t%.2e\n", recon_error);
        }
        
        free_complex_array(x_blue);
        free_complex_array(x_naive);
        free_complex_array(original);
    }
    
    /* Demonstrate with prime-length FFT */
    printf("\n\nPrime-length FFT example (n = 31):\n");
    printf("----------------------------------\n");
    
    int n = 31;
    complex_t* signal = allocate_complex_array(n);
    
    /* Create unit impulse */
    generate_impulse(signal, n);
    printf("Input: Impulse at index 0\n");
    
    /* Compute FFT */
    fft_bluestein(signal, n);
    
    /* Verify all magnitudes are 1 (property of impulse FFT) */
    printf("FFT magnitudes: ");
    int all_ones = 1;
    for (int i = 0; i < 10; i++) {
        double mag = cabs(signal[i]);
        printf("%.3f ", mag);
        if (fabs(mag - 1.0) > 1e-10) all_ones = 0;
    }
    printf("...\n");
    printf("All magnitudes = 1.0? %s\n", all_ones ? "YES (PASS)" : "NO (FAIL)");
    
    /* Demonstrate flexibility with various sizes */
    printf("\n\nFlexibility demonstration:\n");
    printf("-------------------------\n");
    for (int size = 10; size <= 20; size++) {
        complex_t* x = allocate_complex_array(size);
        
        /* Create ramp signal */
        for (int i = 0; i < size; i++) {
            x[i] = i + 0 * I;
        }
        
        fft_bluestein(x, size);
        printf("FFT of size %d: DC component = ", size);
        print_complex(x[0]);
        printf(" (sum = %.0f)\n", creal(x[0]));
        
        free_complex_array(x);
    }
    
    /* Algorithm characteristics */
    printf("\n\nAlgorithm Characteristics:\n");
    printf("=========================\n");
    printf("✓ Works for any input size (prime, composite, power-of-2)\n");
    printf("✓ O(n log n) complexity via FFT-based convolution\n");
    printf("✓ Higher constant factor than specialized algorithms\n");
    printf("✓ Ideal for sizes with large prime factors\n");
    printf("✓ Numerically stable for well-conditioned problems\n");
    
    free_complex_array(signal);
    
    return 0;
}
#endif /* LIB_BUILD */
