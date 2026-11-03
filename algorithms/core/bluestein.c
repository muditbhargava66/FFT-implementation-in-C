#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/*
 * Bluestein's chirp-z transform: an n-point DFT of any length via one
 * power-of-two convolution. With c[m] = exp(dir*pi*i*m^2/n), the identity
 * k*n = (n^2 + k^2 - (k-n)^2)/2 gives
 *
 *   X[k] = c[k] * SUM_j x[j]*c[j]*conj(c[k-j]),
 *
 * which is a linear convolution of x*c with conj(c). Evaluating it with an
 * m-point FFT, m = next_pow2(2n-1), costs O(n log n) for any n, including
 * primes. The price is three m-point FFTs and about 9n complex numbers of
 * scratch, so fft_auto only uses this path when no cheaper factorisation
 * exists.
 */

/*
 * Bluestein FFT. Valid for any positive n. The inverse direction keeps the
 * same construction with the opposite chirp sign, then scales by 1/n.
 */
fft_error_t bluestein_fft(complex_t* x, int n, fft_direction dir) {
    /* Input validation */
    if (!x) return FFT_ERR_NULL;
    if (n <= 0) return FFT_ERR_SIZE;
    /* Guards the 2*n - 1 below: past this, the length would wrap negative and
     * the m-sized scratch arrays would be smaller than the n-sized loops. */
    if (n > (1 << 29)) return FFT_ERR_RANGE;

    /* Find next power of 2 >= 2n - 1 for circular convolution */
    int m = next_power_of_two(2 * n - 1);
    
    /* Allocate working arrays */
    complex_t* a = allocate_complex_array(m);
    complex_t* b = allocate_complex_array(m);
    
    /* 
     * Compute chirp sequence
     * Standard definition: c[k] = exp(-j * pi * k^2 / N) for Forward
     * We use 'dir' to control the sign.
     * If dir == FFT_FORWARD (-1), we want exp(-j...).
     * So phase = dir * PI * k^2 / N.
     */
    complex_t* c = allocate_complex_array(n);
    
    if (!a || !b || !c) {
        free_complex_array(a);
        free_complex_array(b);
        free_complex_array(c);
        return FFT_ERR_ALLOC;
    }
    
    /* Precompute chirp sequence with high precision */
    for (int k = 0; k < n; k++) {
        long double phase = (long double)dir * PI * (long double)k * (long double)k / (long double)n;
        c[k] = cexp(I * (double)phase);
    }
    
    /* Initialize arrays to zero */
    memset(a, 0, m * sizeof(complex_t));
    memset(b, 0, m * sizeof(complex_t));
    
    /* 
     * Step 1: Modulate input
     * a[k] = x[k] * c[k]
     */
    for (int k = 0; k < n; k++) {
        a[k] = x[k] * c[k];
    }
    
    /* 
     * Step 2: Prepare convolution kernel
     * b[k] = conj(c[k])
     * b[m-k] = conj(c[k])
     */
    for (int k = 0; k < n; k++) {
        complex_t conj_c = conj(c[k]);
        b[k] = conj_c;
        if (k > 0) {
            b[m - k] = conj_c;
        }
    }
    
    /*
     * Step 3: Forward FFT of both sequences. m is a power of two by
     * construction, so the kernel's size check cannot fire here.
     */
    radix2_dit_fft(a, m, FFT_FORWARD);
    radix2_dit_fft(b, m, FFT_FORWARD);

    /* Step 4: Pointwise multiplication in frequency domain */
    for (int k = 0; k < m; k++) {
        a[k] *= b[k];
    }
    
    /* Step 5: Inverse FFT to get convolution result. Its 1/m factor is part of
     * the convolution theorem, not the 1/n normalisation applied below. */
    radix2_dit_fft(a, m, FFT_INVERSE);

    /* 
     * Step 6: Final modulation and scaling
     * x[k] = a[k] * c[k]
     * 
     * Base magnitude is correct (1.0).
     * For Inverse FFT, we scale by 1/n.
     */
    double scale_factor = 1.0;
    if (dir == FFT_INVERSE) {
        scale_factor /= n;
    }

    for (int k = 0; k < n; k++) {
        x[k] = a[k] * c[k] * scale_factor;
    }

    free_complex_array(a);
    free_complex_array(b);
    free_complex_array(c);
    return FFT_OK;
}

/**
 * @brief Compute forward FFT using Bluestein's algorithm
 * @param x Input/output array
 * @param n Array length
 */
fft_error_t fft_bluestein(complex_t* x, int n) {
    return bluestein_fft(x, n, FFT_FORWARD);
}

/**
 * @brief Compute inverse FFT using Bluestein's algorithm
 * @param x Input/output array
 * @param n Array length
 */
fft_error_t ifft_bluestein(complex_t* x, int n) {
    return bluestein_fft(x, n, FFT_INVERSE);
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
    printf("- Handles any length, including primes\n");
    printf("- O(n log n) via FFT-based convolution\n");
    printf("- Higher constant factor than the radix kernels\n");
    printf("- Needs about 9n complex numbers of scratch\n");
    printf("- Round-off grows with the padded transform length\n");
    
    free_complex_array(signal);
    
    return 0;
}
#endif /* LIB_BUILD */
