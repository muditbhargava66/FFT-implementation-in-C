#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Cooley-Tukey Radix-2 Decimation-in-Frequency (DIF) FFT
 * 
 * Unlike DIT, DIF performs the twiddle factor multiplication
 * before the recursive calls, dividing based on the first and
 * second halves of the input.
 * 
 * Time Complexity: O(n log n)
 * Space Complexity: O(1) for in-place operation
 */

void radix2_dif_fft(complex_t* x, int n, fft_direction dir) {
    CHECK_POWER_OF_TWO(n);
    
    int log2n = log2_int(n);
    
    // Danielson-Lanczos algorithm (DIF variant)
    for (int stage = log2n; stage >= 1; stage--) {
        int m = 1 << stage;
        int half_m = m >> 1;
        
        complex_t w_m = twiddle_factor(1, m, dir);
        
        for (int k = 0; k < n; k += m) {
            complex_t w = 1.0;
            
            for (int j = 0; j < half_m; j++) {
                int t = k + j;
                int u = t + half_m;
                
                complex_t temp = x[t];
                x[t] = temp + x[u];
                x[u] = (temp - x[u]) * w;
                
                w *= w_m;
            }
        }
    }
    
    // Bit-reversal permutation (done after computation in DIF)
    for (int i = 0; i < n; i++) {
        int j = bit_reverse(i, log2n);
        if (i < j) {
            complex_t temp = x[i];
            x[i] = x[j];
            x[j] = temp;
        }
    }
    
    // Scale for inverse FFT
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

// Wrapper functions
void fft_radix2_dif(complex_t* x, int n) {
    radix2_dif_fft(x, n, FFT_FORWARD);
}

void ifft_radix2_dif(complex_t* x, int n) {
    radix2_dif_fft(x, n, FFT_INVERSE);
}

// Demonstration and comparison with DIT
int main() {
    int n = 16;
    complex_t* signal_dit = allocate_complex_array(n);
    complex_t* signal_dif = allocate_complex_array(n);
    
    // Generate identical test signals
    for (int i = 0; i < n; i++) {
        double value = sin(2 * PI * 3 * i / n) + 0.5 * cos(2 * PI * 7 * i / n);
        signal_dit[i] = value;
        signal_dif[i] = value;
    }
    
    printf("Radix-2 DIF FFT Implementation\n");
    printf("==============================\n\n");
    
    printf("Input signal (first 8 samples): ");
    for (int i = 0; i < 8; i++) {
        printf("%.3f ", creal(signal_dif[i]));
    }
    printf("...\n\n");
    
    // Compare DIT and DIF results
    timer_t timer_dit, timer_dif;
    
    timer_start(&timer_dit);
    radix2_dit_fft(signal_dit, n, FFT_FORWARD);
    timer_stop(&timer_dit);
    
    timer_start(&timer_dif);
    fft_radix2_dif(signal_dif, n);
    timer_stop(&timer_dif);
    
    printf("Timing comparison:\n");
    printf("DIT FFT: %.3f ms\n", timer_dit.elapsed_ms);
    printf("DIF FFT: %.3f ms\n\n", timer_dif.elapsed_ms);
    
    // Verify results match
    printf("Verifying DIT and DIF produce identical results:\n");
    double max_diff = 0.0;
    for (int i = 0; i < n; i++) {
        double diff = cabs(signal_dit[i] - signal_dif[i]);
        if (diff > max_diff) max_diff = diff;
    }
    printf("Maximum difference: %.2e %s\n\n", max_diff, 
           max_diff < 1e-10 ? "(PASS)" : "(FAIL)");
    
    // Show FFT results
    print_complex_array("DIF FFT result", signal_dif, n);
    
    // Test inverse transform
    ifft_radix2_dif(signal_dif, n);
    printf("\nAfter inverse DIF FFT (first 8 samples): ");
    for (int i = 0; i < 8; i++) {
        printf("%.3f ", creal(signal_dif[i]));
    }
    printf("...\n");
    
    // Performance test with larger sizes
    printf("\n\nPerformance comparison for various sizes:\n");
    printf("Size\tDIT (ms)\tDIF (ms)\tRatio\n");
    printf("----------------------------------------\n");
    
    for (int size = 256; size <= 8192; size *= 2) {
        complex_t* test_dit = allocate_complex_array(size);
        complex_t* test_dif = allocate_complex_array(size);
        
        // Generate random data
        for (int i = 0; i < size; i++) {
            double real = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            double imag = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
            test_dit[i] = real + I * imag;
            test_dif[i] = test_dit[i];
        }
        
        timer_start(&timer_dit);
        radix2_dit_fft(test_dit, size, FFT_FORWARD);
        timer_stop(&timer_dit);
        
        timer_start(&timer_dif);
        fft_radix2_dif(test_dif, size);
        timer_stop(&timer_dif);
        
        printf("%d\t%.3f\t\t%.3f\t\t%.2f\n", 
               size, timer_dit.elapsed_ms, timer_dif.elapsed_ms,
               timer_dif.elapsed_ms / timer_dit.elapsed_ms);
        
        free_complex_array(test_dit);
        free_complex_array(test_dif);
    }
    
    free_complex_array(signal_dit);
    free_complex_array(signal_dif);
    
    return 0;
}
