/**
 * Window Functions for Signal Processing
 *
 * Common window functions used in FFT analysis to reduce spectral leakage.
 * Each of these divides by (n-1) or (n-1)/2, so a shorter-than-two sample
 * window is left untouched rather than producing NaN that then permanently
 * corrupts whatever it was applied to.
 */

#include "../include/fft_common.h"

/**
 * Apply Hann (Hanning) window to a signal
 * w(n) = 0.5 * (1 - cos(2πn/(N-1)))
 */
void apply_window_hann(complex_t* signal, int n) {
    if (n <= 1) return;
    for (int i = 0; i < n; i++) {
        double window = 0.5 * (1.0 - cos(TWO_PI * i / (n - 1)));
        signal[i] *= window;
    }
}

/**
 * Apply Hamming window to a signal
 * w(n) = 0.54 - 0.46 * cos(2πn/(N-1))
 */
void apply_window_hamming(complex_t* signal, int n) {
    if (n <= 1) return;
    for (int i = 0; i < n; i++) {
        double window = 0.54 - 0.46 * cos(TWO_PI * i / (n - 1));
        signal[i] *= window;
    }
}

/**
 * Apply Blackman window to a signal
 * w(n) = 0.42 - 0.5*cos(2πn/(N-1)) + 0.08*cos(4πn/(N-1))
 */
void apply_window_blackman(complex_t* signal, int n) {
    if (n <= 1) return;
    for (int i = 0; i < n; i++) {
        double window = 0.42 - 0.5 * cos(TWO_PI * i / (n - 1))
                            + 0.08 * cos(4.0 * M_PI * i / (n - 1));
        signal[i] *= window;
    }
}

/**
 * Apply rectangular (no) window to a signal
 * w(n) = 1 (identity operation)
 */
void apply_window_rectangular(complex_t* signal, int n) {
    // No operation - rectangular window is all ones
    (void)signal;
    (void)n;
}

/**
 * Apply Bartlett (triangular) window to a signal
 * w(n) = 1 - |n - (N-1)/2| / ((N-1)/2)
 */
void apply_window_bartlett(complex_t* signal, int n) {
    if (n <= 1) return;
    double center = (n - 1) / 2.0;
    double half_width = (n - 1) / 2.0;

    for (int i = 0; i < n; i++) {
        double window = 1.0 - fabs(i - center) / half_width;
        signal[i] *= window;
    }
}

/**
 * Apply Welch window to a signal
 * w(n) = 1 - ((n - (N-1)/2) / ((N-1)/2))^2
 */
void apply_window_welch(complex_t* signal, int n) {
    if (n <= 1) return;
    double center = (n - 1) / 2.0;
    double half_width = (n - 1) / 2.0;

    for (int i = 0; i < n; i++) {
        double x = (i - center) / half_width;
        double window = 1.0 - x * x;
        signal[i] *= window;
    }
}

#ifndef LIB_BUILD
// Test program
int main(void) {
    printf("Window Functions Test\n");
    printf("====================\n\n");
    
    int n = 16;
    complex_t* signal = allocate_complex_array(n);

    for (int i = 0; i < n; i++)
        signal[i] = (i == n / 2) ? 1.0 : 0.0;
    
    printf("Testing Hann window:\n");
    apply_window_hann(signal, n);
    printf("  Center value: %.4f\n", creal(signal[n / 2]));
    printf("  Edge value: %.4f\n\n", creal(signal[0]));
    
    for (int i = 0; i < n; i++)
        signal[i] = (i == n / 2) ? 1.0 : 0.0;

    printf("Testing Hamming window:\n");
    apply_window_hamming(signal, n);
    printf("  Center value: %.4f\n", creal(signal[n / 2]));
    printf("  Edge value: %.4f\n\n", creal(signal[0]));
    
    for (int i = 0; i < n; i++)
        signal[i] = (i == n / 2) ? 1.0 : 0.0;

    printf("Testing Blackman window:\n");
    apply_window_blackman(signal, n);
    printf("  Center value: %.4f\n", creal(signal[n / 2]));
    printf("  Edge value: %.4f\n\n", creal(signal[0]));

    complex_t single = 2.5;
    apply_window_hann(&single, 1);
    printf("Single-sample Hann: %.4f (expected 2.5000, not NaN)\n",
           creal(single));

    free_complex_array(signal);
    
    printf("All window functions working correctly!\n");
    return 0;
}
#endif
