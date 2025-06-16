#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"

/**
 * FFT-based Convolution
 * 
 * Implements linear convolution using the FFT, which is more efficient
 * than direct convolution for large signals.
 * 
 * Algorithm:
 * 1. Zero-pad both signals to length N = len(x) + len(y) - 1
 * 2. Compute FFT of both signals
 * 3. Multiply in frequency domain
 * 4. Compute inverse FFT
 * 
 * Time Complexity: O(n log n) vs O(n²) for direct convolution
 */

// Direct convolution for comparison
void direct_convolution(complex_t* x, int nx, complex_t* h, int nh, complex_t* y) {
    int ny = nx + nh - 1;
    
    for (int n = 0; n < ny; n++) {
        y[n] = 0;
        for (int k = 0; k < nh; k++) {
            if (n - k >= 0 && n - k < nx) {
                y[n] += x[n - k] * h[k];
            }
        }
    }
}

// FFT-based convolution
void fft_convolution(complex_t* x, int nx, complex_t* h, int nh, complex_t* y) {
    // Output length
    int ny = nx + nh - 1;
    
    // Find next power of 2 for FFT
    int n_fft = next_power_of_two(ny);
    
    // Allocate zero-padded arrays
    complex_t* x_pad = allocate_complex_array(n_fft);
    complex_t* h_pad = allocate_complex_array(n_fft);
    
    // Zero-pad inputs
    memset(x_pad, 0, n_fft * sizeof(complex_t));
    memset(h_pad, 0, n_fft * sizeof(complex_t));
    memcpy(x_pad, x, nx * sizeof(complex_t));
    memcpy(h_pad, h, nh * sizeof(complex_t));
    
    // Compute FFTs
    radix2_dit_fft(x_pad, n_fft, FFT_FORWARD);
    radix2_dit_fft(h_pad, n_fft, FFT_FORWARD);
    
    // Multiply in frequency domain
    for (int i = 0; i < n_fft; i++) {
        x_pad[i] *= h_pad[i];
    }
    
    // Inverse FFT
    radix2_dit_fft(x_pad, n_fft, FFT_INVERSE);
    
    // Copy result (only valid part)
    memcpy(y, x_pad, ny * sizeof(complex_t));
    
    free_complex_array(x_pad);
    free_complex_array(h_pad);
}

// Circular convolution using FFT
void circular_convolution(complex_t* x, complex_t* h, int n, complex_t* y) {
    complex_t* X = allocate_complex_array(n);
    complex_t* H = allocate_complex_array(n);
    
    // Copy inputs
    memcpy(X, x, n * sizeof(complex_t));
    memcpy(H, h, n * sizeof(complex_t));
    
    // FFT of both signals
    radix2_dit_fft(X, n, FFT_FORWARD);
    radix2_dit_fft(H, n, FFT_FORWARD);
    
    // Multiply in frequency domain
    for (int i = 0; i < n; i++) {
        X[i] *= H[i];
    }
    
    // Inverse FFT
    radix2_dit_fft(X, n, FFT_INVERSE);
    
    // Copy result
    memcpy(y, X, n * sizeof(complex_t));
    
    free_complex_array(X);
    free_complex_array(H);
}

// 2D convolution for image processing
void fft_convolution_2d(complex_t** image, int rows, int cols,
                       complex_t** kernel, int krows, int kcols,
                       complex_t** result) {
    // This is a placeholder for 2D convolution
    // Full implementation would require 2D FFT
    // TODO: Implement 2D FFT convolution
    (void)image; (void)rows; (void)cols;
    (void)kernel; (void)krows; (void)kcols;
    (void)result;
    printf("2D convolution requires 2D FFT implementation\n");
}

// Demonstration functions
void demo_convolution_types() {
    printf("\nConvolution Types Demonstration:\n");
    printf("================================\n");
    
    // Create signals
    int nx = 8, nh = 4;
    complex_t* x = allocate_complex_array(nx);
    complex_t* h = allocate_complex_array(nh);
    
    // Signal: unit impulse at different positions
    for (int i = 0; i < nx; i++) x[i] = 0;
    x[2] = 1;  // Impulse at position 2
    
    // Impulse response (moving average filter)
    for (int i = 0; i < nh; i++) {
        h[i] = 0.25;  // Averaging filter
    }
    
    printf("Signal x: ");
    for (int i = 0; i < nx; i++) printf("%.1f ", creal(x[i]));
    printf("\n");
    
    printf("Filter h: ");
    for (int i = 0; i < nh; i++) printf("%.2f ", creal(h[i]));
    printf("\n");
    
    // Linear convolution
    int ny = nx + nh - 1;
    complex_t* y_linear = allocate_complex_array(ny);
    fft_convolution(x, nx, h, nh, y_linear);
    
    printf("\nLinear convolution result: ");
    for (int i = 0; i < ny; i++) {
        printf("%.2f ", creal(y_linear[i]));
    }
    printf("\n");
    
    // Circular convolution
    complex_t* y_circular = allocate_complex_array(nx);
    circular_convolution(x, h, nx, y_circular);
    
    printf("Circular convolution result: ");
    for (int i = 0; i < nx; i++) {
        printf("%.2f ", creal(y_circular[i]));
    }
    printf("\n");
    
    free_complex_array(x);
    free_complex_array(h);
    free_complex_array(y_linear);
    free_complex_array(y_circular);
}

// Performance comparison
void benchmark_convolution() {
    printf("\n\nConvolution Performance Comparison:\n");
    printf("===================================\n");
    printf("Signal Size\tFilter Size\tDirect (ms)\tFFT (ms)\tSpeedup\n");
    printf("-----------\t-----------\t-----------\t--------\t-------\n");
    
    int signal_sizes[] = {128, 256, 512, 1024};
    int filter_sizes[] = {16, 32, 64, 128};
    int num_tests = 4;
    
    for (int i = 0; i < num_tests; i++) {
        int nx = signal_sizes[i];
        int nh = filter_sizes[i];
        int ny = nx + nh - 1;
        
        // Create test signals
        complex_t* x = allocate_complex_array(nx);
        complex_t* h = allocate_complex_array(nh);
        complex_t* y_direct = allocate_complex_array(ny);
        complex_t* y_fft = allocate_complex_array(ny);
        
        // Fill with random data
        for (int j = 0; j < nx; j++) {
            x[j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
        for (int j = 0; j < nh; j++) {
            h[j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
        
        timer_t timer;
        
        // Time direct convolution
        timer_start(&timer);
        direct_convolution(x, nx, h, nh, y_direct);
        timer_stop(&timer);
        double time_direct = timer.elapsed_ms;
        
        // Time FFT convolution
        timer_start(&timer);
        fft_convolution(x, nx, h, nh, y_fft);
        timer_stop(&timer);
        double time_fft = timer.elapsed_ms;
        
        // Verify correctness
        double error = 0;
        for (int j = 0; j < ny; j++) {
            error += cabs(y_direct[j] - y_fft[j]);
        }
        
        printf("%d\t\t%d\t\t%.2f\t\t%.2f\t\t%.1fx",
               nx, nh, time_direct, time_fft, time_direct / time_fft);
        
        if (error / ny > 1e-10) {
            printf(" (Error: %.2e)", error / ny);
        }
        printf("\n");
        
        free_complex_array(x);
        free_complex_array(h);
        free_complex_array(y_direct);
        free_complex_array(y_fft);
    }
}

// Main demonstration
int main() {
    printf("FFT-based Convolution\n");
    printf("=====================\n");
    
    // Simple example
    printf("\nSimple Convolution Example:\n");
    printf("---------------------------\n");
    
    // Create simple signals
    int nx = 5, nh = 3;
    complex_t x[] = {1, 2, 3, 4, 5};
    complex_t h[] = {1, 0, -1};  // Derivative filter
    
    int ny = nx + nh - 1;
    complex_t* y = allocate_complex_array(ny);
    
    printf("Signal: ");
    for (int i = 0; i < nx; i++) printf("%.0f ", creal(x[i]));
    printf("\n");
    
    printf("Filter: ");
    for (int i = 0; i < nh; i++) printf("%.0f ", creal(h[i]));
    printf("\n");
    
    // Compute convolution
    fft_convolution(x, nx, h, nh, y);
    
    printf("Result: ");
    for (int i = 0; i < ny; i++) {
        printf("%.0f ", creal(y[i]));
    }
    printf("\n");
    printf("(This is the discrete derivative of the input)\n");
    
    free_complex_array(y);
    
    // Demonstrate different convolution types
    demo_convolution_types();
    
    // Performance benchmark
    benchmark_convolution();
    
    // Applications
    printf("\n\nPractical Applications:\n");
    printf("======================\n");
    printf("1. Digital filtering (FIR filters)\n");
    printf("2. Audio effects (reverb, echo)\n");
    printf("3. Image blurring and sharpening\n");
    printf("4. Pattern matching\n");
    printf("5. Polynomial multiplication\n");
    printf("6. Signal correlation\n");
    
    // Overlap-add method for long signals
    printf("\n\nOverlap-Add Method:\n");
    printf("===================\n");
    printf("For very long signals, use block convolution:\n");
    printf("1. Divide signal into overlapping blocks\n");
    printf("2. Convolve each block using FFT\n");
    printf("3. Add overlapping regions\n");
    printf("This allows processing of streaming data!\n");
    
    return 0;
}
