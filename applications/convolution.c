#include "fft_common.h"
#include "fft_algorithms.h"
#include <limits.h>

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
    if (!x || !h || !y || nx <= 0 || nh <= 0)
        return;

    long long ny_ll = (long long)nx + nh - 1;
    if (ny_ll > INT_MAX)
        return;
    int ny = (int)ny_ll;
    
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
    if (!x || !h || !y || nx <= 0 || nh <= 0) return;

    // Output length. The sum is evaluated in 64 bits so two large ints cannot
    // wrap to a small positive length.
    long long ny_ll = (long long)nx + nh - 1;
    if (ny_ll > INT_MAX) return;
    int ny = (int)ny_ll;
    
    // Find next power of 2 for FFT
    int n_fft = next_power_of_two(ny);
    if (n_fft <= 0) return;   /* ny large enough to overflow the FFT size */
    
    // Allocate zero-padded arrays (both allocators zero the memory)
    complex_t* x_pad = allocate_complex_array(n_fft);
    complex_t* h_pad = allocate_complex_array(n_fft);
    if (!x_pad || !h_pad) {
        free_complex_array(x_pad);
        free_complex_array(h_pad);
        return;
    }
    memcpy(x_pad, x, (size_t)nx * sizeof(complex_t));
    memcpy(h_pad, h, (size_t)nh * sizeof(complex_t));
    
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
    memcpy(y, x_pad, (size_t)ny * sizeof(complex_t));
    
    free_complex_array(x_pad);
    free_complex_array(h_pad);
}

// Circular convolution using FFT.
// x, h and y are all n elements long; shorter inputs must be zero-extended by
// the caller, otherwise the copies below read past the end of the array.
void circular_convolution(complex_t* x, complex_t* h, int n, complex_t* y) {
    if (!x || !h || !y || n <= 0) return;

    complex_t* X = allocate_complex_array(n);
    complex_t* H = allocate_complex_array(n);
    if (!X || !H) {
        free_complex_array(X);
        free_complex_array(H);
        return;
    }

    // Copy inputs
    memcpy(X, x, (size_t)n * sizeof(complex_t));
    memcpy(H, h, (size_t)n * sizeof(complex_t));
    
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

// 2D convolution via row-column FFT decomposition
//
// 1. Zero-pad image and kernel to (fft_rows x fft_cols)
// 2. Row-wise FFT, then column-wise FFT on both
// 3. Pointwise multiply in 2-D frequency domain
// 4. Inverse column-wise FFT, then inverse row-wise FFT
// 5. Extract (rows+krows-1) x (cols+kcols-1) valid region
void fft_convolution_2d(complex_t** image, int rows, int cols,
                       complex_t** kernel, int krows, int kcols,
                       complex_t** result) {
    if (!image || !kernel || !result || rows <= 0 || cols <= 0 ||
        krows <= 0 || kcols <= 0)
        return;

    long long out_rows_ll = (long long)rows + krows - 1;
    long long out_cols_ll = (long long)cols + kcols - 1;
    if (out_rows_ll > INT_MAX || out_cols_ll > INT_MAX)
        return;

    int out_rows = (int)out_rows_ll;
    int out_cols = (int)out_cols_ll;
    int fft_rows = next_power_of_two(out_rows);
    int fft_cols = next_power_of_two(out_cols);
    if (fft_rows <= 0 || fft_cols <= 0)
        return;
    if ((long long)fft_rows * fft_cols > INT_MAX)
        return;

    int total = fft_rows * fft_cols;

    /* Flatten into row-major 1-D arrays, zero-padded (the allocator zeroes). */
    complex_t* img_buf = allocate_complex_array(total);
    complex_t* ker_buf = allocate_complex_array(total);
    complex_t* tmp = allocate_complex_array(total);
    if (!img_buf || !ker_buf || !tmp) {
        free_complex_array(img_buf);
        free_complex_array(ker_buf);
        free_complex_array(tmp);
        return;
    }

    for (int r = 0; r < rows; r++)
        for (int c = 0; c < cols; c++)
            img_buf[r * fft_cols + c] = image[r][c];

    for (int r = 0; r < krows; r++)
        for (int c = 0; c < kcols; c++)
            ker_buf[r * fft_cols + c] = kernel[r][c];

    /* Row-wise forward FFT */
    for (int r = 0; r < fft_rows; r++) {
        radix2_dit_fft(&img_buf[r * fft_cols], fft_cols, FFT_FORWARD);
        radix2_dit_fft(&ker_buf[r * fft_cols], fft_cols, FFT_FORWARD);
    }

    /* Transpose both (swap to column-major for column-wise FFT) */
    for (int r = 0; r < fft_rows; r++)
        for (int c = 0; c < fft_cols; c++)
            tmp[c * fft_rows + r] = img_buf[r * fft_cols + c];
    memcpy(img_buf, tmp, (size_t)total * sizeof(complex_t));

    for (int r = 0; r < fft_rows; r++)
        for (int c = 0; c < fft_cols; c++)
            tmp[c * fft_rows + r] = ker_buf[r * fft_cols + c];
    memcpy(ker_buf, tmp, (size_t)total * sizeof(complex_t));

    /* Column-wise forward FFT (now stored as rows after transpose) */
    for (int c = 0; c < fft_cols; c++) {
        radix2_dit_fft(&img_buf[c * fft_rows], fft_rows, FFT_FORWARD);
        radix2_dit_fft(&ker_buf[c * fft_rows], fft_rows, FFT_FORWARD);
    }

    /* Pointwise multiply */
    for (int i = 0; i < total; i++)
        img_buf[i] *= ker_buf[i];

    /* Inverse column-wise FFT */
    for (int c = 0; c < fft_cols; c++)
        radix2_dit_fft(&img_buf[c * fft_rows], fft_rows, FFT_INVERSE);

    /* Transpose back to row-major */
    for (int r = 0; r < fft_rows; r++)
        for (int c = 0; c < fft_cols; c++)
            tmp[r * fft_cols + c] = img_buf[c * fft_rows + r];
    memcpy(img_buf, tmp, (size_t)total * sizeof(complex_t));

    /* Inverse row-wise FFT */
    for (int r = 0; r < fft_rows; r++)
        radix2_dit_fft(&img_buf[r * fft_cols], fft_cols, FFT_INVERSE);

    /* Copy valid region to result */
    for (int r = 0; r < out_rows; r++)
        for (int c = 0; c < out_cols; c++)
            result[r][c] = img_buf[r * fft_cols + c];

    free_complex_array(img_buf);
    free_complex_array(ker_buf);
    free_complex_array(tmp);
}

// Demonstration functions
void demo_convolution_types() {
    printf("\nConvolution Types Demonstration:\n");
    printf("================================\n");
    
    // Create signals
    int nx = 8, nh = 4;
    complex_t* x = allocate_complex_array(nx);
    complex_t* h = allocate_complex_array(nh);
    if (!x || !h) {
        fprintf(stderr, "Error: allocation failed for convolution demo\n");
        free_complex_array(x);
        free_complex_array(h);
        return;
    }
    
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
    if (!y_linear) {
        fprintf(stderr, "Error: allocation failed for linear convolution\n");
        free_complex_array(x);
        free_complex_array(h);
        return;
    }
    fft_convolution(x, nx, h, nh, y_linear);
    
    printf("\nLinear convolution result: ");
    for (int i = 0; i < ny; i++) {
        printf("%.2f ", creal(y_linear[i]));
    }
    printf("\n");
    
    // Circular convolution of length nx needs both operands nx long; h is only
    // nh long, so zero-extend it first.
    complex_t* y_circular = allocate_complex_array(nx);
    complex_t* h_pad = allocate_complex_array(nx);
    if (!y_circular || !h_pad) {
        fprintf(stderr, "Error: allocation failed for circular convolution\n");
        free_complex_array(y_circular);
        free_complex_array(h_pad);
        free_complex_array(x);
        free_complex_array(h);
        free_complex_array(y_linear);
        return;
    }
    memcpy(h_pad, h, (size_t)nh * sizeof(complex_t));
    circular_convolution(x, h_pad, nx, y_circular);
    free_complex_array(h_pad);
    
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
        if (!x || !h || !y_direct || !y_fft) {
            fprintf(stderr, "benchmark_convolution: allocation failed\n");
            free_complex_array(x);
            free_complex_array(h);
            free_complex_array(y_direct);
            free_complex_array(y_fft);
            continue;
        }
        
        // Fill with random data
        for (int j = 0; j < nx; j++) {
            x[j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
        for (int j = 0; j < nh; j++) {
            h[j] = ((double)rand() / RAND_MAX) * 2.0 - 1.0;
        }
        
        fft_timer_t timer;

        /* Batch until clock() resolves the interval; both functions only
         * read x and h, so repeating them is safe. */
        int reps = 1;
        double total = 0.0;
        for (;;) {
            timer_start(&timer);
            for (int r = 0; r < reps; r++)
                direct_convolution(x, nx, h, nh, y_direct);
            timer_stop(&timer);
            total = timer.elapsed_ms;
            if (total >= 10.0 || reps >= (1 << 20)) break;
            reps *= 4;
        }
        double time_direct = total / reps;

        reps = 1;
        for (;;) {
            timer_start(&timer);
            for (int r = 0; r < reps; r++)
                fft_convolution(x, nx, h, nh, y_fft);
            timer_stop(&timer);
            total = timer.elapsed_ms;
            if (total >= 10.0 || reps >= (1 << 20)) break;
            reps *= 4;
        }
        double time_fft = total / reps;
        
        // Verify correctness
        double error = 0;
        for (int j = 0; j < ny; j++) {
            error += cabs(y_direct[j] - y_fft[j]);
        }
        
        printf("%d\t\t%d\t\t%.3f\t\t%.3f\t\t%.1fx",
               nx, nh, time_direct, time_fft,
               time_fft > 0.0 ? time_direct / time_fft : 0.0);
        
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
    if (!y) {
        fprintf(stderr, "convolution: output allocation failed\n");
        return 1;
    }
    
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
    printf("This allows processing of streaming data.\n");
    
    return 0;
}
