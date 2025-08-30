#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"
#include <string.h>

/**
 * 2D FFT for Image Processing
 * 
 * Implements 2D FFT using row-column decomposition method
 * for image processing applications.
 * 
 * Features:
 * - 2D forward and inverse FFT
 * - Image filtering in frequency domain
 * - Frequency domain visualization
 * - Common image processing operations
 */

// 2D complex array allocation
complex_t** allocate_2d_complex(int rows, int cols) {
    complex_t** array = (complex_t**)malloc(rows * sizeof(complex_t*));
    for (int i = 0; i < rows; i++) {
        array[i] = allocate_complex_array(cols);
    }
    return array;
}

void free_2d_complex(complex_t** array, int rows) {
    for (int i = 0; i < rows; i++) {
        free_complex_array(array[i]);
    }
    free(array);
}

// 2D FFT using row-column decomposition
void fft_2d(complex_t** data, int rows, int cols, fft_direction dir) {
    // Check if dimensions are powers of 2
    CHECK_POWER_OF_TWO(rows);
    CHECK_POWER_OF_TWO(cols);
    
    // Transform rows
    for (int i = 0; i < rows; i++) {
        radix2_dit_fft(data[i], cols, dir);
    }
    
    // Transform columns
    complex_t* column = allocate_complex_array(rows);
    for (int j = 0; j < cols; j++) {
        // Extract column
        for (int i = 0; i < rows; i++) {
            column[i] = data[i][j];
        }
        
        // FFT of column
        radix2_dit_fft(column, rows, dir);
        
        // Put back
        for (int i = 0; i < rows; i++) {
            data[i][j] = column[i];
        }
    }
    free_complex_array(column);
    
    // Scale for inverse transform
    if (dir == FFT_INVERSE) {
        double scale = 1.0 / (rows * cols);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                data[i][j] *= scale;
            }
        }
    }
}

// Shift zero frequency to center (fftshift)
void fft_shift_2d(complex_t** data, int rows, int cols) {
    int half_rows = rows / 2;
    int half_cols = cols / 2;
    
    complex_t** temp = allocate_2d_complex(rows, cols);
    
    // Copy to temp with shift
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            int new_i = (i + half_rows) % rows;
            int new_j = (j + half_cols) % cols;
            temp[new_i][new_j] = data[i][j];
        }
    }
    
    // Copy back
    for (int i = 0; i < rows; i++) {
        memcpy(data[i], temp[i], cols * sizeof(complex_t));
    }
    
    free_2d_complex(temp, rows);
}

// Generate 2D test patterns
void generate_2d_sinusoid(complex_t** data, int rows, int cols, 
                         double fx, double fy) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double phase = TWO_PI * (fx * i / rows + fy * j / cols);
            data[i][j] = cos(phase) + I * sin(phase);
        }
    }
}

void generate_2d_gaussian(complex_t** data, int rows, int cols, 
                         double sigma) {
    double cx = rows / 2.0;
    double cy = cols / 2.0;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double dx = i - cx;
            double dy = j - cy;
            double r2 = dx * dx + dy * dy;
            data[i][j] = exp(-r2 / (2 * sigma * sigma));
        }
    }
}

void generate_rectangle(complex_t** data, int rows, int cols,
                       int rect_height, int rect_width) {
    // Clear array
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            data[i][j] = 0;
        }
    }
    
    // Draw rectangle in center
    int start_row = (rows - rect_height) / 2;
    int start_col = (cols - rect_width) / 2;
    
    for (int i = 0; i < rect_height; i++) {
        for (int j = 0; j < rect_width; j++) {
            if (start_row + i < rows && start_col + j < cols) {
                data[start_row + i][start_col + j] = 1.0;
            }
        }
    }
}

// 2D filters
void apply_ideal_lowpass_2d(complex_t** spectrum, int rows, int cols, double cutoff) {
    double cx = rows / 2.0;
    double cy = cols / 2.0;
    double cutoff2 = cutoff * cutoff;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double dx = i - cx;
            double dy = j - cy;
            double r2 = dx * dx + dy * dy;
            
            if (r2 > cutoff2) {
                spectrum[i][j] = 0;
            }
        }
    }
}

void apply_gaussian_filter_2d(complex_t** spectrum, int rows, int cols, double sigma) {
    double cx = rows / 2.0;
    double cy = cols / 2.0;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double dx = i - cx;
            double dy = j - cy;
            double r2 = dx * dx + dy * dy;
            double filter = exp(-r2 / (2 * sigma * sigma));
            spectrum[i][j] *= filter;
        }
    }
}

// Visualization helpers
void display_magnitude_2d(complex_t** data, int rows, int cols, const char* title) {
    printf("\n%s (magnitude, center region):\n", title);
    printf("==================================\n");
    
    // Find max magnitude for scaling
    double max_mag = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double mag = cabs(data[i][j]);
            if (mag > max_mag) max_mag = mag;
        }
    }
    
    // Display center region
    int display_size = 16;
    if (display_size > rows) display_size = rows;
    if (display_size > cols) display_size = cols;
    
    int start_row = (rows - display_size) / 2;
    int start_col = (cols - display_size) / 2;
    
    for (int i = 0; i < display_size; i++) {
        for (int j = 0; j < display_size; j++) {
            double mag = cabs(data[start_row + i][start_col + j]);
            int level = (int)(9 * mag / max_mag);
            if (level > 9) level = 9;
            printf("%c ", '0' + level);
        }
        printf("\n");
    }
}

// Image processing operations
void edge_detection_2d(complex_t** image, int rows, int cols) {
    // High-pass filter in frequency domain
    fft_2d(image, rows, cols, FFT_FORWARD);
    fft_shift_2d(image, rows, cols);
    
    // Apply high-pass filter (subtract low frequencies)
    double cutoff = fmin(rows, cols) * 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double dx = i - rows/2.0;
            double dy = j - cols/2.0;
            double r = sqrt(dx*dx + dy*dy);
            
            if (r < cutoff) {
                image[i][j] = 0;
            }
        }
    }
    
    fft_shift_2d(image, rows, cols);
    fft_2d(image, rows, cols, FFT_INVERSE);
}

// Demonstration
void demo_2d_fft_properties() {
    int rows = 32, cols = 32;
    complex_t** image = allocate_2d_complex(rows, cols);
    complex_t** spectrum = allocate_2d_complex(rows, cols);
    
    printf("\n2D FFT Properties Demonstration:\n");
    printf("================================\n");
    
    // 1. Rectangle -> Sinc pattern
    printf("\n1. Rectangle Transform:\n");
    generate_rectangle(image, rows, cols, 8, 8);
    
    // Copy for FFT
    for (int i = 0; i < rows; i++) {
        memcpy(spectrum[i], image[i], cols * sizeof(complex_t));
    }
    
    fft_2d(spectrum, rows, cols, FFT_FORWARD);
    fft_shift_2d(spectrum, rows, cols);
    
    display_magnitude_2d(image, rows, cols, "Rectangle (spatial)");
    display_magnitude_2d(spectrum, rows, cols, "Sinc pattern (frequency)");
    
    // 2. Gaussian -> Gaussian
    printf("\n2. Gaussian Transform:\n");
    generate_2d_gaussian(image, rows, cols, 4.0);
    
    for (int i = 0; i < rows; i++) {
        memcpy(spectrum[i], image[i], cols * sizeof(complex_t));
    }
    
    fft_2d(spectrum, rows, cols, FFT_FORWARD);
    fft_shift_2d(spectrum, rows, cols);
    
    printf("Gaussian remains Gaussian in frequency domain\n");
    
    // 3. Rotation property
    printf("\n3. Rotation Property:\n");
    printf("Rotating image by θ rotates spectrum by θ\n");
    
    free_2d_complex(image, rows);
    free_2d_complex(spectrum, rows);
}

// Main program
int main() {
    printf("2D FFT for Image Processing\n");
    printf("===========================\n");
    
    // Test basic 2D FFT
    int rows = 16, cols = 16;
    complex_t** test_image = allocate_2d_complex(rows, cols);
    
    // Create simple test pattern
    printf("\nTest 1: Simple 2D sinusoid\n");
    printf("--------------------------\n");
    
    generate_2d_sinusoid(test_image, rows, cols, 2, 3);
    
    // Forward 2D FFT
    complex_t** spectrum = allocate_2d_complex(rows, cols);
    for (int i = 0; i < rows; i++) {
        memcpy(spectrum[i], test_image[i], cols * sizeof(complex_t));
    }
    
    fft_timer_t timer;
    timer_start(&timer);
    fft_2d(spectrum, rows, cols, FFT_FORWARD);
    timer_stop(&timer);
    
    printf("2D FFT (%dx%d) completed in %.3f ms\n", rows, cols, timer.elapsed_ms);
    
    // Check for peaks at expected frequencies
    fft_shift_2d(spectrum, rows, cols);
    printf("\nSpectrum peaks at (row, col):\n");
    
    double threshold = 0.1;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double mag = cabs(spectrum[i][j]);
            if (mag > threshold) {
                printf("(%d, %d): magnitude = %.2f\n", i, j, mag);
            }
        }
    }
    
    // Test inverse transform
    fft_shift_2d(spectrum, rows, cols);
    fft_2d(spectrum, rows, cols, FFT_INVERSE);
    
    // Verify reconstruction
    double error = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            error += cabs(spectrum[i][j] - test_image[i][j]);
        }
    }
    printf("\nReconstruction error: %.2e\n", error / (rows * cols));
    
    // Demonstrate properties
    demo_2d_fft_properties();
    
    // Image filtering example
    printf("\n\nImage Filtering Example:\n");
    printf("========================\n");
    
    // Create test image with noise
    rows = cols = 64;
    complex_t** noisy_image = allocate_2d_complex(rows, cols);
    
    // Add low frequency signal + high frequency noise
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            double signal = sin(TWO_PI * 2 * i / rows) * sin(TWO_PI * 2 * j / cols);
            double noise = 0.5 * ((double)rand() / RAND_MAX - 0.5);
            noisy_image[i][j] = signal + noise;
        }
    }
    
    printf("Created %dx%d image with signal + noise\n", rows, cols);
    
    // Apply low-pass filter
    complex_t** filtered = allocate_2d_complex(rows, cols);
    for (int i = 0; i < rows; i++) {
        memcpy(filtered[i], noisy_image[i], cols * sizeof(complex_t));
    }
    
    fft_2d(filtered, rows, cols, FFT_FORWARD);
    fft_shift_2d(filtered, rows, cols);
    apply_gaussian_filter_2d(filtered, rows, cols, 10.0);
    fft_shift_2d(filtered, rows, cols);
    fft_2d(filtered, rows, cols, FFT_INVERSE);
    
    printf("Applied Gaussian low-pass filter (σ = 10)\n");
    
    // Performance analysis
    printf("\n\nPerformance Analysis:\n");
    printf("====================\n");
    printf("Size\t\tTime (ms)\tPixels/ms\n");
    printf("----\t\t---------\t---------\n");
    
    int test_sizes[] = {32, 64, 128, 256};
    for (int s = 0; s < 4; s++) {
        int size = test_sizes[s];
        complex_t** test = allocate_2d_complex(size, size);
        
        // Fill with random data
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                test[i][j] = (double)rand() / RAND_MAX;
            }
        }
        
        timer_start(&timer);
        fft_2d(test, size, size, FFT_FORWARD);
        timer_stop(&timer);
        
        printf("%dx%d\t\t%.3f\t\t%.0f\n", 
               size, size, timer.elapsed_ms, 
               (size * size) / timer.elapsed_ms);
        
        free_2d_complex(test, size);
    }
    
    // Applications
    printf("\n\nApplications of 2D FFT:\n");
    printf("======================\n");
    printf("1. Image filtering (blur, sharpen, denoise)\n");
    printf("2. Image compression (JPEG uses DCT)\n");
    printf("3. Pattern recognition\n");
    printf("4. Image registration and alignment\n");
    printf("5. Texture analysis\n");
    printf("6. Optical flow computation\n");
    printf("7. Medical imaging (MRI reconstruction)\n");
    
    // Cleanup
    free_2d_complex(test_image, 16);
    free_2d_complex(spectrum, 16);
    free_2d_complex(noisy_image, rows);
    free_2d_complex(filtered, rows);
    
    return 0;
}
