#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Mixed-Radix FFT Implementation
 * 
 * Handles FFT for composite lengths (e.g., 6, 10, 12, 15, 20, etc.)
 * by factoring the length and using smaller DFTs.
 * 
 * Uses the Cooley-Tukey factorization:
 * - If N = N1 * N2, compute N2 DFTs of size N1, then N1 DFTs of size N2
 * - Supports any factorization, optimized for small prime factors
 * 
 * Time Complexity: O(n log n) for highly composite n
 */

// Structure to hold factorization
typedef struct {
    int* factors;
    int num_factors;
} factorization_t;

// Factor a number into prime factors (simple trial division)
static factorization_t factorize(int n) {
    factorization_t result;
    result.factors = (int*)malloc(32 * sizeof(int));  // Max 32 factors
    result.num_factors = 0;
    
    // Factor out 2s
    while (n % 2 == 0) {
        result.factors[result.num_factors++] = 2;
        n /= 2;
    }
    
    // Factor out odd numbers
    for (int i = 3; i * i <= n; i += 2) {
        while (n % i == 0) {
            result.factors[result.num_factors++] = i;
            n /= i;
        }
    }
    
    // If n is prime > 2
    if (n > 2) {
        result.factors[result.num_factors++] = n;
    }
    
    return result;
}

// Small-N DFT implementations
static void dft_2(complex_t* x, int stride, fft_direction dir) {
    complex_t t = x[0];
    x[0] = t + x[stride];
    x[stride] = t - x[stride];
}

static void dft_3(complex_t* x, int stride, fft_direction dir) {
    complex_t w = cexp(dir * I * TWO_PI / 3);
    complex_t w2 = w * w;
    
    complex_t t0 = x[0];
    complex_t t1 = x[stride];
    complex_t t2 = x[2 * stride];
    
    x[0] = t0 + t1 + t2;
    x[stride] = t0 + w * t1 + w2 * t2;
    x[2 * stride] = t0 + w2 * t1 + w * t2;
}

static void dft_5(complex_t* x, int stride, fft_direction dir) {
    complex_t w = cexp(dir * I * TWO_PI / 5);
    complex_t w2 = w * w;
    complex_t w3 = w2 * w;
    complex_t w4 = w3 * w;
    
    complex_t t0 = x[0];
    complex_t t1 = x[stride];
    complex_t t2 = x[2 * stride];
    complex_t t3 = x[3 * stride];
    complex_t t4 = x[4 * stride];
    
    x[0] = t0 + t1 + t2 + t3 + t4;
    x[stride] = t0 + w * t1 + w2 * t2 + w3 * t3 + w4 * t4;
    x[2 * stride] = t0 + w2 * t1 + w4 * t2 + w * t3 + w3 * t4;
    x[3 * stride] = t0 + w3 * t1 + w * t2 + w4 * t3 + w2 * t4;
    x[4 * stride] = t0 + w4 * t1 + w3 * t2 + w2 * t3 + w * t4;
}

// General small-N DFT (for prime factors)
static void dft_general(complex_t* x, int n, int stride, fft_direction dir) {
    complex_t* temp = allocate_complex_array(n);
    
    // Copy strided input
    for (int i = 0; i < n; i++) {
        temp[i] = x[i * stride];
    }
    
    // Compute DFT
    for (int k = 0; k < n; k++) {
        x[k * stride] = 0;
        for (int j = 0; j < n; j++) {
            x[k * stride] += temp[j] * cexp(dir * I * TWO_PI * j * k / n);
        }
    }
    
    free_complex_array(temp);
}

// Cooley-Tukey decomposition for N = N1 * N2
static void cooley_tukey_decompose(complex_t* x, int n, int n1, int n2, 
                                  int stride, fft_direction dir) {
    // Perform N2 DFTs of size N1
    for (int j = 0; j < n2; j++) {
        mixed_radix_fft_recursive(x + j * stride, n1, n2 * stride, dir);
    }
    
    // Multiply by twiddle factors
    for (int k1 = 0; k1 < n1; k1++) {
        for (int k2 = 1; k2 < n2; k2++) {
            int idx = (k1 * n2 + k2) * stride;
            x[idx] *= twiddle_factor(k1 * k2, n, dir);
        }
    }
    
    // Perform N1 DFTs of size N2
    for (int k1 = 0; k1 < n1; k1++) {
        mixed_radix_fft_recursive(x + k1 * n2 * stride, n2, stride, dir);
    }
}

// Recursive mixed-radix FFT
void mixed_radix_fft_recursive(complex_t* x, int n, int stride, fft_direction dir) {
    if (n == 1) return;
    
    // Use optimized small-N DFTs
    switch (n) {
        case 2:
            dft_2(x, stride, dir);
            return;
        case 3:
            dft_3(x, stride, dir);
            return;
        case 4:
            // Use radix-2
            mixed_radix_fft_recursive(x, 2, stride, dir);
            mixed_radix_fft_recursive(x + 2 * stride, 2, stride, dir);
            dft_2(x, 2 * stride, dir);
            dft_2(x + stride, 2 * stride, dir);
            x[stride] *= twiddle_factor(1, 4, dir);
            x[3 * stride] *= twiddle_factor(3, 4, dir);
            return;
        case 5:
            dft_5(x, stride, dir);
            return;
    }
    
    // Factor n and use Cooley-Tukey
    factorization_t factors = factorize(n);
    
    if (factors.num_factors == 1) {
        // n is prime, use general DFT
        dft_general(x, n, stride, dir);
    } else {
        // Find best factorization (prefer balanced factors)
        int n1 = 1, n2 = n;
        
        for (int i = 0; i < factors.num_factors / 2; i++) {
            n1 *= factors.factors[i];
        }
        n2 = n / n1;
        
        // Apply Cooley-Tukey decomposition
        cooley_tukey_decompose(x, n, n1, n2, stride, dir);
    }
    
    free(factors.factors);
}

// Main mixed-radix FFT function
void mixed_radix_fft(complex_t* x, int n, fft_direction dir) {
    // Create working array for reordering
    complex_t* work = allocate_complex_array(n);
    memcpy(work, x, n * sizeof(complex_t));
    
    // Perform mixed-radix FFT
    mixed_radix_fft_recursive(work, n, 1, dir);
    
    // Copy result back
    memcpy(x, work, n * sizeof(complex_t));
    free_complex_array(work);
    
    // Scale for inverse FFT
    if (dir == FFT_INVERSE) {
        for (int i = 0; i < n; i++) {
            x[i] /= n;
        }
    }
}

// Wrapper functions
void fft_mixed_radix(complex_t* x, int n) {
    mixed_radix_fft(x, n, FFT_FORWARD);
}

void ifft_mixed_radix(complex_t* x, int n) {
    mixed_radix_fft(x, n, FFT_INVERSE);
}

// Test and demonstration
int main() {
    printf("Mixed-Radix FFT Implementation\n");
    printf("==============================\n\n");
    
    // Test various composite sizes
    int test_sizes[] = {6, 10, 12, 15, 20, 24, 30, 35, 40, 42, 
                       48, 60, 72, 84, 90, 100, 120, 144};
    int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);
    
    printf("Testing composite-length FFTs:\n");
    printf("Size\tFactors\t\tMixed-Radix\tNaive DFT\tError\n");
    printf("----\t-------\t\t-----------\t---------\t-----\n");
    
    for (int t = 0; t < num_tests; t++) {
        int n = test_sizes[t];
        
        // Show factorization
        printf("%d\t", n);
        factorization_t factors = factorize(n);
        for (int i = 0; i < factors.num_factors; i++) {
            printf("%d", factors.factors[i]);
            if (i < factors.num_factors - 1) printf("×");
        }
        printf("\t\t");
        if (factors.num_factors < 3) printf("\t");
        free(factors.factors);
        
        // Prepare test data
        complex_t* x_mixed = allocate_complex_array(n);
        complex_t* x_naive = allocate_complex_array(n);
        complex_t* original = allocate_complex_array(n);
        
        // Generate test signal
        for (int i = 0; i < n; i++) {
            double value = sin(2 * PI * 2 * i / n) + 
                          0.3 * cos(2 * PI * 5 * i / n);
            x_mixed[i] = x_naive[i] = original[i] = value;
        }
        
        // Time mixed-radix FFT
        timer_t timer;
        timer_start(&timer);
        fft_mixed_radix(x_mixed, n);
        timer_stop(&timer);
        double time_mixed = timer.elapsed_ms;
        
        // Time naive DFT
        timer_start(&timer);
        naive_dft(x_naive, n, FFT_FORWARD);
        timer_stop(&timer);
        double time_naive = timer.elapsed_ms;
        
        // Verify correctness
        double max_error = 0;
        for (int i = 0; i < n; i++) {
            double error = cabs(x_mixed[i] - x_naive[i]);
            if (error > max_error) max_error = error;
        }
        
        printf("%.3f ms\t%.3f ms\t%.2e\n", 
               time_mixed, time_naive, max_error);
        
        free_complex_array(x_mixed);
        free_complex_array(x_naive);
        free_complex_array(original);
    }
    
    // Demonstrate specific examples
    printf("\n\nDetailed example: 12-point FFT (12 = 3 × 4)\n");
    printf("-------------------------------------------\n");
    
    int n = 12;
    complex_t* signal = allocate_complex_array(n);
    
    // Create simple test signal
    for (int i = 0; i < n; i++) {
        signal[i] = i;
    }
    
    printf("Input: ");
    for (int i = 0; i < n; i++) {
        printf("%.0f ", creal(signal[i]));
    }
    printf("\n");
    
    // Compute FFT
    fft_mixed_radix(signal, n);
    
    printf("\nFFT output:\n");
    for (int i = 0; i < n; i++) {
        printf("X[%2d] = ", i);
        print_complex(signal[i]);
        printf("\n");
    }
    
    // Test inverse
    ifft_mixed_radix(signal, n);
    
    printf("\nAfter inverse FFT: ");
    for (int i = 0; i < n; i++) {
        printf("%.0f ", creal(signal[i]));
    }
    printf("\n");
    
    free_complex_array(signal);
    
    return 0;
}
