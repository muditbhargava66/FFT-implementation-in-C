#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"
#include "radix4_internal.h"

/*
 * Radix-4 decimation-in-time. Base-4 digit reversal, one radix-2 stage when
 * log2(n) is odd, then radix-4 stages, so any power-of-two length works (not
 * only powers of four). A radix-4 stage needs three twiddle multiplies per
 * four outputs where two radix-2 stages need four; the trade is an n-element
 * reorder buffer per call.
 */

/*
 * Reference illustration of the radix-4 butterfly for the demo program. Not
 * used by radix4_fft(), which pre-multiplies its inputs instead of applying
 * twiddles to the outputs.
 */
#ifndef LIB_BUILD
/* Only the demonstration main() calls this; compiling it into the library
 * leaves an unused static function behind, which is a Clang warning. */
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
#endif /* LIB_BUILD */

/*
 * A base-4 digit reversal needs a second buffer: the source and destination
 * indices overlap, and in-place digit reversal is not a permutation cycle that
 * a single pass can follow safely. The cost is one n-element allocation and
 * two passes over n elements per transform.
 */
int fft_radix4_reorder_into(complex_t* x, int n, complex_t* temp) {
    if (!x || n <= 1) return 0;
    if (!temp) return -1;

    int log2n = log2_int(n);
    int is_power_of_4 = (log2n % 2 == 0);

    for (int i = 0; i < n; i++) {
        int weight = n;
        int temp_i = i;
        int accumulated_rev = 0;

        /* For N = 2 * 4^k the leading digit is base 2; every later digit is
         * base 4. */
        if (!is_power_of_4) {
            int radix = 2;
            weight /= radix;
            int digit = temp_i % radix;
            temp_i /= radix;
            accumulated_rev += digit * weight;
        }

        while (weight > 1) {
            int radix = 4;
            weight /= radix;
            int digit = temp_i % radix;
            temp_i /= radix;
            accumulated_rev += digit * weight;
        }

        temp[i] = x[accumulated_rev];
    }

    memcpy(x, temp, (size_t)n * sizeof(complex_t));
    return 0;
}

int fft_radix4_reorder(complex_t* x, int n) {
    if (!x || n <= 1) return 0;

    complex_t* temp = (complex_t*)malloc((size_t)n * sizeof(complex_t));
    if (!temp) return -1;

    int rc = fft_radix4_reorder_into(x, n, temp);
    free(temp);
    return rc;
}

/*
 * Radix-4 DIT: base-4 digit reversal, one radix-2 stage when log2(n) is odd,
 * then radix-4 stages. The inverse scale is applied once at the end.
 */
static fft_error_t radix4_fft_core(complex_t* x, int n, fft_direction dir,
                                   complex_t* temp) {
    int log2n = log2_int(n);
    int is_power_of_4 = (log2n % 2 == 0);

    if (fft_radix4_reorder_into(x, n, temp) != 0)
        return FFT_ERR_ALLOC;

    /* log2(n) odd leaves one factor of 2. A length-2 DFT is just sum and
     * difference, with unity twiddles, so one pass over adjacent pairs
     * handles it. */
    int start_length = 4;
    
    if (!is_power_of_4) {
        int m = 2;
        
        for (int k = 0; k < n; k += m) {
            complex_t t = x[k + 1];
            complex_t u = x[k];
            x[k] = u + t;
            x[k + 1] = u - t;
        }
        /* Combine four DFTs of size 2 into size 8 at the next stage. */
        start_length = 8;
    }
    
    /* Remaining stages are radix-4. The loop advances by hand so that the
     * final m *= 4 on a 2^29 transform cannot overflow int. */
    for (int m = start_length; m <= n;) {
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
                /* Upper three branches carry the twiddle for this j; the
                 * combine below is the 4-point DFT with t3 rotated by -j
                 * (forward) or +j (inverse). */
                int i0 = k + j;
                int i1 = i0 + m4;
                int i2 = i0 + m2;
                int i3 = i0 + m3;
                
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

        if (m == n) break;
        m *= 4;
    }
    
    /* Scale for inverse transform. */    if (dir == FFT_INVERSE) {
        double scale = 1.0 / n;
        for (int i = 0; i < n; i++) {
            x[i] *= scale;
        }
    }
    return FFT_OK;
}

fft_error_t radix4_fft_scratch(complex_t* x, int n, fft_direction dir,
                               complex_t* temp) {
    if (!x || !temp) return FFT_ERR_NULL;
    if (!is_power_of_two(n)) return FFT_ERR_SIZE;
    if (n <= 1) return FFT_OK;
    return radix4_fft_core(x, n, dir, temp);
}

fft_error_t radix4_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    if (!is_power_of_two(n)) return FFT_ERR_SIZE;
    if (n <= 1) return FFT_OK;

    complex_t* temp = (complex_t*)malloc((size_t)n * sizeof(complex_t));
    if (!temp) return FFT_ERR_ALLOC;

    fft_error_t err = radix4_fft_core(x, n, dir, temp);
    free(temp);
    return err;
}

/* Convenience wrappers. Any power-of-two length is valid; powers of four skip
 * the extra radix-2 stage. */
fft_error_t fft_radix4(complex_t* x, int n) {
    return radix4_fft(x, n, FFT_FORWARD);
}

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
