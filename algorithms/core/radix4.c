#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * @file radix4.c
 * @brief True Radix-4 Decimation-in-Time (DIT) FFT Implementation
 *
 * @details
 * This file implements a true radix-4 Cooley-Tukey DIT FFT.
 *
 * The radix-4 butterfly combines four inputs into four outputs:
 *   X0 = a + b + c + d
 *   X1 = a - j*b - c + j*d
 *   X2 = a - b + c - d
 *   X3 = a + j*b - c - j*d
 *
 * where b, c, and d already include their twiddle factors.
 *
 * Features:
 * - True radix-4 butterfly implementation
 * - Base-4 digit-reversal permutation
 * - In-place computation
 * - Precomputed twiddle tables via execution plans
 *
 * Limitations:
 * - Input size must be a power of 4: 4, 16, 64, 256, ...
 *
 * Time Complexity: O(n log n)
 * Space Complexity: O(1) for in-place computation, plus O(n) for the plan twiddle table
 */

/* =========================================================
 * Utilities
 * ========================================================= */

/**
 * @brief Check whether n is a power of 4.
 *
 * @param n Input size
 * @return 1 if n is a power of 4, 0 otherwise
 */
static inline int is_power_of_four(int n) {
    if (n < 1) {
        return 0;
    }

    while ((n % 4) == 0) {
        n /= 4;
    }

    return (n == 1);
}

/**
 * @brief Compute integer log base 4.
 *
 * @param n Input size (must be a power of 4)
 * @return log4(n)
 */
static inline int log4_int(int n) {
    int p = 0;

    while (n > 1) {
        n /= 4;
        p++;
    }

    return p;
}

/**
 * @brief Reverse the base-4 digits of an index.
 *
 * @param x Index to reverse
 * @param digits Number of base-4 digits
 * @return Base-4 digit-reversed index
 */
static inline int digit_reverse_base4(int x, int digits) {
    int rev = 0;

    for (int i = 0; i < digits; i++) {
        rev = (rev << 2) | (x & 0x3);
        x >>= 2;
    }

    return rev;
}

/* =========================================================
 * Radix-4 Butterfly
 * ========================================================= */

/**
 * @brief Apply a radix-4 butterfly.
 *
 * @details
 * Inputs b, c, and d are assumed to already include their twiddle factors.
 *
 * Forward FFT:
 *   X0 = a + b + c + d
 *   X1 = a - j*b - c + j*d
 *   X2 = a - b + c - d
 *   X3 = a + j*b - c - j*d
 *
 * Inverse FFT uses the conjugate sign convention.
 *
 * @param a Pointer to first value
 * @param b Pointer to second value
 * @param c Pointer to third value
 * @param d Pointer to fourth value
 * @param dir Transform direction
 */
static inline void radix4_butterfly(complex_t* a, complex_t* b, complex_t* c, complex_t* d,
                                    fft_direction dir) {
    complex_t a0 = *a;
    complex_t b0 = *b;
    complex_t c0 = *c;
    complex_t d0 = *d;

    complex_t t0 = a0 + c0;
    complex_t t1 = a0 - c0;
    complex_t t2 = b0 + d0;

    complex_t j_factor = (dir == FFT_FORWARD) ? (-I) : (I);
    complex_t t3 = (b0 - d0) * j_factor;

    *a = t0 + t2;
    *b = t1 + t3;
    *c = t0 - t2;
    *d = t1 - t3;
}

/* =========================================================
 * Plan
 * ========================================================= */

typedef struct {
    int n;
    fft_direction dir;
    complex_t* W;
} radix4_plan_t;

/**
 * @brief Initialize a radix-4 execution plan.
 *
 * @param plan Plan to initialize
 * @param n FFT size
 * @param dir Transform direction
 * @return 1 on success, 0 on failure
 */
static int radix4_plan_init(radix4_plan_t* plan, int n, fft_direction dir) {
    if (!plan) {
        return 0;
    }

    if (!is_power_of_four(n)) {
        return 0;
    }

    plan->n = n;
    plan->dir = dir;
    plan->W = allocate_complex_array(n);

    if (!plan->W) {
        return 0;
    }

    for (int k = 0; k < n; k++) {
        double angle = dir * TWO_PI * k / n;
        plan->W[k] = cos(angle) + I * sin(angle);
    }

    return 1;
}

/**
 * @brief Free a radix-4 execution plan.
 *
 * @param plan Plan to free
 */
static void radix4_plan_free(radix4_plan_t* plan) {
    if (!plan) {
        return;
    }

    if (plan->W) {
        free_complex_array(plan->W);
        plan->W = NULL;
    }
}

/* =========================================================
 * Core FFT using a precomputed plan
 * ========================================================= */

/**
 * @brief Execute radix-4 FFT using a precomputed plan.
 *
 * @param x Input/output array
 * @param plan Precomputed radix-4 plan
 */
static void radix4_fft_with_plan(complex_t* x, const radix4_plan_t* plan) {
    if (!x || !plan || !plan->W) {
        fprintf(stderr, "Error: Invalid arguments passed to radix4_fft_with_plan\n");
        return;
    }

    int n = plan->n;
    fft_direction dir = plan->dir;
    const complex_t* W = plan->W;

    if (n <= 1) {
        return;
    }

    int log4n = log4_int(n);

    /* Step 1: Base-4 digit-reversal permutation */
    for (int i = 0; i < n; i++) {
        int j = digit_reverse_base4(i, log4n);
        if (i < j) {
            complex_t tmp = x[i];
            x[i] = x[j];
            x[j] = tmp;
        }
    }

    /* Step 2: Radix-4 iterative stages */
    int m = 4;
    for (int stage = 1; stage <= log4n; stage++, m *= 4) {
        int quarter = m >> 2;
        int step = n / m;

        /* Fast path for first stage: all twiddles are 1 */
        if (m == 4) {
            for (int k = 0; k < n; k += 4) {
                complex_t a = x[k];
                complex_t b = x[k + 1];
                complex_t c = x[k + 2];
                complex_t d = x[k + 3];

                radix4_butterfly(&a, &b, &c, &d, dir);

                x[k]     = a;
                x[k + 1] = b;
                x[k + 2] = c;
                x[k + 3] = d;
            }
            continue;
        }

        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < quarter; j++) {
                int i0 = k + j;
                int i1 = i0 + quarter;
                int i2 = i1 + quarter;
                int i3 = i2 + quarter;

                int idx1 = j * step;
                int idx2 = 2 * j * step;
                int idx3 = 3 * j * step;

                complex_t a = x[i0];
                complex_t b = x[i1] * W[idx1];
                complex_t c = x[i2] * W[idx2];
                complex_t d = x[i3] * W[idx3];

                radix4_butterfly(&a, &b, &c, &d, dir);

                x[i0] = a;
                x[i1] = b;
                x[i2] = c;
                x[i3] = d;
            }
        }
    }

    /* Step 3: Scale for inverse FFT */
    if (dir == FFT_INVERSE) {
        double scale = 1.0 / (double)n;
        for (int i = 0; i < n; i++) {
            x[i] *= scale;
        }
    }
}

/* =========================================================
 * Public API
 * ========================================================= */

/**
 * @brief Main radix-4 FFT implementation.
 *
 * @param x Input/output array of complex numbers
 * @param n Array length (must be power of 4)
 * @param dir Transform direction (FFT_FORWARD or FFT_INVERSE)
 */
void radix4_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) {
        fprintf(stderr, "Error: Null pointer passed to radix4_fft\n");
        return;
    }

    if (!is_power_of_four(n)) {
        fprintf(stderr, "Error: Radix-4 FFT requires size to be a power of 4\n");
        return;
    }

    radix4_plan_t plan;
    if (!radix4_plan_init(&plan, n, dir)) {
        fprintf(stderr, "Error: Failed to initialize radix-4 plan\n");
        return;
    }

    radix4_fft_with_plan(x, &plan);
    radix4_plan_free(&plan);
}

/**
 * @brief Compute forward FFT using radix-4.
 *
 * @param x Input/output array
 * @param n Array length (must be power of 4)
 */
void fft_radix4(complex_t* x, int n) {
    radix4_fft(x, n, FFT_FORWARD);
}

/**
 * @brief Compute inverse FFT using radix-4.
 *
 * @param x Input/output array
 * @param n Array length (must be power of 4)
 */
void ifft_radix4(complex_t* x, int n) {
    radix4_fft(x, n, FFT_INVERSE);
}

#ifndef LIB_BUILD

/* =========================================================
 * Demo helpers
 * ========================================================= */

/**
 * @brief Demonstrate radix-4 butterfly behavior.
 */
static void demonstrate_butterfly(void) {
    printf("\nRadix-4 Butterfly Demonstration:\n");
    printf("================================\n");

    complex_t a = 1.0 + 0.0 * I;
    complex_t b = 0.0 + 1.0 * I;
    complex_t c = -1.0 + 0.0 * I;
    complex_t d = 0.0 - 1.0 * I;

    printf("Input values:\n");
    printf("  a = "); print_complex(a); printf("\n");
    printf("  b = "); print_complex(b); printf(" = j\n");
    printf("  c = "); print_complex(c); printf("\n");
    printf("  d = "); print_complex(d); printf(" = -j\n");

    radix4_butterfly(&a, &b, &c, &d, FFT_FORWARD);

    printf("\nOutput values (forward FFT butterfly):\n");
    printf("  A = "); print_complex(a); printf("\n");
    printf("  B = "); print_complex(b); printf("\n");
    printf("  C = "); print_complex(c); printf("\n");
    printf("  D = "); print_complex(d); printf("\n");
}

/**
 * @brief Analyze implementation-level multiplication counts.
 *
 * @details
 * This counts the complex multiplications executed by the FFT core code:
 *
 * Radix-2 implementation (radix2_dit.c):
 * - temp = x[u] * w
 * - w *= w_m
 *
 * Radix-4 implementation (this file, plan-based execution only):
 * - b = x[i1] * W[idx1]
 * - c = x[i2] * W[idx2]
 * - d = x[i3] * W[idx3]
 * - t3 = (b0 - d0) * j_factor
 *
 * It excludes:
 * - plan creation cost
 * - digit-reversal swaps
 * - additions/subtractions
 * - memory traffic
 */
static void analyze_operation_count(void) {
    printf("\n\nImplementation-Level Multiplication Count Analysis:\n");
    printf("===================================================\n");
    printf("N\tRadix-2 code muls\tRadix-4 code muls\tSavings\n");
    printf("--------\t-----------------\t-----------------\t-------\n");

    for (int n = 16; n <= 65536; n *= 4) {
        int log2n = log2_int(n);
        int log4n = log4_int(n);

        /* radix2_dit.c:
         * - x[u] * w : (n/2) * log2(n)
         * - w *= w_m : (n/2) * log2(n)
         */
        long long radix2_code_muls = (long long)n * log2n;

        /* radix4.c core (excluding plan creation):
         * - later stages only: 3 multiplications per butterfly
         *   -> 3 * (n/4) * (log4(n) - 1)
         * - every stage butterfly: (b0 - d0) * j_factor
         *   -> (n/4) * log4(n)
         */
        long long radix4_code_muls =
            3LL * (n / 4) * (log4n - 1) +
            1LL * (n / 4) * log4n;

        double savings =
            (1.0 - (double)radix4_code_muls / (double)radix2_code_muls) * 100.0;

        printf("%d\t\t%lld\t\t\t%lld\t\t\t%.1f%%\n",
               n, radix2_code_muls, radix4_code_muls, savings);
    }

    printf("\nThese counts reflect the FFT core code as written.\n");
    printf("They exclude plan construction, swaps, additions, and memory access costs.\n");
}

/**
 * @brief Choose benchmark repeat count based on size.
 *
 * @param n FFT size
 * @return Number of repeats
 */
static int choose_repeats(int n) {
    if (n <= 64) return 10000;
    if (n <= 256) return 5000;
    if (n <= 1024) return 2000;
    if (n <= 4096) return 500;
    if (n <= 16384) return 100;
    return 20;
}

/**
 * @brief Main test and benchmark program.
 */
int main(void) {
    printf("True Radix-4 FFT Implementation (with plan)\n");
    printf("===========================================\n");

    int test_sizes[] = {16, 64, 256, 1024, 4096, 16384, 65536};
    int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);

    printf("\nPerformance Comparison:\n");
    printf("======================\n");

    for (int t = 0; t < num_tests; t++) {
        int n = test_sizes[t];
        int repeats = choose_repeats(n);

        printf("\nTesting size n = %d (repeats = %d):\n", n, repeats);

        complex_t* signal_r2 = allocate_complex_array(n);
        complex_t* signal_r4 = allocate_complex_array(n);
        complex_t* original  = allocate_complex_array(n);

        if (!signal_r2 || !signal_r4 || !original) {
            fprintf(stderr, "Error: Memory allocation failed for n = %d\n", n);
            free_complex_array(signal_r2);
            free_complex_array(signal_r4);
            free_complex_array(original);
            return 1;
        }

        for (int i = 0; i < n; i++) {
            double value = sin(2 * PI * 10 * i / n)
                         + 0.5 * sin(2 * PI * 25 * i / n)
                         + 0.25 * sin(2 * PI * 40 * i / n);

            signal_r2[i] = value;
            signal_r4[i] = value;
            original[i]  = value;
        }

        fft_timer_t timer_r2, timer_r4;

        radix4_plan_t plan_fwd, plan_inv;
        if (!radix4_plan_init(&plan_fwd, n, FFT_FORWARD)) {
            fprintf(stderr, "Error: Failed to initialize radix-4 forward plan for n = %d\n", n);
            free_complex_array(signal_r2);
            free_complex_array(signal_r4);
            free_complex_array(original);
            return 1;
        }

        if (!radix4_plan_init(&plan_inv, n, FFT_INVERSE)) {
            fprintf(stderr, "Error: Failed to initialize radix-4 inverse plan for n = %d\n", n);
            radix4_plan_free(&plan_fwd);
            free_complex_array(signal_r2);
            free_complex_array(signal_r4);
            free_complex_array(original);
            return 1;
        }

        /* Benchmark radix-4 FFT only */
        timer_start(&timer_r4);
        for (int r = 0; r < repeats; r++) {
            memcpy(signal_r4, original, n * sizeof(complex_t));
            radix4_fft_with_plan(signal_r4, &plan_fwd);
        }
        timer_stop(&timer_r4);

        /* Benchmark radix-2 FFT only */
        timer_start(&timer_r2);
        for (int r = 0; r < repeats; r++) {
            memcpy(signal_r2, original, n * sizeof(complex_t));
            radix2_dit_fft(signal_r2, n, FFT_FORWARD);
        }
        timer_stop(&timer_r2);

        /* Recompute once for correctness comparison */
        memcpy(signal_r2, original, n * sizeof(complex_t));
        memcpy(signal_r4, original, n * sizeof(complex_t));

        radix2_dit_fft(signal_r2, n, FFT_FORWARD);
        radix4_fft_with_plan(signal_r4, &plan_fwd);

        printf("  Radix-2 avg time: %.6f ms\n", timer_r2.elapsed_ms / repeats);
        printf("  Radix-4 avg time: %.6f ms\n", timer_r4.elapsed_ms / repeats);

        if (timer_r2.elapsed_ms > 0.0) {
            printf("  Relative gain: %.1f%%\n",
                   (1.0 - timer_r4.elapsed_ms / timer_r2.elapsed_ms) * 100.0);
        }

        double max_diff = 0.0;
        for (int i = 0; i < n; i++) {
            double diff = cabs(signal_r2[i] - signal_r4[i]);
            if (diff > max_diff) {
                max_diff = diff;
            }
        }

        printf("  Max difference vs radix-2: %.2e\n", max_diff);

        /* Inverse transform correctness */
        radix4_fft_with_plan(signal_r4, &plan_inv);

        double error = 0.0;
        for (int i = 0; i < n; i++) {
            error += cabs(signal_r4[i] - original[i]);
        }

        printf("  Reconstruction error: %.2e\n", error / n);

        radix4_plan_free(&plan_fwd);
        radix4_plan_free(&plan_inv);

        free_complex_array(signal_r2);
        free_complex_array(signal_r4);
        free_complex_array(original);
    }

    printf("\n\nImpulse Test (n = 64):\n");
    printf("======================\n");

    int n = 64;
    complex_t* signal = allocate_complex_array(n);
    if (!signal) {
        fprintf(stderr, "Error: Memory allocation failed for impulse test\n");
        return 1;
    }

    radix4_plan_t impulse_plan;
    if (!radix4_plan_init(&impulse_plan, n, FFT_FORWARD)) {
        fprintf(stderr, "Error: Failed to initialize radix-4 plan for impulse test\n");
        free_complex_array(signal);
        return 1;
    }

    generate_impulse(signal, n);
    radix4_fft_with_plan(signal, &impulse_plan);

    int pass = 1;
    for (int i = 0; i < n; i++) {
        if (fabs(cabs(signal[i]) - 1.0) > 1e-10) {
            pass = 0;
            break;
        }
    }

    printf("Impulse response test: %s\n", pass ? "PASS" : "FAIL");

    radix4_plan_free(&impulse_plan);
    free_complex_array(signal);

    demonstrate_butterfly();
    analyze_operation_count();

    printf("\n\nAlgorithm Characteristics:\n");
    printf("=========================\n");
    printf("✓ True radix-4 butterfly implementation\n");
    printf("✓ Requires FFT size to be a power of 4\n");
    printf("✓ Twiddle tables precomputed once per plan\n");
    printf("✓ FFT timing excludes plan construction\n");
    printf("✓ Implementation-level multiplication analysis included\n");
    printf("✓ Fewer stages than radix-2\n");
    printf("✓ In-place computation\n");

    return 0;
}

#endif /* LIB_BUILD */