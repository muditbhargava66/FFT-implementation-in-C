#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"

/**
 * Mixed-Radix FFT Implementation
 *
 * Handles any composite length (6, 10, 12, 15, 20, 24, 360, 1000, ...) by
 * recursive Cooley-Tukey factorisation, so no length that factors into small
 * primes is quadratic any more.
 *
 * The decomposition uses n = N2*n1 + n2 and k = k1 + N1*k2 with n = N1*N2:
 *
 *   X[k1 + N1*k2] = SUM_{n2} w_N^{n2*k1} * w_N2^{n2*k2}
 *                     * ( SUM_{n1} w_N1^{n1*k1} * x[N2*n1 + n2] )
 *
 * i.e. N2 transforms of length N1, then a twiddle multiply, then N1 transforms
 * of length N2. Because k = k1 + N1*k2 enumerates 0..N-1 exactly once as k1 and
 * k2 range, the output is already in natural order: there is no digit-reversal
 * step to get wrong.
 *
 * Prime lengths have no non-trivial factorisation and take a direct O(p^2)
 * DFT; fft_auto routes prime sizes to Bluestein instead, so this only affects
 * sizes whose largest prime factor is large.
 *
 * Time Complexity: O(n log n) whenever n is not prime
 */

/* Smallest prime factor of n, or 0 when n < 4 or n is prime. */
static int smallest_factor(int n) {
    if (n < 4) return 0;
    if (n % 2 == 0) return 2;
    for (int i = 3; (long)i * i <= (long)n; i += 2) {
        if (n % i == 0) return i;
    }
    return 0;
}

/* Unnormalised direct DFT; used for prime-length leaves. */
static void dft_direct(complex_t* out, const complex_t* in, int n,
                       fft_direction dir) {
    for (int k = 0; k < n; k++) {
        complex_t acc = 0.0;
        for (int j = 0; j < n; j++) {
            /* Reduce the angle argument modulo n before scaling: forming j*k
             * outright loses precision and can overflow for larger primes. */
            double angle = (double)dir * TWO_PI *
                           (double)((long)j * k % (long)n) / (double)n;
            acc += in[j] * (cos(angle) + I * sin(angle));
        }
        out[k] = acc;
    }
}

/* Unnormalised power-of-two transform, for either direction.
 * radix2_dit_fft normalises its inverse pass, so the inverse is obtained via
 * the conjugation identity IDFT_unnorm(x) = conj(DFT(conj(x))) to keep every
 * path in this file unnormalised and let the caller scale once. */
static void dft_pow2_unnormalised(complex_t* out, const complex_t* in, int n,
                                  fft_direction dir) {
    if (dir == FFT_FORWARD) {
        memcpy(out, in, (size_t)n * sizeof(complex_t));
        radix2_dit_fft(out, n, FFT_FORWARD);
        return;
    }
    for (int i = 0; i < n; i++) out[i] = conj(in[i]);
    radix2_dit_fft(out, n, FFT_FORWARD);
    for (int i = 0; i < n; i++) out[i] = conj(out[i]);
}

/*
 * Recursive mixed-radix kernel. out and in must not overlap; in is read fully
 * before out is written at each level, but the sub-transforms write into out
 * while later stages still read from in, so aliasing would corrupt the result.
 */
static int ct_recursion(complex_t* out, const complex_t* in, int n,
                        fft_direction dir) {
    if (n == 1) {
        out[0] = in[0];
        return 0;
    }

    if (is_power_of_two(n)) {
        dft_pow2_unnormalised(out, in, n, dir);
        return 0;
    }

    int n1 = smallest_factor(n);
    if (n1 == 0) {              /* prime length */
        dft_direct(out, in, n, dir);
        return 0;
    }
    int n2 = n / n1;

    complex_t* e = (complex_t*)malloc((size_t)n * sizeof(complex_t));
    complex_t* row = (complex_t*)malloc((size_t)n1 * sizeof(complex_t));
    complex_t* col = (complex_t*)malloc((size_t)n2 * sizeof(complex_t));
    complex_t* res = (complex_t*)malloc((size_t)n2 * sizeof(complex_t));
    if (!e || !row || !col || !res) {
        free(e); free(row); free(col); free(res);
        return -1;
    }

    /* Stage A: n2 transforms of length n1 over strided input rows. */
    for (int j2 = 0; j2 < n2; j2++) {
        for (int j1 = 0; j1 < n1; j1++) row[j1] = in[j2 + n2 * j1];
        if (ct_recursion(e + (size_t)j2 * n1, row, n1, dir) != 0) {
            free(e); free(row); free(col); free(res);
            return -1;
        }
    }

    /* Twiddle: e[j2*n1 + k1] *= w_n^(j2*k1). w_n^1 is computed once and the
     * per-row rotation is accumulated, so this whole stage costs a single
     * transcendental call rather than one per butterfly. */
    complex_t w_n = twiddle_factor(1, n, dir);
    complex_t w_j2 = 1.0;
    for (int j2 = 0; j2 < n2; j2++) {
        complex_t w = 1.0;
        complex_t* erow = e + (size_t)j2 * n1;
        for (int k1 = 0; k1 < n1; k1++) {
            erow[k1] *= w;
            w *= w_j2;
        }
        w_j2 *= w_n;
    }

    /* Stage B: n1 transforms of length n2 down the columns, scattered into the
     * natural output position k = k1 + n1*k2. */
    for (int k1 = 0; k1 < n1; k1++) {
        for (int j2 = 0; j2 < n2; j2++) col[j2] = e[(size_t)j2 * n1 + k1];
        if (ct_recursion(res, col, n2, dir) != 0) {
            free(e); free(row); free(col); free(res);
            return -1;
        }
        for (int k2 = 0; k2 < n2; k2++) out[k1 + n1 * k2] = res[k2];
    }

    free(e); free(row); free(col); free(res);
    return 0;
}

// General small-N DFT (for any factors), used by the strided entry point.
static void dft_general(complex_t* x, int n, int stride, fft_direction dir) {
    if (n <= 0) return;
    complex_t* temp = allocate_complex_array(n);
    if (!temp) return;

    for (int i = 0; i < n; i++) {
        temp[i] = x[i * stride];
    }

    for (int k = 0; k < n; k++) {
        x[k * stride] = 0;
        for (int j = 0; j < n; j++) {
            x[k * stride] += temp[j] * cexp(dir * I * TWO_PI * j * k / n);
        }
    }

    /* Scale for inverse DFT */
    if (dir == FFT_INVERSE) {
        for (int k = 0; k < n; k++) {
            x[k * stride] /= n;
        }
    }

    free_complex_array(temp);
}

// Main mixed-radix FFT function, in place, scaled like the other transforms.
fft_error_t mixed_radix_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    if (n <= 0) return FFT_ERR_SIZE;

    if (n == 1) return FFT_OK;

    if (is_power_of_two(n)) return radix2_dit_fft(x, n, dir);

    complex_t* in = (complex_t*)malloc((size_t)n * sizeof(complex_t));
    if (!in) return FFT_ERR_ALLOC;
    memcpy(in, x, (size_t)n * sizeof(complex_t));

    int exhausted = ct_recursion(x, in, n, dir);
    free(in);
    if (exhausted) return FFT_ERR_ALLOC;

    if (dir == FFT_INVERSE) {
        double scale = 1.0 / (double)n;
        for (int i = 0; i < n; i++) x[i] *= scale;
    }
    return FFT_OK;
}

// Strided entry point from fft_algorithms.h. A strided transform must read and
// write x[i*stride]; ignoring stride would transform the wrong elements and
// overwrite x[0..n-1] regardless of stride. The kernel above is contiguous, so
// strided sub-transforms take the direct DFT, which does honour the stride.
void mixed_radix_fft_recursive(complex_t* x, int n, int stride, fft_direction dir) {
    if (!x || n <= 0 || stride <= 0) return;

    if (stride == 1) {
        mixed_radix_fft(x, n, dir);
        return;
    }

    dft_general(x, n, stride, dir);
}

// Wrapper functions
fft_error_t fft_mixed_radix(complex_t* x, int n) {
    return mixed_radix_fft(x, n, FFT_FORWARD);
}

fft_error_t ifft_mixed_radix(complex_t* x, int n) {
    return mixed_radix_fft(x, n, FFT_INVERSE);
}

#ifndef LIB_BUILD

/* Factorization — only needed by the demo program below */
typedef struct {
  int *factors;
  int num_factors;
} factorization_t;

static factorization_t factorize(int n) {
  factorization_t result;
  result.factors = (int *)malloc(32 * sizeof(int));
  result.num_factors = 0;
  if (!result.factors)
    return result;

  while (n % 2 == 0) {
    result.factors[result.num_factors++] = 2;
    n /= 2;
  }
  for (int i = 3; i * i <= n; i += 2) {
    while (n % i == 0) {
      result.factors[result.num_factors++] = i;
      n /= i;
    }
  }
  if (n > 2)
    result.factors[result.num_factors++] = n;

  return result;
}

int main() {
    printf("Mixed-Radix FFT Implementation\n");
    printf("==============================\n\n");

    int test_sizes[] = {6, 10, 12, 15, 20, 24, 30, 35, 40, 42,
                        48, 60, 72, 84, 90, 100, 120, 144};
    int num_tests = sizeof(test_sizes) / sizeof(test_sizes[0]);

    printf("Testing composite-length FFTs:\n");
    printf("Size\tFactors\t\tMixed-Radix\tNaive DFT\tError\n");
    printf("----\t-------\t\t-----------\t---------\t-----\n");

    for (int t = 0; t < num_tests; t++) {
        int n = test_sizes[t];

        printf("%d\t", n);
        factorization_t factors = factorize(n);
        for (int i = 0; i < factors.num_factors; i++) {
            printf("%d", factors.factors[i]);
            if (i < factors.num_factors - 1) printf("x");
        }
        printf("\t\t");
        if (factors.num_factors < 3) printf("\t");
        free(factors.factors);

        complex_t* x_mixed = allocate_complex_array(n);
        complex_t* x_naive = allocate_complex_array(n);

        for (int i = 0; i < n; i++) {
            double value = sin(2 * PI * 2 * i / n) +
                           0.3 * cos(2 * PI * 5 * i / n);
            x_mixed[i] = x_naive[i] = value;
        }

        fft_timer_t timer;
        timer_start(&timer);
        fft_mixed_radix(x_mixed, n);
        timer_stop(&timer);
        double time_mixed = timer.elapsed_ms;

        timer_start(&timer);
        naive_dft(x_naive, n, FFT_FORWARD);
        timer_stop(&timer);
        double time_naive = timer.elapsed_ms;

        double max_error = 0;
        for (int i = 0; i < n; i++) {
            double error = cabs(x_mixed[i] - x_naive[i]);
            if (error > max_error) max_error = error;
        }

        printf("%.3f ms\t%.3f ms\t%.2e\n",
               time_mixed, time_naive, max_error);

        free_complex_array(x_mixed);
        free_complex_array(x_naive);
    }

    printf("\n\nDetailed example: 12-point FFT (12 = 3 x 4)\n");
    printf("--------------------------------------------\n");

    int n = 12;
    complex_t* signal = allocate_complex_array(n);

    for (int i = 0; i < n; i++)
        signal[i] = i;

    printf("Input: ");
    for (int i = 0; i < n; i++)
        printf("%.0f ", creal(signal[i]));
    printf("\n");

    fft_mixed_radix(signal, n);

    printf("\nFFT output:\n");
    for (int i = 0; i < n; i++) {
        printf("X[%2d] = ", i);
        print_complex(signal[i]);
        printf("\n");
    }

    ifft_mixed_radix(signal, n);

    printf("\nAfter inverse FFT: ");
    for (int i = 0; i < n; i++)
        printf("%.0f ", creal(signal[i]));
    printf("\n");

    free_complex_array(signal);

    return 0;
}

#endif /* LIB_BUILD */

