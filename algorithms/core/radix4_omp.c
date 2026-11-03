#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"
#include "radix4_internal.h"

#ifdef _OPENMP

#include <omp.h>

/* Below this length the parallel region costs more than the stage it
 * parallelizes. The value is a starting point; benchmark_all reports the
 * crossover on the host. */
#define OMP_MIN_N (1 << 19)

int radix4_omp_available(void) { return omp_get_max_threads() > 1; }

fft_error_t radix4_omp_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    if (!is_power_of_two(n)) return FFT_ERR_SIZE;
    if (n <= 1) return FFT_OK;
    if (n < OMP_MIN_N || omp_get_max_threads() <= 1)
        return radix4_fft(x, n, dir);

    const int log2n = log2_int(n);
    const int is_power_of_4 = (log2n % 2 == 0);

    if (fft_radix4_reorder(x, n) != 0)
        return FFT_ERR_ALLOC;

    int start_length = 4;
    if (!is_power_of_4) {
#pragma omp parallel for schedule(static)
        for (int k = 0; k < n; k += 2) {
            const complex_t t = x[k + 1];
            const complex_t u = x[k];
            x[k] = u + t;
            x[k + 1] = u - t;
        }
        start_length = 8;
    }

    for (int m = start_length; m <= n;) {
        const int m4 = m / 4;
        complex_t* tw = (complex_t*)malloc((size_t)m4 * sizeof(complex_t));
        if (!tw) return FFT_ERR_ALLOC;

        /* Build one table and derive W^2 and W^3 in the butterfly. Two extra
         * complex multiplies are cheaper than reading two more tables from
         * memory at every butterfly. */
        const complex_t step = twiddle_factor(1, m, dir);
        complex_t w = 1.0;
        for (int j = 0; j < m4; j++) {
            tw[j] = w;
            w *= step;
        }

        /* Every (block, j) pair is independent once the stage twiddles are in
         * tw[]. Collapsing both loops keeps the final stages, which have only
         * a few blocks, spread across threads. */
#pragma omp parallel for collapse(2) schedule(static)
        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < m4; j++) {
                const int i0 = k + j;
                const int i1 = i0 + m4;
                const int i2 = i0 + 2 * m4;
                const int i3 = i0 + 3 * m4;

                const complex_t w1 = tw[j];
                const complex_t w2 = w1 * w1;
                const complex_t w3 = w2 * w1;

                const complex_t x0 = x[i0];
                const complex_t x1 = x[i1] * w1;
                const complex_t x2 = x[i2] * w2;
                const complex_t x3 = x[i3] * w3;

                const complex_t t0 = x0 + x2;
                const complex_t t1 = x0 - x2;
                const complex_t t2 = x1 + x3;
                const complex_t t3 = (dir == FFT_FORWARD) ? (x1 - x3) * (-I)
                                                          : (x1 - x3) * I;

                x[i0] = t0 + t2;
                x[i1] = t1 + t3;
                x[i2] = t0 - t2;
                x[i3] = t1 - t3;
            }
        }

        free(tw);
        if (m == n) break;
        m *= 4;
    }

    if (dir == FFT_INVERSE) {
        const double scale = 1.0 / (double)n;
#pragma omp parallel for schedule(static)
        for (int i = 0; i < n; i++) x[i] *= scale;
    }
    return FFT_OK;
}

#else /* !_OPENMP */

int radix4_omp_available(void) { return 0; }

fft_error_t radix4_omp_fft(complex_t* x, int n, fft_direction dir) {
    (void)x;
    (void)n;
    (void)dir;
    return FFT_ERR_UNSUPPORTED;
}

#endif /* _OPENMP */
