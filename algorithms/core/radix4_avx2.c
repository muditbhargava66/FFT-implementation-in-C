#include "../../include/fft_common.h"
#include "../../include/fft_algorithms.h"
#include "radix4_internal.h"

/*
 * AVX2 radix-4 DIT kernel.
 *
 * The vector code assumes the GCC/Clang representation of `double complex`:
 * two consecutive doubles, real then imaginary. That is the ABI on every
 * target this file is compiled for; MSVC cannot build the library anyway.
 * Two complex numbers fit in one 256-bit register.
 */

#if (defined(__x86_64__) || defined(__i386__)) && \
    (defined(__GNUC__) || defined(__clang__))
#define FFT_HAVE_AVX2_BUILD 1
#include <immintrin.h>
#else
#define FFT_HAVE_AVX2_BUILD 0
#endif

/* The scalar radix-4 kernel still wins below this size because the AVX2
 * kernel builds a twiddle table per stage. Measured on the verification host;
 * benchmark_all reports the local crossover. */
#define AVX2_MIN_N (1 << 12)

int radix4_avx2_available(void) {
#if FFT_HAVE_AVX2_BUILD
    __builtin_cpu_init();
    return __builtin_cpu_supports("avx2") != 0;
#else
    return 0;
#endif
}

#if FFT_HAVE_AVX2_BUILD

#define AVX2_TARGET __attribute__((target("avx2")))

AVX2_TARGET
static inline __m256d cmul(__m256d a, __m256d b) {
    /* [ar,ai] * [br,bi] = [ar*br - ai*bi, ar*bi + ai*br], for both complex
     * values in the register. addsub subtracts on even lanes and adds on odd
     * lanes, which is exactly the real/imaginary split. */
    const __m256d ar = _mm256_movedup_pd(a);
    const __m256d ai = _mm256_permute_pd(a, 0xF);
    const __m256d bs = _mm256_permute_pd(b, 0x5);
    return _mm256_addsub_pd(_mm256_mul_pd(ar, b), _mm256_mul_pd(ai, bs));
}

AVX2_TARGET
static inline __m256d mul_minus_i(__m256d z) {
    /* (re + i*im) * (-i) = im - i*re */
    const __m256d swap = _mm256_permute_pd(z, 0x5);
    const __m256d sign = _mm256_set_pd(-1.0, 1.0, -1.0, 1.0);
    return _mm256_mul_pd(swap, sign);
}

AVX2_TARGET
static inline __m256d mul_plus_i(__m256d z) {
    /* (re + i*im) * i = -im + i*re */
    const __m256d swap = _mm256_permute_pd(z, 0x5);
    const __m256d sign = _mm256_set_pd(1.0, -1.0, 1.0, -1.0);
    return _mm256_mul_pd(swap, sign);
}

AVX2_TARGET
static void radix4_stage_avx2(complex_t* x, int n, int m, int m4,
                              fft_direction dir, const complex_t* tw) {
    const int m2 = 2 * m4;
    const int m3 = 3 * m4;

    for (int k = 0; k < n; k += m) {
        int j = 0;

        /* m4 is 1 for the m=4 stage and even for every later stage, so the
         * vector loop covers all but a possible single trailing point. */
        for (; j + 1 < m4; j += 2) {
            const int i0 = k + j;
            const int i1 = i0 + m4;
            const int i2 = i0 + m2;
            const int i3 = i0 + m3;

            const __m256d w1 = _mm256_loadu_pd((const double*)&tw[j]);
            const __m256d w2 = cmul(w1, w1);
            const __m256d w3 = cmul(w2, w1);

            const __m256d x0 = _mm256_loadu_pd((const double*)&x[i0]);
            const __m256d x1 = cmul(_mm256_loadu_pd((const double*)&x[i1]), w1);
            const __m256d x2 = cmul(_mm256_loadu_pd((const double*)&x[i2]), w2);
            const __m256d x3 = cmul(_mm256_loadu_pd((const double*)&x[i3]), w3);

            const __m256d t0 = _mm256_add_pd(x0, x2);
            const __m256d t1 = _mm256_sub_pd(x0, x2);
            const __m256d t2 = _mm256_add_pd(x1, x3);
            const __m256d d = _mm256_sub_pd(x1, x3);
            const __m256d t3 =
                (dir == FFT_FORWARD) ? mul_minus_i(d) : mul_plus_i(d);

            _mm256_storeu_pd((double*)&x[i0], _mm256_add_pd(t0, t2));
            _mm256_storeu_pd((double*)&x[i1], _mm256_add_pd(t1, t3));
            _mm256_storeu_pd((double*)&x[i2], _mm256_sub_pd(t0, t2));
            _mm256_storeu_pd((double*)&x[i3], _mm256_sub_pd(t1, t3));
        }

        if (j < m4) {
            const int i0 = k + j;
            const int i1 = i0 + m4;
            const int i2 = i0 + m2;
            const int i3 = i0 + m3;

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
            const complex_t t3 =
                (dir == FFT_FORWARD) ? (x1 - x3) * (-I) : (x1 - x3) * I;

            x[i0] = t0 + t2;
            x[i1] = t1 + t3;
            x[i2] = t0 - t2;
            x[i3] = t1 - t3;
        }
    }
}

fft_error_t radix4_avx2_fft(complex_t* x, int n, fft_direction dir) {
    if (!x) return FFT_ERR_NULL;
    if (!is_power_of_two(n)) return FFT_ERR_SIZE;
    if (n <= 1) return FFT_OK;
    if (!radix4_avx2_available()) return FFT_ERR_UNSUPPORTED;
    if (n < AVX2_MIN_N) return radix4_fft(x, n, dir);

    const int log2n = log2_int(n);
    const int is_power_of_4 = (log2n % 2 == 0);

    if (fft_radix4_reorder(x, n) != 0)
        return FFT_ERR_ALLOC;

    int start_length = 4;
    if (!is_power_of_4) {
        for (int k = 0; k < n; k += 2) {
            const complex_t t = x[k + 1];
            const complex_t u = x[k];
            x[k] = u + t;
            x[k + 1] = u - t;
        }
        start_length = 8;
    }

    complex_t* tw = (complex_t*)malloc(((size_t)n / 4) * sizeof(complex_t));
    if (!tw) return FFT_ERR_ALLOC;

    for (int m = start_length; m <= n;) {
        const int m4 = m / 4;
        const complex_t step = twiddle_factor(1, m, dir);
        complex_t w = 1.0;
        for (int j = 0; j < m4; j++) {
            tw[j] = w;
            w *= step;
        }

        radix4_stage_avx2(x, n, m, m4, dir, tw);

        if (m == n) break;
        m *= 4;
    }

    free(tw);

    if (dir == FFT_INVERSE) {
        const double scale = 1.0 / (double)n;
        for (int i = 0; i < n; i++) x[i] *= scale;
    }
    return FFT_OK;
}

#else /* !FFT_HAVE_AVX2_BUILD */

fft_error_t radix4_avx2_fft(complex_t* x, int n, fft_direction dir) {
    (void)x;
    (void)n;
    (void)dir;
    return FFT_ERR_UNSUPPORTED;
}

#endif /* FFT_HAVE_AVX2_BUILD */
