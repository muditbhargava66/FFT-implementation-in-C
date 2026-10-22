#ifndef FFT_COMMON_H
#define FFT_COMMON_H

/* Feature-test macros must be set before the first system header is included.
 * Strict -std=c11/c99 makes glibc hide posix_memalign/sysconf, and GCC 14+
 * turns the resulting implicit declarations into hard errors. Defining it here
 * makes the fix independent of compiler flags; keep the CMake definition too,
 * because a translation unit that includes <unistd.h> before this header is
 * already past the point where it can take effect. */
#ifndef _WIN32
#  ifndef _GNU_SOURCE
#    define _GNU_SOURCE 1
#  endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <assert.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifdef _WIN32
    #include <malloc.h>  // For _aligned_malloc on Windows
#endif

/* Every transform here is written against C11 complex arithmetic: `x[i] =
 * x[j] * w`, `a += b`, `x /= n`. That needs the _Complex keyword and the
 * operator overloading the standard defines for it.
 *
 * The Microsoft C compiler has a complex.h (cexp, cabs, _Dcomplex, _Cbuild) but
 * no _Complex keyword, and Microsoft's own documentation is explicit that "the
 * complex math library functions are NOT overloaded with the +, -, *, and /
 * operators" because complex numbers are not a native type there. Supporting
 * cl.exe in C mode would therefore mean rewriting every kernel against
 * _Cmulcc/_Caddcc-style helpers.
 *
 * The branch that used to sit here did neither: it included <complex>, which is
 * a C++ header, and typedef'd _Dcomplex, so it could not compile either. It is
 * replaced by this diagnostic rather than a half-written emulation that no one
 * can test, because MSVC builds are only reachable from a machine with Visual
 * Studio installed. On Windows, build with MinGW-w64 gcc or clang, both of
 * which the CMake configuration supports. */
#if defined(_MSC_VER) && !defined(__clang__)
#error This library requires C11 _Complex arithmetic, which the MSVC C compiler does not provide. Use MinGW-w64 gcc or clang on Windows.
#endif

#include <complex.h>
typedef double complex complex_t;

// Compiler optimization hints
#ifdef __GNUC__
    #define LIKELY(x)   __builtin_expect(!!(x), 1)
    #define UNLIKELY(x) __builtin_expect(!!(x), 0)
    #define FORCE_INLINE __attribute__((always_inline)) inline
#else
    #define LIKELY(x)   (x)
    #define UNLIKELY(x) (x)
    #define FORCE_INLINE inline
#endif

// Constants
#define PI 3.14159265358979323846
#define TWO_PI (2.0 * PI)

// FFT direction
typedef enum {
    FFT_FORWARD = -1,
    FFT_INVERSE = 1
} fft_direction;

/*
 * Error codes returned by every transform entry point. A library must not
 * terminate the process it is linked into, so the transforms report instead of
 * calling exit(); FFT_ERR_SIZE replaces the old CHECK_POWER_OF_TWO macro, which
 * aborted the host program on a bad length.
 */
typedef enum {
    FFT_OK = 0,             /* success */
    FFT_ERR_NULL = -1,      /* null buffer */
    FFT_ERR_SIZE = -2,      /* length invalid for this algorithm */
    FFT_ERR_ALLOC = -3,     /* allocation failed */
    FFT_ERR_RANGE = -4      /* argument out of supported range */
} fft_error_t;

/* Validate a required power-of-two length from inside a function returning
 * fft_error_t. */
#define FFT_REQUIRE_POW2(n) \
    do { if (!is_power_of_two(n)) return FFT_ERR_SIZE; } while (0)

// Common utility functions
static FORCE_INLINE int is_power_of_two(int n) {
    return LIKELY(n > 0) && ((n & (n - 1)) == 0);
}

static inline int next_power_of_two(int n) {
    if (n <= 1) return 1;
    unsigned int v = (unsigned int)(n - 1);
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    return (int)(v + 1u);
}

/* Returns 0 for n <= 1. Shifting a signed negative right never reaches zero,
 * so the previous form spun forever on negative input. */
static inline int log2_int(int n) {
    int log = 0;
    while (n > 1) { n >>= 1; log++; }
    return log;
}

// Bit reversal function with optimization for small sizes
static inline unsigned int bit_reverse(unsigned int x, int log2n) {
    if (log2n <= 0) return 0;

    if (log2n <= 8) {
        /* All four swap stages must run: skipping the byte swap leaves a
         * partial reversal and misplaces every index for log2n <= 4. */
        x = ((x & 0xAAAA) >> 1) | ((x & 0x5555) << 1);
        x = ((x & 0xCCCC) >> 2) | ((x & 0x3333) << 2);
        x = ((x & 0xF0F0) >> 4) | ((x & 0x0F0F) << 4);
        x = ((x & 0xFF00) >> 8) | ((x & 0x00FF) << 8);
        return x >> (16 - log2n);
    }

    /* General case for log2n > 8 */
    unsigned int r = 0;
    for (int i = 0; i < log2n; i++) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

// Memory allocation helpers with alignment for SIMD
static inline complex_t* allocate_complex_array(int n) {
    if (n <= 0) return NULL;
#ifdef _WIN32
    // Use aligned allocation on Windows for better SIMD performance
    void* ptr = _aligned_malloc((size_t)n * sizeof(complex_t), 32);
    if (ptr) {
        memset(ptr, 0, (size_t)n * sizeof(complex_t));
    }
    return (complex_t*)ptr;
#else
    return (complex_t*)calloc((size_t)n, sizeof(complex_t));
#endif
}

static inline void free_complex_array(complex_t* arr) {
    if (arr) {
#ifdef _WIN32
        _aligned_free(arr);
#else
        free(arr);
#endif
    }
}

// Twiddle factor computation with optimization for common cases
static inline complex_t twiddle_factor(int k, int n, fft_direction dir) {
    // Optimize for common cases
    if (k == 0) return 1.0;
    if (k * 4 == n) return (dir == FFT_FORWARD) ? -I : I;  // ±j
    if (k * 2 == n) return -1.0;  // -1
    if (k * 4 == 3 * n) return (dir == FFT_FORWARD) ? I : -I;  // ∓j
    
    double angle = dir * TWO_PI * k / n;
    return cexp(I * angle);
}

// Performance timing (renamed to avoid conflict with system timer_t)
typedef struct {
    clock_t start;
    clock_t end;
    double elapsed_ms;
} fft_timer_t;

static inline void timer_start(fft_timer_t* timer) {
    timer->start = clock();
}

static inline void timer_stop(fft_timer_t* timer) {
    timer->end = clock();
    timer->elapsed_ms = ((double)(timer->end - timer->start) / CLOCKS_PER_SEC) * 1000.0;
}

/*
 * Abort-on-error helper for STANDALONE PROGRAMS only. It calls exit(), which is
 * never acceptable inside the library; library code returns fft_error_t via
 * FFT_REQUIRE_POW2 instead. Its one remaining user is applications/image_fft.c.
 * The matching CHECK_NULL was removed once the library helpers it guarded began
 * returning NULL/fft_error_t: a header macro that kills the host process is not
 * something a library should export.
 */
#define CHECK_POWER_OF_TWO(n) \
    if (!is_power_of_two(n)) { \
        fprintf(stderr, "Error: Size %d is not a power of two\n", n); \
        exit(EXIT_FAILURE); \
    }

// Complex number utilities
static inline void print_complex(complex_t c) {
    double real = creal(c);
    double imag = cimag(c);
    if (fabs(real) < 1e-10) real = 0;
    if (fabs(imag) < 1e-10) imag = 0;
    printf("(%.3f, %.3fi)", real, imag);
}

static inline void print_complex_array(const char* label, complex_t* arr, int n) {
    printf("%s: ", label);
    for (int i = 0; i < n; i++) {
        print_complex(arr[i]);
        printf(" ");
    }
    printf("\n");
}

// Signal generation utilities
static inline void generate_sine_wave(complex_t* signal, int n, double freq, double fs) {
    for (int i = 0; i < n; i++) {
        signal[i] = sin(TWO_PI * freq * i / fs) + 0 * I;
    }
}

static inline void generate_square_wave(complex_t* signal, int n, double freq, double fs) {
    int period = (int)(fs / freq);
    for (int i = 0; i < n; i++) {
        signal[i] = (i % period < period / 2) ? 1.0 : -1.0;
    }
}

static inline void generate_impulse(complex_t* signal, int n) {
    memset(signal, 0, n * sizeof(complex_t));
    signal[0] = 1.0;
}

/* The compute_* helpers return NULL on failure. They are called from library
 * code and from applications, so they must not exit() the process. */
static inline double* compute_magnitude(complex_t* fft_result, int n) {
    if (!fft_result || n <= 0) return NULL;
    double* mag = (double*)malloc((size_t)n * sizeof(double));
    if (!mag) return NULL;

    for (int i = 0; i < n; i++) {
        mag[i] = cabs(fft_result[i]);
    }
    return mag;
}

static inline double* compute_phase(complex_t* fft_result, int n) {
    if (!fft_result || n <= 0) return NULL;
    double* phase = (double*)malloc((size_t)n * sizeof(double));
    if (!phase) return NULL;

    for (int i = 0; i < n; i++) {
        phase[i] = carg(fft_result[i]);
    }
    return phase;
}

static inline double* compute_power_spectrum(complex_t* fft_result, int n) {
    if (!fft_result || n <= 0) return NULL;
    double* power = (double*)malloc((size_t)n * sizeof(double));
    if (!power) return NULL;

    for (int i = 0; i < n; i++) {
        double mag = cabs(fft_result[i]);
        power[i] = mag * mag / n;
    }
    return power;
}

// Window functions for signal processing
void apply_window_hann(complex_t* signal, int n);
void apply_window_hamming(complex_t* signal, int n);
void apply_window_blackman(complex_t* signal, int n);
void apply_window_rectangular(complex_t* signal, int n);
void apply_window_bartlett(complex_t* signal, int n);
void apply_window_welch(complex_t* signal, int n);

#endif // FFT_COMMON_H
