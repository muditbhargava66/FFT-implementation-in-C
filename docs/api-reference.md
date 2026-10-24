# API Reference

The public headers are `include/fft_common.h`, `include/fft_algorithms.h`,
`include/fft_auto.h` and `include/fft_gpu.h`.

## Error handling

Every transform returns `fft_error_t`. On failure the caller's buffer is left
unchanged.

| Value | Meaning |
| --- | --- |
| `FFT_OK` (0) | The transform ran and the buffer holds the result. |
| `FFT_ERR_NULL` (-1) | A required pointer was null. |
| `FFT_ERR_SIZE` (-2) | The length is invalid for this algorithm, e.g. a non-power-of-two for radix-4. |
| `FFT_ERR_ALLOC` (-3) | Scratch memory could not be obtained. |
| `FFT_ERR_RANGE` (-4) | An argument is outside the supported range, e.g. Bluestein above 2^29. |

Discarding the return value is legal C, but then a failed transform is
indistinguishable from a successful one. `CHECK_POWER_OF_TWO(n)` is an
abort-on-error macro for standalone programs: it prints a message and calls
`exit()`. Library code uses `FFT_REQUIRE_POW2(n)`, which returns `FFT_ERR_SIZE`
from the enclosing function.

## Core transforms

```c
fft_error_t radix2_dit_fft(complex_t *x, int n, fft_direction dir);
fft_error_t radix2_dif_fft(complex_t *x, int n, fft_direction dir);
fft_error_t radix4_fft(complex_t *x, int n, fft_direction dir);
fft_error_t split_radix_fft(complex_t *x, int n, fft_direction dir);
fft_error_t bluestein_fft(complex_t *x, int n, fft_direction dir);
fft_error_t mixed_radix_fft(complex_t *x, int n, fft_direction dir);
fft_error_t fft_recursive(complex_t *x, int n, fft_direction dir);
fft_error_t naive_dft(complex_t *x, int n, fft_direction dir);
fft_error_t optimized_dft(complex_t *x, int n, fft_direction dir);
```

All operate in place on `x`. `dir` is `FFT_FORWARD` (-1) or `FFT_INVERSE` (+1);
the inverse divides by `n`.

Size requirements:

| Function | Length |
| --- | --- |
| `radix2_dit_fft`, `radix2_dif_fft`, `fft_recursive` | power of two |
| `radix4_fft` | power of two (one leading radix-2 stage when `log2(n)` is odd) |
| `split_radix_fft` | power of two |
| `bluestein_fft` | any length up to 2^29 |
| `mixed_radix_fft` | any length; factors 2, 3, 5 and 7 are used, and a prime length falls back to a direct DFT |
| `naive_dft`, `optimized_dft` | any length; O(n²) reference implementations |

The sixteen convenience wrappers (`fft_radix2_dit`, `ifft_radix2_dit`, and so
on) take `(complex_t *x, int n)` and return `fft_error_t`.

`mixed_radix_fft_recursive(x, n, stride, dir)` is exposed for callers that
already own a strided layout; the stride is the distance between the elements
belonging to one transform.

## Planning API

### `fft_plan_dft_1d`

```c
fft_plan_t fft_plan_dft_1d(int n, complex_t *in, complex_t *out,
                           int sign, unsigned flags);
```

Creates a plan for a complex 1-D transform. `in` and `out` may be the same
buffer. `sign` is -1 for forward and +1 for inverse.

Returns `NULL` if `n <= 0`, either buffer is null, allocation fails, or
`FFT_WISDOM_ONLY` was set and no usable choice is cached. Check for `NULL`
before executing.

### Planning flags

| Flag | Value | Effect |
| --- | --- | --- |
| `FFT_ESTIMATE` | `0` | Choose by transform size. No timing. |
| `FFT_MEASURE` | `1<<0` | Time the candidates that accept `n`, keep the fastest, and cache it. |
| `FFT_WISDOM_ONLY` | `1<<3` | Use a cached choice or fail. Never measures. |
| `FFT_REAL_INPUT` | `1<<5` | Set automatically by `fft_plan_r2c_1d`. |
| `FFT_REAL_OUTPUT` | `1<<6` | Set automatically by `fft_plan_c2r_1d`. |
| `FFT_CONSERVE_MEMORY` | `1<<8` | Prefer in-place and low-scratch kernels. Overrides `FFT_MEASURE`. |
| `FFT_PREFER_GPU` | `1<<9` | Select a GPU algorithm when one is available. A default build has none, so this falls back to the CPU. |

`FFT_WISDOM_ONLY` is masked out of the cache key. A choice recorded under
`FFT_MEASURE` is therefore found by a later call passing
`FFT_MEASURE | FFT_WISDOM_ONLY`, which is the intended use: measure once, reuse
without measuring again.

`FFT_PATIENT`, `FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and `FFT_THREADED` were removed
because nothing implemented them. Their bit positions are left unused so
exported wisdom keeps its meaning.

### `fft_execute`

```c
fft_error_t fft_execute(fft_plan_t plan);
```

Executes the plan on the buffers it was created with. Returns `FFT_ERR_NULL` for
a null plan, otherwise the selected algorithm's error code.

### `fft_execute_dft`

```c
fft_error_t fft_execute_dft(fft_plan_t plan, complex_t *in, complex_t *out);
```

Executes a complex 1-D or 2-D plan against different buffers. Real-to-complex
and complex-to-real plans bind a `double *` at planning time and have no complex
buffer pair to substitute, so this returns `FFT_ERR_SIZE` for them.

### `fft_destroy_plan`

```c
void fft_destroy_plan(fft_plan_t plan);
```

Frees the plan, including the internal buffer used by r2c/c2r plans. Passing
`NULL` is allowed.

### `fft_plan_algorithm_name`

```c
const char *fft_plan_algorithm_name(fft_plan_t plan);
```

Returns a static string such as `"Radix-4"`, or `"Unknown"` for a null plan.
Use it instead of duplicating the selection rules in display code.

### `fft_auto`

```c
fft_error_t fft_auto(complex_t *in, complex_t *out, int n, int sign);
```

Plans with `FFT_ESTIMATE`, executes, and destroys the plan. Useful for one-shot
transforms; a reusable plan avoids the planning cost on repeated calls.

### Real-valued plans

```c
fft_plan_t fft_plan_r2c_1d(int n, double *in, complex_t *out, unsigned flags);
fft_plan_t fft_plan_c2r_1d(int n, complex_t *in, double *out, unsigned flags);
```

`r2c` takes `n` real samples and writes `n/2+1` complex bins. `c2r` takes
`n/2+1` complex bins and writes `n` real samples. Both are unnormalised on the
forward direction and divide by `n` on the inverse, so a forward/inverse pair
returns the original signal.

### 2-D plan

```c
fft_plan_t fft_plan_dft_2d(int rows, int cols, complex_t *in,
                           complex_t *out, int sign, unsigned flags);
```

Row-major input. Both dimensions must be powers of two; the planner returns
`NULL` otherwise. The transform runs a row pass, transposes, runs a column pass,
and transposes back.

## Wisdom cache

The cache maps `(transform size, planning flags)` to the algorithm chosen for
that pair. Plans consult it before deciding. It is process-global and unlocked;
do not plan concurrently.

### `fft_export_wisdom_to_string`

```c
char *fft_export_wisdom_to_string(void);
```

Returns a `malloc`'d string that the caller must `free`, or `NULL` on allocation
failure. The format is a `# FFT Wisdom` header followed by `<n> <flags> <algo>`
lines.

### `fft_import_wisdom_from_string`

```c
int fft_import_wisdom_from_string(const char *wisdom);
```

Replaces the cache with the entries in the string and returns the number
imported. Malformed lines, non-positive sizes and unknown algorithm ids are
skipped. Entries that name a GPU algorithm are kept even on a machine without
one; they are checked again when used. Entries pairing an algorithm with a size
it cannot transform are skipped.

```c
char *wisdom = fft_export_wisdom_to_string();
if (wisdom) {
    FILE *f = fopen("fft_wisdom.txt", "w");
    if (f) { fputs(wisdom, f); fclose(f); }
    free(wisdom);
}
```

## Hardware capabilities

```c
unsigned fft_get_hardware_capabilities(void);
```

Returns a bitmask of `FFT_HW_CPU_SSE`, `FFT_HW_CPU_AVX`, `FFT_HW_CPU_AVX2`,
`FFT_HW_CPU_NEON`, `FFT_HW_GPU_CUDA`, `FFT_HW_GPU_MPS` and
`FFT_HW_GPU_OPENCL`.

The CPU bits report what the processor supports. No kernel dispatches on them:
there is no SIMD or NEON implementation in this library. On x86 the AVX bits
also require the OS to preserve YMM state, which is checked through `XGETBV`.

## Memory allocation

```c
complex_t *fft_alloc_complex(size_t n);
double *fft_alloc_real(size_t n);
void fft_free(void *p);
```

`fft_alloc_complex` and `fft_alloc_real` return 64-byte aligned buffers, or
`NULL`. Free them with `fft_free`, not with `free`.

The inline helpers `allocate_complex_array` and `free_complex_array` are
available for small examples. They zero the buffer; on Windows they use
32-byte-aligned allocation, and elsewhere `calloc`.

## Utility functions

Declared as `static inline` in `fft_common.h`:

```c
int is_power_of_two(int n);
int next_power_of_two(int n);
int log2_int(int n);
unsigned int bit_reverse(unsigned int x, int log2n);
complex_t twiddle_factor(int k, int n, fft_direction dir);
void print_complex(complex_t c);
void print_complex_array(const char *label, complex_t *arr, int n);
void generate_sine_wave(complex_t *signal, int n, double freq, double fs);
void generate_square_wave(complex_t *signal, int n, double freq, double fs);
void generate_impulse(complex_t *signal, int n);
double *compute_magnitude(complex_t *fft_result, int n);
double *compute_phase(complex_t *fft_result, int n);
double *compute_power_spectrum(complex_t *fft_result, int n);
```

The three `compute_*` functions return a `malloc`'d array of `n` doubles that
the caller frees.

Window functions (implemented in `utils/window_functions.c`):

```c
void apply_window_rectangular(complex_t *signal, int n);
void apply_window_hann(complex_t *signal, int n);
void apply_window_hamming(complex_t *signal, int n);
void apply_window_blackman(complex_t *signal, int n);
void apply_window_bartlett(complex_t *signal, int n);
void apply_window_welch(complex_t *signal, int n);
```

A window shorter than two samples is left unchanged.

## Types and constants

```c
typedef double complex complex_t;

typedef enum {
    FFT_FORWARD = -1,
    FFT_INVERSE = 1
} fft_direction;

#define PI 3.14159265358979323846
#define TWO_PI (2.0 * PI)

typedef struct {
    clock_t start;
    clock_t end;
    double elapsed_ms;
} fft_timer_t;
```

`timer_start(&timer)` and `timer_stop(&timer)` fill `elapsed_ms`. `clock()`
resolution differs between platforms, so time a batch of transforms rather than
one small transform. The benchmark and planner do this internally.

`fft_version()` returns a static version string.

## Complete example

```c
#include <stdio.h>
#include <stdlib.h>
#include "fft_common.h"
#include "fft_algorithms.h"

int main(void) {
    int n = 1024;
    double fs = 44100.0;

    complex_t *signal = allocate_complex_array(n);
    if (!signal)
        return 1;

    generate_sine_wave(signal, n, 1000.0, fs);
    apply_window_hann(signal, n);
    radix2_dit_fft(signal, n, FFT_FORWARD);

    double *mag = compute_magnitude(signal, n);
    if (!mag) {
        free_complex_array(signal);
        return 1;
    }

    int peak_bin = 0;
    for (int i = 1; i < n / 2; i++)
        if (mag[i] > mag[peak_bin])
            peak_bin = i;

    printf("peak: %.1f Hz\n", peak_bin * fs / n);

    free(mag);
    free_complex_array(signal);
    return 0;
}
```

Compile with:

```bash
gcc -I./include example.c -L./build/lib -lfft_lib_c -lm -o example
```
