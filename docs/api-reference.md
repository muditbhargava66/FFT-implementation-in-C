# API Reference

Complete reference for all functions in the FFT library.

## Platform Compatibility

Built and tested on Windows 11 with MinGW-w64 GCC 13.2.0 (1361/1361 checks, no
compiler warnings). Supported compilers are GCC and Clang; MSVC's `cl.exe`
cannot build this library because it has no C11 `_Complex` type — see the
*Platform support* table in the README.

## Core FFT Functions

### Transform Functions

Every transform below returns `fft_error_t` (declared in `fft_common.h`) rather
than `void`, because a library must not terminate the program that links it. A
bad length or a failed allocation is reported to the caller and the buffer is
left untouched:

| Value                | Meaning                                                                |
| -------------------- | ---------------------------------------------------------------------- |
| `FFT_OK` (0)         | the transform ran and the buffer holds the result                      |
| `FFT_ERR_NULL` (-1)  | `x` was null                                                           |
| `FFT_ERR_SIZE` (-2)  | length invalid for this algorithm, e.g. a non-power-of-two for radix-4 |
| `FFT_ERR_ALLOC` (-3) | scratch memory could not be obtained                                   |
| `FFT_ERR_RANGE` (-4) | argument outside the supported range (Bluestein above 2^29)            |

Callers that do not care can discard the value; C does not require using it.
Discarding is what the examples below do.

#### `radix2_dit_fft`
```c
fft_error_t radix2_dit_fft(complex_t* x, int n, fft_direction dir);
```
Performs Radix-2 Decimation-in-Time FFT.

**Parameters:**
- `x`: Input/output array of complex numbers
- `n`: Array length (must be power of 2)
- `dir`: Transform direction (`FFT_FORWARD` or `FFT_INVERSE`)

**Example:**
```c
complex_t* signal = allocate_complex_array(1024);
// Fill signal...
radix2_dit_fft(signal, 1024, FFT_FORWARD);
```

---

#### `radix2_dif_fft`
```c
fft_error_t radix2_dif_fft(complex_t* x, int n, fft_direction dir);
```
Performs Radix-2 Decimation-in-Frequency FFT.

**Parameters:** Same as `radix2_dit_fft`

---

#### `radix4_fft`
```c
fft_error_t radix4_fft(complex_t* x, int n, fft_direction dir);
```
Performs a radix-4 FFT.

**Features:**
- Works for any power of 2 (not just power of 4), with one leading radix-2
  stage when `log2(n)` is odd
- True radix-4 butterflies; twiddles advance multiplicatively per stage
- Allocates an O(n) digit-reversal buffer; `FFT_CONSERVE_MEMORY` avoids it by
  selecting in-place radix-2 DIT
- Includes input validation and error checking

---

#### `split_radix_fft`
```c
fft_error_t split_radix_fft(complex_t* x, int n, fft_direction dir);
```
Performs Split-Radix FFT with minimal operation count.

---

#### `bluestein_fft`
```c
fft_error_t bluestein_fft(complex_t* x, int n, fft_direction dir);
```
Performs FFT for arbitrary size using Bluestein's algorithm.

**Note:** Works for any size, not just powers of 2.

---

#### `mixed_radix_fft`
```c
fft_error_t mixed_radix_fft(complex_t* x, int n, fft_direction dir);
```
Performs Mixed-Radix FFT for composite sizes.

## Planning API (v3.0.0)

### `fft_plan_dft_1d`
```c
fft_plan_t fft_plan_dft_1d(int n, complex_t* in, complex_t* out,
                           int sign, unsigned flags);
```
Creates an optimized 1D complex FFT plan with automatic algorithm selection.

**Parameters:**
- `n`: Transform size
- `in`, `out`: Input and output arrays (may be the same for in-place)
- `sign`: `-1` for forward, `+1` for inverse
- `flags`: OR-combination of the planning flags below

**Returns:** a plan, or `NULL` if `n <= 0`, `in` or `out` is `NULL`, allocation
fails, or `FFT_WISDOM_ONLY` was set and no cached choice exists for this
`(n, flags)`. Callers must check for `NULL` before executing.

| Flag                  | Value  | Effect                                                                                                                                                                                                                                                                                                                                    |
| --------------------- | ------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `FFT_ESTIMATE`        | `0`    | Pick by the size heuristic. No timing, no measurement.                                                                                                                                                                                                                                                                                    |
| `FFT_MEASURE`         | `1<<0` | Time the candidate algorithms that accept `n`, keep the fastest, and cache the decision.                                                                                                                                                                                                                                                  |
| `FFT_WISDOM_ONLY`     | `1<<3` | Use a cached decision or fail. Never measures. See the note below.                                                                                                                                                                                                                                                                        |
| `FFT_REAL_INPUT`      | `1<<5` | Input is real-valued. Set automatically by `fft_plan_r2c_1d`.                                                                                                                                                                                                                                                                             |
| `FFT_REAL_OUTPUT`     | `1<<6` | Output is real-valued.                                                                                                                                                                                                                                                                                                                    |
| `FFT_CONSERVE_MEMORY` | `1<<8` | Prefer the smaller-footprint algorithm. Power-of-two lengths use in-place radix-2 DIT instead of radix-4's O(n) scratch buffer, and highly composite non-power-of-two lengths use mixed-radix instead of Bluestein. It is a constraint, so it also suppresses `FFT_MEASURE`, which would otherwise trade the saved memory back for speed. |
| `FFT_PREFER_GPU`      | `1<<9` | Select a GPU algorithm when one is available. In a default build none is, so this has no effect.                                                                                                                                                                                                                                          |

`FFT_WISDOM_ONLY` is masked out of the cache key, so a decision made under
`FFT_MEASURE` is found by a later call passing `FFT_MEASURE | FFT_WISDOM_ONLY`.
That combination is the intended use: measure once at startup, then plan without
ever paying for a timing sweep again.

`FFT_PATIENT`, `FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and `FFT_THREADED` were removed.
They were accepted and silently ignored, and the threading one had no
implementation to select. Their bit positions are deliberately left unused so
that exported wisdom keeps its meaning.

---

### `fft_plan_r2c_1d`
```c
fft_plan_t fft_plan_r2c_1d(int n, double* in, complex_t* out, fft_flags_t flags);
```
Creates a real-to-complex FFT plan. Produces `n/2+1` Hermitian-symmetric complex bins.

---

### `fft_plan_c2r_1d`
```c
fft_plan_t fft_plan_c2r_1d(int n, complex_t* in, double* out, fft_flags_t flags);
```
Creates a complex-to-real inverse FFT plan. Takes `n/2+1` complex bins and produces `n` real values.

---

### `fft_plan_dft_2d`
```c
fft_plan_t fft_plan_dft_2d(int rows, int cols, complex_t* in, complex_t* out,
                           int direction, fft_flags_t flags);
```
Creates a 2D FFT plan using row-column decomposition with internal transpose.

---

### `fft_execute`
```c
fft_error_t fft_execute(fft_plan_t plan);
```
Executes a previously created plan. Can be called repeatedly with different
data. Returns `FFT_ERR_NULL` for a null plan, otherwise the code from the
selected algorithm; nothing is written to the buffers when it fails.

---

### `fft_destroy_plan`
```c
void fft_destroy_plan(fft_plan_t plan);
```
Frees all resources associated with a plan.

---

### Wisdom Cache

The cache maps a `(transform size, planning flags)` pair to the algorithm chosen
for it. `fft_plan_dft_1d`, `fft_plan_dft_2d` and the real-valued planners all
consult it before deciding, so planning the same size twice is free, and a plan
built with `FFT_MEASURE` pays for its timing sweep once per process. The two
functions below are the only way it crosses a process boundary.

`FFT_WISDOM_ONLY` is masked out of the key before the lookup, so an entry stored
for `FFT_MEASURE` is found by a later plan asking for
`FFT_MEASURE | FFT_WISDOM_ONLY`. Every other flag bit is part of the key, so
`FFT_ESTIMATE` and `FFT_CONSERVE_MEMORY` entries are distinct — which is correct,
since they can legitimately select different algorithms.

Note that the cache is process-global and unlocked. Do not plan from more than
one thread.

#### `fft_export_wisdom_to_string`
```c
char* fft_export_wisdom_to_string(void);
```
Exports the current algorithm-selection wisdom cache as a string. Returns a `malloc`'d string that the caller must `free`, or `NULL` if it cannot be allocated. Format: a `# FFT Wisdom` header line, then one `<n> <flags> <algo>` line per entry.

---

#### `fft_import_wisdom_from_string`
```c
int fft_import_wisdom_from_string(const char* wisdom);
```
Replaces the cache with the entries in a previously exported string. Lines that
are malformed, name a non-positive size, carry an unknown algorithm id, or pair
an algorithm with a size it cannot transform (such as radix-4 with a prime
length) are skipped. Entries naming a GPU algorithm are kept even on a machine
without one, since the file is portable and the hardware is checked again when
the entry is used. Returns the number of entries imported.

**Example:**
```c
// Save wisdom after planning
char* wisdom = fft_export_wisdom_to_string();
FILE* f = fopen("fft_wisdom.txt", "w");
fputs(wisdom, f);
fclose(f);
free(wisdom);

// Load wisdom in future runs
// ... read file into string ...
fft_import_wisdom_from_string(loaded_wisdom);
```

---

#### `fft_plan_algorithm_name`
```c
const char* fft_plan_algorithm_name(fft_plan_t plan);
```
Returns the name of the algorithm the plan actually selected, such as
`"Radix-4"`, or `"Unknown"` for a null plan. Use this instead of restating the
selection rules where they are displayed or logged; a hand-written copy of them
in `demo_features.c` had already drifted from the library it described.

---

### One-Shot Interface

#### `fft_auto`
```c
fft_error_t fft_auto(complex_t* in, complex_t* out, int n, int direction);
```
Performs FFT with automatic optimization in a single call.
Internally creates a plan, executes, and destroys it.

---

## Utility Functions

### Memory Management

#### `allocate_complex_array`
```c
complex_t* allocate_complex_array(int n);
```
Allocates memory for complex array.

**Parameters:**
- `n`: Number of complex elements

**Returns:** Pointer to allocated array or NULL on failure

**Example:**
```c
complex_t* data = allocate_complex_array(1024);
if (!data) {
    fprintf(stderr, "Memory allocation failed\n");
    exit(1);
}
```

---

#### `free_complex_array`
```c
void free_complex_array(complex_t* arr);
```
Frees memory allocated by `allocate_complex_array`.

### Signal Generation

#### `generate_sine_wave`
```c
void generate_sine_wave(complex_t* signal, int n, double freq, double fs);
```
Generates a sine wave signal.

**Parameters:**
- `signal`: Output array
- `n`: Number of samples
- `freq`: Frequency in Hz
- `fs`: Sampling frequency in Hz

**Example:**
```c
generate_sine_wave(signal, 1024, 440.0, 44100.0);  // 440 Hz at 44.1 kHz
```

---

#### `generate_square_wave`
```c
void generate_square_wave(complex_t* signal, int n, double freq, double fs);
```
Generates a square wave signal.

---

#### `generate_impulse`
```c
void generate_impulse(complex_t* signal, int n);
```
Generates an impulse signal (delta function).

### Spectrum Analysis

#### `compute_magnitude`
```c
double* compute_magnitude(complex_t* fft_result, int n);
```
Computes magnitude spectrum from FFT result.

**Returns:** Array of magnitudes (must be freed by caller)

**Example:**
```c
double* mag = compute_magnitude(fft_result, n);
for (int i = 0; i < n/2; i++) {
    printf("Bin %d: %.3f\n", i, mag[i]);
}
free(mag);
```

---

#### `compute_phase`
```c
double* compute_phase(complex_t* fft_result, int n);
```
Computes phase spectrum in radians.

---

#### `compute_power_spectrum`
```c
double* compute_power_spectrum(complex_t* fft_result, int n);
```
Computes power spectrum (magnitude squared).

### Window Functions

#### `apply_window_hann`
```c
void apply_window_hann(complex_t* signal, int n);
```
Applies Hann (Hanning) window to signal.

---

#### `apply_window_hamming`
```c
void apply_window_hamming(complex_t* signal, int n);
```
Applies Hamming window to signal.

---

#### `apply_window_blackman`
```c
void apply_window_blackman(complex_t* signal, int n);
```
Applies Blackman window to signal.

## Helper Functions

### Size Validation

#### `is_power_of_two`
```c
int is_power_of_two(int n);
```
Checks if a number is a power of 2.

**Returns:** 1 if power of 2, 0 otherwise

---

#### `next_power_of_two`
```c
int next_power_of_two(int n);
```
Finds the next power of 2 greater than or equal to n.

---

#### `log2_int`
```c
int log2_int(int n);
```
Computes integer log base 2.

### Complex Number Operations

#### `print_complex`
```c
void print_complex(complex_t c);
```
Prints a complex number in format "(real, imag i)".

---

#### `print_complex_array`
```c
void print_complex_array(const char* label, complex_t* arr, int n);
```
Prints an array of complex numbers with a label.

### Performance Measurement

#### Timer Structure
```c
typedef struct {
    clock_t start;
    clock_t end;
    double elapsed_ms;
} fft_timer_t;
```

#### `timer_start`
```c
void timer_start(fft_timer_t* timer);
```
Starts timing measurement.

---

#### `timer_stop`
```c
void timer_stop(fft_timer_t* timer);
```
Stops timing and calculates elapsed time in milliseconds.

**Example:**
```c
fft_timer_t timer;
timer_start(&timer);
radix2_dit_fft(signal, n, FFT_FORWARD);
timer_stop(&timer);
printf("FFT took %.3f ms\n", timer.elapsed_ms);
```

## Constants and Types

### Constants
```c
#define PI 3.14159265358979323846
#define TWO_PI (2.0 * PI)
```

### Types
```c
typedef double complex complex_t;  // Complex number type

typedef enum {
    FFT_FORWARD = -1,
    FFT_INVERSE = 1
} fft_direction;
```

## Error Handling

### Macros

#### `CHECK_NULL`
```c
CHECK_NULL(ptr, msg);
```
Checks for NULL pointer and exits with error message.

---

#### `CHECK_POWER_OF_TWO`
```c
CHECK_POWER_OF_TWO(n);
```
Validates that size is a power of 2.

## Complete Example

```c
#include "fft_common.h"
#include "fft_algorithms.h"

int main() {
    // Setup
    int n = 1024;
    double fs = 44100.0;
    
    // Allocate memory
    complex_t* signal = allocate_complex_array(n);
    CHECK_NULL(signal, "Failed to allocate memory");
    
    // Generate test signal
    generate_sine_wave(signal, n, 1000.0, fs);  // 1 kHz tone
    
    // Apply window
    apply_window_hann(signal, n);
    
    // Perform FFT
    fft_timer_t timer;
    timer_start(&timer);
    radix2_dit_fft(signal, n, FFT_FORWARD);
    timer_stop(&timer);
    
    printf("FFT completed in %.3f ms\n", timer.elapsed_ms);
    
    // Analyze results
    double* magnitude = compute_magnitude(signal, n);
    
    // Find peak
    int peak_bin = 0;
    double peak_mag = 0.0;
    for (int i = 0; i < n/2; i++) {
        if (magnitude[i] > peak_mag) {
            peak_mag = magnitude[i];
            peak_bin = i;
        }
    }
    
    double peak_freq = peak_bin * fs / n;
    printf("Peak at %.1f Hz\n", peak_freq);
    
    // Cleanup
    free(magnitude);
    free_complex_array(signal);
    
    return 0;
}
```
