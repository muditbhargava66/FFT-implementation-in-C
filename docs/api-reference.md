# API Reference

Complete reference for all functions in the FFT library.

## Core FFT Functions

### Transform Functions

#### `radix2_dit_fft`
```c
void radix2_dit_fft(complex_t* x, int n, fft_direction dir);
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
void radix2_dif_fft(complex_t* x, int n, fft_direction dir);
```
Performs Radix-2 Decimation-in-Frequency FFT.

**Parameters:** Same as `radix2_dit_fft`

---

#### `radix4_fft`
```c
void radix4_fft(complex_t* x, int n, fft_direction dir);
```
Performs Radix-4 FFT (size must be power of 4).

---

#### `split_radix_fft`
```c
void split_radix_fft(complex_t* x, int n, fft_direction dir);
```
Performs Split-Radix FFT with minimal operation count.

---

#### `bluestein_fft`
```c
void bluestein_fft(complex_t* x, int n, fft_direction dir);
```
Performs FFT for arbitrary size using Bluestein's algorithm.

**Note:** Works for any size, not just powers of 2.

---

#### `mixed_radix_fft`
```c
void mixed_radix_fft(complex_t* x, int n, fft_direction dir);
```
Performs Mixed-Radix FFT for composite sizes.

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
} timer_t;
```

#### `timer_start`
```c
void timer_start(timer_t* timer);
```
Starts timing measurement.

---

#### `timer_stop`
```c
void timer_stop(timer_t* timer);
```
Stops timing and calculates elapsed time in milliseconds.

**Example:**
```c
timer_t timer;
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
    timer_t timer;
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
