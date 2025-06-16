# Quick Start Guide

Get up and running with the FFT library in 5 minutes!

## Installation

### Using the Quick Start Script (Recommended)

```bash
# Clone the repository
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

# Run quick start
chmod +x quickstart.sh
./quickstart.sh
```

Choose option 1 to build everything for first-time users.

### Manual Installation

```bash
# Build the library and all examples
make all

# Or build only what you need
make algorithms    # Core FFT algorithms
make applications  # Example applications
make tests        # Test suite
```

## Your First FFT

### Simple Example

Create a file `test_fft.c`:

```c
#include <stdio.h>
#include "fft_common.h"
#include "fft_algorithms.h"

int main() {
    // 1. Create a signal
    int n = 16;
    complex_t* signal = allocate_complex_array(n);
    
    // 2. Generate a sine wave (5 Hz at 16 Hz sample rate)
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 5 * i / 16.0);
    }
    
    // 3. Perform FFT
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // 4. Display magnitude spectrum
    printf("Frequency Bin | Magnitude\n");
    printf("-------------|----------\n");
    for (int i = 0; i < n/2; i++) {
        double magnitude = cabs(signal[i]);
        printf("     %2d      | %.3f\n", i, magnitude);
    }
    
    // 5. Clean up
    free_complex_array(signal);
    return 0;
}
```

Compile and run:

```bash
gcc -I./include test_fft.c -L./lib -lfft -lm -o test_fft
./test_fft
```

### Expected Output

```
Frequency Bin | Magnitude
-------------|----------
      0      | 0.000
      1      | 0.000
      2      | 0.000
      3      | 0.000
      4      | 0.000
      5      | 8.000    <- Peak at 5 Hz
      6      | 0.000
      7      | 0.000
```

## Common Use Cases

### 1. Audio Spectrum Analysis

```c
// Analyze audio signal
int sample_rate = 44100;
int fft_size = 4096;
complex_t* audio = allocate_complex_array(fft_size);

// Load or generate audio data...
generate_sine_wave(audio, fft_size, 440.0, sample_rate);  // A4 note

// Apply window to reduce spectral leakage
apply_window_hann(audio, fft_size);

// Compute FFT
radix2_dit_fft(audio, fft_size, FFT_FORWARD);

// Get magnitude spectrum
double* magnitude = compute_magnitude(audio, fft_size);

// Find frequency of each bin
for (int i = 0; i < fft_size/2; i++) {
    double freq = i * sample_rate / (double)fft_size;
    printf("%.1f Hz: %.3f\n", freq, magnitude[i]);
}
```

### 2. Signal Filtering

```c
// Low-pass filter using FFT
void lowpass_filter(complex_t* signal, int n, double cutoff_freq, double sample_rate) {
    // Transform to frequency domain
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Zero out high frequencies
    int cutoff_bin = (int)(cutoff_freq * n / sample_rate);
    for (int i = cutoff_bin; i < n - cutoff_bin; i++) {
        signal[i] = 0;
    }
    
    // Transform back
    radix2_dit_fft(signal, n, FFT_INVERSE);
}
```

### 3. Fast Convolution

```c
// Convolve two signals using FFT
complex_t* fft_convolve(complex_t* x, int nx, complex_t* h, int nh) {
    int n = next_power_of_two(nx + nh - 1);
    
    complex_t* X = allocate_complex_array(n);
    complex_t* H = allocate_complex_array(n);
    
    // Zero-pad and copy
    memset(X, 0, n * sizeof(complex_t));
    memset(H, 0, n * sizeof(complex_t));
    memcpy(X, x, nx * sizeof(complex_t));
    memcpy(H, h, nh * sizeof(complex_t));
    
    // FFT both signals
    radix2_dit_fft(X, n, FFT_FORWARD);
    radix2_dit_fft(H, n, FFT_FORWARD);
    
    // Multiply in frequency domain
    for (int i = 0; i < n; i++) {
        X[i] *= H[i];
    }
    
    // Inverse FFT
    radix2_dit_fft(X, n, FFT_INVERSE);
    
    free_complex_array(H);
    return X;
}
```

## Try the Examples

### Run Pre-built Examples

```bash
# Basic FFT demo
./bin/radix2_dit

# Audio spectrum analyzer
./bin/audio_spectrum

# Convolution demo
./bin/convolution

# Run all benchmarks
./bin/benchmark_all
```

### Explore Different Algorithms

```c
// For power-of-2 sizes
radix2_dit_fft(signal, 1024, FFT_FORWARD);    // Standard
radix4_fft(signal, 1024, FFT_FORWARD);        // Fewer operations
split_radix_fft(signal, 1024, FFT_FORWARD);   // Optimal

// For arbitrary sizes
bluestein_fft(signal, 1000, FFT_FORWARD);     // Any size
mixed_radix_fft(signal, 360, FFT_FORWARD);    // Composite sizes
```

## Performance Tips

1. **Use power-of-2 sizes** when possible (faster algorithms)
2. **Apply window functions** for spectral analysis
3. **Reuse memory** for multiple transforms
4. **Consider real-FFT** for real-valued signals
5. **Use appropriate algorithm** for your size

## Next Steps

- 📖 Read the [Algorithm Guide](algorithms.md) for deep understanding
- 🔬 Check [Applications](applications.md) for more examples  
- 📊 Review [Performance Guide](performance.md) for optimization
- 🛠️ See [API Reference](api-reference.md) for all functions

## Troubleshooting

### "Size is not a power of two"
Use Bluestein or Mixed-Radix algorithms for arbitrary sizes:
```c
bluestein_fft(signal, 1000, FFT_FORWARD);
```

### Compilation Errors
Ensure you're linking with the math library:
```bash
gcc your_file.c -lfft -lm  # -lm is important!
```

### Wrong Results
- Check if you need to apply a window function
- Verify correct normalization for inverse FFT
- Ensure proper complex number handling

## Getting Help

- 📚 Check the documentation
- 🐛 Report issues on GitHub
- 💬 Ask questions in discussions
- 📧 Contact maintainers
