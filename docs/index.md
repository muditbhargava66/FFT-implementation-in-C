# FFT Implementation in C Documentation

Welcome to the comprehensive documentation for the FFT Implementation in C project. This repository provides educational and production-ready implementations of various Fast Fourier Transform algorithms.

## Overview

The Fast Fourier Transform (FFT) is one of the most important algorithms in digital signal processing, enabling efficient computation of the Discrete Fourier Transform (DFT). This project provides:

- **Multiple Algorithm Implementations**: From basic to advanced FFT algorithms
- **Real-World Applications**: Audio processing, filtering, and more
- **Performance Optimizations**: SIMD, parallel, and embedded implementations
- **Educational Resources**: Detailed explanations and visualizations

## Features

### 🎯 Core Algorithms
- **Radix-2 (DIT/DIF)**: The classic Cooley-Tukey algorithms
- **Radix-4**: Higher radix for better performance
- **Split-Radix**: Optimal operation count
- **Bluestein**: Arbitrary size FFTs
- **Mixed-Radix**: Composite size optimization

### 🚀 Optimizations
- **SIMD Implementation**: Using SSE/AVX instructions
- **Parallel Processing**: OpenMP multi-core support
- **Fixed-Point**: For embedded systems
- **Cache-Optimized**: Memory-efficient implementations

### 📊 Applications
- Audio spectrum analysis
- Digital signal filtering
- Fast convolution
- Image processing
- Power spectrum computation

## Quick Navigation

- **[Getting Started](getting-started.md)**: Installation and first steps
- **[Algorithm Guide](algorithms.md)**: Detailed algorithm explanations
- **[API Reference](api-reference.md)**: Complete function documentation
- **[Applications](applications.md)**: Real-world usage examples
- **[Performance Guide](performance.md)**: Benchmarks and optimization tips

## Mathematical Background

The Discrete Fourier Transform of a sequence $x[n]$ is defined as:

$$X[k] = \sum_{n=0}^{N-1} x[n] \cdot e^{-j2\pi kn/N}$$

The FFT algorithms reduce the computational complexity from $O(N^2)$ to $O(N \log N)$ by exploiting symmetries in the twiddle factors.

## Code Example

```c
#include "fft_common.h"
#include "fft_algorithms.h"

int main() {
    int n = 1024;
    complex_t* signal = allocate_complex_array(n);
    
    // Generate a test signal
    double fs = 1000.0;  // Sampling frequency
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 50 * i / fs);  // 50 Hz sine wave
    }
    
    // Apply FFT
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Process spectrum...
    double* magnitude = compute_magnitude(signal, n);
    
    // Clean up
    free(magnitude);
    free_complex_array(signal);
    return 0;
}
```

## Project Structure

```
FFT-implementation-in-C/
├── include/          # Header files
├── algorithms/       # Core FFT implementations
├── applications/     # Real-world applications
├── optimizations/    # Performance-optimized versions
├── benchmarks/       # Performance testing
├── tests/           # Unit tests
├── docs/            # Documentation
└── examples/        # Example programs
```

## Why This Project?

This project serves multiple purposes:

1. **Educational Resource**: Learn FFT algorithms with clear, documented code
2. **Reference Implementation**: Production-quality code with proper error handling
3. **Performance Comparison**: Benchmark different algorithms
4. **Practical Examples**: See FFT used in real applications

## License

This project is licensed under the MIT License. See [LICENSE](https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/LICENSE) for details.
