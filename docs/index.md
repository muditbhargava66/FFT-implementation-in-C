# FFT Implementation in C Documentation

Welcome to the comprehensive documentation for the FFT Implementation in C project. This repository provides educational and production-ready implementations of various Fast Fourier Transform algorithms.

## Overview

The Fast Fourier Transform (FFT) is one of the most important algorithms in digital signal processing, enabling efficient computation of the Discrete Fourier Transform (DFT). This project provides:

- **Multiple Algorithm Implementations**: From basic to advanced FFT algorithms
- **Real-World Applications**: Audio processing, filtering, and more
- **Measured Algorithm Selection**: kernel choice backed by recorded cost ratios,
  with an optional timing sweep and a reusable wisdom cache
- **Educational Resources**: Detailed explanations and diagrams

## Features

### Core Algorithms
- **Radix-2 (DIT/DIF)**: The classic Cooley-Tukey algorithms
- **Radix-4**: Higher radix for better performance
- **Split-Radix**: Optimal operation count
- **Bluestein**: Arbitrary size FFTs
- **Mixed-Radix**: Composite size optimization

### Planning
- **`FFT_ESTIMATE`**: pick by the size heuristic, no timing
- **`FFT_MEASURE`**: time the candidates that accept this length and keep the
  fastest, then cache the decision
- **`FFT_WISDOM_ONLY`**: use a cached decision or fail to plan; never measures.
  Useful for latency-bounded callers that would rather get `NULL` than pay for a
  timing sweep
- **`fft_plan_algorithm_name()`**: report what a plan actually selected
- Wisdom exports to and imports from a text string, so a decision survives a
  process restart

### What this library does not do
There is no SIMD, threaded, fixed-point or GPU execution path. Earlier revisions
of these docs listed SSE2/AVX kernels, a pthreads/OpenMP stage splitter and a
Q15 path; those files were compiled into the library but declared in no header
and called by nothing, so no plan could reach them. They have been deleted rather
than left as unreachable objects. See the Changelog.

### Applications
- Audio spectrum analysis with beautiful Unicode displays
- Digital signal filtering
- Fast convolution
- Image processing
- Power spectrum computation

### Platform status
- **1361/1361 checks passing**, verified against an independent direct DFT
- **Supported compilers**: GCC and Clang, on Linux, macOS and Windows
  (MinGW-w64 or clang). MSVC's `cl.exe` cannot build this library: it needs C11
  `_Complex` arithmetic, which it does not provide
- **CMake build** — simple `cmake -S . -B build && cmake --build build`
- **Not thread-safe**: selection and the wisdom cache use shared globals without
  locking
- **GPU code does not execute**: in a default build every `fft_gpu_*` symbol
  resolves to `utils/fft_gpu_stubs.c`, where `fft_gpu_available()` returns 0.
  Configuring with `-DFFT_ENABLE_GPU=ON` compiles the CUDA objects without
  connecting them to the planner. No GPU result in this repository has been
  measured on hardware

## Quick Navigation

- **[Getting Started](getting-started.md)**: Installation and first steps
- **[Algorithm Guide](algorithms.md)**: Detailed algorithm explanations
- **[API Reference](api-reference.md)**: Complete function documentation
- **[Applications](applications.md)**: Real-world usage examples
- **[Performance Guide](performance.md)**: Benchmarks and optimization tips
- **[Optimizations](optimizations.md)**: v3.0.0 performance improvements and techniques
- **[Security Policy](https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/SECURITY.md)**:
  private vulnerability reporting and scope

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
├── include/          # Public headers
├── algorithms/       # Core FFT implementations (core/, dft/, auto/)
├── utils/            # Shared helpers, window functions, GPU stubs
├── applications/     # Real-world applications
├── benchmarks/       # Performance testing
├── tests/            # Test suite (test_all.c)
├── gpu/              # CUDA and Metal sources, outside the default build
├── cmake/            # CMake package configuration
├── docs/             # Documentation
├── .github/          # CI, Dependabot, issue and pull-request templates
├── SECURITY.md       # Private vulnerability reporting policy
└── CONTRIBUTING.md   # Contribution and review workflow
```

There is no `optimizations/` directory and no `examples/` directory. The former
held unreachable SIMD, threading and fixed-point kernels and was deleted; the
example programs live in `applications/`.

## Why This Project?

This project serves multiple purposes:

1. **Educational Resource**: Learn FFT algorithms with clear, documented code
2. **Reference Implementation**: Production-quality code with proper error handling
3. **Performance Comparison**: Benchmark different algorithms
4. **Practical Examples**: See FFT used in real applications

## License

This project is licensed under the Apache License 2.0. See [LICENSE](https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/LICENSE) for details.
