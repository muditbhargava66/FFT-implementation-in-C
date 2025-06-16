# FFT Study Repository

A comprehensive collection of Fast Fourier Transform (FFT) algorithms and implementations in C, designed for educational purposes and performance comparison.

## 📚 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Repository Structure](#repository-structure)
- [Algorithms Implemented](#algorithms-implemented)
- [Getting Started](#getting-started)
- [Building and Running](#building-and-running)
- [Code Quality](#code-quality)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Overview

This repository serves as a comprehensive study resource for understanding FFT algorithms, their implementations, optimizations, and applications. It includes various FFT algorithms from basic to advanced, real-world applications, and performance optimizations.

### Key Features

- **Educational Focus**: Each implementation includes detailed documentation explaining the algorithm, its complexity, and use cases
- **Production Quality**: Professional coding standards with proper error handling, memory management, and testing
- **Performance Analysis**: Comprehensive benchmarking suite for comparing different implementations
- **Real Applications**: Practical examples showing how FFT is used in audio processing, filtering, and more
- **Multiple Optimizations**: SIMD, parallel, and fixed-point implementations for different platforms

## Repository Structure

```
FFT-implementation-in-C/
├── include/               # Header files
│   ├── fft_common.h      # Common utilities and definitions
│   └── fft_algorithms.h  # FFT function declarations
├── algorithms/           # Core algorithm implementations
│   ├── core/            # FFT algorithms
│   └── dft/             # DFT implementations
├── applications/         # Real-world applications
├── optimizations/        # Performance optimizations
├── benchmarks/          # Performance testing
├── tests/               # Test suites
├── examples/            # Example programs
├── utils/               # Utility functions
└── lib/                 # Compiled library files
```

## Algorithms Implemented

### ✅ Core FFT Algorithms

| Algorithm | Description | Complexity | Status |
|-----------|-------------|------------|--------|
| Radix-2 DIT | Decimation-in-Time FFT | O(n log n) | ✅ |
| Radix-2 DIF | Decimation-in-Frequency FFT | O(n log n) | ✅ |
| Radix-4 | Base-4 FFT algorithm | O(n log n) | ✅ |
| Split-Radix | Optimal operation count | O(n log n) | ✅ |
| Bluestein | Arbitrary size FFT | O(n log n) | ✅ |
| Mixed-Radix | Composite size FFT | O(n log n) | ✅ |
| Recursive | Educational implementation | O(n log n) | ✅ |
| Iterative | Cache-efficient | O(n log n) | ✅ |

### 🔧 Optimizations

| Optimization | Target Platform | Speedup | Status |
|--------------|----------------|---------|--------|
| SIMD (SSE/AVX) | x86/x64 | 2-4x | ✅ |
| OpenMP Parallel | Multi-core | ~Nx cores | ✅ |
| Fixed-Point | Embedded/DSP | Platform-specific | ✅ |

### 📊 Applications

- **Audio Spectrum Analyzer**: Real-time frequency analysis
- **Digital Filtering**: Low-pass, high-pass, band-pass filters
- **Convolution**: Fast convolution using FFT
- **Power Spectrum**: Signal power analysis
- **2D Image FFT**: Image processing applications
- **Pitch Detection**: Musical note detection

## Getting Started

### Prerequisites

- GCC compiler (or compatible C compiler)
- Make build system
- POSIX-compliant system (Linux, macOS, WSL)
- Optional: OpenMP support for parallel implementations

### Quick Start

1. Clone the repository:
   ```bash
   git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
   cd FFT-implementation-in-C
   ```

2. Run the quick start script:
   ```bash
   chmod +x quickstart.sh
   ./quickstart.sh
   ```

3. Or build manually:
   ```bash
   make all
   ```

## Building and Running

### Build Commands

```bash
# Build everything
make all

# Build specific components
make algorithms      # Core FFT algorithms
make applications   # Application examples
make optimizations  # Optimized versions
make benchmarks     # Benchmark suite
make tests          # Test suite

# Build with debug symbols
make debug

# Build with profiling
make profile

# Clean build artifacts
make clean
```

### Running Examples

```bash
# Run a specific algorithm
./bin/radix2_dit

# Run benchmarks
./bin/benchmark_all

# Run tests
./bin/test_all

# Run demonstrations
make demo
```

### Library Usage

The project builds a static library `libfft.a` that can be linked with your programs:

```c
#include "fft_common.h"
#include "fft_algorithms.h"

int main() {
    int n = 1024;
    complex_t* signal = allocate_complex_array(n);
    
    // Generate signal...
    
    // Apply FFT
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Process spectrum...
    
    free_complex_array(signal);
    return 0;
}
```

Compile with:
```bash
gcc -I/path/to/include myprogram.c -L/path/to/lib -lfft -lm
```

## Code Quality

### Coding Standards

- **C99 Standard**: Modern C features for better code clarity
- **Consistent Style**: 4-space indentation, clear naming conventions
- **Comprehensive Documentation**: Doxygen-style comments for all functions
- **Error Handling**: Proper validation and error reporting
- **Memory Safety**: Careful allocation/deallocation, no memory leaks

### Algorithm Implementation

Each algorithm implementation includes:

1. **Detailed Header Documentation**: Mathematical background, complexity analysis, references
2. **Step-by-Step Comments**: Clear explanation of each algorithmic step
3. **Input Validation**: Checking for valid sizes, null pointers, etc.
4. **Test Cases**: Built-in testing in main() function
5. **Performance Measurement**: Timing and complexity verification

Example structure:
```c
/**
 * @brief Algorithm name and brief description
 * 
 * @details
 * Detailed explanation of the algorithm, including:
 * - Mathematical foundation
 * - Step-by-step process
 * - Data structures used
 * 
 * @param x Input/output array
 * @param n Array length
 * @param dir Transform direction
 * 
 * Time Complexity: O(n log n)
 * Space Complexity: O(1)
 */
void algorithm_fft(complex_t* x, int n, fft_direction dir) {
    /* Input validation */
    if (!x || n <= 0) {
        fprintf(stderr, "Error: Invalid input\n");
        return;
    }
    
    /* Algorithm implementation with clear comments */
    // Step 1: ...
    // Step 2: ...
}
```

## Documentation

### Algorithm Documentation

Each algorithm includes:
- Mathematical derivation
- Complexity analysis
- Implementation notes
- Usage examples
- Performance characteristics

### API Documentation

See `include/fft_common.h` and `include/fft_algorithms.h` for the complete API reference.

Key functions:
- `radix2_dit_fft()` - Standard FFT for power-of-2 sizes
- `bluestein_fft()` - FFT for arbitrary sizes
- `allocate_complex_array()` - Memory allocation helper
- `compute_magnitude()` - Convert to magnitude spectrum

## Performance

Benchmark results on typical hardware (Intel i7, 3.6GHz):

| Algorithm | N=1024 | N=4096 | N=16384 |
|-----------|--------|--------|---------|
| Radix-2 DIT | 0.08ms | 0.4ms | 1.8ms |
| Split-Radix | 0.06ms | 0.3ms | 1.4ms |
| SIMD Optimized | 0.03ms | 0.15ms | 0.7ms |
| Parallel (4 cores) | 0.02ms | 0.1ms | 0.5ms |

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

Areas for contribution:
- Additional algorithms (Prime Factor, Winograd)
- GPU implementations (CUDA/OpenCL)
- More applications
- Performance optimizations
- Documentation improvements

## References

1. Cooley, J. W., & Tukey, J. W. (1965). "An algorithm for the machine calculation of complex Fourier series"
2. Duhamel, P., & Vetterli, M. (1990). "Fast Fourier transforms: a tutorial review"
3. Frigo, M., & Johnson, S. G. (2005). "The design and implementation of FFTW3"
4. Van Loan, C. (1992). "Computational frameworks for the fast Fourier transform"

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

**Note**: This is an educational repository. For production use, consider established libraries like FFTW, Intel MKL, or cuFFT.
