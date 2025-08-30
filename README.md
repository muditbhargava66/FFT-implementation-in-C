<div align="center">

# FFT Implementation in C

[![Documentation](https://readthedocs.org/projects/fft-implementation-in-c/badge/?version=latest)](https://fft-implementation-in-c.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue)]()
[![GPU Support](https://img.shields.io/badge/GPU-CUDA%20%7C%20Metal-green)]()

**A comprehensive, production-ready Fast Fourier Transform (FFT) library with automatic algorithm selection, GPU acceleration, and cross-platform support.**

</div>

## 🚀 What's New in v2.0.0

- **Automatic Algorithm Selection**: New intelligent API that chooses the best algorithm
- **GPU Acceleration**: CUDA support for NVIDIA GPUs and Metal Performance Shaders for Apple Silicon
- **Redesigned API**: Simplified interface with `fft_auto()` for ease of use
- **Cross-Platform**: Full compatibility across Linux, macOS, and Windows
- **Bug Fixes**: Critical Bluestein algorithm fix for prime-sized transforms

## 📚 Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Performance](#performance)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## ✨ Features

### Core Algorithms
- **Radix-2 DIT/DIF**: Classic Cooley-Tukey implementations
- **Radix-4**: Higher radix for 25% fewer operations  
- **Split-Radix**: Optimal operation count
- **Bluestein**: Arbitrary size FFTs (with v2.0.0 bug fixes)
- **Mixed-Radix**: Efficient for composite sizes
- **Automatic Selection**: Let the library choose the best algorithm

### GPU Acceleration (New in v2.0.0)
- **NVIDIA CUDA**: Massive speedups on NVIDIA GPUs
- **Apple Metal**: Optimized for M1/M2/M3 processors
- **Automatic GPU Detection**: Falls back to CPU if GPU unavailable

### Applications
- Audio spectrum analysis with windowing
- Digital filtering (low-pass, high-pass, band-pass)
- Fast convolution
- Power spectrum estimation
- 2D image FFT processing

### Optimizations
- SIMD vectorization (SSE, AVX, AVX-512, NEON)
- Multi-threaded execution with OpenMP
- Cache-optimized memory access
- Fixed-point arithmetic for embedded systems

## 🏃 Quick Start

### Simple FFT (v2.0.0 API)

```c
#include <fft_auto.h>

int main() {
    int n = 1024;
    complex_t* signal = fft_alloc_complex(n);
    
    // Generate signal
    for (int i = 0; i < n; i++) {
        signal[i] = sin(2 * PI * 50 * i / 1000.0);  // 50 Hz
    }
    
    // Automatic FFT - chooses best algorithm and uses GPU if available
    fft_auto(signal, signal, n, -1);  // -1 for forward FFT
    
    // Process results...
    
    fft_free(signal);
    return 0;
}
```

### With Planning (Advanced)

```c
// Create optimized plan
fft_plan_t plan = fft_plan_dft_1d(n, signal, signal, -1, 
                                   FFT_MEASURE | FFT_PREFER_GPU);

// Execute multiple times with same plan
for (int i = 0; i < 1000; i++) {
    generate_signal(signal, n);
    fft_execute(plan);
    process_spectrum(signal, n);
}

fft_destroy_plan(plan);
```

## 💾 Installation

### Prerequisites

- C compiler with C99 support (GCC 4.8+, Clang 3.4+, MSVC 2015+)
- Optional: CUDA Toolkit 11.0+ for NVIDIA GPU support
- Optional: Xcode 12+ for Metal support on macOS

### Build from Source

```bash
# Clone repository
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

# Quick build
./quickstart.sh

# Or manual build
make all              # Build everything
make gpu-demo         # Build GPU demonstrations
make install          # Install system-wide
```

### Platform-Specific Notes

#### macOS
```bash
# For OpenMP support
brew install gcc
export CC=gcc-13

# Metal support is automatic on Apple Silicon
```

#### Linux
```bash
# For CUDA support
# Install CUDA Toolkit from NVIDIA

# Build with GPU support
make all
```

#### Windows
```bash
# Use WSL or MinGW
# Visual Studio project coming soon
```

## 📖 Usage

### Basic Usage

```c
#include <fft_auto.h>

// Allocate aligned memory
complex_t* data = fft_alloc_complex(1024);

// Perform FFT with automatic optimization
fft_auto(data, data, 1024, FFT_FORWARD);

// For inverse FFT
fft_auto(data, data, 1024, FFT_INVERSE);

// Free memory
fft_free(data);
```

### GPU Acceleration

```c
// Check GPU availability
if (fft_gpu_available()) {
    printf("GPU: %s\n", fft_gpu_get_device_name());
}

// Force GPU usage
fft_plan_t plan = fft_plan_dft_1d(n, in, out, -1, FFT_PREFER_GPU);
```

### Real-valued FFT

```c
double* real_signal = fft_alloc_real(1024);
complex_t* spectrum = fft_alloc_complex(513);  // n/2 + 1

fft_plan_t plan = fft_plan_r2c_1d(1024, real_signal, spectrum, FFT_ESTIMATE);
fft_execute(plan);
```

## 📊 Performance

### Benchmark Results (Intel i9-12900K + RTX 3090)

| Size | CPU (AVX2) | GPU (CUDA) | Speedup |
|------|------------|------------|---------|
| 1K   | 0.08 ms    | 0.02 ms    | 4x      |
| 16K  | 1.8 ms     | 0.15 ms    | 12x     |
| 256K | 35 ms      | 1.2 ms     | 29x     |
| 1M   | 150 ms     | 4.5 ms     | 33x     |

### Apple M2 Max Performance

| Size | CPU (NEON) | GPU (Metal) | Speedup |
|------|------------|-------------|---------|
| 1K   | 0.06 ms    | 0.03 ms     | 2x      |
| 16K  | 1.2 ms     | 0.20 ms     | 6x      |
| 256K | 28 ms      | 2.1 ms      | 13x     |

## 📚 Documentation

Full documentation is available at: https://fft-implementation-in-c.readthedocs.io/

- [Getting Started Guide](https://fft-implementation-in-c.readthedocs.io/en/latest/getting-started/)
- [API Reference](https://fft-implementation-in-c.readthedocs.io/en/latest/api-reference/)
- [Algorithm Details](https://fft-implementation-in-c.readthedocs.io/en/latest/algorithms/)
- [GPU Programming Guide](https://fft-implementation-in-c.readthedocs.io/en/latest/gpu-guide/)

## 🤝 Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.

### Areas for Contribution
- Additional GPU backends (OpenCL, ROCm)
- More algorithms (Prime Factor, Winograd)
- Language bindings (Python, Julia, Rust)
- Performance optimizations

## 🔄 Migration from v1.x

The v2.0.0 API is mostly backward compatible. Key changes:

```c
// Old API (v1.x)
radix2_dit_fft(signal, n, FFT_FORWARD);

// New API (v2.0) - automatic optimization
fft_auto(signal, signal, n, -1);

// Or use planning for repeated transforms
fft_plan_t plan = fft_plan_dft_1d(n, signal, signal, -1, FFT_MEASURE);
fft_execute(plan);
```

## 🙏 Acknowledgments

- Original FFT algorithm by Cooley and Tukey
- Inspired by FFTW's planning approach
- GPU implementations based on cuFFT and Metal Performance Shaders
- Community contributors and testers

---

**Note**: This is a high-performance library suitable for production use. For educational purposes, explore the `algorithms/core/` directory for well-documented implementations.





<div align="center">

## Star History

<a href="https://www.star-history.com/#muditbhargava66/FFT-implementation-in-C&Date">
 <picture>
   <source media="(prefers-color-scheme: dark)" srcset="https://api.star-history.com/svg?repos=muditbhargava66/FFT-implementation-in-C&type=Date&theme=dark" />
   <source media="(prefers-color-scheme: light)" srcset="https://api.star-history.com/svg?repos=muditbhargava66/FFT-implementation-in-C&type=Date" />
   <img alt="Star History Chart" src="https://api.star-history.com/svg?repos=muditbhargava66/FFT-implementation-in-C&type=Date" />
 </picture>
</a>

---  

⭐️ Star the repo and consider contributing!  
  
📫 **Contact**: [@muditbhargava66](https://github.com/muditbhargava66)
🐛 **Report Issues**: [Issue Tracker](https://github.com/muditbhargava66/FFT-implementation-in-C/issues)
📚 [Documentation](https://fft-implementation-in-c.readthedocs.io/)
💬 [Discussions](https://github.com/muditbhargava66/FFT-implementation-in-C/discussions)
  
© 2025 Mudit Bhargava. [MIT License](LICENSE)  
<!-- Copyright symbol using HTML entity for better compatibility -->
</div>