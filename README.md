<div align="center">

# FFT Implementation in C

[![Documentation](https://readthedocs.org/projects/fft-implementation-in-c/badge/?version=latest)](https://fft-implementation-in-c.readthedocs.io/)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue)]()
[![GPU Support](https://img.shields.io/badge/GPU-source--only-lightgrey)]()
[![Release](https://img.shields.io/badge/release-3.0.0%20(2026--10--23)-success)]()

**A Fast Fourier Transform library in C11: eight algorithms, automatic
algorithm selection backed by measurement, a planning API, and cross-platform
builds for GCC and Clang.**

</div>

## What's New in v3.0.0

- **Real-to-Complex & Complex-to-Real FFT**: Proper `fft_plan_r2c_1d` / `fft_plan_c2r_1d` with Hermitian-symmetric output
- **2D FFT Planning**: `fft_plan_dft_2d` via row-column decomposition with internal transpose
- **Wisdom Cache**: Save and reload algorithm selections with `fft_export_wisdom_to_string` / `fft_import_wisdom_from_string`
- **2D FFT Convolution**: Full `fft_convolution_2d` using row-column FFT decomposition
- **Cross-PSD Coherence**: Proper magnitude-squared coherence using Welch's method
- **Measured Algorithm Selection**: The size rules and the `FFT_MEASURE` timing sweep are both backed by recorded numbers, and the chosen algorithm is readable with `fft_plan_algorithm_name()`
- **1361/1361 Checks Passing**: 100% pass rate across all 8 algorithm implementations, verified against an independent direct DFT
- **Validated Planning Wisdom**: Cached and imported algorithm choices are checked against the transform size and the machine's hardware; `FFT_CONSERVE_MEMORY` selects an in-place kernel and is no longer overridden by `FFT_MEASURE`
- **Error codes instead of `exit()`**: Every transform returns `fft_error_t`; linking the library no longer puts your process at risk of being terminated by it
- **CI on Linux and Windows**: GCC and Clang lanes on Linux, MinGW-w64 on Windows, each failing on compiler warnings and checking the benchmark's pass/fail line
- **License**: Moved from MIT to Apache 2.0 for patent protection

## Table of Contents

- [Features](#features)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
- [Performance](#performance)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [License](#license)

## Features

### Core Algorithms
- **Radix-2 DIT/DIF**: Classic Cooley-Tukey implementations
- **Radix-4**: Higher radix for 25% fewer operations
- **Split-Radix**: Fewest multiplications of the group, though 1.6-2.9x slower in wall time here because the recursion moves about 2.5n elements per level
- **Bluestein**: Arbitrary size FFTs
- **Mixed-Radix**: True Cooley-Tukey for composite sizes; efficient when memory matters more than speed
- **Automatic Selection**: Let the library choose the best algorithm

### Planning & Execution API
- **1D Complex FFT**: `fft_plan_dft_1d` with automatic algorithm selection
- **Real-to-Complex FFT**: `fft_plan_r2c_1d` for real-valued signals
- **Complex-to-Real FFT**: `fft_plan_c2r_1d` for inverse real-valued transforms
- **2D Complex FFT**: `fft_plan_dft_2d` via row-column decomposition
- **Wisdom Cache**: Save/restore algorithm selections across sessions

### GPU Acceleration (source present, not functional)
- **NVIDIA CUDA**: `gpu/fft_cuda.cu`. `-DFFT_ENABLE_GPU=ON` compiles the CUDA
  objects, but does not connect them to the C planner.
- **Apple Metal**: `gpu/fft_metal.m`, outside the default build
- **Stubbed in the default build**: `utils/fft_gpu_stubs.c` supplies every
  `fft_gpu_*` symbol, so `fft_gpu_available()` returns 0 and `FFT_PREFER_GPU`
  falls back to a CPU algorithm

The GPU dispatch calls in `fft_auto.c` are compiled only for Apple or
CUDA-compiler builds, and no backend has been run against hardware since it left
the default build. This repository therefore publishes no GPU timings.

### Applications
- Audio spectrum analysis with windowing
- Digital filtering (low-pass, high-pass, band-pass)
- Fast 1D and 2D convolution
- Power spectrum estimation with coherence analysis
- 2D image FFT processing

### Optimizations

What the library actually does to go fast, and what it does not:

- **Algorithm selection**: `select_algorithm()` picks a kernel from measured
  cost ratios rather than from folklore, and `FFT_MEASURE` times the candidates
  on the host and caches the winner.
- **Wisdom cache**: `FFT_WISDOM_ONLY` reuses a cached choice and refuses to plan
  if there is none, so planning never measures. Exportable and importable as
  text via `fft_export_wisdom_to_string` / `fft_import_wisdom_from_string`.
- **Twiddle accumulation**: the radix kernels build twiddles multiplicatively
  per stage instead of calling `cexp()` per butterfly.
- **Aligned allocation**: 64-byte aligned buffers from `fft_alloc_complex`.
- **No SIMD, threading, fixed-point or GPU execution.** Earlier revisions of
  this README listed SSE2/AVX kernels, a pthreads/OpenMP stage splitter and a
  Q15 fixed-point path. Those files compiled into `libfft_lib_c.a` but were
  declared in no header and called by nothing, so no plan could ever reach
  them; they have been deleted rather than left as unreachable objects. The GPU
  sources under `gpu/` are outside the default build; `FFT_ENABLE_GPU=ON`
  compiles the CUDA objects without wiring them to the planner, so a default
  build is CPU-only. No GPU number in this repository has been measured.
- **Not thread-safe**: algorithm selection and the wisdom cache share global
  state with no locking, so do not plan concurrently.

## Quick Start

### Simple FFT (One-Shot API)

```c
#include <fft_auto.h>

int main() {
    int n = 1024;
    complex_t* signal = fft_alloc_complex(n);
    
    // Generate signal
    for (int i = 0; i < n; i++)
        signal[i] = sin(2 * PI * 50 * i / 1000.0);  // 50 Hz
    
    // Automatic FFT — chooses the best CPU algorithm for this machine
    fft_auto(signal, signal, n, -1);  // -1 for forward FFT
    
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

### Real-to-Complex FFT (New in v3.0.0)

```c
double* real_signal = fft_alloc_real(1024);
complex_t* spectrum = fft_alloc_complex(513);  // n/2 + 1

fft_plan_t plan = fft_plan_r2c_1d(1024, real_signal, spectrum, FFT_ESTIMATE);
fft_execute(plan);
fft_destroy_plan(plan);
```

### Wisdom Cache (New in v3.0.0)

```c
// Export wisdom after planning
char* wisdom = fft_export_wisdom_to_string();
save_to_file("fft_wisdom.txt", wisdom);
free(wisdom);

// Import wisdom in future runs
char* saved = load_from_file("fft_wisdom.txt");
fft_import_wisdom_from_string(saved);
free(saved);
```

## Installation

### Prerequisites

- C compiler with C11 `complex.h` and `_Complex` support: GCC 4.8+, Clang 3.4+.
  `cl.exe` is not supported — see *Platform support* under Installation
- CMake 3.10+
- Optional: CUDA Toolkit 11.0+ for NVIDIA GPU support
- Optional: Xcode 12+ for Metal support on macOS

### Build from Source

```bash
# Clone repository
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

# Configure and build
cmake -S . -B build
cmake --build build --config Release

# Run tests
./build/bin/test_all    # 1361/1361 checks expected
```

#### macOS
```bash
cmake -S . -B build
cmake --build build
./build/bin/test_all
```

#### Platform support

| Platform | Compiler           | Status                                                            |
| -------- | ------------------ | ----------------------------------------------------------------- |
| Windows  | MinGW-w64 GCC      | Verified: 1361/1361 checks, no warnings                           |
| Windows  | clang / clang-cl   | Expected to work, not verified                                    |
| Linux    | GCC, Clang         | Target platform, not verified in this repository's recent history |
| macOS    | Clang (Apple), GCC | Target platform, not verified                                     |
| Windows  | MSVC `cl.exe`      | **Not supported**                                                 |

Only the MinGW-w64 line above was executed while these docs were last updated.
The others describe what the code targets — it is standard C11 plus POSIX where
guarded — and are not claims about measured behaviour. If you build on one of
them, `ctest` is the check.

`cl.exe` cannot build this library. Every kernel is written against C11 complex
arithmetic (`x[i] = x[j] * w`, `a += b`, `x /= n`), which needs the `_Complex`
keyword and the operators the standard defines for it. Microsoft's C compiler
has neither: its `complex.h` offers structure types such as `_Dcomplex` and
states that the arithmetic operators are not defined for them. Supporting it
would mean rewriting every kernel against `_Cmulcc`-style helpers, so CMake now
stops with that explanation instead of failing partway through the build.

```cmd
cmake -S . -B build
cmake --build build
.\build\bin\test_all.exe
```

## Usage

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

### 2D FFT (New in v3.0.0)

```c
int rows = 256, cols = 256;
complex_t* image = fft_alloc_complex(rows * cols);

fft_plan_t plan = fft_plan_dft_2d(rows, cols, image, image, -1, FFT_ESTIMATE);
fft_execute(plan);
fft_destroy_plan(plan);
```

## Performance

### Test Results — 1361/1361 checks passing

Every implementation is compared against a direct O(n^2) DFT computed inside the
test file with no shared code, across 30 sizes from 1 to 16384 including primes,
`2 * 4^k` and sizes above 1024. See [docs/performance.md](docs/performance.md)
for the full measurement method.

### Release 3.0.0 verification (2026-10-23)

| Check                                                  | Result                                   |
| ------------------------------------------------------ | ---------------------------------------- |
| Clean build, MinGW-w64 GCC 13.2.0, `-Wall -Wextra -O3` | 0 warnings, 0 errors                     |
| Correctness suite (`test_all`)                         | 1361/1361 checks                         |
| `ctest`                                                | 1/1 test passed                          |
| `benchmark_all`                                        | 52 passed, 0 failed, 2 skipped in 17.1 s |
| Example programs (9 executables)                       | all exit 0                               |
| `realtime_analyzer` launches                           | 5/5, 0 crashes                           |

The skipped benchmark rows are the two implementations that cannot run at the
largest size in that harness; they are counted, not silently dropped.

### Relative cost, all 8 kernels

Each cell is that kernel's time divided by **radix-2 DIT at the same n**, so
below 1.00 is faster. Measured on an Intel Core i7-9750H with MinGW-w64 GCC
13.2.0 and the project's default `-O3` (no `-march=native`; that is opt-in via
`-DFFT_NATIVE_ARCH=ON`).

| n     | Radix-2 DIF | Radix-4  | Split-Radix | Mixed-Radix | Bluestein |
| ----- | ----------- | -------- | ----------- | ----------- | --------- |
| 16    | 1.06        | 1.41     | 1.64        | 0.99        | 11.11     |
| 64    | 1.06        | 1.21     | 2.44        | 1.02        | 11.30     |
| 256   | 1.00        | 1.10     | 2.38        | 0.99        | 11.90     |
| 512   | 0.95        | **0.83** | 2.00        | 1.01        | 10.34     |
| 1024  | 0.97        | 0.97     | 2.30        | 1.04        | 11.41     |
| 4096  | 1.02        | **0.89** | 1.99        | 1.03        | 9.54      |
| 16384 | 1.00        | **0.93** | 2.86        | 1.01        | 9.68      |

Ratios rather than milliseconds, because absolute times are not reproducible on
this machine: three runs of one probe put radix-2 DIT at n=16384 at 0.501, 1.372
and 1.433 ms — a 2.9x spread for byte-identical code. See
[docs/performance.md](docs/performance.md) for the method and the per-cell ranges.

What the numbers support, and what they do not:

- **The Mixed-Radix column is the control.** `mixed_radix_fft()` forwards
  power-of-two lengths straight to `radix2_dit_fft()`, so it has to read 1.00 —
  and it lands at 0.99-1.04. That is the resolution of the whole method, so
  differences smaller than it are not measurements.
- **Radix-2 DIT and DIF are indistinguishable**, 0.95-1.07 at every size.
- **Radix-4 crosses over between 256 and 512**: 1.41x and 1.21x slower at 16 and
  64, then 0.83x, 0.89x and 0.93x faster from 512 up. That is where
  `select_algorithm()` puts the threshold.
- Split-radix is reproducibly 1.6-2.9x slower despite the lowest multiplication
  count: its recursion moves about 2.5n elements per level.
- Bluestein pads the convolution to `next_power_of_two(2n-1)`, so using it where
  a radix kernel applies costs about 10x — the most stable column here.

### Non-power-of-two: Bluestein vs mixed-radix

| n               | Bluestein | Mixed-Radix | Bluestein / Mixed |
| --------------- | --------- | ----------- | ----------------- |
| 360 (2^3*3^2*5) | 0.24069   | 0.37200     | 0.65x             |
| 1000 (2^3*5^3)  | 0.46667   | 1.40107     | 0.33x             |
| 4095 (3*5*7*13) | 2.22222   | 9.50000     | 0.23x             |
| 97 (prime)      | 0.05233   | 0.84127     | 0.06x             |
| 1021 (prime)    | 0.46200   | 81.00000    | 0.01x             |

Bluestein wins at every length measured. Mixed-radix has no factor to divide by
at a prime and falls back to a direct O(n^2) DFT, which is where the two-order-
of-magnitude gaps come from. It stays reachable for callers that pass
`FFT_CONSERVE_MEMORY`, because Bluestein holds roughly 9n complex numbers of
scratch against mixed-radix's n. The same flag keeps power-of-two lengths on
in-place radix-2 DIT rather than radix-4, which permutes through an O(n) buffer.

## Documentation

Full documentation is available at: https://fft-implementation-in-c.readthedocs.io/

- [Getting Started Guide](https://fft-implementation-in-c.readthedocs.io/en/latest/getting-started/)
- [API Reference](https://fft-implementation-in-c.readthedocs.io/en/latest/api-reference/)
- [Algorithm Details](https://fft-implementation-in-c.readthedocs.io/en/latest/algorithms/)
- [Performance Guide](https://fft-implementation-in-c.readthedocs.io/en/latest/performance/)

## Contributing

We welcome contributions! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines.
Security reports should follow [SECURITY.md](SECURITY.md) and use private
reporting rather than a public issue.

### Areas for Contribution
- Additional GPU backends (OpenCL, ROCm)
- More algorithms (Prime Factor, Winograd)
- Language bindings (Python, Julia, Rust)
- Performance optimizations

## Migration from v2.x

The v3.0.0 API is fully backward compatible with v2.0.0. New additions:

```c
// v2.0.0 API (still works)
fft_auto(signal, signal, n, -1);
fft_plan_t plan = fft_plan_dft_1d(n, signal, signal, -1, FFT_MEASURE);

// New v3.0.0 additions
fft_plan_t r2c = fft_plan_r2c_1d(n, real_in, complex_out, FFT_ESTIMATE);
fft_plan_t c2r = fft_plan_c2r_1d(n, complex_in, real_out, FFT_ESTIMATE);
fft_plan_t fft2 = fft_plan_dft_2d(rows, cols, in, out, -1, FFT_ESTIMATE);

char* wisdom = fft_export_wisdom_to_string();
fft_import_wisdom_from_string(wisdom);
```

## Acknowledgments

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
  
**Contact**: [@muditbhargava66](https://github.com/muditbhargava66)
**Report Issues**: [Issue Tracker](https://github.com/muditbhargava66/FFT-implementation-in-C/issues)
[Documentation](https://fft-implementation-in-c.readthedocs.io/)
[Discussions](https://github.com/muditbhargava66/FFT-implementation-in-C/discussions)
  
© 2024-2026 Mudit Bhargava. [Apache License 2.0](LICENSE)  
</div>