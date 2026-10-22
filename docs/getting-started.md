# Getting Started

This guide will help you get started with the FFT Implementation in C project.

## Prerequisites

Before building the project, ensure you have the following installed:

- **C Compiler**: GCC 4.8+ or Clang 3.4+, with C11 `complex.h` and `_Complex`.
  MSVC's `cl.exe` is not supported; CMake stops with an explanation if it is
  selected. On Windows use MinGW-w64 or clang.
- **CMake**: 3.10+
- **Math Library**: Standard C math library (libm)

### Platform-Specific Notes

=== "Linux"
    ```bash
    # Install build essentials
    sudo apt-get update
    sudo apt-get install build-essential
    ```

=== "macOS"
    ```bash
    # Install Xcode Command Line Tools
    xcode-select --install
    ```

=== "Windows"
    ```cmd
    # CMake Build
    cmake -S . -B build
    cmake --build build
    ```

    Tested with MinGW-w64 GCC 13.2.0: 1361/1361 checks pass with no compiler
    warnings.

    **Prerequisites**:
    - MinGW-w64 (or clang). Visual Studio's `cl.exe` cannot build this library —
      it requires C11 `_Complex` arithmetic that MSVC does not provide, and CMake
      stops with that explanation if `cl.exe` is the selected C compiler
    - CMake 3.10+
    - Git for Windows
    
=== "Windows (WSL)"
    ```bash
    # In WSL terminal
    sudo apt-get update
    sudo apt-get install build-essential cmake
    ```

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C
```

### 2. Build with CMake

The project uses CMake for cross-platform builds:

```bash
# Configure
cmake -S . -B build

# Build
cmake --build build --config Release
```

### 3. Run Tests

```bash
cd build/bin
./test_all
```

## First Program

Create a simple program to test the FFT:

```c
// my_first_fft.c
#include <stdio.h>
#include "fft_common.h"
#include "fft_algorithms.h"

int main() {
    // Create a simple signal
    int n = 8;
    complex_t* signal = allocate_complex_array(n);
    
    // Initialize with a simple pattern
    for (int i = 0; i < n; i++) {
        signal[i] = (i < n/2) ? 1.0 : -1.0;
    }
    
    printf("Input signal:\n");
    print_complex_array("x", signal, n);
    
    // Perform FFT
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    printf("\nFFT result:\n");
    print_complex_array("X", signal, n);
    
    // Compute magnitude spectrum
    double* mag = compute_magnitude(signal, n);
    printf("\nMagnitude spectrum:\n");
    for (int i = 0; i < n; i++) {
        printf("|X[%d]| = %.3f\n", i, mag[i]);
    }
    
    // Clean up
    free(mag);
    free_complex_array(signal);
    
    return 0;
}
```

Compile and run:

```bash
gcc -I./include my_first_fft.c -L./build/lib -lfft_lib_c -lm -o my_first_fft
./my_first_fft
```

## Using the Library

### Linking with Your Project

To use the FFT library in your own project:

1. Include the headers:
   ```c
   #include "fft_common.h"
   #include "fft_algorithms.h"
   ```

2. Link with the library:
   ```bash
   gcc -I/path/to/include your_program.c -L/path/to/build/lib -lfft_lib_c -lm
   ```

### Basic Usage Pattern

```c
// 1. Allocate memory
complex_t* data = allocate_complex_array(size);

// 2. Fill with your data
generate_signal(data, size);

// 3. Apply FFT
radix2_dit_fft(data, size, FFT_FORWARD);

// 4. Process results
process_spectrum(data, size);

// 5. Optional: Inverse transform
radix2_dit_fft(data, size, FFT_INVERSE);

// 6. Clean up
free_complex_array(data);
```

## Directory Structure

After building, you'll find:

- `build/lib/`: the static library (`libfft_lib_c.a` with GCC/Clang)
- `build/bin/`: `test_all`, `benchmark_all`, and the `applications/` demos

## Common Issues and Solutions

### Issue: Power of 2 Requirement

Most FFT algorithms require input sizes to be powers of 2. Use:

```c
// Check if size is power of 2
if (!is_power_of_two(n)) {
    // Use Bluestein or Mixed-Radix for arbitrary sizes
    bluestein_fft(data, n, FFT_FORWARD);
}
```

### Issue: Compilation Warnings

The build compiles with `-Wall -Wextra` and is warning-clean, so any warning is
a real finding rather than expected noise. CI fails the build on the first one.

## IDE Integration

### Visual Studio Code

1. Install extensions:
   - C/C++
   - CMake Tools
   - CUDA (if using GPU)

2. Configure `.vscode/settings.json`:
   ```json
   {
       "C_Cpp.default.compilerPath": "C:/msys64/mingw64/bin/gcc.exe",
       "C_Cpp.default.includePath": [
           "${workspaceFolder}/include",
           "C:/Program Files/NVIDIA GPU Computing Toolkit/CUDA/v12.4/include"
       ]
   }
   ```

### Visual Studio Code (Windows)

Point the C/C++ extension at a MinGW-w64 or clang compiler as above. Visual
Studio's own `cl.exe` toolset cannot compile this library — see
[Prerequisites](#prerequisites).

## Performance Optimization

### Compiler Optimizations

The default build is already `-O3`. Two things worth knowing before adding flags:

```bash
# Tune for the CPU you are building on. Off by default: the binaries will not
# run on older CPUs, and MinGW-w64 GCC emits misaligned 32-byte AVX stack stores
# under it that crash some programs intermittently. See the comment above the
# FFT_NATIVE_ARCH option in CMakeLists.txt.
cmake -S . -B build -DFFT_NATIVE_ARCH=ON
cmake --build build --config Release
```

Do **not** add `-ffast-math`. It permits the compiler to reassociate the
butterfly sums and to assume no value is ever NaN or Inf, both of which are
load-bearing assumptions for an FFT; it trades correctness for a speedup this
library is not memory-bound enough to need.

### Windows-Specific Optimizations

1. **Disable Windows Defender real-time scanning** for the build directory
2. **Set high performance power plan**
3. **Use SSD** for faster compilation
4. **Increase RAM** if building large projects

## Next Steps

- Explore the [Algorithm Guide](algorithms.md) to understand different FFT implementations
- Check out [Applications](applications.md) for real-world examples
- Review the [API Reference](api-reference.md) for detailed function documentation
- Run benchmarks to compare performance: `./bin/benchmark_all`
- For Windows users: See [Windows Support](windows.md) for detailed Windows-specific information
