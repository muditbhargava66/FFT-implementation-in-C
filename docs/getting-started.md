# Getting Started

This guide will help you get started with the FFT Implementation in C project.

## Prerequisites

Before building the project, ensure you have the following installed:

- **C Compiler**: GCC 4.8+ or Clang 3.4+ (with C99 support)
- **Make**: GNU Make 3.81+
- **Math Library**: Standard C math library (libm)
- **Optional**: OpenMP support for parallel implementations

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
    
    # For OpenMP support, install GCC via Homebrew
    brew install gcc
    ```

=== "Windows (WSL)"
    ```bash
    # In WSL terminal
    sudo apt-get update
    sudo apt-get install build-essential
    ```

## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C
```

### 2. Quick Build

The easiest way to get started is using the quickstart script:

```bash
chmod +x quickstart.sh
./quickstart.sh
```

Select option 1 to build everything for the first time.

### 3. Manual Build

Alternatively, you can use Make directly:

```bash
# Build everything
make all

# Build specific components
make algorithms     # Core FFT algorithms only
make applications  # Applications only
make optimizations # Optimized versions only
make tests        # Test suite only
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
gcc -I./include my_first_fft.c -L./lib -lfft -lm -o my_first_fft
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
   gcc -I/path/to/include your_program.c -L/path/to/lib -lfft -lm
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

- `bin/`: Executable demonstrations and examples
- `lib/`: Static library file (`libfft.a`)
- `obj/`: Object files (intermediate build artifacts)

## Common Issues and Solutions

### Issue: OpenMP Not Found (macOS)

On macOS, the default Clang compiler doesn't support OpenMP. Solutions:

1. Install GCC via Homebrew:
   ```bash
   brew install gcc
   export CC=gcc-13  # Use appropriate version
   ```

2. Or build without OpenMP:
   ```bash
   make OPENMP_FLAGS=""
   ```

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

Some warnings are expected due to optimization flags. Critical warnings are errors and will prevent compilation.

## Next Steps

- Explore the [Algorithm Guide](algorithms.md) to understand different FFT implementations
- Check out [Applications](applications.md) for real-world examples
- Review the [API Reference](api-reference.md) for detailed function documentation
- Run benchmarks to compare performance: `./bin/benchmark_all`
