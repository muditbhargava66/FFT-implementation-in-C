# Getting Started

## Prerequisites

- A C11 compiler with `_Complex` support: GCC 4.8+ or Clang 3.4+. On Windows,
  use MinGW-w64 or clang; MSVC's `cl.exe` is rejected at configure time because
  it has no C11 `_Complex` arithmetic.
- CMake 3.20 or newer.
- The C math library (`libm` on Linux; macOS resolves these functions from its
  C runtime).

=== "Linux"

    ```bash
    sudo apt-get update
    sudo apt-get install build-essential cmake
    ```

=== "macOS"

    ```bash
    xcode-select --install
    ```

    CMake is available from Homebrew (`brew install cmake`) or from cmake.org.

=== "Windows"

    Install MinGW-w64 and CMake, and run the build from an MSYS2/MinGW shell or
    with `-DCMAKE_C_COMPILER` pointing at `gcc.exe`. See
    [Windows Support](windows.md).

=== "Windows (WSL)"

    ```bash
    sudo apt-get update
    sudo apt-get install build-essential cmake
    ```

## Build and test

```bash
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

`ctest` runs `build/bin/test_all`, which prints a summary line. A successful run
ends with `1361/1361 checks passed`.

To choose a compiler, pass it at configure time:

```bash
cmake -S . -B build -DCMAKE_C_COMPILER=clang
cmake --build build
```

The default build uses `-O3`. `-march=native` is opt-in with
`-DFFT_NATIVE_ARCH=ON`; see [Optimizations](optimizations.md) for when that is a
good idea.

## First program

```c
/* my_first_fft.c */
#include <stdio.h>
#include "fft_common.h"
#include "fft_algorithms.h"

int main(void) {
    int n = 8;
    complex_t *signal = allocate_complex_array(n);
    if (!signal) return 1;

    for (int i = 0; i < n; i++)
        signal[i] = (i < n / 2) ? 1.0 : -1.0;

    print_complex_array("x", signal, n);
    radix2_dit_fft(signal, n, FFT_FORWARD);
    print_complex_array("X", signal, n);

    double *mag = compute_magnitude(signal, n);
    for (int i = 0; i < n; i++)
        printf("|X[%d]| = %.3f\n", i, mag[i]);

    free(mag);
    free_complex_array(signal);
    return 0;
}
```

Compile it against the built library:

```bash
gcc -I./include my_first_fft.c -L./build/lib -lfft_lib_c -lm -o my_first_fft
./my_first_fft
```

On Windows the library file is `libfft_lib_c.a` as well; use the same command
from a MinGW shell.

## Linking from your own project

```bash
gcc -I/path/to/include your_program.c \
    -L/path/to/build/lib -lfft_lib_c -lm
```

If you are using CMake, link the `fft_lib` target. There is no installed CMake
package configuration; add the source directory with `add_subdirectory()` or
link the static library and include directory directly.

## Basic pattern

```c
complex_t *data = allocate_complex_array(n);
/* fill data */
if (radix2_dit_fft(data, n, FFT_FORWARD) != FFT_OK) {
    /* handle a null buffer or an invalid length */
}
/* read the spectrum */
free_complex_array(data);
```

Every transform returns `fft_error_t`. If it fails after it has started, the
buffer may be partly modified, so do not use the result. Discarding the return
value is allowed in C, but checking it is the only way to tell a failed
transform from a successful one.

For repeated transforms of the same size, the planning API avoids re-deciding
the algorithm on every call:

```c
fft_plan_t plan = fft_plan_dft_1d(n, data, data, -1, FFT_ESTIMATE);
for (int i = 0; i < iterations; i++) {
    fill_signal(data, n);
    fft_execute(plan);
    consume(data, n);
}
fft_destroy_plan(plan);
```

## Directory layout after building

- `build/lib/libfft_lib_c.a`: the static library
- `build/bin/`: `test_all`, `benchmark_all`, `gpu_benchmark`, and the
  `applications/` demos

## Common problems

**Size rejected by radix-2 or radix-4.** Those kernels require a power of two.
Use Bluestein for an arbitrary length, or let `fft_auto()` choose.

**Undefined references to `cexp`, `cabs` or `sincos`.** Add `-lm` on Linux.
The CMake targets already do this through the library target.

**A warning appears during the build.** The project compiles with `-Wall
-Wextra` and is warning-clean on the compilers listed in the README. Treat a new
warning as something to investigate rather than expected noise.

## Next steps

- [Quick Start](quickstart.md) for a shorter tour
- [Algorithm Guide](algorithms.md) for the decompositions and their trade-offs
- [Applications](applications.md) for the runnable demos
- [API Reference](api-reference.md) for the full surface
