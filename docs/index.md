# FFT Implementation in C

A C11 FFT library with eight algorithm implementations, a planning API, and
algorithm selection based on measured cost. It builds with GCC and Clang on
Linux, macOS and Windows (MinGW-w64). MSVC is not supported.

```c
#include <fft_auto.h>

int main(void) {
    int n = 1024;
    complex_t *signal = fft_alloc_complex(n);
    for (int i = 0; i < n; i++)
        signal[i] = sin(2 * PI * 50 * i / 1000.0);

    fft_auto(signal, signal, n, FFT_FORWARD);
    fft_free(signal);
    return 0;
}
```

## What the library provides

**Kernels:** radix-2 DIT, radix-2 DIF, radix-4, split-radix, Bluestein,
mixed-radix, recursive radix-2, and naive/optimized DFT references.

**Planning:** one-shot calls and reusable plans for 1-D complex, real-to-complex,
complex-to-real and 2-D transforms. The planner chooses by transform size, or by
timing the candidates when `FFT_MEASURE` is set. A wisdom cache can save and
restore those choices.

**Errors:** every transform returns `fft_error_t` and leaves the caller's buffer
untouched when it fails. The library never exits the process.

## What it does not do

- No SIMD, threaded or fixed-point kernels.
- No GPU execution. The CUDA and Metal sources are outside the default build;
  `fft_gpu_*` resolves to CPU stubs and `fft_gpu_available()` returns 0.
- Not thread-safe. Algorithm selection and the wisdom cache use global state.
- 2-D transforms require power-of-two dimensions.

## Build

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The library is `build/lib/libfft_lib_c.a`; the test and benchmark binaries are
in `build/bin`. See [Getting Started](getting-started.md) for platform-specific
notes and [Windows Support](windows.md) for MinGW details.

## Verification

The suite compares every implementation against a direct DFT computed
independently inside the test file, over 30 sizes from 1 to 16384. Release
3.1.0 was run locally on Ubuntu 24.04 (GCC 13.3.0, Clang 18.1.3) and macOS 27
arm64 (Apple Clang 21.0.0, Homebrew GCC 16.2.0): 1361/1361 checks, no compiler
warnings, and 52 passed / 0 failed / 2 skipped benchmark rows where the
benchmark was run. The 3.0.0 release was measured on Windows 11 with MinGW-w64
GCC 13.2.0 with the same counts. The 3.1.0 CI run passed on Linux GCC/Clang,
macOS Apple Clang, Windows MinGW-w64 GCC 16.2.0 and clang 23.1.3, including the
benchmark pass/fail check and the application demos. Details are in
[Performance](performance.md).

## Documentation

- [Getting Started](getting-started.md)
- [Quick Start](quickstart.md)
- [Algorithm Guide](algorithms.md)
- [Applications](applications.md)
- [API Reference](api-reference.md)
- [Performance Guide](performance.md)
- [Optimizations](optimizations.md)
- [Security Policy](https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/SECURITY.md)

## Project structure

```
include/        public headers
algorithms/     kernels and planner (core/, dft/, auto/)
utils/          window functions, helpers, GPU stubs
applications/   runnable demos
benchmarks/     benchmark_all and gpu_benchmark
tests/          test_all.c
gpu/            CUDA and Metal sources, not in the default build
docs/           this site
```

## Background

The DFT is

$$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-j2\pi kn/N},$$

and the FFT factors it into $O(N\log N)$ operations by reusing twiddle
factors. The [Algorithm Guide](algorithms.md) describes the decompositions used
here and when each one is worth choosing.

## License

Apache License 2.0. See [LICENSE](https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/LICENSE).
