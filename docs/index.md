# FFT Implementation in C

A C11 FFT library with seven FFT kernels, two DFT references, a planning API,
and algorithm selection based on measured cost. GCC and Clang on Linux, macOS
and Windows (MinGW-w64); MSVC's `cl.exe` is not supported because it has no C11
`_Complex` arithmetic.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The library is `build/lib/libfft_lib_c.a`. `test_all`, `test_utils`,
`benchmark_all` and the demos are in `build/bin`. A `-DFFT_ENABLE_METAL=ON`
build on macOS 14+ also produces `test_gpu_metal`.

## Start here

- [Getting Started](getting-started.md): requirements and a first build
- [Quick Start](quickstart.md): plans, wisdom, r2c/c2r and 2-D transforms
- [Algorithm Guide](algorithms.md): what each kernel does and how the planner
  chooses
- [API Reference](api-reference.md): public functions and the error contract
- [Performance](performance.md): measured ratios and how they were taken

## Verification

Release 3.2.0 was run locally on macOS 27 arm64 (Apple Clang 21.0.0, Homebrew
GCC 16.2.0) and Ubuntu 24.04 arm64 (GCC 13.3.0, Clang 18.1.3): 1361/1361 checks
in `test_all`, 40/40 in `test_utils`, no compiler warnings with `-Werror`, and
52 passed / 0 failed / 2 skipped benchmark rows. The Metal build on the same M1
Pro passed 70/70 GPU checks, also with `-Werror` for C and Objective-C. CI runs
both test binaries, the benchmark pass/fail check and the application demos on
Linux, macOS and Windows MinGW-w64, and builds and runs the Metal backend on
macOS. The performance ratios come from a Windows 11 machine running the 3.0.0
kernels; the kernels have not changed since.

## Background

The DFT is

$$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-j2\pi kn/N},$$

and the FFT factors it into $O(N\log N)$ operations by reusing twiddle
factors. The [Algorithm Guide](algorithms.md) describes the decompositions used
here and when each one is worth choosing.
