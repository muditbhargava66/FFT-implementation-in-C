# FFT Implementation in C

A C11 FFT library with seven FFT kernels, two DFT references, a planning API,
algorithm selection based on measured cost, and optional CUDA (Linux/WSL) and
Metal (macOS) backends. GCC and Clang on Linux, macOS and Windows (MinGW-w64);
MSVC's `cl.exe` is not supported because it has no C11 `_Complex` arithmetic.

```bash
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

The library is `build/lib/libfft_lib_c.a`. `test_all`, `test_utils`,
`benchmark_all` and the demos are in `build/bin`. A `-DFFT_ENABLE_METAL=ON`
build on macOS 14+ also produces `test_gpu_metal`; a `-DFFT_ENABLE_CUDA=ON`
build on Linux or WSL produces `test_gpu_cuda`.

## Start here

- [Getting Started](getting-started.md): requirements and a first build
- [Quick Start](quickstart.md): plans, wisdom, r2c/c2r and 2-D transforms
- [Algorithm Guide](algorithms.md): what each kernel does and how the planner
  chooses
- [API Reference](api-reference.md): public functions and the error contract
- [Performance](performance.md): measured ratios and how they were taken

## Verification

Release 3.3.0 was run locally on Windows 11 with MinGW-w64 GCC 13.2.0
(1362/1362 `test_all` checks, 54/54 utility checks, 0 warnings), and on Ubuntu
22.04 under WSL with GCC 11.4.0 and CUDA 11.5.119 on a GTX 1050
(1362/1362 + 54/54 + 74/74 CUDA checks, CTest 3/3, benchmarks 52 passed /
0 failed / 2 skipped, 0 warnings). The eight application demos and
`gpu_benchmark` pass under both configurations; the Windows demos ran from a
Release build with `-g` because application control blocked three optimized
binaries by hash. The 3.2.0 macOS and Linux
results are retained in the README history; its Metal build passed 70/70 GPU
checks on an M1 Pro. CI runs both CPU test binaries, the benchmark pass/fail
check and the demos on Linux, macOS and Windows MinGW-w64; 3.3.0 adds a CUDA
compile-and-skip lane because hosted runners have no GPU. The CPU performance
ratios come from a Windows 11 machine running the 3.0.0 kernels; the kernels
have not changed since.

## Background

The DFT is

$$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-j2\pi kn/N},$$

and the FFT factors it into $O(N\log N)$ operations by reusing twiddle
factors. The [Algorithm Guide](algorithms.md) describes the decompositions used
here and when each one is worth choosing.
