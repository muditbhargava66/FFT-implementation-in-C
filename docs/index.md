# FFT Implementation in C

A C11 FFT library with ten FFT kernels, two DFT references, a planning API,
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

Release 3.5.0 changes documentation, code comments and benchmark reporting; no
kernel or planner behavior changed. The 3.4.0 kernel set was verified on
Windows 11 with MinGW-w64 GCC 13.2.0 and on Ubuntu 22.04 under WSL with GCC
11.4.0, covering the default, OpenMP, `-Werror`, ASan/UBSan and CUDA builds. The
exact check and benchmark counts are in the README table. The strict MkDocs
build and the nine-job CI workflow pass; the macOS lane builds and tests the
Metal backend, and the CUDA lane compiles without a device. CPU performance
ratios come from the WSL 3.4.0 build and the CUDA numbers from the GTX 1050.

## Background

The DFT is

$$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-j2\pi kn/N},$$

and the FFT factors it into $O(N\log N)$ operations by reusing twiddle
factors. The [Algorithm Guide](algorithms.md) describes the decompositions used
here and when each one is worth choosing.
