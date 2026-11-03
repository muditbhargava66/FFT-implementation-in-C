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

Release 3.4.0 was verified locally on Windows 11 with MinGW-w64 GCC 13.2.0:
1406/1406 `test_all` checks and 54/54 utility checks in a default build,
1415/1415 and 17712 threaded checks with `-DFFT_ENABLE_OPENMP=ON`, both
warning-free with `-Werror`; CTest 2/2 and 3/3; all nine demos pass;
`benchmark_all` reports 66 passed / 0 failed / 54 skipped (68/0/52 with OpenMP).
On Ubuntu 22.04 under WSL, GCC 11.4.0 passed the default, OpenMP, `-Werror` and
ASan/UBSan builds, and the CUDA 11.5.119 build on a GTX 1050 passed 74/74 GPU
checks with CTest 3/3. The strict MkDocs build passes. The 3.4.0 CI workflow
adds an OpenMP lane to the previous eight jobs, for nine total; Linux GCC and
Clang, macOS 15 arm64 with Metal, ASan/UBSan, Windows MinGW-w64 GCC and Clang,
the CUDA 12.6 compile-and-skip container lane, and strict docs are unchanged.
The first push run (`cc15dca`, run 41640913349) passed all nine jobs.
The macOS job is pinned to `macos-15` because `macos-latest` now resolves to
macOS 26 Arm64 and queued during the macOS 14 retirement migration. The 3.2.0
macOS and Linux results are retained in the README history; its Metal build
passed 70/70 GPU checks on an M1 Pro. The CPU performance ratios come from the
WSL 3.4.0 build; the CUDA numbers come from the GTX 1050.

## Background

The DFT is

$$X[k] = \sum_{n=0}^{N-1} x[n]\, e^{-j2\pi kn/N},$$

and the FFT factors it into $O(N\log N)$ operations by reusing twiddle
factors. The [Algorithm Guide](algorithms.md) describes the decompositions used
here and when each one is worth choosing.
