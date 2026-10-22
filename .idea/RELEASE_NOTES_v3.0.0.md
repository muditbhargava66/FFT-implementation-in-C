# FFT Implementation in C v3.0.0

Released **2026-10-23**. This is a major release: the planning API is now the
supported interface, the kernels that carry an algorithm's name actually
implement it, and the documentation no longer claims features the build cannot
execute.

## Highlights

- **A planner that can be trusted.** Plans are built from validated kernel
  choices; an imported wisdom entry that names an impossible algorithm or size
  is rejected instead of producing a plan that fails every execution.
  `FFT_WISDOM_ONLY` is implemented and fails closed.
- **Real split-radix and real mixed-radix kernels.** Split-radix previously ran
  radix-2 DIT under a different name. Mixed-radix previously fell through to an
  O(n^2) DFT for non-power-of-two lengths; at N = 4095 that went from 1513 ms to
  5.90 ms.
- **Correct small transforms.** A byte-swap guard in `bit_reverse()` corrupted
  N = 4, 8 and 16 for radix-2 DIT/DIF and everything built on them.
- **Memory safety throughout.** The c2r path no longer reads past the caller's
  half-spectrum, allocations are checked, and the transform entry points return
  error codes instead of terminating the host process.
- **Honest performance numbers.** Timings are published as cost relative to
  radix-2 DIT at the same length, measured within one run. Absolute
  milliseconds on the reference machine vary by up to 2.9x between runs of the
  same binary.
- **Warning-clean CI** for Linux GCC, Linux Clang, Windows MinGW-w64, and a
  strict MkDocs build, plus Dependabot, issue forms, a PR template and a
  security policy.

## Verification

Windows 11, MinGW-w64 GCC 13.2.0, `-Wall -Wextra -O3`:

| Check | Result |
| --- | --- |
| Build | 0 warnings, 0 errors |
| `test_all` | 1361/1361 checks against an independent direct DFT |
| `ctest` | 1/1 passed |
| `benchmark_all` | 52 passed, 0 failed, 2 skipped in 17.1 s |
| Example programs | 9 built, all exit 0 |
| `realtime_analyzer` | 5/5 launches without a crash |
| `mkdocs build --strict` | passed |
| CI (Linux GCC, Linux Clang, Windows MinGW-w64, docs) | all four lanes passed |

The two skipped benchmark rows are reported by the harness, not dropped. The
Linux GCC suite was also run directly under WSL Ubuntu 22.04 with the same
1361/1361 result. A strict C99 build with GCC 15.3 in a container — the
toolchain generation from issue #4 — also passed all 1361 checks with implicit
function declarations promoted to errors. The Clang and hosted-runner lanes are
exercised by CI.

## Algorithm and planning notes

- Radix-2 DIT is the baseline; DIF is indistinguishable from it, and radix-4 is
  0.83x-0.97x its cost from N = 512 up.
- Split-radix has the lower operation count but is 1.6x-2.9x slower here,
  because its recursion moves about 2.5n elements per level.
- Bluestein remains the choice for non-power-of-two lengths, measured 2.1x to
  241x faster than mixed-radix.
- `FFT_MEASURE` times the candidates that accept the length, requires a 15%
  improvement over the heuristic plus a repeated paired confirmation, and
  caches the result in exported wisdom.

See `docs/performance.md` for the table with per-cell ranges and the method.

## Breaking changes

- **Error codes:** every transform returns `fft_error_t`. Ignoring the value
  still compiles, but callers that relied on `exit()` behaviour must handle
  `FFT_ERR_NULL`, `FFT_ERR_SIZE`, `FFT_ERR_ALLOC` and `FFT_ERR_RANGE`.
- **Removed flags:** `FFT_PATIENT`, `FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and
  `FFT_THREADED` were accepted and ignored. Their bit positions are left unused
  so existing exported wisdom and compiled callers keep their meaning.
- **MSVC is refused** at configure time. The library needs C11 `_Complex`
  arithmetic, which `cl.exe` does not provide. GCC, Clang and clang-cl build it.
- **License:** Apache License 2.0 replaces MIT, adding an explicit patent grant.
- **Build:** the legacy Makefile/batch/PowerShell scripts are gone; use CMake.
  `-march=native` is no longer a default; enable `-DFFT_NATIVE_ARCH=ON`
  deliberately.

## Known limitations

- GPU sources are shipped but no plan reaches them in a default build. Every
  `fft_gpu_*` symbol resolves to a stub; no GPU number in this repository was
  measured on hardware.
- There is no SIMD, threaded or fixed-point execution path. `optimizations/`
  and the duplicate iterative radix-2 kernel were removed rather than left
  unreachable.
- Planning and the wisdom cache use shared globals without locking; the library
  is not thread-safe.
- MSVC is unsupported.
- Issue #5 (Meson build support) is not part of this release.

## Links

- Documentation: https://fft-implementation-in-c.readthedocs.io/
- Full changelog: https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/CHANGELOG.md
- Compare v2.0.0...v3.0.0: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v2.0.0...v3.0.0
- Security policy: https://github.com/muditbhargava66/FFT-implementation-in-C/blob/main/SECURITY.md

**Full Changelog**: https://github.com/muditbhargava66/FFT-implementation-in-C/compare/v2.0.0...v3.0.0
