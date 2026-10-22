## Summary

Release **v3.0.0 (2026-10-23)** of FFT Implementation in C. This branch
consolidates the v3 planning API with the correctness, performance, build and
documentation work required to ship it. It fixes the planner bugs where a
non-NULL plan could never execute, replaces two kernels that did not implement
the algorithm in their name, removes code and claims that the build could not
reach, and puts the repository on a warning-clean CI footing.

Head: `56b9351` plus the audit/docs cleanup commit on `release/v3.0.0` (use the
branch tip). CI run `39836257536` passed all four lanes; the audit commit only
changes documentation. Two earlier runs failed and produced fixes on the
branch: the Linux targets now link `libm`, and the radix-4 demo helper is no
longer compiled into the library.

## What changed

**Correctness**

- `bit_reverse()` returned a wrong permutation for `log2n <= 4`, corrupting
  radix-2 DIT/DIF at N = 4, 8 and 16 (and Bluestein/mixed-radix at those sizes).
- `test_all` returned 0 even when checks failed, and had no independent DFT
  reference to compare against. Both are fixed; the suite now runs 1361 checks,
  including a direct-summation reference over a checked bin set.
- Fixed memory-safety defects: c2r over-read of the caller's half-spectrum,
  unchecked allocations in several transforms, `mixed_radix_fft_recursive()`
  ignoring its `stride`, `optimized_dft_real()` leaving a bin uninitialised for
  odd `n`, `goertzel_single_bin()` phase error, and `frequency_to_note()`
  indexing below C0.
- `FFT_CONSERVE_MEMORY` is honoured at power-of-two lengths and is no longer
  overridden by `FFT_MEASURE`.
- Wisdom entries are validated on import and on lookup. An unusable entry is
  skipped, and `FFT_WISDOM_ONLY` fails closed with `NULL` instead of returning a
  dead plan.

**Algorithms and performance**

- `split_radix_fft()` is now a real split-radix FFT (one `n/2` plus two `n/4`
  decomposition), not radix-2 DIT under another name.
- `mixed_radix_fft()` is now a real Cooley-Tukey mixed-radix FFT; N = 4095 went
  from 1513 ms to 5.90 ms.
- Radix-4 no longer calls `cexp()` per butterfly; it is now 0.83x-0.97x the
  cost of radix-2 DIT from N = 512 up.
- The planner threshold moved from 64 to 256 based on measured ratios, and
  `FFT_MEASURE` now performs an interleaved, batch-calibrated sweep.
- Published timings are now relative to radix-2 DIT within a single run,
  because absolute milliseconds on the test machine vary by up to 2.9x between
  runs of the same binary.

**API**

- Every transform returns `fft_error_t` instead of calling `exit()`.
- Added `fft_plan_r2c_1d`, `fft_plan_c2r_1d`, `fft_plan_dft_2d`,
  `fft_convolution_2d`, `compute_coherence`, and
  `fft_plan_algorithm_name()`.
- Removed accepted-but-ignored flags: `FFT_PATIENT`, `FFT_EXHAUSTIVE`,
  `FFT_UNALIGNED`, `FFT_THREADED`. `FFT_WISDOM_ONLY` is now implemented.
- License changed from MIT to Apache License 2.0.

**Build, CI and repository hygiene**

- Modern CMake build; MSVC is refused at configure time with an explanation
  because it cannot compile C11 `_Complex` arithmetic. The library target now
  propagates `libm` on Unix, which the first CI run proved was missing.
- Removed the unreachable `optimizations/` directory, `iterative_fft.c`, the
  legacy build scripts and the duplicate radix-2 kernel.
- The radix-4 demonstration helper is compiled only with its demo `main()`,
  clearing the Clang `-Wunused-function` warning that failed the first run.
- `-march=native` is opt-in (`FFT_NATIVE_ARCH`) after it exposed a MinGW-w64
  misaligned-AVX-store crash in `realtime_analyzer`.
- CI now covers Linux GCC, Linux Clang, Windows MinGW-w64 and a strict MkDocs
  build, runs on `release/**` pushes, and fails on compiler warnings.
- Added `.github/dependabot.yml`, issue forms, a PR template and `SECURITY.md`.
- `.gitignore` now covers MkDocs output, in-source CMake artifacts, Python
  caches, local `.env` files and crash dumps.

**Documentation**

- README, CHANGELOG and all docs pages now match the code: no SIMD, threading,
  fixed-point or GPU execution claims; no stale build scripts or options.
- Fixed broken MkDocs navigation and removed missing assets and the
  polyfill.io reference.
- Read the Docs builds on Ubuntu 24.04 with Python 3.12 (`mkdocs build --strict`
  passes).

## Verification

Verified on Windows with MinGW-w64 GCC 13.2.0, `-Wall -Wextra -O3`:

| Check | Result |
| --- | --- |
| Clean build | 0 warnings, 0 errors |
| Correctness suite (`test_all`) | 1361/1361 checks |
| `ctest` | 1/1 test passed |
| `benchmark_all` | 52 passed, 0 failed, 2 skipped in 17.1 s |
| Example programs | 9 executables, all exit 0 |
| `realtime_analyzer` | 5/5 launches, no crash |
| MkDocs strict build | passes in 1.36 s |
| Strict C99, GCC 15.3 container | 1361/1361 checks |

The two skipped benchmark rows are documented as such; they are counted, not
silently dropped. The Clang lane and the GitHub-hosted runners are validated by
CI, not on this machine (no local Clang or Linux host). GPU code remains
source-only and is not executed by any plan.

## Breaking changes and upgrade notes

- Transforms now return `fft_error_t`; ignoring the return value still compiles.
- `FFT_PATIENT`, `FFT_EXHAUSTIVE`, `FFT_UNALIGNED` and `FFT_THREADED` are gone.
  Their bit positions are deliberately not reused, so imported wisdom remains
  interpretable.
- MSVC is unsupported and now fails at configure time.
- License is Apache 2.0, not MIT.

## Related issues

- Fixes #3 (radix-4 was not radix-4; stale documentation contradicted the code)
- Fixes #4 (Debian compilation errors under strict C11/C99 feature-test macros)
- Does **not** implement #5 (Meson build support); that remains open.
- Supersedes PR #7 (base-4-only rewrite, one commit behind), whose implementation
  is replaced by the branch's radix-4 kernel; do not merge both.

## Checklist

- [x] Build is warning-clean with the documented toolchain
- [x] Full test suite and CTest pass locally
- [x] Benchmarks run with `Failed: 0`
- [x] All example programs build and exit 0
- [x] Documentation builds with `--strict`
- [x] README, CHANGELOG, docs and repository metadata updated
