# GitHub Issue and PR Reply Drafts - v3.0.0

Prepared for the v3.0.0 release branch, 2026-10-23. Post these after the release
PR is merged so the links and the version on `main` are correct.

## Issue #3: radix4_fft is not radix4

Hi @Stenzel, you were right, and thank you for the report.

The v2 implementation did use radix-2 butterflies for every stage; the comment
in the file said so explicitly. v3.0.0 replaces it with a real radix-4 DIT
kernel:

- base-4 digit reversal, with a leading radix-2 stage when `log2(n)` is odd, so
  every power-of-two length is accepted (not only powers of four);
- `m *= 4` stages with the four-point butterfly and per-block twiddles;
- error codes instead of terminating the caller;
- an independent direct-DFT comparison in the test suite, including n = 4, 8,
  16, 32, 64, 256, 512, 1024, 4096 and 16384. The suite reports 1361/1361
  checks and passes under Linux GCC, Linux Clang, Windows MinGW-w64 and GCC 15
  in strict C99 mode.

The planner selects radix-4 above 256 points; below that radix-2 DIT is faster
because radix-4 pays for its O(n) digit-reversal buffer. Fixing this in v3.0.0,
so I am closing the issue.

## Issue #4: Compilation errors on Debian

Hi @x-logLT, thanks for the detailed report. All three points are addressed in
v3.0.0:

1. **`strdup` and `posix_memalign` implicit declarations.** The feature-test
   macros are now defined in `include/fft_common.h` before any system header is
   included, not only on the CMake command line. A strict `-std=c99` build with
   implicit function declarations promoted to errors now compiles without
   diagnostics. There are no `strdup` calls left.
2. **GCC 14+ hard errors.** The same change covers the stricter implicit
   declaration behaviour; CI builds warning-free under GCC and Clang.
3. **`<sys/sysctl.h>`.** There is no `sysctl` usage left. The old macOS
   core-count helper went away with the unreachable parallel code, so no
   replacement header or platform branch is required.

Verification: GCC 15.3 in a container, strict C99, `-Werror=implicit-function-declaration`,
reported 1361/1361 checks. Closing this with the v3.0.0 release.

## Issue #5: Adding meson build system support

Hi @antoniovazquezblanco, thanks again for the suggestion. v3.0.0 does not add a
Meson build; it ships the CMake build only. I am leaving this open because a
`meson.build` and a WrapDB package are still welcome, either from you or as
future work. The library is C11, has no threading dependency, and the GPU
sources are not wired into the planner, so a Meson port should not need to
reproduce optional CUDA or pthread options.

## PR #7: Fix Issue #3: implement true radix-4 FFT (DIT)

Hi @noro-ovo, thank you for taking the time to implement this and for linking it
to issue #3.

v3.0.0 now ships its own radix-4 rewrite, so this PR is superseded and I am
closing it without merging. Two concrete reasons:

- The release implementation supports every power-of-two length, including
  `2 * 4^k` such as n = 8, 32 and 128, via a leading radix-2 stage. The PR
  accepts powers of four only.
- The PR changes `radix4_fft()` back to `void` and reports failures on stderr.
  v3.0.0 has moved the whole API to `fft_error_t`, and merging the PR would
  regress that contract.

The n = 16 mismatch recorded in the PR description is not in the PR butterfly:
the base branch at 15c49c4 had a buggy `bit_reverse()` for `log2(n) <= 4`, so
radix-2 returned wrong spectra at n = 4, 8 and 16, and the PR's comparison demo
starts at 16. v3.0.0 fixes that permutation and compares radix-4 against an
independent direct DFT at n = 4, 8, 16 and every other tested size; the suite
passes 1361/1361.

If you would still like to contribute, issue #5 (Meson support) is open and
independent of this kernel work.
