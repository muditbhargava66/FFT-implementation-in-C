# Security Policy

## Supported versions

| Version | Supported |
| ------- | --------- |
| 3.0.x   | Yes       |
| 2.x     | No        |
| 1.x     | No        |

The default build is CPU-only. The CUDA and Metal sources under `gpu/` are not
compiled into the default library, and no GPU backend has been run against
hardware. Reports about GPU code are still welcome if you state the
configuration used.

## Reporting a vulnerability

Please do not open a public issue for a suspected vulnerability. Use GitHub's
private vulnerability reporting:

1. Open the repository's **Security** tab.
2. Select **Report a vulnerability**.
3. Provide the details listed below.

If private reporting is unavailable, contact the maintainer through the GitHub
profile linked in the repository rather than disclosing the issue publicly.

## What to include

- Affected version or commit, platform, compiler, and build flags.
- A minimal reproducer: transform size, planning flags, buffer layout, and the
  input that triggers the problem.
- Expected and actual behaviour, and whether the impact is a crash, wrong
  result, or resource exhaustion.
- Sanitizer output if available (ASan, UBSan, or Valgrind).
- Whether the issue reproduces through the public planning API
  (`fft_plan_dft_1d`, `fft_plan_r2c_1d`, `fft_plan_c2r_1d`,
  `fft_plan_dft_2d`) or requires calling a kernel directly.

## Scope

In scope:

- Memory-safety defects: out-of-bounds reads or writes, use-after-free, double
  free, and uninitialised reads.
- Integer overflow or underflow that leads to an incorrect allocation or index.
- Malformed `fft_import_wisdom_from_string()` input that produces an unusable,
  unsafe, or silently non-executing plan.
- Supply-chain or build issues in the sources, CMake files, workflows, or
  dependency configuration shipped in this repository.

Out of scope:

- Numerical differences within floating-point round-off. The test suite's
  tolerance is documented in `tests/test_all.c`.
- Deliberate denial of service through enormous transform sizes or allocation
  requests; the API reports allocation failure rather than aborting.
- Vulnerabilities in third-party dependencies. Report those upstream; this
  repository updates them through Dependabot.
- The unbuilt GPU backends, unless the report also affects a supported build.

## Disclosure

The maintainer will acknowledge a report as soon as practical, confirm the
affected versions, and prepare a fix. Please allow time for a release before
public disclosure. Credit will be given in the changelog unless you prefer
otherwise.
