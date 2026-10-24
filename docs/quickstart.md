# Quick Start

This is a short tour of the library for version 3.0.0. The full setup notes are
in [Getting Started](getting-started.md).

## Build

```bash
git clone https://github.com/muditbhargava66/FFT-implementation-in-C.git
cd FFT-implementation-in-C

cmake -S . -B build
cmake --build build
./build/bin/test_all       # ends with 1361/1361 checks passed
```

## A first FFT

```c
#include <stdio.h>
#include "fft_common.h"
#include "fft_algorithms.h"

int main(void) {
    int n = 16;
    complex_t *signal = allocate_complex_array(n);

    /* 5 cycles across the window: 5 Hz at a 16 Hz sample rate */
    for (int i = 0; i < n; i++)
        signal[i] = sin(2 * PI * 5 * i / 16.0);

    if (radix2_dit_fft(signal, n, FFT_FORWARD) != FFT_OK)
        return 1;

    printf("Frequency Bin | Magnitude\n");
    printf("-------------|----------\n");
    for (int i = 0; i < n / 2; i++)
        printf("     %2d      | %.3f\n", i, cabs(signal[i]));

    free_complex_array(signal);
    return 0;
}
```

Compile and run:

```bash
gcc -I./include test_fft.c -L./build/lib -lfft_lib_c -lm -o test_fft
./test_fft
```

The output is:

```
Frequency Bin | Magnitude
-------------|----------
      0      | 0.000
      1      | 0.000
      2      | 0.000
      3      | 0.000
      4      | 0.000
      5      | 8.000
      6      | 0.000
      7      | 0.000
```

The peak at bin 5 has magnitude 8 because a real sine of amplitude 1 is split
between the positive and negative frequency bins, and the transform is
unnormalised in the forward direction.

## Let the planner choose

```c
#include <fft_auto.h>

complex_t *data = fft_alloc_complex(1024);
/* fill data */
fft_auto(data, data, 1024, FFT_FORWARD);
```

`fft_auto()` plans, executes and destroys the plan in one call. If you transform
the same size repeatedly, keep the plan instead:

```c
fft_plan_t plan = fft_plan_dft_1d(1024, data, data, -1, FFT_MEASURE);
if (!plan) return 1;

printf("chosen: %s\n", fft_plan_algorithm_name(plan));
fft_execute(plan);
fft_destroy_plan(plan);
```

`FFT_MEASURE` times the candidates once and caches the result. Without it the
planner uses the size heuristic and does no timing.

## Real input

```c
double *audio = fft_alloc_real(4096);
complex_t *spectrum = fft_alloc_complex(2049);   /* n/2 + 1 */

fft_plan_t plan = fft_plan_r2c_1d(4096, audio, spectrum, FFT_ESTIMATE);
fft_execute(plan);
fft_destroy_plan(plan);
```

The inverse is `fft_plan_c2r_1d()`: give it `n/2+1` complex bins and it writes
`n` real samples.

## 2-D input

```c
int rows = 64, cols = 64;
complex_t *image = fft_alloc_complex(rows * cols);

fft_plan_t plan = fft_plan_dft_2d(rows, cols, image, image, -1, FFT_ESTIMATE);
fft_execute(plan);
fft_destroy_plan(plan);
```

Both dimensions must be powers of two.

## Choosing a kernel directly

The individual kernels are public. For power-of-two lengths:

```c
radix2_dit_fft(signal, 1024, FFT_FORWARD);   /* in-place, no scratch */
radix4_fft(signal, 1024, FFT_FORWARD);       /* faster from n=512 up */
split_radix_fft(signal, 1024, FFT_FORWARD);  /* fewest multiplications */
```

For arbitrary lengths:

```c
bluestein_fft(signal, 1000, FFT_FORWARD);    /* any size */
mixed_radix_fft(signal, 360, FFT_FORWARD);   /* composite sizes */
```

The measured trade-offs are in the [Performance Guide](performance.md).
Split-radix has the lowest operation count but is the slowest of the three
power-of-two kernels in wall time, so "fewest operations" is not the same as
"fastest".

## Run the demos

```bash
./build/bin/demo_features      # planning, selection and wisdom
./build/bin/convolution        # FFT convolution, 1-D and 2-D
./build/bin/audio_spectrum     # windowed spectrum analysis
./build/bin/pitch_detection    # peak and harmonic product spectrum
./build/bin/benchmark_all      # timings for every implementation
```

The per-algorithm teaching programs are not separate build targets. Their
`main()` functions are excluded from the library by `LIB_BUILD`, and the demo
files depend on each other, so they are not compiled as executables. Use
`test_all`, `benchmark_all`, and the demos under `applications/` instead.

## Troubleshooting

**"Size is not a power of two"** — use `bluestein_fft()`, `mixed_radix_fft()`,
or let the planner choose. Radix-2 and radix-4 require powers of two.

**Undefined references to math functions** — add `-lm` on Linux.

**Unexpected spectrum** — check whether you applied a window, and remember that
the forward transform is scaled by n while the inverse divides by n.
