# Applications

The programs in `applications/` are runnable demonstrations. Each one is a
single C file with its own helpers; those helpers are not part of the library
API and are not declared in `include/`. They show how the kernels can be used,
not a stable interface to build against.

Build them all with CMake and run them from `build/bin/`.

## audio_spectrum

`./build/bin/audio_spectrum`

Generates a synthetic signal with three tones (440, 554 and 659 Hz), then
analyses it with rectangular, Hann, Hamming and Blackman windows. It prints an
ASCII spectrum, the strongest peaks, and the nearest musical notes for a C major
chord. The window comparison is the useful part: it shows the leakage and
resolution differences without needing an audio file.

## realtime_analyzer

`./build/bin/realtime_analyzer`

Simulates a streaming analyzer: a 2048-sample window, 512-sample hop, Hann
window and four-frame averaging. It waits for Enter, then processes ten seconds
of generated audio and prints peak tracks. The program uses a sliding buffer and
`fft_plan_dft_1d()` rather than transforming each hop with a fresh one-shot call.

## fft_filtering

`./build/bin/fft_filtering`

Builds low-pass, high-pass and band-pass filters in the frequency domain, plots
their responses in ASCII, and filters a multi-tone signal. It also designs a
windowed-sinc FIR filter with `design_fir_filter()` and compares the
frequency-domain and time-domain approaches. All three filter types are
implemented in this file; the library itself only supplies the FFT and window
functions.

## convolution

`./build/bin/convolution`

Shows linear convolution (zero-padded FFT), circular convolution, 2-D
convolution with zero padding, and a direct-versus-FFT timing comparison. The
2-D routine is `fft_convolution_2d()`, implemented in this file:
zero-pad both operands, forward transform, multiply pointwise, inverse
transform, and copy the valid region.

For a linear convolution of lengths `nx` and `nh`, the FFT length is
`next_power_of_two(nx + nh - 1)` so the circular wrap does not corrupt the
result.

## power_spectrum

`./build/bin/power_spectrum`

Compares spectral estimators on a two-tone signal with added noise:

- periodogram (`compute_periodogram`)
- Welch's method with 50% overlap and a Hann window (`welch_psd`)
- autocorrelation via FFT
- cross-correlation via FFT
- magnitude-squared coherence for two signals (`compute_coherence`)

It also computes spectral centroid, spread, skewness, kurtosis, rolloff and
flatness from a PSD. These are demo-local functions, not library exports.

## pitch_detection

`./build/bin/pitch_detection`

Runs three pitch estimators on a pure tone, a harmonic-rich note, and a C major
chord:

- peak picking with quadratic interpolation
- harmonic product spectrum
- autocorrelation

`detect_pitch_with_confidence()` combines them and reports the note name and
cents offset. The chord case is intentionally hard; its lower confidence is part
of the demonstration.

## image_fft

`./build/bin/image_fft`

A 2-D tour using `complex_t **` arrays: forward and inverse transforms, FFT
shift, ideal and Gaussian low-pass filters, edge detection, and a small
reconstruction-error check. The 2-D helpers here (`fft_2d`, `fft_shift_2d`)
predate the planner's `fft_plan_dft_2d()` and operate on row pointers; both
produce the same row-column decomposition.

## demo_features

`./build/bin/demo_features`

Exercises the planning API: hardware capability bits, algorithm selection at
several sizes, `fft_plan_algorithm_name()`, the wisdom cache, r2c/c2r and the
2-D plan. It prints which algorithm the planner actually chose rather than a
copy of the selection rules.

## GPU benchmark

`./build/bin/gpu_benchmark`

Reports CPU timings and checks `fft_gpu_available()`. In a default build no GPU
backend is compiled in, so it prints that the GPU is not available. It is built
so the GPU API cannot drift unnoticed, not because a GPU path is usable.

## Reusing the code

If you want one of these routines in your own program, copy the function and
its declarations into your source, or use the library API for the transforms
themselves. The application files are compiled as standalone executables; they
are not linked into `libfft_lib_c.a`.

Two details worth keeping if you adapt the code:

- The forward transform is unnormalised, so a forward/inverse pair returns the
  original signal only after the inverse divides by `n` (which it does).
- Windowed analysis needs the window applied before the transform, and the
  magnitude scaling should account for the window's coherent gain if absolute
  amplitudes matter.
