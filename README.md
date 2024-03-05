# FFT Implementation in C

This repository contains two implementations of the Fast Fourier Transform (FFT) algorithm in C: a sequential version and a parallel version using OpenMP.

## Folder Structure

- `sequential-fft`: Contains the sequential FFT implementation targeting the ARM Cortex-A9 processor on the Zynq-7000 All Programmable SoC.
- `parallel-fft`: Contains the parallel FFT implementation optimized using OpenMP for improved performance.

## Sequential FFT

- Implements the Cooley-Tukey radix-2 decimation-in-time (DIT) FFT algorithm.
- Targeted for the ARM Cortex-A9 processor on the Zynq-7000 All Programmable SoC.
- Build and run instructions are provided in the `sequential-fft` folder.

## Parallel FFT using OpenMP

- Implements the Cooley-Tukey radix-2 DIT FFT algorithm with OpenMP parallelization.
- Supports both forward and inverse FFT computations.
- Requires a C compiler with OpenMP support (e.g., GCC) and C99 standard or later.
- Usage instructions are provided in the `parallel-fft` folder.