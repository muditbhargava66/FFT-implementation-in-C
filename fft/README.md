# FFT Implementation in C 

This repository contains an implementation of the Cooley-Tukey radix-2 decimation-in-time (DIT) Fast Fourier Transform (FFT) algorithm in C. The code is targeted for the ARM Cortex-A9 processor on the Zynq-7000 All Programmable SoC.

## Table of Contents

- [Introduction](#introduction)
- [Algorithm](#algorithm)
- [Code Structure](#code-structure)
- [Build Instructions](#build-instructions)
- [Running the FFT](#running-the-fft)
- [Performance Analysis](#performance-analysis)
- [References](#references)

## Introduction

The Fast Fourier Transform (FFT) is a highly efficient algorithm for computing the Discrete Fourier Transform (DFT) of a sequence. It reduces the complexity from O(N^2) to O(NlogN). This implementation uses the Cooley-Tukey radix-2 decimation-in-time (DIT) algorithm which recursively divides the DFT into smaller DFTs of even and odd indexed terms.

## Algorithm

The radix-2 DIT FFT algorithm can be summarized in the following steps:

1. Divide the input sequence x[n] into two parts: f1[n] with even-indexed samples and f2[n] with odd-indexed samples.
2. Recursively compute the FFT of f1[n] and f2[n].
3. Combine the results to produce the FFT of x[n] using twiddle factors.

The mathematical formulation is explained in detail in the code comments.

## Code Structure

- `fft.c`: Contains the main FFT function `fft()` which implements the radix-2 DIT algorithm, along with utility functions for bit reversal and printing the output. The `main()` function demonstrates the usage with an example 8-point input.

## Build Instructions

1. Install Xilinx SDK on your machine.
2. Create a new C project in Xilinx SDK.
3. Add `fft.c` to your project source files. 
4. Select the target processor as "ARM Cortex-A9".
5. Build the project.

## Running the FFT

1. Connect the Zedboard to your machine.
2. Power on the Zedboard and ensure it is detected by Xilinx SDK.
3. Program the FPGA with the generated bitstream.
4. Run the compiled FFT application on the ARM processor.
5. View the output in the SDK console. It will display the input sequence, the computed FFT output, and the execution time.

## Performance Analysis

The execution time of the FFT is measured using the `clock()` function from the `time.h` library. The code is run on the ARM Cortex-A9 processor at 667 MHz. 

For an example 8-point input sequence, the measured execution time is printed in the output, demonstrating the real-time performance of the implementation.

## References

1. [Cooley, James W., and John W. Tukey. "An algorithm for the machine calculation of complex Fourier series." Mathematics of computation 19.90 (1965): 297-301.](https://community.ams.org/journals/mcom/1965-19-090/S0025-5718-1965-0178586-1/S0025-5718-1965-0178586-1.pdf)

2. Zynq-7000 SoC Technical Reference Manual, Xilinx Inc.

3. Xilinx SDK Documentation, Xilinx Inc.