Certainly! Here's a detailed README for the FFT implementation using OpenMP:

# Fast Fourier Transform (FFT) Implementation using OpenMP

This project implements the Fast Fourier Transform (FFT) algorithm using the Cooley-Tukey radix-2 decimation-in-time (DIT) method. The implementation is optimized for performance using OpenMP parallelization.

## Table of Contents

1. [Introduction](#introduction)
2. [Features](#features)
3. [Requirements](#requirements)
4. [Usage](#usage)
5. [Algorithm](#algorithm)
6. [Performance Analysis](#performance-analysis)
7. [Acknowledgements](#acknowledgements)
8. [References](#references)

## Introduction

The Fast Fourier Transform (FFT) is a widely used algorithm for efficiently computing the Discrete Fourier Transform (DFT) of a sequence. It has numerous applications in signal processing, image processing, and other domains. This project provides an optimized implementation of the FFT using the Cooley-Tukey radix-2 DIT method and leverages OpenMP for parallelization to improve performance.

## Features

- Implements the Cooley-Tukey radix-2 DIT FFT algorithm
- Supports both forward and inverse FFT computations
- Utilizes OpenMP for parallelization to accelerate the computation
- Provides a simple and intuitive interface for easy integration
- Includes performance analysis to measure execution time

## Requirements

- C compiler with OpenMP support (e.g., GCC)
- OpenMP library
- C99 standard or later

## Usage

1. Clone the repository:
   ```
   git clone https://github.com/your-username/fft-openmp.git
   ```

2. Navigate to the project directory:
   ```
   cd ./FFT-implementation-in-C/fft-openmp
   ```

3. Compile the code using a C compiler with OpenMP support:
   ```
   gcc-13 -fopenmp fft_openmp.c -o fft_openmp -lm
   ```

4. Run the executable:
   ```
   ./fft_openmp
   ```

5. The program will output the input sequence, the computed FFT, the inverse FFT, and the execution time.

## Algorithm

The FFT algorithm implemented in this project follows the Cooley-Tukey radix-2 DIT method. The main steps of the algorithm are as follows:

1. Perform bit-reversal permutation on the input sequence.
2. Recursively divide the input sequence into two halves until reaching the base case of size 2.
3. Compute the FFT of the even-indexed and odd-indexed subsequences separately.
4. Combine the results using the twiddle factors to obtain the FFT of the original sequence.
5. Repeat steps 3-4 until the entire FFT is computed.

The inverse FFT is computed by performing the forward FFT with the conjugate of the twiddle factors and scaling the result by the size of the input sequence.

## Performance Analysis

The project includes performance analysis to measure the execution time of the FFT computation. The execution time is measured using the `clock()` function from the `time.h` library. The time taken for both the forward FFT and the inverse FFT is displayed in seconds.

OpenMP parallelization is used to accelerate the computation by distributing the workload across multiple threads. The bit-reversal permutation, the outer loop of the Danielson-Lanczos algorithm, and the scaling loop for the inverse FFT are parallelized using OpenMP directives.

## Acknowledgements

This project was inspired by the need for efficient FFT implementations in various domains. We would like to acknowledge the contributions of the open-source community and the developers of the OpenMP library for their valuable work.

## References

[1] [Cooley, James W., and John W. Tukey. "An algorithm for the machine calculation of complex Fourier series." Mathematics of computation 19.90 (1965): 297-301.](https://community.ams.org/journals/mcom/1965-19-090/S0025-5718-1965-0178586-1/S0025-5718-1965-0178586-1.pdf)

[2] OpenMP Architecture Review Board. (2018). OpenMP Application Programming Interface Version 5.0. Retrieved from https://www.openmp.org/wp-content/uploads/OpenMP-API-Specification-5.0.pdf

[3] [Press, W. H., Teukolsky, S. A., Vetterling, W. T., & Flannery, B. P. (2007). Numerical recipes: The art of scientific computing (3rd ed.). Cambridge University Press.](https://books.google.com/books?hl=en&lr=&id=1aAOdzK3FegC&oi=fnd&pg=PA33&dq=Numerical+recipes:+The+art+of+scientific+computing&ots=3mMpJcEooc&sig=MnOVAz381TobdYa6Krsj0Bp4-_E)