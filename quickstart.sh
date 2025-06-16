#!/bin/bash

# Quick Start Script for FFT Study Repository
# This script helps users quickly build and test the FFT implementations

echo "FFT Study Repository - Quick Start"
echo "=================================="
echo

# Check for required tools
echo "Checking prerequisites..."
command -v gcc >/dev/null 2>&1 || { echo "Error: gcc is required but not installed."; exit 1; }
command -v make >/dev/null 2>&1 || { echo "Error: make is required but not installed."; exit 1; }
echo "✓ Prerequisites satisfied"
echo

# Build options
echo "Select build option:"
echo "1. Build everything (recommended for first time)"
echo "2. Build core algorithms only"
echo "3. Build applications only"
echo "4. Build with optimizations"
echo "5. Run demos"
echo "6. Run benchmarks"
echo "7. Run tests"
echo "8. Clean build"
echo
read -p "Enter choice (1-8): " choice

case $choice in
    1)
        echo "Building all components..."
        make clean
        make all
        echo
        echo "Build complete! All binaries are in the bin/ directory."
        echo
        echo "Try these examples:"
        echo "  ./bin/radix2_dit     - Basic FFT algorithm"
        echo "  ./bin/audio_spectrum - Audio spectrum analyzer"
        echo "  ./bin/convolution    - FFT-based convolution"
        ;;
    
    2)
        echo "Building core algorithms..."
        make algorithms
        echo
        echo "Core algorithms built! Try:"
        echo "  ./bin/radix2_dit"
        echo "  ./bin/split_radix"
        echo "  ./bin/bluestein"
        ;;
    
    3)
        echo "Building applications..."
        make applications
        echo
        echo "Applications built! Try:"
        echo "  ./bin/audio_spectrum"
        echo "  ./bin/fft_filtering"
        echo "  ./bin/power_spectrum"
        ;;
    
    4)
        echo "Building optimized implementations..."
        make optimizations
        echo
        echo "Optimized versions built! Try:"
        echo "  ./bin/simd_fft"
        echo "  ./bin/parallel_fft"
        echo "  ./bin/fixed_point_fft"
        ;;
    
    5)
        echo "Running demonstrations..."
        make demo
        ;;
    
    6)
        echo "Building and running benchmarks..."
        make benchmarks
        ./bin/benchmark_all
        ;;
    
    7)
        echo "Building and running tests..."
        make tests
        ./bin/test_all
        ;;
    
    8)
        echo "Cleaning build artifacts..."
        make clean
        echo "Clean complete!"
        ;;
    
    *)
        echo "Invalid choice. Please run the script again."
        exit 1
        ;;
esac

echo
echo "For more information:"
echo "  - Read the main README.md"
echo "  - Check algorithm-specific READMEs in subdirectories"
echo "  - Run 'make help' for all build options"
echo
echo "Happy FFT studying! 🎵"
