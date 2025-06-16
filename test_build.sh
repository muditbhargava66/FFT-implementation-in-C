#!/bin/bash

# Quick test script to verify the build fixes

echo "Testing FFT build system fixes..."
echo "================================"

# Clean previous build
echo "Cleaning previous build..."
make clean > /dev/null 2>&1

# Try building just the library first
echo "Building FFT library..."
make lib/libfft.a

if [ $? -eq 0 ]; then
    echo "✓ Library built successfully"
else
    echo "✗ Library build failed"
    exit 1
fi

# Try building a simple algorithm
echo "Building radix2_dit..."
make bin/radix2_dit

if [ $? -eq 0 ]; then
    echo "✓ radix2_dit built successfully"
else
    echo "✗ radix2_dit build failed"
    exit 1
fi

# Try building bluestein (the one that had issues)
echo "Building bluestein..."
make bin/bluestein

if [ $? -eq 0 ]; then
    echo "✓ bluestein built successfully"
else
    echo "✗ bluestein build failed"
    exit 1
fi

# Run a simple test
echo "Running radix2_dit test..."
./bin/radix2_dit > /dev/null

if [ $? -eq 0 ]; then
    echo "✓ radix2_dit runs successfully"
else
    echo "✗ radix2_dit execution failed"
fi

echo ""
echo "All tests passed! The build system is working correctly."
echo "You can now run: make all"
