#!/bin/bash

# Cleanup script for FFT v2.0.0 release
# Removes unnecessary files and prepares for release

echo "FFT v2.0.0 Release Preparation"
echo "=============================="

# Files and directories to remove
REMOVE_FILES=(
    # Temporary files
    "*.tmp"
    "*.temp"
    "*.bak"
    "*.backup"
    "*.swp"
    "*.swo"
    "*~"
    ".DS_Store"
    
    # Build artifacts that shouldn't be in repo
    "*.o"
    "*.a"
    "*.so"
    "*.dylib"
    "*.exe"
    "*.out"
    
    # Output files from examples
    "*.dat"
    "*.txt"
    "*.png"
    "*.wav"
    "spectrum.png"
    "fft_data.txt"
    "plot_spectrum.gnu"
    
    # Profiling data
    "gmon.out"
    "callgrind.out.*"
    "cachegrind.out.*"
    "perf.data"
    "perf.data.old"
    
    # Coverage files
    "*.gcda"
    "*.gcno"
    "*.gcov"
)

REMOVE_DIRS=(
    # Build directories
    "bin/"
    "obj/"
    "lib/"
    "build/"
    "dist/"
    
    # IDE directories  
    ".idea/"
    "*.dSYM/"
    
    # Documentation build
    "docs/_build/"
    "docs/.doctrees/"
    
    # Python cache
    "__pycache__/"
    "*.egg-info/"
    ".pytest_cache/"
    
    # Test outputs
    "test_output/"
    "benchmark_results/"
)

# Function to remove files
remove_files() {
    echo "Removing temporary and build files..."
    for pattern in "${REMOVE_FILES[@]}"; do
        find . -name "$pattern" -type f -exec rm -f {} \; 2>/dev/null
    done
}

# Function to remove directories
remove_directories() {
    echo "Removing build directories..."
    for dir in "${REMOVE_DIRS[@]}"; do
        find . -name "$dir" -type d -exec rm -rf {} \; 2>/dev/null
    done
}

# Function to fix line endings
fix_line_endings() {
    echo "Fixing line endings..."
    find . -name "*.c" -o -name "*.h" -o -name "*.md" | xargs dos2unix 2>/dev/null || true
}

# Function to remove empty directories
remove_empty_dirs() {
    echo "Removing empty directories..."
    find . -type d -empty -delete 2>/dev/null
}

# Function to update version strings
update_version() {
    echo "Updating version strings to v2.0.0..."
    
    # Update version in source files
    find . -name "*.c" -o -name "*.h" | xargs sed -i.bak 's/v1\.[0-9]\.[0-9]/v2.0.0/g'
    find . -name "*.c" -o -name "*.h" | xargs sed -i.bak 's/Version 1\.[0-9]\.[0-9]/Version 2.0.0/g'
    
    # Remove backup files
    find . -name "*.bak" -type f -exec rm -f {} \;
}

# Function to check for large files
check_large_files() {
    echo "Checking for large files..."
    find . -type f -size +1M -exec ls -lh {} \; | grep -v ".git"
}

# Main cleanup process
main() {
    # Confirm before proceeding
    echo "This will clean up the repository for v2.0.0 release."
    echo "Make sure you have committed all important changes!"
    read -p "Continue? (y/n) " -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cleanup cancelled."
        exit 1
    fi
    
    # Start cleanup
    echo
    echo "Starting cleanup..."
    
    remove_files
    remove_directories
    fix_line_endings
    remove_empty_dirs
    update_version
    
    echo
    echo "Checking for issues..."
    check_large_files
    
    echo
    echo "Repository cleanup complete!"
    echo
    echo "Next steps:"
    echo "1. Review changes: git status"
    echo "2. Add changes: git add -A"
    echo "3. Commit: git commit -m 'chore: Cleanup for v2.0.0 release'"
    echo "4. Tag release: git tag -a v2.0.0"
}

# Run main function
main
