# FFT Implementation Makefile v2.0.0
# Now with GPU support and automatic algorithm selection

CC = gcc
CXX = g++
NVCC = nvcc
CFLAGS = -Wall -Wextra -O3 -march=native -ffast-math -std=c99
CXXFLAGS = -Wall -Wextra -O3 -march=native -ffast-math -std=c++11
NVCCFLAGS = -O3 -arch=sm_70
LDFLAGS = -lm -lpthread

# Platform detection
UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# Platform-specific settings
ifeq ($(UNAME_S),Darwin)
    # macOS
    LDFLAGS += -framework Metal -framework MetalPerformanceShaders -framework Foundation
    GPU_BACKEND = METAL
    
    # Check for Homebrew GCC for OpenMP
    ifneq ($(shell which gcc-13),)
        CC = gcc-13
        OPENMP_FLAGS = -fopenmp
    else
        OPENMP_FLAGS = 
    endif
else ifeq ($(UNAME_S),Linux)
    # Linux
    OPENMP_FLAGS = -fopenmp
    
    # Check for CUDA
    ifneq ($(shell which nvcc),)
        GPU_BACKEND = CUDA
        LDFLAGS += -lcudart -lcufft
    endif
endif

# Directories
INC_DIR = include
ALGO_DIR = algorithms
APP_DIR = applications
OPT_DIR = optimizations
GPU_DIR = gpu
BENCH_DIR = benchmarks
TEST_DIR = tests
BIN_DIR = bin
OBJ_DIR = obj
LIB_DIR = lib
DOC_DIR = docs

# Create directories
$(shell mkdir -p $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR))
$(shell mkdir -p $(OBJ_DIR)/algorithms/core $(OBJ_DIR)/algorithms/dft $(OBJ_DIR)/algorithms/auto)
$(shell mkdir -p $(OBJ_DIR)/applications $(OBJ_DIR)/optimizations $(OBJ_DIR)/gpu)
$(shell mkdir -p $(ALGO_DIR)/auto)

# Include paths
INCLUDES = -I$(INC_DIR)

# Version info
VERSION = 2.0.0
CFLAGS += -DFFT_VERSION=\"$(VERSION)\"

# Source files
CORE_SRCS = $(wildcard $(ALGO_DIR)/core/*.c)
DFT_SRCS = $(wildcard $(ALGO_DIR)/dft/*.c)
AUTO_SRCS = $(wildcard $(ALGO_DIR)/auto/*.c)
APP_SRCS = $(wildcard $(APP_DIR)/*.c)
OPT_SRCS = $(wildcard $(OPT_DIR)/*.c)
BENCH_SRCS = $(wildcard $(BENCH_DIR)/*.c)
TEST_SRCS = $(wildcard $(TEST_DIR)/*.c)

# GPU sources
GPU_COMMON_SRCS = $(GPU_DIR)/fft_gpu.c
GPU_COMMON_OBJS = $(OBJ_DIR)/gpu/fft_gpu.o

ifeq ($(GPU_BACKEND),CUDA)
    GPU_SRCS = $(GPU_DIR)/fft_cuda.cu $(GPU_COMMON_SRCS)
    GPU_OBJS = $(OBJ_DIR)/gpu/fft_cuda.o $(GPU_COMMON_OBJS)
else ifeq ($(GPU_BACKEND),METAL)
    GPU_SRCS = $(GPU_DIR)/fft_metal.m $(GPU_COMMON_SRCS)
    GPU_OBJS = $(OBJ_DIR)/gpu/fft_metal.o $(GPU_COMMON_OBJS)
else
    GPU_SRCS = $(GPU_COMMON_SRCS)
    GPU_OBJS = $(GPU_COMMON_OBJS)
endif

# Object files for library
LIB_OBJS = $(patsubst $(ALGO_DIR)/core/%.c,$(OBJ_DIR)/algorithms/core/%_lib.o,$(CORE_SRCS))
LIB_OBJS += $(patsubst $(ALGO_DIR)/dft/%.c,$(OBJ_DIR)/algorithms/dft/%_lib.o,$(DFT_SRCS))
LIB_OBJS += $(patsubst $(ALGO_DIR)/auto/%.c,$(OBJ_DIR)/algorithms/auto/%_lib.o,$(AUTO_SRCS))
LIB_OBJS += $(GPU_OBJS)

# Static library
STATIC_LIB = $(LIB_DIR)/libfft.a

# Default target
all: info $(STATIC_LIB) algorithms applications optimizations benchmarks tests
	@echo "Build complete for v$(VERSION)!"

# Version info
info:
	@echo "========================================="
	@echo "FFT Implementation v$(VERSION)"
	@echo "Platform: $(UNAME_S) $(UNAME_M)"
	@echo "Compiler: $(CC)"
	@echo "GPU Backend: $(GPU_BACKEND)"
	@echo "========================================="

# Build static library
$(STATIC_LIB): $(LIB_OBJS)
	@echo "Creating library..."
	ar rcs $@ $^

# Compile library objects
$(OBJ_DIR)/algorithms/core/%_lib.o: $(ALGO_DIR)/core/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -DLIB_BUILD -c $< -o $@

$(OBJ_DIR)/algorithms/dft/%_lib.o: $(ALGO_DIR)/dft/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -DLIB_BUILD -c $< -o $@

$(OBJ_DIR)/algorithms/auto/%_lib.o: $(ALGO_DIR)/auto/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -DLIB_BUILD -c $< -o $@

# GPU common compilation
$(OBJ_DIR)/gpu/fft_gpu.o: $(GPU_DIR)/fft_gpu.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) -DLIB_BUILD -c $< -o $@

# GPU compilation
ifeq ($(GPU_BACKEND),CUDA)
$(OBJ_DIR)/gpu/fft_cuda.o: $(GPU_DIR)/fft_cuda.cu
	@mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(INCLUDES) -c $< -o $@
else ifeq ($(GPU_BACKEND),METAL)
$(OBJ_DIR)/gpu/fft_metal.o: $(GPU_DIR)/fft_metal.m
	@mkdir -p $(dir $@)
	$(CC) -ObjC $(CFLAGS) $(INCLUDES) -c $< -o $@
endif

# Targets
algorithms: $(STATIC_LIB) core dft auto

core: $(patsubst $(ALGO_DIR)/core/%.c,$(BIN_DIR)/%,$(CORE_SRCS))
dft: $(patsubst $(ALGO_DIR)/dft/%.c,$(BIN_DIR)/%,$(DFT_SRCS))
auto: $(patsubst $(ALGO_DIR)/auto/%.c,$(BIN_DIR)/%,$(AUTO_SRCS))
applications: $(STATIC_LIB) $(patsubst $(APP_DIR)/%.c,$(BIN_DIR)/%,$(APP_SRCS))
optimizations: $(STATIC_LIB) $(patsubst $(OPT_DIR)/%.c,$(BIN_DIR)/%,$(OPT_SRCS))
benchmarks: $(STATIC_LIB) $(BIN_DIR)/benchmark_all
tests: $(STATIC_LIB) $(BIN_DIR)/test_all

# Build executables
$(BIN_DIR)/%: $(ALGO_DIR)/core/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(ALGO_DIR)/dft/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(ALGO_DIR)/auto/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(APP_DIR)/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/%: $(OPT_DIR)/%.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# Special rules for parallel implementations
$(BIN_DIR)/parallel_fft: $(OPT_DIR)/parallel_fft.c $(STATIC_LIB)
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# Benchmark and test compilation
$(BIN_DIR)/benchmark_all: $(BENCH_DIR)/benchmark_all.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

$(BIN_DIR)/test_all: $(TEST_DIR)/test_all.c $(STATIC_LIB)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INCLUDES) $< -L$(LIB_DIR) -lfft -o $@ $(LDFLAGS)

# GPU demos
gpu-demo: $(STATIC_LIB)
	@echo "Building GPU demos..."
	$(CC) $(CFLAGS) $(INCLUDES) examples/gpu_demo.c -L$(LIB_DIR) -lfft -o $(BIN_DIR)/gpu_demo $(LDFLAGS)

# Documentation
docs:
	@echo "Building documentation..."
	cd $(DOC_DIR) && mkdocs build

docs-serve:
	cd $(DOC_DIR) && mkdocs serve

# Installation
install: all
	@echo "Installing FFT library v$(VERSION)..."
	mkdir -p /usr/local/include/fft
	cp $(INC_DIR)/*.h /usr/local/include/fft/
	cp $(STATIC_LIB) /usr/local/lib/
	@echo "Installation complete!"

uninstall:
	rm -rf /usr/local/include/fft
	rm -f /usr/local/lib/libfft.a

# Cleaning
clean:
	rm -rf $(BIN_DIR) $(OBJ_DIR) $(LIB_DIR)
	rm -f gmon.out *.dat *.png *.txt

distclean: clean
	rm -rf $(DOC_DIR)/_build
	find . -name "*.o" -delete
	find . -name "*.pyc" -delete
	find . -name "__pycache__" -delete

# Testing
test: tests
	@echo "Running tests..."
	./$(BIN_DIR)/test_all

benchmark: benchmarks
	@echo "Running benchmarks..."
	./$(BIN_DIR)/benchmark_all

# Development helpers
format:
	@echo "Formatting code..."
	find . -name "*.c" -o -name "*.h" | xargs clang-format -i

lint:
	@echo "Running static analysis..."
	cppcheck --enable=all --suppress=missingIncludeSystem -I$(INC_DIR) .

# Package release
release: clean all docs
	@echo "Creating release package v$(VERSION)..."
	mkdir -p fft-v$(VERSION)
	cp -r $(INC_DIR) $(ALGO_DIR) $(APP_DIR) $(OPT_DIR) $(GPU_DIR) fft-v$(VERSION)/
	cp -r $(DOC_DIR) $(TEST_DIR) $(BENCH_DIR) fft-v$(VERSION)/
	cp README.md LICENSE Makefile CHANGELOG.md fft-v$(VERSION)/
	tar -czf fft-v$(VERSION).tar.gz fft-v$(VERSION)
	rm -rf fft-v$(VERSION)
	@echo "Release package created: fft-v$(VERSION).tar.gz"

# Help
help:
	@echo "FFT Implementation Build System v$(VERSION)"
	@echo "========================================="
	@echo "Targets:"
	@echo "  all          - Build everything"
	@echo "  algorithms   - Build core FFT algorithms"
	@echo "  applications - Build FFT applications"
	@echo "  optimizations- Build optimized implementations"
	@echo "  benchmarks   - Build benchmark suite"
	@echo "  tests        - Build test suite"
	@echo "  gpu-demo     - Build GPU demonstration"
	@echo "  docs         - Generate documentation"
	@echo "  install      - Install library system-wide"
	@echo "  test         - Run all tests"
	@echo "  benchmark    - Run all benchmarks"
	@echo "  clean        - Remove build artifacts"
	@echo "  release      - Create release package"
	@echo ""
	@echo "Platform: $(UNAME_S) with $(GPU_BACKEND) GPU support"

.PHONY: all info algorithms core dft auto applications optimizations benchmarks \
        tests gpu-demo docs docs-serve install uninstall clean distclean test \
        benchmark format lint release help
