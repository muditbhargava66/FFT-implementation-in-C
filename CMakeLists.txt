cmake_minimum_required(VERSION 3.10)
project(FFT_Implementation_C VERSION 3.1.0 LANGUAGES C)

# Optional: Enable CUDA if available
include(CheckLanguage)
check_language(CUDA)
if(CMAKE_CUDA_COMPILER)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 14)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
    # Set CUDA architectures (modify based on your GPU)
    # GTX 1050 detected: compute capability 6.1 (Pascal)
    # Format: 61 = sm_61, 70 = sm_70, etc.
    set(CMAKE_CUDA_ARCHITECTURES "61;70;75;80" CACHE STRING "CUDA architectures")
    
    # CUDA compiler flags for MSVC compatibility
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda --expt-relaxed-constexpr")
    if(MSVC)
        set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=/wd4819,/wd4996")
    endif()
    
    message(STATUS "CUDA support enabled with architectures: ${CMAKE_CUDA_ARCHITECTURES}")
else()
    message(STATUS "CUDA not found - GPU acceleration disabled")
endif()

# Set C standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# The kernels use C11 _Complex arithmetic, which the MSVC C compiler does not
# provide. Fail at configure time with the reason instead of partway through
# the build. Gate on the compiler id because clang-cl sets MSVC and does
# support _Complex.
if(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    message(FATAL_ERROR
        "MSVC's C compiler cannot build this library: it requires C11 _Complex "
        "arithmetic, which cl.exe does not provide. Configure with MinGW-w64 "
        "(for example `cmake -G \"MinGW Makefiles\" ..`) or with clang.")
endif()

# Output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directories
include_directories(include)

# -march=native is opt-in: it produces binaries that only run on the build CPU,
# and MinGW-w64 GCC 13.2.0 miscompiles 32-byte AVX stack stores at some
# alignment phases, which crashed realtime_analyzer intermittently.
option(FFT_NATIVE_ARCH "Tune for the build host CPU (-march=native); non-portable" OFF)

# Compiler options
add_compile_options(-Wall -Wextra -O3)
if(FFT_NATIVE_ARCH)
    message(WARNING
        "FFT_NATIVE_ARCH is ON: binaries will only run on CPUs compatible "
        "with the build host, and MinGW-w64 GCC may emit misaligned 32-byte "
        "AVX stack stores that crash intermittently.")
    add_compile_options(-march=native)
endif()
# Expose POSIX functions (posix_memalign, sysconf) in strict C11 mode
add_compile_definitions(_GNU_SOURCE _POSIX_C_SOURCE=200809L)

# CONFIGURE_DEPENDS makes the glob re-run when files are added or removed.
file(GLOB_RECURSE FFT_C_SOURCES CONFIGURE_DEPENDS
    "algorithms/core/*.c"
    "algorithms/dft/*.c"
    "algorithms/auto/*.c"
    "utils/*.c"
)

# The GPU sources are outside the default build and are not wired to the
# planner. Enabling this compiles the CUDA objects; fft_gpu_* still resolves to
# utils/fft_gpu_stubs.c, so execution remains on the CPU.
option(FFT_ENABLE_GPU "Build the CUDA GPU objects (unmaintained; dispatch not wired up)" OFF)

if(FFT_ENABLE_GPU AND NOT CMAKE_CUDA_COMPILER)
    message(WARNING "FFT_ENABLE_GPU is ON but no CUDA compiler was found; building CPU-only.")
endif()
if(FFT_ENABLE_GPU AND CMAKE_CUDA_COMPILER)
    message(STATUS "GPU objects will be compiled, but fft_gpu_* still resolves to the CPU-only stubs (see CMakeLists.txt)")
else()
    message(STATUS "GPU backend disabled: CPU-only build")
endif()

add_library(fft_lib_c STATIC ${FFT_C_SOURCES})
target_compile_definitions(fft_lib_c PRIVATE LIB_BUILD)

# libm is separate on Linux. PUBLIC so tests, benchmarks and applications
# inherit it; MinGW and macOS resolve these symbols from the C runtime.
if(UNIX AND NOT APPLE)
    target_link_libraries(fft_lib_c PUBLIC m)
endif()

if(FFT_ENABLE_GPU AND CMAKE_CUDA_COMPILER)
    file(GLOB GPU_CUDA_SOURCES "gpu/*.cu")

    # Create CUDA library (only .cu files)
    add_library(fft_lib_cuda STATIC ${GPU_CUDA_SOURCES})
    target_compile_definitions(fft_lib_cuda PRIVATE LIB_BUILD __CUDACC__)
    target_link_libraries(fft_lib_cuda PUBLIC cudart cufft)
    set_target_properties(fft_lib_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_RESOLVE_DEVICE_SYMBOLS ON
    )

    # Combine into single library
    add_library(fft_lib INTERFACE)
    target_link_libraries(fft_lib INTERFACE fft_lib_c fft_lib_cuda)
else()
    # No CUDA - just use C library
    add_library(fft_lib ALIAS fft_lib_c)
endif()

# Tests
enable_testing()

# Helper to add test executables
function(add_fft_test name source)
    add_executable(${name} ${source})
    target_link_libraries(${name} PRIVATE fft_lib)
    add_test(NAME ${name} COMMAND ${name})
endfunction()

# Add test_all
add_executable(test_all tests/test_all.c)
target_link_libraries(test_all PRIVATE fft_lib)
add_test(NAME test_all COMMAND test_all)

# Add benchmarks
add_executable(benchmark_all benchmarks/benchmark_all.c)
target_link_libraries(benchmark_all PRIVATE fft_lib)

# Built as a target so it cannot silently drift from the GPU API again.
add_executable(gpu_benchmark benchmarks/gpu_benchmark.c)
target_link_libraries(gpu_benchmark PRIVATE fft_lib)

# Add applications
file(GLOB APP_SOURCES "applications/*.c")
foreach(APP_SRC ${APP_SOURCES})
    get_filename_component(APP_NAME ${APP_SRC} NAME_WE)
    add_executable(${APP_NAME} ${APP_SRC})
    target_link_libraries(${APP_NAME} PRIVATE fft_lib)
endforeach()

# CLI tools (individual algorithm demos) are disabled in the CMake build.
# Each algorithm .c file defines both the algorithm functions and a main()
# guarded by #ifndef LIB_BUILD. Building standalone CLIs would cause duplicate
# symbols against the library, and the demos depend on each other. Use test_all
# or benchmark_all instead.
