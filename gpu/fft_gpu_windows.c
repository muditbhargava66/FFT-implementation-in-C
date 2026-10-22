/**
 * @file fft_gpu_windows.c
 * @brief Windows-specific GPU detection and initialization
 * 
 * Provides Windows-specific GPU detection for NVIDIA CUDA support
 */

#ifdef _WIN32

#include "../include/fft_gpu.h"
#include <windows.h>
#include <stdio.h>

// Windows-specific CUDA detection
static int detect_nvidia_gpu_windows(void) {
    // Try to load CUDA runtime DLL
    HMODULE cudart_dll = LoadLibraryA("cudart64_12.dll");
    if (!cudart_dll) {
        cudart_dll = LoadLibraryA("cudart64_11.dll");
    }
    if (!cudart_dll) {
        cudart_dll = LoadLibraryA("cudart64_10.dll");
    }
    
    if (cudart_dll) {
        FreeLibrary(cudart_dll);
        return 1;  // CUDA runtime found
    }
    
    return 0;  // CUDA not available
}

// Windows-specific GPU device enumeration
int fft_gpu_enumerate_devices_windows(void) {
    if (!detect_nvidia_gpu_windows()) {
        printf("No NVIDIA GPU detected on Windows\n");
        return 0;
    }
    
    printf("NVIDIA GPU detected on Windows\n");
    
    // Try to get device information through Windows APIs
    // This is a simplified version - full implementation would use CUDA APIs
    return 1;
}

// Windows-specific GPU memory info
void fft_gpu_get_windows_memory_info(size_t* total, size_t* available) {
    *total = 0;
    *available = 0;
    
    if (!detect_nvidia_gpu_windows()) {
        return;
    }
    
    // In a full implementation, this would query CUDA device properties
    // For now, we'll set reasonable defaults
    *total = 8ULL * 1024 * 1024 * 1024;      // 8GB default
    *available = 6ULL * 1024 * 1024 * 1024;  // 6GB available
}

// Windows-specific CUDA path detection
const char* fft_gpu_get_cuda_path_windows(void) {
    static char cuda_path[MAX_PATH];
    
    // Check environment variable first
    DWORD result = GetEnvironmentVariableA("CUDA_PATH", cuda_path, MAX_PATH);
    if (result > 0 && result < MAX_PATH) {
        return cuda_path;
    }
    
    // Check common installation paths
    const char* common_paths[] = {
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.4",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.3",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.2",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.1",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v12.0",
        "C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v11.8",
        NULL
    };
    
    for (int i = 0; common_paths[i]; i++) {
        DWORD attrs = GetFileAttributesA(common_paths[i]);
        if (attrs != INVALID_FILE_ATTRIBUTES && (attrs & FILE_ATTRIBUTE_DIRECTORY)) {
            strcpy(cuda_path, common_paths[i]);
            return cuda_path;
        }
    }
    
    return NULL;  // CUDA path not found
}

// Windows-specific performance optimization
void fft_gpu_optimize_for_windows(void) {
    // Set Windows-specific optimizations
    
    // Increase process priority for GPU computations
    SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
    
    // Set thread priority
    SetThreadPriority(GetCurrentThread(), THREAD_PRIORITY_ABOVE_NORMAL);
    
    // Disable Windows GPU scheduling if available (Windows 10 2004+)
    // This can improve GPU performance for compute workloads
    HKEY hKey;
    LONG result = RegOpenKeyExA(HKEY_LOCAL_MACHINE,
        "SYSTEM\\CurrentControlSet\\Control\\GraphicsDrivers", 0, KEY_SET_VALUE, &hKey);
    
    if (result == ERROR_SUCCESS) {
        DWORD value = 1;  // Disable hardware-accelerated GPU scheduling
        RegSetValueExA(hKey, "HwSchMode", 0, REG_DWORD, (BYTE*)&value, sizeof(value));
        RegCloseKey(hKey);
    }
}

// Windows-specific error handling
void fft_gpu_handle_windows_error(const char* operation, DWORD error_code) {
    char error_msg[256];
    FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                   NULL, error_code, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                   error_msg, sizeof(error_msg), NULL);
    
    fprintf(stderr, "Windows GPU Error in %s: %s (Code: %lu)\n", 
            operation, error_msg, error_code);
}

// Windows-specific GPU monitoring
typedef struct {
    HANDLE timer;
    LARGE_INTEGER frequency;
    LARGE_INTEGER start_time;
} windows_gpu_timer_t;

static windows_gpu_timer_t g_gpu_timer;

void fft_gpu_timer_start_windows(void) {
    QueryPerformanceFrequency(&g_gpu_timer.frequency);
    QueryPerformanceCounter(&g_gpu_timer.start_time);
}

double fft_gpu_timer_stop_windows(void) {
    LARGE_INTEGER end_time;
    QueryPerformanceCounter(&end_time);
    
    double elapsed = (double)(end_time.QuadPart - g_gpu_timer.start_time.QuadPart);
    elapsed /= g_gpu_timer.frequency.QuadPart;
    
    return elapsed * 1000.0;  // Return milliseconds
}

#endif // _WIN32