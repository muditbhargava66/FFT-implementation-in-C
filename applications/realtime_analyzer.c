#include "fft_common.h"
#include "fft_algorithms.h"
#include <unistd.h>

/**
 * Real-time FFT Spectrum Analyzer
 * 
 * Simulates a real-time spectrum analyzer that processes
 * streaming audio data using FFT with overlapping windows.
 * 
 * Features:
 * - Sliding window analysis
 * - Real-time peak detection
 * - Frequency tracking
 * - Moving average smoothing
 */

// Analyzer configuration
typedef struct {
    int fft_size;
    int hop_size;
    double sample_rate;
    const char* window_type;
    int averaging_frames;
} analyzer_config_t;

// Analyzer state
typedef struct {
    analyzer_config_t config;
    complex_t* input_buffer;
    complex_t* window_buffer;
    double* magnitude_avg;
    int buffer_pos;
    int frame_count;
} analyzer_t;

// Create analyzer
analyzer_t* create_analyzer(analyzer_config_t config) {
    /* radix2_dit_fft needs a power-of-two length, the hop cannot exceed the
     * window, and the moving average divides by averaging_frames. */
    if (config.fft_size <= 0 || !is_power_of_two(config.fft_size) ||
        config.hop_size <= 0 || config.hop_size > config.fft_size ||
        config.averaging_frames < 1 || !config.window_type)
        return NULL;

    analyzer_t* analyzer = (analyzer_t*)malloc(sizeof(analyzer_t));
    if (!analyzer)
        return NULL;

    analyzer->config = config;
    analyzer->input_buffer = allocate_complex_array(config.fft_size);
    analyzer->window_buffer = allocate_complex_array(config.fft_size);
    analyzer->magnitude_avg = (double*)calloc(config.fft_size/2 + 1, sizeof(double));
    analyzer->buffer_pos = 0;
    analyzer->frame_count = 0;

    if (!analyzer->input_buffer || !analyzer->window_buffer ||
        !analyzer->magnitude_avg) {
        free_complex_array(analyzer->input_buffer);
        free_complex_array(analyzer->window_buffer);
        free(analyzer->magnitude_avg);
        free(analyzer);
        return NULL;
    }

    return analyzer;
}

void destroy_analyzer(analyzer_t* analyzer) {
    if (!analyzer)
        return;
    free_complex_array(analyzer->input_buffer);
    free_complex_array(analyzer->window_buffer);
    free(analyzer->magnitude_avg);
    free(analyzer);
}

// Process audio frame
void process_frame(analyzer_t* analyzer, double* audio_samples, int num_samples) {
    if (!analyzer || !audio_samples || num_samples <= 0)
        return;

    // Add samples to circular buffer
    for (int i = 0; i < num_samples; i++) {
        analyzer->input_buffer[analyzer->buffer_pos] = audio_samples[i];
        analyzer->buffer_pos = (analyzer->buffer_pos + 1) % analyzer->config.fft_size;
    }
    
    // Check if we have enough samples for FFT
    if (analyzer->buffer_pos % analyzer->config.hop_size == 0) {
        // Copy and reorder buffer for FFT
        int start = analyzer->buffer_pos;
        for (int i = 0; i < analyzer->config.fft_size; i++) {
            int idx = (start + i) % analyzer->config.fft_size;
            analyzer->window_buffer[i] = analyzer->input_buffer[idx];
        }
        
        // Apply window
        if (strcmp(analyzer->config.window_type, "hann") == 0) {
            apply_window_hann(analyzer->window_buffer, analyzer->config.fft_size);
        } else if (strcmp(analyzer->config.window_type, "hamming") == 0) {
            apply_window_hamming(analyzer->window_buffer, analyzer->config.fft_size);
        }
        
        // Compute FFT
        radix2_dit_fft(analyzer->window_buffer, analyzer->config.fft_size, FFT_FORWARD);
        
        // Update magnitude average
        double alpha = 1.0 / analyzer->config.averaging_frames;
        for (int i = 0; i <= analyzer->config.fft_size/2; i++) {
            double mag = cabs(analyzer->window_buffer[i]);
            analyzer->magnitude_avg[i] = (1 - alpha) * analyzer->magnitude_avg[i] + alpha * mag;
        }
        
        analyzer->frame_count++;
    }
}

// Get current spectrum
void get_spectrum(analyzer_t* analyzer, double* magnitude, double* frequencies) {
    for (int i = 0; i <= analyzer->config.fft_size/2; i++) {
        magnitude[i] = analyzer->magnitude_avg[i];
        frequencies[i] = i * analyzer->config.sample_rate / analyzer->config.fft_size;
    }
}

// Display spectrum with ASCII art
void display_spectrum(analyzer_t* analyzer) {
    const int display_bins = 32;
    const int bar_width = 50;
    
    double* magnitude = (double*)calloc(analyzer->config.fft_size/2 + 1, sizeof(double));
    double* frequencies = (double*)calloc(analyzer->config.fft_size/2 + 1, sizeof(double));
    if (!magnitude || !frequencies) {
        free(magnitude);
        free(frequencies);
        return;
    }
    
    get_spectrum(analyzer, magnitude, frequencies);
    
    // Clear screen (Unix/Linux)
    printf("\033[2J\033[H");
    
    printf("Real-time FFT Spectrum Analyzer\n");
    printf("===============================\n");
    printf("Sample Rate: %.0f Hz | FFT Size: %d | Frame: %d\n\n",
           analyzer->config.sample_rate, analyzer->config.fft_size, analyzer->frame_count);
    
    // Find max magnitude for scaling
    double max_mag = 0;
    for (int i = 1; i < analyzer->config.fft_size/2; i++) {
        if (magnitude[i] > max_mag) max_mag = magnitude[i];
    }
    
    // Display spectrum bars
    int bin_step = (analyzer->config.fft_size/2) / display_bins;
    if (bin_step < 1) bin_step = 1;
    int bars = analyzer->config.fft_size/2;
    if (bars > display_bins) bars = display_bins;
    
    for (int i = 0; i < bars; i++) {
        int bin = i * bin_step;
        double freq = frequencies[bin];
        /* Before the first hop has completed, every average is still zero and
         * dividing by max_mag printed NaN for the whole first display. */
        double mag_db = (max_mag > 0.0)
                            ? 20 * log10(magnitude[bin] / max_mag + 1e-10)
                            : -60.0;
        int bar_length = (int)((mag_db + 60) / 60 * bar_width);
        if (bar_length < 0) bar_length = 0;
        if (bar_length > bar_width) bar_length = bar_width;
        
        printf("%5.0f Hz |", freq);
        for (int j = 0; j < bar_length; j++) printf("█");
        for (int j = bar_length; j < bar_width; j++) printf(" ");
        printf("| %6.1f dB\n", mag_db);
    }
    
    free(magnitude);
    free(frequencies);
}

// Simulate audio stream
void simulate_audio_stream(analyzer_t* analyzer, int duration_seconds) {
    int samples_per_frame = 256;
    double* audio_frame = (double*)malloc(samples_per_frame * sizeof(double));
    if (!audio_frame)
        return;
    
    int total_frames = (duration_seconds * analyzer->config.sample_rate) / samples_per_frame;
    
    for (int frame = 0; frame < total_frames; frame++) {
        // Generate test signal with time-varying frequency
        double base_freq = 100 + 50 * sin(2 * PI * 0.5 * frame * samples_per_frame / analyzer->config.sample_rate);
        
        for (int i = 0; i < samples_per_frame; i++) {
            double t = (frame * samples_per_frame + i) / analyzer->config.sample_rate;
            audio_frame[i] = 0.5 * sin(2 * PI * base_freq * t) +
                           0.3 * sin(2 * PI * 200 * t) +
                           0.2 * sin(2 * PI * 500 * t) +
                           0.1 * ((double)rand() / RAND_MAX - 0.5);
        }
        
        // Process frame
        process_frame(analyzer, audio_frame, samples_per_frame);
        
        // Display every 10th frame
        if (frame % 10 == 0) {
            display_spectrum(analyzer);
            usleep(50000);  // 50ms delay for visualization
        }
    }
    
    free(audio_frame);
}

/* Time one forward transform, growing the batch until clock() resolves it.
 * The analyzer reports this number, so it must come from the machine rather
 * than from a constant in the demo. */
static double measure_fft_ms(int n) {
    complex_t* work = allocate_complex_array(n);
    complex_t* seed = allocate_complex_array(n);
    if (!work || !seed) {
        free_complex_array(work);
        free_complex_array(seed);
        return 0.0;
    }
    for (int i = 0; i < n; i++)
        seed[i] = sin(TWO_PI * 440.0 * i / 44100.0);

    fft_timer_t timer;
    int reps = 1;
    double total = 0.0;
    for (;;) {
        timer_start(&timer);
        for (int r = 0; r < reps; r++) {
            memcpy(work, seed, (size_t)n * sizeof(complex_t));
            radix2_dit_fft(work, n, FFT_FORWARD);
        }
        timer_stop(&timer);
        total = timer.elapsed_ms;
        if (total >= 10.0 || reps >= (1 << 16)) break;
        reps *= 4;
    }

    free_complex_array(work);
    free_complex_array(seed);
    return total / reps;
}

// Frequency tracker
typedef struct {
    double frequency;
    double magnitude;
    double phase;
    int active;
} tracked_peak_t;

void track_peaks(analyzer_t* analyzer, tracked_peak_t* peaks, int max_peaks) {
    double* magnitude = (double*)calloc(analyzer->config.fft_size/2 + 1, sizeof(double));
    double* frequencies = (double*)calloc(analyzer->config.fft_size/2 + 1, sizeof(double));

    // Reset peaks first so a failed allocation leaves the caller with a
    // complete, inactive array rather than stale entries.
    for (int i = 0; i < max_peaks; i++) {
        peaks[i].active = 0;
    }
    if (!magnitude || !frequencies) {
        free(magnitude);
        free(frequencies);
        return;
    }

    get_spectrum(analyzer, magnitude, frequencies);
    
    // Find local maxima
    int peak_count = 0;
    for (int i = 2; i < analyzer->config.fft_size/2 - 2 && peak_count < max_peaks; i++) {
        if (magnitude[i] > magnitude[i-1] && magnitude[i] > magnitude[i+1] &&
            magnitude[i] > magnitude[i-2] && magnitude[i] > magnitude[i+2]) {
            
            // Quadratic interpolation for precise frequency
            double y1 = magnitude[i-1];
            double y2 = magnitude[i];
            double y3 = magnitude[i+1];
            double delta = 0.5 * (y1 - y3) / (y1 - 2*y2 + y3);
            
            peaks[peak_count].frequency = (i + delta) * analyzer->config.sample_rate / analyzer->config.fft_size;
            peaks[peak_count].magnitude = y2;
            peaks[peak_count].phase = carg(analyzer->window_buffer[i]);
            peaks[peak_count].active = 1;
            peak_count++;
        }
    }
    
    free(magnitude);
    free(frequencies);
}

// Main demonstration
int main() {
    printf("Real-time FFT Spectrum Analyzer Demo\n");
    printf("====================================\n\n");
    
    // Configuration
    analyzer_config_t config = {
        .fft_size = 2048,
        .hop_size = 512,
        .sample_rate = 44100,
        .window_type = "hann",
        .averaging_frames = 4
    };
    
    // Create analyzer
    analyzer_t* analyzer = create_analyzer(config);
    if (!analyzer) {
        fprintf(stderr, "Failed to create analyzer\n");
        return 1;
    }
    
    printf("Configuration:\n");
    printf("- FFT Size: %d\n", config.fft_size);
    printf("- Hop Size: %d (%.1f%% overlap)\n", 
           config.hop_size, 100.0 * (1.0 - (double)config.hop_size / config.fft_size));
    printf("- Sample Rate: %.0f Hz\n", config.sample_rate);
    printf("- Window: %s\n", config.window_type);
    printf("- Frequency Resolution: %.2f Hz\n", config.sample_rate / config.fft_size);
    printf("- Time Resolution: %.2f ms\n\n", 1000.0 * config.hop_size / config.sample_rate);
    
    printf("Press Enter to start real-time analysis...");
    getchar();
    
    // Run real-time simulation
    printf("\nSimulating real-time audio stream...\n");
    simulate_audio_stream(analyzer, 10);  // 10 seconds
    
    // Peak tracking demo
    printf("\n\nPeak Tracking Demo:\n");
    printf("==================\n");
    
    tracked_peak_t peaks[10];
    
    // Generate test signal with multiple tones
    double* test_signal = (double*)malloc(config.fft_size * sizeof(double));
    if (!test_signal) {
        destroy_analyzer(analyzer);
        return 1;
    }
    for (int i = 0; i < config.fft_size; i++) {
        double t = i / config.sample_rate;
        test_signal[i] = sin(2 * PI * 440 * t) +     // A4
                        0.7 * sin(2 * PI * 554 * t) + // C#5
                        0.5 * sin(2 * PI * 659 * t);  // E5
    }
    
    // Process and track peaks
    process_frame(analyzer, test_signal, config.fft_size);
    track_peaks(analyzer, peaks, 10);
    
    printf("\nDetected peaks:\n");
    for (int i = 0; i < 10; i++) {
        if (peaks[i].active) {
            printf("Peak %d: %.1f Hz (magnitude: %.2f, phase: %.2f rad)\n",
                   i+1, peaks[i].frequency, peaks[i].magnitude, peaks[i].phase);
        }
    }
    
    // Performance metrics
    printf("\n\nPerformance Metrics:\n");
    printf("===================\n");
    
    double fft_time = measure_fft_ms(config.fft_size);
    if (fft_time > 0.0) {
        /* One FFT is due every hop, so the budget per transform is the hop
         * duration, not fft_size / hop_size times it. */
        double hop_ms = 1000.0 * config.hop_size / config.sample_rate;
        double max_throughput = hop_ms / fft_time;
        printf("- FFT computation time: %.4f ms (measured)\n", fft_time);
        printf("- Max real-time factor: %.1fx\n", max_throughput);
    } else {
        printf("- FFT computation time: unavailable (allocation failed)\n");
    }
    printf("- Latency: %.1f ms\n", 1000.0 * config.fft_size / config.sample_rate);
    
    // Tips for real implementation
    printf("\n\nImplementation Tips:\n");
    printf("===================\n");
    printf("1. Use a ring buffer for the sliding window\n");
    printf("2. Pre-compute window functions\n");
    printf("3. Reuse an fft_plan_t instead of planning per frame\n");
    printf("4. Implement double-buffering for capture and processing\n");
    printf("5. Size the FFT to the hop budget on the target machine\n");
    
    // Cleanup
    free(test_signal);
    destroy_analyzer(analyzer);
    
    return 0;
}
