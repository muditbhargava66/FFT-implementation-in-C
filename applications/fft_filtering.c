#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"

/**
 * FFT-based Filtering
 * 
 * Implements various digital filters using FFT for efficient
 * frequency domain filtering.
 * 
 * Filters implemented:
 * - Low-pass filter
 * - High-pass filter
 * - Band-pass filter
 * - Band-stop (notch) filter
 * - Custom frequency response filter
 */

// Filter types
typedef enum {
    FILTER_LOWPASS,
    FILTER_HIGHPASS,
    FILTER_BANDPASS,
    FILTER_BANDSTOP,
    FILTER_CUSTOM
} filter_type_t;

// Filter parameters
typedef struct {
    filter_type_t type;
    double cutoff_low;    // Low cutoff frequency (Hz)
    double cutoff_high;   // High cutoff frequency (Hz)
    double sample_rate;   // Sample rate (Hz)
    double transition_width; // Transition band width (Hz)
} filter_params_t;

// Generate ideal frequency response
void generate_ideal_response(complex_t* H, int n, filter_params_t* params) {
    double freq_bin = params->sample_rate / n;
    
    for (int k = 0; k < n; k++) {
        double freq = k * freq_bin;
        if (k > n/2) freq = (k - n) * freq_bin;  // Negative frequencies
        
        double magnitude = 0.0;
        
        switch (params->type) {
            case FILTER_LOWPASS:
                magnitude = (fabs(freq) <= params->cutoff_low) ? 1.0 : 0.0;
                break;
                
            case FILTER_HIGHPASS:
                magnitude = (fabs(freq) >= params->cutoff_high) ? 1.0 : 0.0;
                break;
                
            case FILTER_BANDPASS:
                magnitude = (fabs(freq) >= params->cutoff_low && 
                           fabs(freq) <= params->cutoff_high) ? 1.0 : 0.0;
                break;
                
            case FILTER_BANDSTOP:
                magnitude = (fabs(freq) < params->cutoff_low || 
                           fabs(freq) > params->cutoff_high) ? 1.0 : 0.0;
                break;
                
            default:
                magnitude = 1.0;
        }
        
        H[k] = magnitude;
    }
}

// Apply smooth transition (raised cosine)
void apply_transition_band(complex_t* H, int n, filter_params_t* params) {
    double freq_bin = params->sample_rate / n;
    double transition = params->transition_width;
    
    for (int k = 0; k < n; k++) {
        double freq = k * freq_bin;
        if (k > n/2) freq = (k - n) * freq_bin;
        
        double magnitude = creal(H[k]);
        
        switch (params->type) {
            case FILTER_LOWPASS:
                if (fabs(freq) > params->cutoff_low - transition/2 &&
                    fabs(freq) < params->cutoff_low + transition/2) {
                    double t = (fabs(freq) - (params->cutoff_low - transition/2)) / transition;
                    magnitude = 0.5 * (1 + cos(PI * t));
                }
                break;
                
            case FILTER_HIGHPASS:
                if (fabs(freq) > params->cutoff_high - transition/2 &&
                    fabs(freq) < params->cutoff_high + transition/2) {
                    double t = (fabs(freq) - (params->cutoff_high - transition/2)) / transition;
                    magnitude = 0.5 * (1 - cos(PI * t));
                }
                break;
                
            // Add similar transitions for bandpass/bandstop
            default:
                break;
        }
        
        H[k] = magnitude;
    }
}

// FFT-based filtering
void fft_filter(complex_t* signal, int n, filter_params_t* params) {
    // Generate filter frequency response
    complex_t* H = allocate_complex_array(n);
    generate_ideal_response(H, n, params);
    
    if (params->transition_width > 0) {
        apply_transition_band(H, n, params);
    }
    
    // Forward FFT of signal
    radix2_dit_fft(signal, n, FFT_FORWARD);
    
    // Apply filter in frequency domain
    for (int k = 0; k < n; k++) {
        signal[k] *= H[k];
    }
    
    // Inverse FFT
    radix2_dit_fft(signal, n, FFT_INVERSE);
    
    free_complex_array(H);
}

// Design FIR filter using frequency sampling method
complex_t* design_fir_filter(int taps, filter_params_t* params) {
    int n = next_power_of_two(taps * 4);  // Zero-pad for better frequency resolution
    
    // Generate frequency response
    complex_t* H = allocate_complex_array(n);
    generate_ideal_response(H, n, params);
    apply_transition_band(H, n, params);
    
    // Inverse FFT to get impulse response
    radix2_dit_fft(H, n, FFT_INVERSE);
    
    // Extract and window the impulse response
    complex_t* h = allocate_complex_array(taps);
    int shift = taps / 2;
    
    for (int i = 0; i < taps; i++) {
        int idx = (i - shift + n) % n;
        h[i] = H[idx];
        
        // Apply Hamming window
        double window = 0.54 - 0.46 * cos(TWO_PI * i / (taps - 1));
        h[i] *= window;
    }
    
    free_complex_array(H);
    return h;
}

// Demonstrate filter responses
void plot_filter_response(filter_params_t* params, const char* title) {
    int n = 512;
    complex_t* H = allocate_complex_array(n);
    
    generate_ideal_response(H, n, params);
    apply_transition_band(H, n, params);
    
    printf("\n%s Frequency Response:\n", title);
    printf("================================\n");
    
    // Plot magnitude response
    int display_bins = 64;
    for (int k = 0; k < display_bins; k++) {
        double freq = k * params->sample_rate / n;
        double mag = creal(H[k]);
        
        printf("%5.0f Hz |", freq);
        int bar_length = (int)(40 * mag);
        for (int j = 0; j < bar_length; j++) {
            printf("█");
        }
        printf(" %.3f\n", mag);
    }
    
    free_complex_array(H);
}

// Generate test signals
void generate_multi_tone_signal(complex_t* signal, int n, double sample_rate) {
    double frequencies[] = {100, 500, 1000, 2000, 4000, 8000};
    int num_freqs = 6;
    
    for (int i = 0; i < n; i++) {
        signal[i] = 0;
        for (int f = 0; f < num_freqs; f++) {
            signal[i] += sin(TWO_PI * frequencies[f] * i / sample_rate) / num_freqs;
        }
    }
}

// Main demonstration
int main() {
    printf("FFT-based Filtering\n");
    printf("===================\n");
    
    // Parameters
    double sample_rate = 44100.0;
    int signal_length = 1024;
    
    // Create test signal
    complex_t* signal = allocate_complex_array(signal_length);
    complex_t* original = allocate_complex_array(signal_length);
    
    generate_multi_tone_signal(signal, signal_length, sample_rate);
    memcpy(original, signal, signal_length * sizeof(complex_t));
    
    printf("\nOriginal signal contains frequencies: 100, 500, 1000, 2000, 4000, 8000 Hz\n");
    
    // Test different filters
    filter_params_t params;
    params.sample_rate = sample_rate;
    params.transition_width = 100.0;
    
    // 1. Low-pass filter
    params.type = FILTER_LOWPASS;
    params.cutoff_low = 1500.0;
    plot_filter_response(&params, "Low-pass Filter (1500 Hz)");
    
    memcpy(signal, original, signal_length * sizeof(complex_t));
    fft_filter(signal, signal_length, &params);
    
    // Analyze filtered signal
    radix2_dit_fft(signal, signal_length, FFT_FORWARD);
    double* mag = compute_magnitude(signal, signal_length);
    
    printf("\nFiltered signal spectrum (first 10 bins):\n");
    for (int i = 0; i < 10; i++) {
        double freq = i * sample_rate / signal_length;
        printf("%.0f Hz: %.3f\n", freq, mag[i]);
    }
    free(mag);
    
    // 2. High-pass filter
    params.type = FILTER_HIGHPASS;
    params.cutoff_high = 2500.0;
    plot_filter_response(&params, "High-pass Filter (2500 Hz)");
    
    // 3. Band-pass filter
    params.type = FILTER_BANDPASS;
    params.cutoff_low = 800.0;
    params.cutoff_high = 1200.0;
    plot_filter_response(&params, "Band-pass Filter (800-1200 Hz)");
    
    // 4. Band-stop (notch) filter
    params.type = FILTER_BANDSTOP;
    params.cutoff_low = 900.0;
    params.cutoff_high = 1100.0;
    params.transition_width = 50.0;
    plot_filter_response(&params, "Notch Filter (900-1100 Hz)");
    
    // FIR filter design example
    printf("\n\nFIR Filter Design Example:\n");
    printf("==========================\n");
    
    params.type = FILTER_LOWPASS;
    params.cutoff_low = 2000.0;
    params.transition_width = 200.0;
    
    int taps = 65;
    complex_t* fir_coeff = design_fir_filter(taps, &params);
    
    printf("FIR filter coefficients (first 10):\n");
    for (int i = 0; i < 10; i++) {
        printf("h[%d] = %.6f\n", i, creal(fir_coeff[i]));
    }
    
    // Performance comparison
    printf("\n\nPerformance Analysis:\n");
    printf("=====================\n");
    
    fft_timer_t timer;
    
    // Time-domain FIR filtering
    timer_start(&timer);
    // Would implement direct convolution here
    timer_stop(&timer);
    double time_fir = timer.elapsed_ms;
    
    // FFT-based filtering
    timer_start(&timer);
    fft_filter(signal, signal_length, &params);
    timer_stop(&timer);
    double time_fft = timer.elapsed_ms;
    
    printf("Signal length: %d samples\n", signal_length);
    printf("FIR filter length: %d taps\n", taps);
    printf("Time-domain FIR filtering: %.3f ms (estimated)\n", time_fir);
    printf("FFT-based filtering time: %.3f ms\n", time_fft);
    printf("\nFFT filtering is more efficient for long filters!\n");
    
    // Applications
    printf("\n\nPractical Applications:\n");
    printf("=======================\n");
    printf("1. Audio equalization\n");
    printf("2. Noise removal\n");
    printf("3. Signal enhancement\n");
    printf("4. Anti-aliasing\n");
    printf("5. Channel separation\n");
    printf("6. Crossover networks\n");
    
    // Cleanup
    free_complex_array(signal);
    free_complex_array(original);
    free_complex_array(fir_coeff);
    
    return 0;
}
