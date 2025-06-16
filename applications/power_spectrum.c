#include "../include/fft_common.h"
#include "../include/fft_algorithms.h"

// Window function implementations
void apply_window_hann(complex_t* signal, int n) {
    for (int i = 0; i < n; i++) {
        double window = 0.5 * (1.0 - cos(TWO_PI * i / (n - 1)));
        signal[i] *= window;
    }
}

void apply_window_hamming(complex_t* signal, int n) {
    for (int i = 0; i < n; i++) {
        double window = 0.54 - 0.46 * cos(TWO_PI * i / (n - 1));
        signal[i] *= window;
    }
}

void apply_window_blackman(complex_t* signal, int n) {
    for (int i = 0; i < n; i++) {
        double window = 0.42 - 0.5 * cos(TWO_PI * i / (n - 1)) 
                      + 0.08 * cos(4 * PI * i / (n - 1));
        signal[i] *= window;
    }
}

/**
 * Power Spectrum Analysis
 * 
 * Implements various spectral analysis techniques for signal processing
 * and communications applications.
 * 
 * Features:
 * - Power Spectral Density (PSD) estimation
 * - Welch's method for improved estimates
 * - Autocorrelation via FFT
 * - Cross-correlation and coherence
 * - Spectral statistics
 */

// PSD estimation methods
typedef enum {
    PSD_PERIODOGRAM,
    PSD_WELCH,
    PSD_BARTLETT,
    PSD_BLACKMAN_TUKEY
} psd_method_t;

// Welch's method parameters
typedef struct {
    int window_size;
    int overlap;
    const char* window_type;
    int num_averages;
} welch_params_t;

// Compute periodogram (simple PSD estimate)
double* compute_periodogram(complex_t* signal, int n, double sample_rate) {
    double* psd = (double*)malloc((n/2 + 1) * sizeof(double));
    
    // Copy signal for FFT
    complex_t* x = allocate_complex_array(n);
    memcpy(x, signal, n * sizeof(complex_t));
    
    // Apply window (Hann by default)
    apply_window_hann(x, n);
    
    // Compute FFT
    radix2_dit_fft(x, n, FFT_FORWARD);
    
    // Compute power spectral density
    double window_power = 0.375 * n;  // Hann window power
    double scale = 1.0 / (sample_rate * window_power);
    
    for (int k = 0; k <= n/2; k++) {
        double power = cabs(x[k]) * cabs(x[k]);
        psd[k] = power * scale;
        if (k > 0 && k < n/2) {
            psd[k] *= 2.0;  // Account for negative frequencies
        }
    }
    
    free_complex_array(x);
    return psd;
}

// Welch's method for improved PSD estimation
double* welch_psd(complex_t* signal, int signal_len, double sample_rate, 
                  welch_params_t* params) {
    int window_size = params->window_size;
    int overlap = params->overlap;
    int step = window_size - overlap;
    int num_windows = (signal_len - overlap) / step;
    
    // Allocate averaged PSD
    double* psd_avg = (double*)calloc(window_size/2 + 1, sizeof(double));
    complex_t* window = allocate_complex_array(window_size);
    
    // Process each window
    for (int w = 0; w < num_windows; w++) {
        int start = w * step;
        
        // Extract window
        for (int i = 0; i < window_size; i++) {
            if (start + i < signal_len) {
                window[i] = signal[start + i];
            } else {
                window[i] = 0;
            }
        }
        
        // Compute periodogram for this window
        double* psd_window = compute_periodogram(window, window_size, sample_rate);
        
        // Accumulate
        for (int k = 0; k <= window_size/2; k++) {
            psd_avg[k] += psd_window[k];
        }
        
        free(psd_window);
    }
    
    // Average
    for (int k = 0; k <= window_size/2; k++) {
        psd_avg[k] /= num_windows;
    }
    
    free_complex_array(window);
    return psd_avg;
}

// Compute autocorrelation using FFT
complex_t* autocorrelation_fft(complex_t* signal, int n) {
    // Zero-pad to avoid circular correlation
    int n_fft = next_power_of_two(2 * n);
    complex_t* x = allocate_complex_array(n_fft);
    
    // Copy and zero-pad signal
    memset(x, 0, n_fft * sizeof(complex_t));
    memcpy(x, signal, n * sizeof(complex_t));
    
    // FFT
    radix2_dit_fft(x, n_fft, FFT_FORWARD);
    
    // Compute power spectrum
    for (int k = 0; k < n_fft; k++) {
        x[k] = x[k] * conj(x[k]);
    }
    
    // Inverse FFT
    radix2_dit_fft(x, n_fft, FFT_INVERSE);
    
    // Extract autocorrelation (first n samples)
    complex_t* acf = allocate_complex_array(n);
    memcpy(acf, x, n * sizeof(complex_t));
    
    free_complex_array(x);
    return acf;
}

// Cross-correlation via FFT
complex_t* cross_correlation_fft(complex_t* x, complex_t* y, int n) {
    int n_fft = next_power_of_two(2 * n);
    complex_t* X = allocate_complex_array(n_fft);
    complex_t* Y = allocate_complex_array(n_fft);
    
    // Zero-pad signals
    memset(X, 0, n_fft * sizeof(complex_t));
    memset(Y, 0, n_fft * sizeof(complex_t));
    memcpy(X, x, n * sizeof(complex_t));
    memcpy(Y, y, n * sizeof(complex_t));
    
    // FFT both signals
    radix2_dit_fft(X, n_fft, FFT_FORWARD);
    radix2_dit_fft(Y, n_fft, FFT_FORWARD);
    
    // Cross-spectrum: X* × Y
    for (int k = 0; k < n_fft; k++) {
        X[k] = conj(X[k]) * Y[k];
    }
    
    // Inverse FFT
    radix2_dit_fft(X, n_fft, FFT_INVERSE);
    
    // Extract result
    complex_t* ccf = allocate_complex_array(n);
    memcpy(ccf, X, n * sizeof(complex_t));
    
    free_complex_array(X);
    free_complex_array(Y);
    return ccf;
}

// Compute coherence between two signals
double* compute_coherence(complex_t* x, complex_t* y, int n, int window_size) {
    double* coherence = (double*)malloc((window_size/2 + 1) * sizeof(double));
    
    // Compute cross-spectral density and auto-spectral densities
    welch_params_t params = {
        .window_size = window_size,
        .overlap = window_size / 2,
        .window_type = "hann",
        .num_averages = 0
    };
    
    double* psd_xx = welch_psd(x, n, 1.0, &params);
    double* psd_yy = welch_psd(y, n, 1.0, &params);
    
    // For simplicity, compute magnitude squared coherence
    // |Cxy|² = |Pxy|² / (Pxx × Pyy)
    
    for (int k = 0; k <= window_size/2; k++) {
        if (psd_xx[k] > 0 && psd_yy[k] > 0) {
            // Simplified - would need cross-PSD for full implementation
            coherence[k] = 1.0;  // Placeholder
        } else {
            coherence[k] = 0.0;
        }
    }
    
    free(psd_xx);
    free(psd_yy);
    return coherence;
}

// Spectral statistics
typedef struct {
    double mean_frequency;
    double rms_bandwidth;
    double spectral_centroid;
    double spectral_flux;
    double spectral_rolloff;
    double total_power;
} spectral_stats_t;

spectral_stats_t compute_spectral_statistics(double* psd, int n, double sample_rate) {
    spectral_stats_t stats = {0};
    
    double freq_bin = sample_rate / (2 * n);
    double total_power = 0;
    double weighted_freq_sum = 0;
    
    // Compute total power and weighted frequency sum
    for (int k = 0; k < n; k++) {
        double freq = k * freq_bin;
        total_power += psd[k];
        weighted_freq_sum += freq * psd[k];
    }
    
    stats.total_power = total_power;
    
    // Mean frequency (spectral centroid)
    if (total_power > 0) {
        stats.mean_frequency = weighted_freq_sum / total_power;
        stats.spectral_centroid = stats.mean_frequency;
    }
    
    // RMS bandwidth
    double variance = 0;
    for (int k = 0; k < n; k++) {
        double freq = k * freq_bin;
        double diff = freq - stats.mean_frequency;
        variance += diff * diff * psd[k];
    }
    
    if (total_power > 0) {
        stats.rms_bandwidth = sqrt(variance / total_power);
    }
    
    // Spectral rolloff (95% of power)
    double cumulative_power = 0;
    double rolloff_threshold = 0.95 * total_power;
    
    for (int k = 0; k < n; k++) {
        cumulative_power += psd[k];
        if (cumulative_power >= rolloff_threshold) {
            stats.spectral_rolloff = k * freq_bin;
            break;
        }
    }
    
    return stats;
}

// Demonstration functions
void demo_signal_analysis() {
    printf("\nSignal Analysis Demonstrations:\n");
    printf("===============================\n");
    
    int n = 1024;
    double sample_rate = 1000.0;
    complex_t* signal = allocate_complex_array(n);
    
    // Generate test signal with known properties
    double f1 = 50.0, f2 = 120.0;
    double snr_db = 10.0;
    double noise_power = pow(10, -snr_db/10);
    
    for (int i = 0; i < n; i++) {
        double t = i / sample_rate;
        signal[i] = sin(TWO_PI * f1 * t) + 0.5 * sin(TWO_PI * f2 * t);
        
        // Add white noise
        signal[i] += sqrt(noise_power) * ((double)rand() / RAND_MAX - 0.5);
    }
    
    printf("Test signal: 50 Hz + 120 Hz with %.0f dB SNR\n", snr_db);
    
    // 1. Periodogram
    printf("\n1. Periodogram Analysis:\n");
    double* psd_simple = compute_periodogram(signal, n, sample_rate);
    
    printf("Frequency (Hz)\tPSD (dB/Hz)\n");
    for (int k = 0; k < 10; k++) {
        double freq = k * sample_rate / n;
        double psd_db = 10 * log10(psd_simple[k] + 1e-10);
        printf("%.1f\t\t%.1f\n", freq, psd_db);
    }
    
    // 2. Welch's method
    printf("\n2. Welch's Method (reduced variance):\n");
    welch_params_t welch_params = {
        .window_size = 256,
        .overlap = 128,
        .window_type = "hann"
    };
    
    double* psd_welch = welch_psd(signal, n, sample_rate, &welch_params);
    
    printf("Improved estimates at signal frequencies:\n");
    int bin_50hz = (int)(50.0 * welch_params.window_size / sample_rate);
    int bin_120hz = (int)(120.0 * welch_params.window_size / sample_rate);
    
    printf("50 Hz: %.1f dB/Hz\n", 10 * log10(psd_welch[bin_50hz]));
    printf("120 Hz: %.1f dB/Hz\n", 10 * log10(psd_welch[bin_120hz]));
    
    // 3. Autocorrelation
    printf("\n3. Autocorrelation Analysis:\n");
    complex_t* acf = autocorrelation_fft(signal, n);
    
    printf("Lag\tAutocorrelation\n");
    for (int lag = 0; lag < 10; lag++) {
        printf("%d\t%.3f\n", lag, creal(acf[lag]) / creal(acf[0]));
    }
    
    // 4. Spectral statistics
    printf("\n4. Spectral Statistics:\n");
    spectral_stats_t stats = compute_spectral_statistics(psd_simple, n/2 + 1, sample_rate);
    
    printf("Mean frequency: %.1f Hz\n", stats.mean_frequency);
    printf("RMS bandwidth: %.1f Hz\n", stats.rms_bandwidth);
    printf("Spectral rolloff: %.1f Hz\n", stats.spectral_rolloff);
    printf("Total power: %.3f\n", stats.total_power);
    
    free(psd_simple);
    free(psd_welch);
    free_complex_array(acf);
    free_complex_array(signal);
}

// Main program
int main() {
    printf("Power Spectrum Analysis\n");
    printf("=======================\n");
    
    // Run demonstrations
    demo_signal_analysis();
    
    // Performance comparison
    printf("\n\nPerformance Comparison:\n");
    printf("======================\n");
    
    int sizes[] = {1024, 4096, 16384};
    for (int s = 0; s < 3; s++) {
        int n = sizes[s];
        complex_t* signal = allocate_complex_array(n);
        
        // Generate random signal
        for (int i = 0; i < n; i++) {
            signal[i] = ((double)rand() / RAND_MAX - 0.5);
        }
        
        timer_t timer;
        
        // Time periodogram
        timer_start(&timer);
        double* psd = compute_periodogram(signal, n, 1000.0);
        timer_stop(&timer);
        
        printf("N=%d: Periodogram time = %.3f ms\n", n, timer.elapsed_ms);
        
        free(psd);
        free_complex_array(signal);
    }
    
    // Applications
    printf("\n\nPractical Applications:\n");
    printf("=======================\n");
    printf("1. Vibration analysis in mechanical systems\n");
    printf("2. Communications signal analysis\n");
    printf("3. Biomedical signal processing (EEG, ECG)\n");
    printf("4. Speech and audio analysis\n");
    printf("5. Radar and sonar processing\n");
    printf("6. Seismic data analysis\n");
    printf("7. Financial time series analysis\n");
    
    // Best practices
    printf("\n\nBest Practices:\n");
    printf("===============\n");
    printf("1. Use Welch's method for noisy signals\n");
    printf("2. Choose window size based on frequency resolution needs\n");
    printf("3. Use 50%% overlap for good variance reduction\n");
    printf("4. Apply appropriate windowing to reduce leakage\n");
    printf("5. Consider zero-padding for smoother spectra\n");
    printf("6. Verify stationarity assumption for long signals\n");
    
    return 0;
}
