#include "fft_common.h"
#include "fft_algorithms.h"

/**
 * Power Spectrum Analysis
 * 
 * PSD estimation (periodogram, Welch), autocorrelation, cross-correlation
 * and coherence, with a summary of spectral statistics.
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

// Compute autocorrelation using FFT.
// Returns a heap array the caller must release with free_complex_array().
complex_t* autocorrelation_fft(complex_t* signal, int n) {
    if (!signal || n <= 0) return NULL;

    // Zero-pad to avoid circular correlation
    int n_fft = next_power_of_two(2 * n);
    if (n_fft <= 0) return NULL;

    /* allocate_complex_array() already zeroes the buffer. */
    complex_t* x = allocate_complex_array(n_fft);
    complex_t* acf = allocate_complex_array(n);
    if (!x || !acf) {
        free_complex_array(x);
        free_complex_array(acf);
        return NULL;
    }
    memcpy(x, signal, (size_t)n * sizeof(complex_t));
    
    // FFT
    radix2_dit_fft(x, n_fft, FFT_FORWARD);
    
    // Compute power spectrum
    for (int k = 0; k < n_fft; k++) {
        x[k] = x[k] * conj(x[k]);
    }
    
    // Inverse FFT
    radix2_dit_fft(x, n_fft, FFT_INVERSE);
    
    // Extract autocorrelation (first n samples)
    memcpy(acf, x, (size_t)n * sizeof(complex_t));
    
    free_complex_array(x);
    return acf;
}

// Cross-correlation via FFT.
// Returns a heap array the caller must release with free_complex_array().
complex_t* cross_correlation_fft(complex_t* x, complex_t* y, int n) {
    if (!x || !y || n <= 0) return NULL;

    int n_fft = next_power_of_two(2 * n);
    if (n_fft <= 0) return NULL;

    complex_t* X = allocate_complex_array(n_fft);
    complex_t* Y = allocate_complex_array(n_fft);
    if (!X || !Y) {
        free_complex_array(X);
        free_complex_array(Y);
        return NULL;
    }
    memcpy(X, x, (size_t)n * sizeof(complex_t));
    memcpy(Y, y, (size_t)n * sizeof(complex_t));
    
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
    if (!ccf) {
        free_complex_array(X);
        free_complex_array(Y);
        return NULL;
    }
    memcpy(ccf, X, (size_t)n * sizeof(complex_t));
    
    free_complex_array(X);
    free_complex_array(Y);
    return ccf;
}

// Compute magnitude-squared coherence between two signals
//   |Cxy(f)|^2 = |Pxy(f)|^2 / (Pxx(f) * Pyy(f))
double* compute_coherence(complex_t* x, complex_t* y, int n, int window_size) {
    int half = window_size / 2;
    int num_bins = half + 1;
    double* coherence = (double*)malloc(num_bins * sizeof(double));

    int overlap = half;                       /* 50 % overlap */
    int step = window_size - overlap;
    int num_windows = (n - overlap) / step;
    if (num_windows < 1) num_windows = 1;

    /* Accumulate Pxx, Pyy, and complex Pxy per bin */
    double* psd_xx = (double*)calloc(num_bins, sizeof(double));
    double* psd_yy = (double*)calloc(num_bins, sizeof(double));
    complex_t* psd_xy = allocate_complex_array(num_bins);
    for (int k = 0; k < num_bins; k++)
        psd_xy[k] = 0.0;

    complex_t* wx = allocate_complex_array(window_size);
    complex_t* wy = allocate_complex_array(window_size);

    for (int w = 0; w < num_windows; w++) {
        int start = w * step;

        /* Extract & window both channels */
        for (int i = 0; i < window_size; i++) {
            int idx = start + i;
            double hann = 0.5 * (1.0 - cos(TWO_PI * i / (window_size - 1)));
            wx[i] = (idx < n ? x[idx] : 0.0) * hann;
            wy[i] = (idx < n ? y[idx] : 0.0) * hann;
        }

        radix2_dit_fft(wx, window_size, FFT_FORWARD);
        radix2_dit_fft(wy, window_size, FFT_FORWARD);

        for (int k = 0; k < num_bins; k++) {
            psd_xx[k] += cabs(wx[k]) * cabs(wx[k]);
            psd_yy[k] += cabs(wy[k]) * cabs(wy[k]);
            psd_xy[k] += conj(wx[k]) * wy[k];
        }
    }

    /* Coherence: |Pxy|^2 / (Pxx * Pyy) */
    for (int k = 0; k < num_bins; k++) {
        double denom = psd_xx[k] * psd_yy[k];
        if (denom > 0)
            coherence[k] = cabs(psd_xy[k]) * cabs(psd_xy[k]) / denom;
        else
            coherence[k] = 0.0;
    }

    free(psd_xx);
    free(psd_yy);
    free_complex_array(psd_xy);
    free_complex_array(wx);
    free_complex_array(wy);
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
    if (acf) {
        double lag0 = creal(acf[0]);
        for (int lag = 0; lag < 10; lag++) {
            printf("%d\t%.3f\n", lag, lag0 != 0.0 ? creal(acf[lag]) / lag0 : 0.0);
        }
    } else {
        printf("  unavailable (allocation failed)\n");
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
        
        fft_timer_t timer;

        /* Batch until clock() resolves the interval. One 1024-point
         * periodogram is well under its ~1 ms tick and printed 0.000 ms. */
        int reps = 1;
        double total = 0.0;
        double* psd = NULL;
        for (;;) {
            timer_start(&timer);
            for (int r = 0; r < reps; r++) {
                free(psd);
                psd = compute_periodogram(signal, n, 1000.0);
            }
            timer_stop(&timer);
            total = timer.elapsed_ms;
            if (total >= 10.0 || reps >= (1 << 16)) break;
            reps *= 4;
        }
        
        printf("N=%d: Periodogram time = %.4f ms (avg of %d)\n",
               n, total / reps, reps);
        
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
