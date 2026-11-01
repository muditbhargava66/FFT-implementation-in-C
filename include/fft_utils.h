#ifndef FFT_UTILS_H
#define FFT_UTILS_H

#include "fft_common.h"

/* Signal generation, extra windows, file I/O and spectrum helpers. These
 * lived in utils/fft_utils.c without a declaration, so only that translation
 * unit could call them. */

/* Signal generation */
void generate_chirp_signal(complex_t *signal, int n, double f_start,
                           double f_end, double sample_rate);
void generate_noise(complex_t *signal, int n, double amplitude,
                    unsigned int seed);
void generate_multi_tone(complex_t *signal, int n, double *frequencies,
                         double *amplitudes, int num_tones,
                         double sample_rate);

/* Windows that are not in fft_common.h. A window shorter than two samples is
 * left unchanged, matching the others. */
void apply_window_kaiser(complex_t *signal, int n, double beta);
void apply_window_tukey(complex_t *signal, int n, double alpha);

/* File I/O. The format is the one save_complex_array writes:
 * comment lines, a "# Size: N" header, then "index real imag magnitude
 * phase" rows. save returns 0 on success. load allocates the array on
 * success and returns 0; on failure it returns -1, sets *data to NULL and
 * *n to 0. Free the array with free_complex_array(). */
int save_complex_array(const char *filename, complex_t *data, int n);
int load_complex_array(const char *filename, complex_t **data, int *n);

/* Spectrum helpers */
void find_spectral_peaks(complex_t *spectrum, int n, double sample_rate,
                         double *peak_freqs, double *peak_mags, int *num_peaks,
                         int max_peaks);
double compute_snr(complex_t *signal, complex_t *noise, int n);
void print_spectrum_plot(complex_t *spectrum, int n, double sample_rate,
                         int width);

/* Plot script, zero padding and modulation */
void save_gnuplot_script(const char *script_name, const char *data_file,
                         const char *output_image, const char *title);
complex_t *zero_pad(complex_t *signal, int original_size, int padded_size);
void frequency_shift(complex_t *signal, int n, double shift_freq,
                     double sample_rate);

#endif /* FFT_UTILS_H */
