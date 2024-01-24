#' Synthetic signal.
#'
#' Generates a synthetic signal with peaks with normal distribution.
#'
#' @param peaks Can be an scalar with the number of peaks wanted in the signal
#' or a vector with the starting points for each peak.
#' @param length_out The final length of the desired signal (in indexes).
#' @param intensity The maximum intensity of the peaks, if is a scalar all
#' peaks will have the same intensity, if a vector each value will correspond to
#' a peak.
#' @param alpha For the shape of each peak, if a vector each value will
#' correspond to a peak, changes the shape to the first half of each peak
#' @param xi For the shape of each peak, if a vector each value will
#' correspond to a peak, changes the center of the peak
#' @param length_peaks The length of each peak (in indexes), if a scalar all
#' peaks will have the same length, if a vector each value will correspond to
#' a peak.
#' @param baseline A vector that will be added to the base signal, if != NULL
#' length_out will be ignored and calculated with the length of this vector.
#' @param noise An scalar with the intensity of the noise desired.
#'
#' @return A vector with the designed signal.
#' @export
#'
#' @examples
#' synthetic_signal()
#' synthetic_signal(peaks = 5)
#' synthetic_signal(intensity = c(100, 600, 1000),
#'                  baseline = seq(0, 150, by = 0.1))
#' synthetic_signal(peaks = 10, length_out = 1500, noise = 10)

synthetic_signal <- function(peaks = 2,
                             length_out = 1000,
                             intensity = 1000,
                             alpha= 0,
                             xi = 0,
                             length_peaks = NULL,
                             baseline = NULL,
                             noise = NULL){
  # Check dimensions of input arguments
  if (length(peaks) == 1){
    n_peaks <- max(peaks,
                   length(intensity),
                   length(alpha),
                   length(xi),
                   length(length_peaks))
  } else {
    n_peaks <- max(length(peaks),
                   length(intensity),
                   length(alpha),
                   length(xi),
                   length(length_peaks))
  }

  if (!((length(peaks) == 1 | length(peaks) == n_peaks) &
        (length(intensity) == 1 | length(intensity) == n_peaks) &
        (length(alpha) == 1 | length(alpha) == n_peaks) &
        (length(xi) == 1 | length(xi) == n_peaks) &
        (is.null(length_peaks) | length(length_peaks) == n_peaks))){
    stop("One or more arguments do not coincide with the number of peaks")
  }
  # Create the base signal with or without baseline
  signal <- rep(0, length_out)
  if (is.vector(baseline)){
    signal <- baseline
    length_out <- length(baseline)
  } else if (!is.null(baseline)){
    warning("Baseline must be a vector, no baseline will be used")
  }
  # Create length_peaks vector
  if (length(length_peaks) == 1){
    length_peaks <- rep(length_peaks, n_peaks)
  } else if (is.null(length_peaks)){
    length_peaks <- rep(floor(length_out / n_peaks), n_peaks)
  }
  # Create starts vector
  if (length(peaks) == 1){
    starts <- rep(0, n_peaks)
    starts[1] <- 1
    if (n_peaks >= 2){
      for (i in 2 : n_peaks){
        starts[i] <- starts[i - 1] + length_peaks[i - 1]
      }
    }
  } else {
    starts <- peaks
  }
  # Create alpha vector
  if (length(alpha) == 1){
    alpha <- rep(alpha, n_peaks)
  }
  # Create xi vector
  if (length(xi) == 1){
    xi <- rep(xi, n_peaks)
  }
  # Create intensity vector
  if (length(intensity) == 1){
    intensity <- rep(intensity, n_peaks)
  }
  # Add peaks to signal
  for (i in 1:n_peaks){
    peak <- generate_peak(intensity = intensity[i],
                          length_peak = length_peaks[i],
                          alpha = alpha[i],
                          xi = xi[i])
    signal <- add_peak(signal = signal, peak = peak, start = starts[i])
  }
  signal <- signal[1:length_out]
  # Add noise
  if (!is.null(noise)){
    if (length(noise) == 1){
      signal <- signal + sample(seq(-noise, noise, by = 0.1),
                                length_out,
                                replace = TRUE)
    } else {
      warning("Noise must be scalar, no noise will be added")
    }
  }
  return(signal)
}

#' Synthetic chromatogram.
#'
#' Creates a synthetic chromatogram.
#'
#' @param n_peaks Number of peaks to be added
#' @param length_out length of the signal to be created
#' @param intensity average intensity of the peaks
#' @param random_intensity random maximum value to add or subtract to the
#' intensity
#' @param mov_peaks randomness of the peaks, the smaller it is the most each
#' peak can move
#' @param alpha for skewed normal distribution, modifies the first half of the
#' peak
#' @param xi for skewed normal distribution, moves the center of the
#' distribution
#'
#' @return a vector that contains a synthetic chromatogram
#' @export
#'
#' @examples

synthetic_chromatogram<- function(n_peaks,
                                  length_out = 2000,
                                  intensity = 1000,
                                  random_intensity = NULL,
                                  mov_peaks = NULL,
                                  alpha = 3,
                                  xi = -2){
  max_peak <- 104 * exp(0.4 * seq(0,n_peaks)) - 104 ##FIXME seq con n_peaks no ideal
  max_peak <- utils::tail(which(max_peak < (length_out - floor(length_out / n_peaks))), 1) - 1
  x <- seq(0.1, max_peak, length.out = n_peaks)
  if (!is.null(mov_peaks)){
    mov_peaks <- (x[2]-x[1])/mov_peaks
    mov_peaks <- sample(seq(- mov_peaks, mov_peaks, by =0.01),
                        n_peaks,
                        replace = TRUE)
    x <- x + mov_peaks
  }
  peaks <- 104 * exp(0.4  * x) - 104
  if (peaks[1] < 1 ){
    peaks[1] <- 1
  }
    if (!is.null(random_intensity)){
      intensity <- sample(seq((intensity - random_intensity),
                              (intensity + random_intensity)),
                          n_peaks,
                          replace = TRUE)
    }
  chromatogram <- synthetic_signal(peaks = peaks,
                                   length_out = length_out,
                                   intensity = intensity,
                                   alpha = alpha,
                                   xi= xi)
  return(chromatogram)
}

#' Create Synthetic dataset.
#'
#' @param n_samples number of chromatograms to be generated
#' @param n_peaks number of peaks in each chromatogram
#' @param length_out length of chromatograms
#' @param mov_peaks random movement of the position of the peaks
#' @param intensity intensity of the peaks
#' @param random_intensity random variation for the intensity of the peaks
#'
#' @return synthetic dataset.
#' @export
#'
#' @examples
#' create_synthetic_dataset(10,10,1000,5)
create_synthetic_dataset <- function(n_samples,
                                     n_peaks,
                                     length_out,
                                     mov_peaks,
                                     intensity=1000,
                                     random_intensity=NULL){
  dataset <- matrix(nrow = n_samples, ncol = length_out)
  for (i in 1:n_samples){
    dataset[i,] <- synthetic_chromatogram(n_peaks = n_peaks,
                                          length_out = length_out,
                                          intensity = intensity,
                                          random_intensity = random_intensity,
                                          mov_peaks = mov_peaks)
  }
  return(dataset)
}



#' Peak generator.
#'
#' Generates a peak with normal distribution.
#'
#' @param intensity Maximum intensity of the peak.
#' @param length_peak Length of the peak in indexes.
#' @param alpha for skewed normal distribution, modifies the first half of the
#' peak
#' @param xi for skewed normal distribution, moves the center of the
#' distribution
#'
#' @return A vector with a peak.
#' @keywords internal


generate_peak <- function(intensity = 1, length_peak = 100, alpha = 3, xi = -2){
  p <- sn::dsn(seq(-3, 3, length.out = length_peak), alpha = alpha, xi = xi)
  peak <- ((p - min(p)) / (max(p) - min(p))) * intensity
  return(peak)
}

#' Add peak.
#'
#'Adds a peak to a signal in the desired place.
#'
#' @param signal The original signal.
#' @param peak A vector containing the peak that wants o be added.
#' @param start The position (index) in the signal where the peak will begin.
#'
#' @return The signal with the added peak.
#' @keywords internal

add_peak <- function(signal, peak, start){
  length_peak <- length(peak)
  signal[(start) : (start + length_peak - 1)] <- signal[(start) : (start + length_peak - 1)] + peak
  return(signal)
}
