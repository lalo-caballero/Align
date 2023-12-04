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
#' @param width The width of the peaks, if a scalar the same width will be used
#' for all peaks, if a vector each value will correspond to a peak. The parameter
#' will be passed to function 'dnorm()' where sd = length_peak / width.
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
#' synthetic_signal(intensity = c(100, 600, 1000), baseline = seq(0, 150, by = 0.1))
#' synthetic_signal(peaks = 10, length_out = 1500, noise = 10)

synthetic_signal <- function(peaks = 2,
                             length_out = 1000,
                             intensity = 1000,
                             width= 5,
                             length_peaks = NULL,
                             baseline = NULL,
                             noise = NULL){
  # Check dimensions of input arguments
  if (length(peaks) == 1){
    n_peaks <- max(peaks, length(intensity), length(width), length(length_peaks))
  } else {
    n_peaks <- max(length(peaks), length(intensity), length(width), length(length_peaks))
  }

  if (!((length(peaks) == 1 | length(peaks) == n_peaks) &
        (length(intensity) == 1 | length(intensity) == n_peaks) &
        (length(width) == 1 | length(width) == n_peaks) &
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
  # Create width vector
  if (length(width) == 1){
    width <- rep(width, n_peaks)
  }
  # Create intensity vector
  if (length(intensity) == 1){
    intensity <- rep(intensity, n_peaks)
  }
  # Add peaks to signal
  for (i in 1:n_peaks){
    peak <- generate_peak(intensity = intensity[i], length_peak = length_peaks[i], width = width[i])
    signal <- add_peak(signal = signal, peak = peak, start = starts[i])
  }
  signal <- signal[1:length_out]
  # Add noise
  if (!is.null(noise)){
    if (length(noise) == 1){
      signal <- signal + sample(seq(-noise, noise, by = 0.1), length_out, replace = TRUE)
    } else {
      warning("Noise must be scalar, no noise will be added")
    }
  }
  return(signal)
}

#' Peak generator.
#'
#' Generates a peak with normal distribution.
#'
#' @param intensity Maximum intensity of the peak.
#' @param length_peak Length of the peak in indexes.
#' @param width Width of the peak, this parameter will be passed to 'dnorm()'
#' and used as sd = length_peak / width.
#'
#' @return A vector with a peak.
#' @export
#'
#' @examples
#' generate_peak()
#' generate_peak(intensity = 500)
#' generate_peak(length_peak = 150, width = 20)

generate_peak <- function(intensity = 1, length_peak = 100, width = 5){
  p <- stats::dnorm(seq(( - length_peak / 2) + 1, length_peak / 2), sd = length_peak / width)
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
#' @export
#'
#' @examples
#' signal <- rep(0, 150)
#' peak <- generate_peak()
#' add_peak(signal = signal, peak = peak, start = 25)

add_peak <- function(signal, peak, start){
  length_peak <- length(peak)
  signal[(start) : (start + length_peak - 1)] <- signal[(start) : (start + length_peak - 1)] + peak
  return(signal)
}



synthetic_chromatogram<- function(){
  peaks <- 104 * exp(0.4  * seq(0.5,7,by=0.5)) - 104
  width <- seq(20,3,length.out=length(peaks))
  chromatogram <- synthetic_signal(peaks = peaks,
                                   length_out = 2000,
                                   width = width,
                                   baseline = NULL,
                                   noise = NULL)
  return(chromatogram)
}

