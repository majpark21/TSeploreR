#' Extract Peak
#'
#' A wrapper for different peak extraction methods. In addition, has some options
#' to detect peak position to give a more precise estimate of peak position.
#' @param y numeric vector. The series from which to extract peak positions.
#' @param method Peak extraction method(s). One or several of ("zscore",
#'   "spike", "bump", "wavelet", "extrema").
#' @param all.winlen numeric. Set a common window length for all methods relying
#'   on rolling mean. Will overwrite any xxx.winlen argument.
#' @param z.winlen numeric. Window length for rolling mean in "zscore" method.
#' @param z.nsd numeric. Number of standard deviation away from rolling mean for
#'   "zscore" method.
#' @param z.influence numeric. How much should peaks affect rolling mean? 0
#'   means that they are not taken into account for computing rolling mean. 1
#'   means that the rolling mean is computed without regards of points being
#'   peaks or not. Used in "zscore" method.
#' @param z.robust logical. If TRUE rolling median and MAD replace rolling mean
#'   and standard deviation. Used in "zscore" method.
#' @param sp.winlen numeric. Window length for rolling mean in "spike" method.
#' @param sp.nsd numeric. Number of standard deviation away from rolling mean
#'   for "spike" method.
#' @param sp.roi numeric vector of length 2. Where to look for spikes in "spike"
#'   method.
#' @param bp.winlen numeric. Window length for rolling mean in "bump" method.
#' @param bp.mind numeric. Minimum distance between peaks for "bump" method. See
#'   ?peakpick::peakpick - neighlim.
#' @param bp.maxderiv numeric. Upper limit for the estimatied derivative for a
#'   point to be considered for a peak call. See ?peakpick::peakpick -
#'   deriv.lim.
#' @param bp.nsd numeric. Number of standard deviation away from rolling mean
#'   for "bump" method.
#' @param bp.npos integer. Peak standard deviations and means will be estimated
#'   plus/minus npos positions from peak. See ?peakpick::peakpick - peak.npos.
#' @param wt.widths numeric vector. widths to use for calculating the CWT
#'   matrix. In general, this range should cover the expected width of peaks of
#'   interest. for "wavelet" method. See scipy.signal.find_peaks_cwt
#'   documentation.
#' @param wt.snr Minimum SNR ratio. See scipy.signal.find_peaks_cwt
#'   documentation.
#' @param wt.centile_noise When calculating the noise floor, percentile of data
#'   points examined below which to consider noise. See
#'   scipy.signal.find_peaks_cwt documentation.
#' @param ext.winlen Window length for looking for local extrema. For "extrema"
#'   method.
#' @param ext.what One of "mini" or "maxi". For "extrema" method.
#' @param show.warning Logical, disable warning for the current execution. Useful when providing
#' all.winlen and applying the function multiple times.
#' @param recenter numeric or NULL. If numeric, attempt to reposition peak position to a more precise
#' location. Specifically, starting from peak estimated position, look for maximum in a window centered
#' around it. The length of the window corresponds to the parameter value. Highly recommended for "spike"
#' and "bump" method which tend to give an imprecise estimation of peak position.
#' @param filter.thresh numeric or NULL. If numeric, attempt at filtering peaks on the base of a 
#' minimal amplitude. To do so, TS are detrended by subtracting a rolling mean and
#' amplitude to baseline is computed. Peak is kept if this amplitude is above the "filter" value.
#' If NULL, no filter is applied. It is highly recommended to recenter the peaks if this option is used.
#' @param filter.winlen size of the window for moving average used to detrend prior to filtering. It is
#' highly recommended to check the detrended signal (part of the output if filter is specified): in the
#' detrended signal, area between peaks should oscillate around 0. If not, try to adjust the value.
#' @param detect.early numeric or NULL. If numeric, extrapolate the time series by the corresponding
#' number of points at its start. This aims at improving peak detection at the beginning of the series.
#' Extrapolation is performed by ARIMA model (see ?forecast::auto.arima)
#'
#' @details Zscore method is based on
#'   https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data#22640362.
#'
#'
#'   Spike and Bump methods are wrapper for the function of the peakPick
#'   package: https://cran.r-project.org/web/packages/peakPick/peakPick.pdf. For
#'   the bump method, a rolling mean is first applied to the signal as
#'   recommended by the package authors.
#'
#'   The wavelet method relies on scipy (python module) peak finder
#'   https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks_cwt.html.
#'
#' @return A data table containing the signal along with results all of selected
#'   methods. Each method result is stored in a different column. If filter is specified,
#'   also returns the detrended signal.
#' @export
#' 
extract_peak <-
  function(y,
           method = c("spike", "bump"),
           all.winlen = NULL,
           z.winlen = NULL,
           z.nsd = 3,
           z.influence = 0.5,
           z.robust = TRUE,
           sp.winlen = NULL,
           sp.nsd = 3,
           sp.roi = NULL,
           bp.winlen = NULL,
           bp.mind = NULL,
           bp.maxderiv = 0.04,
           bp.nsd = 0.5,
           bp.npos = 10L,
           wt.widths = NULL,
           wt.snr = 1,
           wt.centile_noise = 10,
           ext.winlen = NULL,
           ext.what = "maxi",
           recenter = 2,
           filter.thresh = NULL,
           filter.winlen = all.winlen*2,
           detect.early = NULL,
           show.warning = TRUE) {
    require(data.table)
    require(peakPick)
    require(reticulate)
    
    # Save vector with call for detect.early recursive call
    argList <- c(as.list(environment()))
    
    vmethod = c("zscore", "spike", "bump", "wavelet", "extrema")
    if (any(!method %in% vmethod)) {
      stop(
        "method must be one (or several) of c(\"zscore\", \"spike\", \"bump\", \"wavelet\", \"extrema\")"
      )
    }
    
    # Default values
    if (!is.null(all.winlen)) {
      z.winlen <- all.winlen
      sp.winlen <- all.winlen
      bp.winlen <- all.winlen
      ext.winlen <- all.winlen
      if(show.warning) warning("When all.winlen is provided, it overwrites all other xxx.winlen arguments.")
    }
    if ("spike" %in% method & is.null(sp.roi)) {
      # Look for spikes everywhere in the signal
      sp.roi = c(sp.winlen + 1, length(y) - sp.winlen)
    }
    
    # Initialize table of results
    dt.out <- data.table(signal = y)
    
    if ("zscore" %in% method) {
      # From https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data#22640362
      # Detect peaks by defining a rolling mean+sd enveloppe
      peak_smoothZ <-
        function(y, lag, threshold, influence, robust = T) {
          signals <- rep(0, length(y))
          filteredY <- y[0:lag]
          avgFilter <- NULL
          stdFilter <- NULL
          if (robust) {
            avgFilter[lag] <- median(y[0:lag])
            stdFilter[lag] <- mad(y[0:lag])
          } else {
            avgFilter[lag] <- mean(y[0:lag])
            stdFilter[lag] <- sd(y[0:lag])
          }
          
          for (i in (lag + 1):length(y)) {
            if (abs(y[i] - avgFilter[i - 1]) > threshold * stdFilter[i - 1]) {
              if (y[i] > avgFilter[i - 1]) {
                signals[i] <- 1
                
              } else {
                signals[i] <- -1
                
              }
              filteredY[i] <-
                influence * y[i] + (1 - influence) * filteredY[i - 1]
            } else {
              signals[i] <- 0
              filteredY[i] <- y[i]
            }
            if (robust) {
              avgFilter[i] <- median(filteredY[(i - lag):i])
              stdFilter[i] <- mad(filteredY[(i - lag):i])
            } else {
              avgFilter[i] <- mean(filteredY[(i - lag):i])
              stdFilter[i] <- sd(filteredY[(i - lag):i])
            }
          }
          return(list(
            "signals" = signals,
            "avgFilter" = avgFilter,
            "stdFilter" = stdFilter
          ))
        }
      dt.out[, zscore := peak_smoothZ(y, z.winlen, z.nsd, z.influence, z.robust)]
    }
    
    if ("spike" %in% method) {
      # Spike detection with peakPick package, add the trajectory twice to a matrix to avoid bug
      mat.y <- matrix(rep(y, 2), ncol = 2)
      dt.out[, spike := detect.spikes(mat.y, sp.roi, sp.winlen, sp.nsd, verbose = F)[, 1]]
    }
    
    if ("bump" %in% method) {
      # Bump detection with peakPick package, add the trajectory twice to a matrix to avoid bug
      mat.y <- matrix(rep(y, 2), ncol = 2)
      # First smooth the signal with rolling mean, padded at extremeties with linear extrapolation
      dt.out[, bump := peakpick(
        apply(mat.y, 2, rollex, k = bp.winlen),
        neighlim = bp.mind,
        deriv.lim = bp.maxderiv,
        peak.min.sd = bp.nsd,
        peak.npos = bp.npos
      )[, 1]]
    }
    
    if ("wavelet" %in% method) {
      scipy.signal <- import("scipy.signal")
      pos.peaks <-
        scipy.signal$find_peaks_cwt(y,
                                    widths = wt.widths,
                                    min_snr = wt.snr,
                                    noise_perc = wt.centile_noise)
      v.wavelet <- rep(FALSE, length(y))
      v.wavelet[pos.peaks] <- TRUE
      dt.out[, wavelet := v.wavelet]
    }
    
    if ("extrema" %in% method) {
      pos.peaks <- TSexploreR::detect.peak(y, ext.winlen, ext.what)
      v.ext <- rep(FALSE, length(y))
      v.ext[pos.peaks] <- TRUE
      dt.out[, extrema := v.ext]
    }
    
    # Extrapolation for early peaks
    if(!is.null(detect.early)){
      # Rerun all peak extraction, with extrapolation
      extrapolate_extract_peak <- function(y, h, args){
        require(forecast)
        # Reverse because peaks are missed in the beginning of series
        # Forecast with arima model
        yr <- rev(y)
        fit <- auto.arima(yr)
        y_extend <- rev(c(yr, forecast(fit, h = h)$mean))
        # Extract peaks and trim back to initial size
        args[["y"]] <- y_extend
        args[["detect.early"]] <- NULL # Prevent infinite recursion
        args[["show.warning"]] <- FALSE
        args[["recenter"]] <- NULL
        args[["filter.thresh"]] <- NULL
        print(args)
        peaks <- do.call(extract_peak, args)
        return(peaks[(h+1):nrow(peaks)])
      }
      dt.out <- extrapolate_extract_peak(y, detect.early, argList)
      print(dt.out)
    }
    
    # Centering (repositioning to peak maximum)
    if(!is.null(recenter)){
      recenter_peaks <- function(y, positions, tol.winlen){
        #y: signal with potentially several single peaks
        #positions boolean indicating max_pos of peaks
        #recenter peaks at maximum around mask in tolerance window
        pos_recenter <- rep(FALSE, length(y))
        for(pos in positions){
          left <- ifelse((pos-tol.winlen) < 1, 1, pos-tol.winlen)
          right <- ifelse((pos+tol.winlen) > length(y), length(y), pos+tol.winlen)
          shift <- which.max(y[left:right]) - 1
          pos_recenter[left + shift] <- TRUE
        }
        return(pos_recenter)
      }
      # Apply to each method
      for(j in 2:ncol(dt.out)){
        dt.out[, (colnames(dt.out)[j]) := recenter_peaks(signal, which(dt.out[[j]]), recenter)]
      }
    }
    
    # Filtering based on amplitude (in detrended data)
    if(!is.null(filter.thresh)){
      filter_in_window <- function(y, positions, thresh){
        y_filtered <- rep(FALSE, length(y))
        for(pos in positions){
          if(y[pos] >= thresh) y_filtered[pos] <- TRUE
        }
        return(y_filtered)
      }
      # Add detrended signal to output and reorder
      dt.out[, detrended.signal := signal - rollex(signal, filter.winlen)]
      setcolorder(dt.out, c("signal", "detrended.signal", method))
      for(j in 3:ncol(dt.out)){
        dt.out[, (colnames(dt.out)[j]) := filter_in_window(detrended.signal, which(dt.out[[j]]), filter.thresh)]
      }
    }
    return(dt.out)
  }
