#' Extract Peak
#'
#' A wrapper for different peak extraction methods. In addition, has some
#' options to detect peak position to give a more precise estimate of peak
#' position.
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
#' @param show.warning Logical, disable warning for the current execution.
#'   Useful when providing all.winlen and applying the function multiple times.
#' @param recenter numeric or NULL. If numeric, attempt to reposition peak
#'   position to a more precise location. Specifically, starting from peak
#'   estimated position, look for maximum in a window centered around it. The
#'   length of the window corresponds to the parameter value. Highly recommended
#'   for "spike" and "bump" method which tend to give an imprecise estimation of
#'   peak position.
#' @param filter.thresh numeric or NULL. If numeric, attempt at filtering peaks
#'   on the base of a minimal amplitude. To do so, TS are detrended by
#'   subtracting a rolling mean and amplitude to baseline is computed. Peak is
#'   kept if this amplitude is above the "filter" value. If NULL, no filter is
#'   applied. It is highly recommended to recenter the peaks if this option is
#'   used.
#' @param filter.winlen size of the window for moving average used to detrend
#'   prior to filtering. It is highly recommended to check the detrended signal
#'   (part of the output if filter is specified): in the detrended signal, area
#'   between peaks should oscillate around 0. If not, try to adjust the value.
#' @param detect.early numeric or NULL. If numeric, extrapolate the time series
#'   by the corresponding number of points at its start. This aims at improving
#'   peak detection at the beginning of the series. Extrapolation is performed
#'   by ARIMA model (see ?forecast::auto.arima)
#'
#' @details Zscore method is based on
#'   https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data#22640362.
#'
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
#'
#' @return A data table containing the signal along with results all of selected
#'   methods. Each method result is stored in a different column. If filter is
#'   specified, also returns the detrended signal.
#'
#' @seealso isolate_peak
#' @examples
#' # Curve coming from real-data:
#' y <- c(1.114, 1.146, 1.106, 1.049, 1.031, 1.01, 1, 0.992, 0.987, 0.981, 1.009, 1.072, 1.053, 1.018, 0.988,
#'  0.979, 0.965, 0.958, 0.956, 0.942, 0.939, 0.932, 0.926, 0.922, 0.919, 0.919, 0.916, 0.915, 0.918, 0.916,
#'  0.913, 0.911, 0.907, 0.908, 0.904, 0.905, 0.9, 0.899, 0.899, 0.905, 0.905, 0.908, 0.913, 0.919, 0.919,
#'  0.911, 0.922, 0.918, 0.916, 0.915, 0.92, 0.917, 0.915, 0.926, 0.964, 1.099, 1.075, 1.018, 0.999, 0.989,
#'  0.977, 0.975, 0.964, 0.96, 0.958, 0.945, 0.946, 0.944, 0.944, 0.938, 0.934, 0.936, 0.933, 0.926, 0.921,
#'  0.923, 0.924, 0.926, 0.929, 0.932, 0.926, 0.932, 0.939, 0.941, 0.976, 1.116, 1.09, 1.053, 1.034, 1.019,
#'  1.038, 1.104, 1.089, 1.059, 1.027, 1.012, 1.009, 1, 0.991, 0.986, 0.973, 0.969, 0.975, 0.966, 0.961,
#'  0.96, 0.949, 0.946, 0.943, 0.938, 0.942, 0.937, 0.934, 0.936, 0.941, 0.939, 0.97, 0.932, 0.931, 0.941,
#'  0.944, 0.949, 0.969, 1.048, 1.052, 1.059, 1.103, 1.083, 1.05, 1.036, 1.015, 1.01, 0.997, 0.995, 0.988,
#'  0.977, 0.979, 0.977, 0.981, 0.982, 0.976, 0.971, 0.97, 0.965, 0.966, 0.967, 0.967, 0.972, 0.904, 0.966,
#'  0.978, 0.983, 0.979, 0.975, 0.98, 0.973, 0.968, 0.97, 0.962, 0.966, 0.962, 0.968, 0.964, 0.956, 0.963,
#'  0.954, 0.958, 0.956, 0.958, 0.952, 0.95, 0.952, 0.946, 0.953, 0.952, 0.952, 0.953, 0.952, 0.954, 0.956)
#'
#' # Extract peaks with "vanilla" bump method. First smooth signal with a moving average of width 5, peaks must be
#' # separated with a minimal distance of 2 points.
#' peaks_vanilla <- extract_peak(y, method = "bump", bp.winlen = 5, bp.mind = 1, recenter = NULL, filter.thresh = NULL)
#' # Add: Recenter peaks to the maximum in a neigbourhood of 2 points.
#' peaks_center <- extract_peak(y, method = "bump", bp.winlen = 5, bp.mind = 1, recenter = 2, filter.thresh = NULL)
#' # Add: detrend the signal with a moving average of width 4, and filter peaks which amplitude is below 0.05 in the detrended signal.
#' peaks_center_filter <- extract_peak(y, method = "bump", bp.winlen = 5, bp.mind = 1, recenter = 2, filter.thresh = 0.05, filter.winlen = 10)
#' # Add: detection of early peak by extrapolating 10 points
#' peaks_final <- extract_peak(y, method = "bump", bp.winlen = 5, bp.mind = 1, recenter = 2, filter.thresh = 0.05, filter.winlen = 10, detect.early = 10)
#'
#' par(mfrow=c(3,2))
#' plot(y, type = "b", main = "Vanilla bump detection")
#' abline(v=which(peaks_vanilla$bump), col = "red", lty = "dashed")
#' plot(y, type = "b", main = "Peak recentering")
#' abline(v=which(peaks_center$bump), col = "red", lty = "dashed")
#' plot(y, type = "b", main = "Peak center + filter - Raw data")
#' abline(v=which(peaks_center_filter$bump), col = "red", lty = "dashed")
#' plot(peaks_center_filter$detrended.signal, type = "b", main = "Peak center + filter - Detrended data\nCheck that baseline is around 0")
#' abline(v=which(peaks_center_filter$bump), col = "red", lty = "dashed")
#' plot(y, type = "b", main = "Peak center + filter + early detection")
#' abline(v=which(peaks_final$bump), col = "red", lty = "dashed")
#' # Plot arima prediction
#' y_extend <- rev(c(rev(y), forecast(auto.arima(rev(y)), h = 10)$mean))
#' plot(y_extend, main = "Extended Series for early peak detection.\nBlue dots are predicted by ARIMA.", type = "b")
#' points(1:10, y_extend[1:10], pch=20, col="blue")
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
    
    # If detect.early is provided, calls the function recursively --> avoid double call
    if(is.null(detect.early)){
      
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
    # End of detect.early is null
    }
    
    # Extrapolation for early peaks
    else if(!is.null(detect.early)){
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
        #args[["show.warning"]] <- FALSE
        peaks <- do.call(extract_peak, args)
        return(peaks[(h+1):nrow(peaks)])
      }
      dt.out <- extrapolate_extract_peak(y, detect.early, argList)
    }
    
    return(dt.out)
  }




#' Isolate peaks
#'
#' Isolate peaks from a signal based on first (and second derivative). Starting
#' from the tip of the peak, walk left (then right) until derivatives overcome
#' thresholds. When this happens, mark a border of the peak. Individual peaks
#' are delimited by these borders. The rational of using the 2nd derivative is
#' to stop whenever the signal is steady (low 2nd derivative) and evolving slowly
#' (low 1st derivative).
#' @param y A numeric vector, from which to isolate peaks.
#' @param positions A logical of same length as y. Indicate the positions of the
#'   tip of the peaks. Can be obtained by running TSexploreR::extract_peak()
#' @param thresh.order1 Threshold for 1st derivative. Whenever the walk
#'   encounters a 1st derivative that is bigger than the threshold, stop and
#'   define a peak border. Default to 0, i.e. whenever the signal increases
#'   again, stop and define a border.
#' @param use.second logical, use 2nd derivative?
#' @param thresh.order22 If use.second = TRUE. If absolute second derivative is
#'   below this threshold, consider the signal as steady. If thresh.order12 is
#'   also below its threshold at the same time point, stop and define a border.
#' @param thresh.order12 If use.second = TRUE. If absolute first derivative is
#'   below this threshold, consider the signal as evolving very slowly. If
#'   thresh.order22 is also below its threshold at the same time point, stop and
#'   define a border.
#' @param One of c("data.table","list", "vector"). See Return.
#'
#' @return If out.type="data.table", returns a data.table with 2 columns. The
#'   first one is the signal, the second is a vector indicating to which peak
#'   the portion of the signal belongs to. In this vector, 0 means that the
#'   portion of the signal does not belong to any peak, while a number > 0
#'   indicate that it belongs to a peak. 1 indicates the first peak referred in
#'   'positions', 2 the second and so on. If for one of the peaks cannot be
#'   properly isolated, its number is skipped. If out.type="vector", returns
#'   only the second column of the data.table.
#'   
#'   Finally, if out.type="list",
#'   returns a list of same length as positions. Each element of the list
#'   posseses 4 elements: $start: starting point of the peak; $center: tip of the
#'   peak (max amplitude); $end: ending point of the peak; $y: part of the
#'   signal comprised between $start and $end.
#' @seealso extract_peak
#' @export
#'
#' @examples
#' # Curve coming from real-data:
#' y <- c(1.114, 1.146, 1.106, 1.049, 1.031, 1.01, 1, 0.992, 0.987, 0.981, 1.009, 1.072, 1.053, 1.018, 0.988,
#'  0.979, 0.965, 0.958, 0.956, 0.942, 0.939, 0.932, 0.926, 0.922, 0.919, 0.919, 0.916, 0.915, 0.918, 0.916,
#'  0.913, 0.911, 0.907, 0.908, 0.904, 0.905, 0.9, 0.899, 0.899, 0.905, 0.905, 0.908, 0.913, 0.919, 0.919,
#'  0.911, 0.922, 0.918, 0.916, 0.915, 0.92, 0.917, 0.915, 0.926, 0.964, 1.099, 1.075, 1.018, 0.999, 0.989,
#'  0.977, 0.975, 0.964, 0.96, 0.958, 0.945, 0.946, 0.944, 0.944, 0.938, 0.934, 0.936, 0.933, 0.926, 0.921,
#'  0.923, 0.924, 0.926, 0.929, 0.932, 0.926, 0.932, 0.939, 0.941, 0.976, 1.116, 1.09, 1.053, 1.034, 1.019,
#'  1.038, 1.104, 1.089, 1.059, 1.027, 1.012, 1.009, 1, 0.991, 0.986, 0.973, 0.969, 0.975, 0.966, 0.961,
#'  0.96, 0.949, 0.946, 0.943, 0.938, 0.942, 0.937, 0.934, 0.936, 0.941, 0.939, 0.97, 0.932, 0.931, 0.941,
#'  0.944, 0.949, 0.969, 1.048, 1.052, 1.059, 1.103, 1.083, 1.05, 1.036, 1.015, 1.01, 0.997, 0.995, 0.988,
#'  0.977, 0.979, 0.977, 0.981, 0.982, 0.976, 0.971, 0.97, 0.965, 0.966, 0.967, 0.967, 0.972, 0.904, 0.966,
#'  0.978, 0.983, 0.979, 0.975, 0.98, 0.973, 0.968, 0.97, 0.962, 0.966, 0.962, 0.968, 0.964, 0.956, 0.963,
#'  0.954, 0.958, 0.956, 0.958, 0.952, 0.95, 0.952, 0.946, 0.953, 0.952, 0.952, 0.953, 0.952, 0.954, 0.956)
#'  
#'  peaks_pos <- extract_peak(y, method = "bump", bp.winlen = 5, bp.mind = 1, recenter = 2, filter.thresh = 0.05, filter.winlen = 10, detect.early = 10)
#'  peaks <- isolate_peak(y=y, positions=which(peaks_pos$bump), thresh.order1=0, use.second=FALSE, out.type="list")
#'  plot(y, type = "b", main = "Points are colored according to their respective peaks.")
#'  abline(v=which(peaks_pos$bump), lty="dashed", col="red")
#'  col_ind <- 1
#'  for(peak in peaks){
#'    if(!any(is.na(c(peak$start, peak$end)))){
#'      points(peak$start:peak$end, y[peak$start:peak$end], pch=20, col=palette()[col_ind])
#'      col_ind <- col_ind+1
#'    }
#'  }
#'
isolate_peak <- function(y, positions, thresh.order1=0, use.second = FALSE, thresh.order22=NULL, thresh.order12=NULL, out.type="data.table"){
  # First and second derivatives, padded with NAs
  # Derivatives according to index, not time!
  d1 <- c(NA, diff(y, differences = 1))
  if(use.second){
    d2 <- abs(c(NA, NA, diff(y, differences = 2)))
    if(thresh.order22 < 0 | thresh.order12 < 0){
      thresh.order22 <- abs(thresh.order22)
      thresh.order12 <- abs(thresh.order12)
      warning("thresh.order22 and/or thresh.order12 was provided as a negative value,
              automatically turn it into its value is set to its absolute value")
    }
    }
  
  isolated_peaks <- list()
  j <- 1
  for(pos in positions){
    # Right walk: stop when d1 is positive (or above thresh) or when both acceleration and speed are low
    right <- NA
    for(i in (pos+1):length(y)){
      if(d1[i] >= thresh.order1){
        #print("stop right d1")
        right <- i - 1
        break
      }
      if(use.second){
        if((i >= pos+2) & (d2[i] <= thresh.order22) & (abs(d1[i]) <= thresh.order12)){
          #print("stop right d2")
          right <- i
          break
        }
      }
    }
    
    # Left walk: same stopping criteria as right walk
    left <- NA
    for(i in seq(pos, 1, -1)){
      # If peak is close to start of TS break or get error 
      if(is.na(d1[i])) break
      # Less or equal order 1 because reverse time
      if(d1[i] <= thresh.order1){
        #print("stop left d1")
        #print(pos)
        left <- i
        break
      }
      if(use.second){
        if((i <= pos-2) & (d2[i] <= thresh.order22) & (abs(d1[i]) <= thresh.order12)){
          #print("stop left d2")
          #print(pos)
          left <- i
          break
        }
      }
    }
    # If one of the borders cannot be defined, return NA
    if(is.na(left) | is.na(right)){
      isolated_peaks[[j]] <- list(start = left, center = pos, end = right, y = NA)
    }
    else{
      isolated_peaks[[j]] <- list(start = left, center = pos, end = right, y = y[left:right])
    }
    j <- j+1
  }
  
  # Output list, vector or data.table
  if(out.type == "list"){
    return(isolated_peaks)
  }
  else if(out.type == "vector" | out.type == "data.table"){
    # 0: no peak, X: belong to peak X
    vout <- rep(0, length(y))
    starts <- unlist(isolated_peaks)[names(unlist(isolated_peaks)) == "start"]
    ends <- unlist(isolated_peaks)[names(unlist(isolated_peaks)) == "end"]
    # if one of the borders are missing, don't consider and skip i
    for(i in 1:length(starts)){
      if(!any(is.na(c(starts[i], ends[i])))){
        vout[starts[i]:ends[i]] <- i
      }
    }
    if(out.type == "data.table"){
      return(data.table(signal=y, peak=vout))
    }
    else if(out.type == "vector") return(vout)
  }
}