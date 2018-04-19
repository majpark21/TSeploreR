###########################################################################
# Functions to extract features from single peak signal                   #
# All functions are build to work on a simple numerical vector, named y   #
###########################################################################


#### Extract All Features ####

#' FeatAllFeat
#'
#' Extract all single-peak features
#' @param y a numerical vector from which features are extracted.
#' @param basal a normalization constant to bring base of the peak around 0.
#' @param start.lag.grow Index of y where ascending phase starts. See
#'   FeatLagGrow, argument 'start'.
#' @param end.exp.dec Index of y where descending phase ends. See FeatExpDec,
#'   argument 'end'.
#' @param ... Additional arguments for FeatFWHM.
#'
#' @return A list of 13 features: \itemize{ \item mini: minimum of y \item maxi:
#'   maximum of y \item diff.min.max: difference maxi-mini \item time.min: index
#'   of y where it is minimized \item time.max: index of y where it is maximized
#'   \item max.amp: amplitude of peak. Obtained by subtracting constant basal.
#'   See FeatMaxAmplitude. \item FWHM: Full-Width at Half maximum. Difference
#'   right-left. See FeatFWHM. \item left: time where half maximum of the peak
#'   is reached on its left flank. \item right: time where half maximum of the
#'   peak is reached on its right flank. \item grow.half,grow.lag: growth rate
#'   estimated by the slope of a linear regression of the ascending phase of the
#'   peak. 2 different methods are used to isolate this phase. First, it is
#'   defined as the time between which the signal has reached half of its
#'   maximum (left) and its maximum; see FeatHalfMaxGrow. Second, as the time
#'   between a hard-encoded starting time point (start.lag.grow) and peak
#'   maximum; see FeatLagGrow. \item dec.half,dec.exp: decay rate estimated by
#'   the slope of a linear regression of the descending phase of the peak. 2
#'   different methods are used to isolate this phase. First, it is defined as
#'   the time between which the signal has reached its maximum and half of its
#'   maximum (right); see FeatHalfMaxDec. Second, as the time between peak
#'   maximum and a hard-encoded endtime point (start.lag.grow); see FeatExpDec.
#'   In this case the signal is modelled by EXPONENTIAL DECAY. }
#'
#' @export
#'
#' @examples
#' signal <- dnorm(x = seq(-4,4, length.out = 500), mean = 0, sd = 1)
#' feats <- FeatAllFeat(y = signal, basal = 0, start.lag.grow = 100, end.exp.dec = 400)
#' plot(signal)
#' # Half-maximum points
#' abline(h = max(signal)/2, col = "green")
#' abline(v = c(feats$left, feats$right), col = c("blue", "red"))
#'
#' # Growth and Decay slopes according to half-max methods
#' growth0 = c(feats$left, max(signal)/2)
#' growth1 = c(which.max(signal), max(signal)/2 + feats$grow.half * (which.max(signal)-feats$left))
#' segments(x0 = growth0[1], y0 = growth0[2], x1 = growth1[1], y1 = growth1[2], col = "blue", lwd = 4, lty = "dashed")
#'
#' decay0 = c(which.max(signal), max(signal))
#' decay1 = c(feats$right, max(signal) + feats$dec.half * (feats$right-which.max(signal)))
#' segments(x0 = growth0[1], y0 = growth0[2], x1 = growth1[1], y1 = growth1[2], col = "blue", lwd = 4, lty = "dashed")
#'
FeatAllFeat <- function(y, basal, start.lag.grow, end.exp.dec, ...){
  minmax <- FeatDiffMinMax(y)
  amplitude <- FeatMaxAmplitude(y, basal = basal)
  fwhm <- FeatFWHM(y, basal = basal, ...)
  grow.half.max <- FeatHalfMaxGrow(y, basal = basal, ...)
  grow.lag <- FeatLagGrow(y, start = start.lag.grow)
  dec.half.max <- FeatHalfMaxDec(y, basal = basal, ...)
  dec.exp <- FeatExpDec(y, end = end.exp.dec)
  return(list(mini = minmax$min,
              maxi = minmax$max,
              diff.min.max = minmax$diff.min.max,
              time.min = minmax$time.min,
              time.max = minmax$time.max,
              max.amp = amplitude$max,
              FWHM = fwhm$fwhm,
              left = fwhm$left,
              right = fwhm$right,
              grow.half = grow.half.max,
              grow.lag = grow.lag,
              dec.half = dec.half.max,
              dec.exp = dec.exp))
}


#### Amplitude Features ####

#' FeatDiffMinMax
#'
#' Difference between minimum and maximum value.
#' @param y a numerical vector
#'
#' @return List of 5:
#' \itemize{
#'   \item min: minimum of y
#'   \item max: maximum of y
#'   \item diff.min.max: $max - $min
#'   \item time.min: first index at which $min is reached in y
#'   \item time.max: first index at which $max is reached in y
#' }
#' @export
#'
FeatDiffMinMax <- function(y){
  mini <- min(y)
  maxi <- max(y)
  tmini <- which.min(y)
  tmaxi <- which.max(y)
  return(list(min = mini, max = maxi, diff.min.max = maxi - mini, time.min = tmini, time.max = tmaxi))
}


#' FeatMaxAmplitude
#'
#' Subtracts a constant value to a vector and return the maximum value of the
#' difference.
#' @param y a numerical vector
#' @param basal value to be subtracted for the signal to obtain height of the
#'   peak as amplitude instead of absolute value
#'
#' @details Interest over FetDiffMinMax, is that min(y) might be an outlier, or
#'   not really at the base of the peak. If basal = min(y); the result is
#'   strictly equivalent to FeatDiffMinMax.
#' @return list of 2: \itemize{ \item max: maximum value of y after subtracting
#'   basal from it. \item time.max: index at which $max is reached. }
#' @export
#'
FeatMaxAmplitude <- function(y, basal = min(y)){
  # Shift the data to have basal at 0, this is important so that maximum peak is
  # measured in amplitude instead of absolute value
  y <- y - basal
  return(list(max = max(y), time.max = which.max(y)))
}


#### FWHM Features ####

#' Find Full Width at Half maximum by spline interpolation
#'
#' @param x numeric, typically time vector
#' @param y numeric, typically a signal recorded over time at x.
#' @param n numeric, number of points returned by interpolation. Increase
#'   improves accuracy of FWHM estimation.
#' @param method one of "minimum" or "walk". Walk method is highly recommended
#'   as it avoids many pitfalls of the minimum method. Minimum simply pick the
#'   points which values are the closest to the half of the max value. Walk
#'   method walks left and right from the maximum point, and stops whenever half
#'   max value is crossed.
#' @param basal Normalization constant to bring base of the signal at 0.
#'
#' @return A list of 3: left and right represent the value of interpolated x at
#'   which half of the maximum in y is reached. fwhm is the difference between
#'   right and left.
#'
FeatFWHM <- function(y, x = seq_along(y), n = 30*length(x), method = "walk", basal = min(y)){
  if(!method %in% c("walk", "minimum")) stop("Method must be one of c('walk','minimum')")
  # Shift the data to have basal at 0, this is important so that maximum peak is
  # measured in amplitude instead of absolute value
  y <- y - basal
  # Do NOT use interpolated data as maximum estimation
  maxi <- max(y)
  tmaxi <- x[which.max(y)]
  # Interpolate data by spline
  fit <- mySpline(x, y, n)
  if(method == "walk"){
    thresh <- maxi/2
    start <- which(fit$x == tmaxi)[1]
    # Right walk
    right <- NA
    for(i in start:length(fit$x)){
      if(fit$y[i] < thresh){
        # Once threshold is passed, choose the value between i and i-1 that is the closest to it
        right <- ifelse(abs(thresh - fit$y[i-1]) < abs(thresh - fit$y[i]), fit$x[i-1], fit$x[i])
        break
      }
    }
    # Left walk
    left <- NA
    for(i in seq(start, 1, -1)){
      if(fit$y[i] < thresh){
        left <- ifelse(abs(thresh - fit$y[i+1]) < abs(thresh - fit$y[i]), fit$x[i+1], fit$x[i])
        break
      }
    }
  } else if(method == "minimum"){
    left <- fit$x[which.min(abs(fit$y[fit$x < tmaxi] - (maxi/2)))]
    right <- fit$x[which.max(which(fit$x <= tmaxi)) + which.min(abs(fit$y[fit$x > tmaxi] - (maxi/2)))]
  }
  if(length(right) == 0){right <- NA; fwhm <- NA}
  if(length(left) == 0){left <- NA; fwhm <- NA}
  return(list(fwhm = right - left, left = left, right = right))
}


#### Decay Features ####

#' FeatHalfMaxDec
#'
#' Estimate decay rate of the peak.
#' Fit a linear model on descending phase of a peak between (possibly interpolated) time at which maximum of the peak is reached
#' and half-max of the peak is reached. The fit is enforced to pass through the peak maximum.
#' @param y numeric, typically measured variable.
#' @param ... extra arguments for FeatFWHM. Notably 'basal' for setting base of the peak around 0.
#'
#' @return Slope of the linear regression, estimator of the decay rate.
#' @seealso FeatFWHM, FeatExpDec
#' @export
#'
FeatHalfMaxDec <- function(y, ...){
  # Enforces fit to pass through maximum
  right <- FeatFWHM(y = y, ...)$right
  if(length(right) == 0) return(NA)
  if(is.na(right)) return(NA)
  tmaxi <- which.max(y)
  # Use the point closest to half-max as additional point for regression (if it is interpolated then it is not an integer)
  if(right %% 1 !=0){
    splinetraj <- mySpline(x = seq_along(y), y, 30*length(y))
    newpoint <- splinetraj$y[splinetraj$x == right]
    ytrim <- c(y[tmaxi:right], newpoint) - max(y)
    t.after.max <- seq_along(ytrim) - 1
    # Replace last time by the time of interpolation
    t.after.max[length(t.after.max)] <- right - tmaxi
  } else {
    ytrim <- y[tmaxi:right] - max(y)
    t.after.max <- seq_along(ytrim) - 1
  }
  fit <- lm(ytrim ~ t.after.max + 0)
  return(coef(fit)[1])
}


#' FeatExpDec
#'
#' Estimate decay rate of the peak. Fit a linear model on the descending phase,
#' using exponential decay model. The fit is performed betwenn peak maximum and
#' til the index specified by end and is enforced to pass through peak maximum.
#' @param y numeric, typically measured variable!
#' @param end int, at which index should the regression end.
#'
#' @return slope of the linear regression, estimator of the decay rate.
#' @seealso FeatHalfMaxDec
#' @export
#'
FeatExpDec <- function(y,end){
  logy <- log(y) - log(max(y))
  ytrim <- logy[which.max(y):end]
  time.reg <- seq_along(ytrim) - 1
  fit <- lm(ytrim ~ time.reg + 0)
  return(coef(fit)[1])
}


#### Growth Features ####

#' FeatLagGrow
#'
#' Estimate growth rate of the peak. Fit a linear model on the growing phase,
#' starting at specified time, til max of the peak.
#' @param y numeric, typically measured variable
#' @param start int, at which index should the regression start.
#'
#' @return Slope of the linear regression, estimator of the growth rate.
#' @export
#'
FeatLagGrow <- function(y, start){
  tmaxi <- which.max(y)
  ytrim <- y[start:tmaxi]
  time.reg <- seq_along(ytrim) -1
  fit <- lm(ytrim ~ time.reg)
  return(coef(fit)[2])
}


#' FeatHalfMaxGrow
#'
#' Estimate growth rate of the peak. Fit a linear model on growing phase of a
#' peak between (possibly interpolated) time at which half-max of the peak is
#' reached and maximum of the peak
#' @param y numeric, typically measured variable.
#' @param ... extra arguments for FeatFWHM. Notably 'basal' for setting base of
#'   the peak around 0.
#'
#' @return Slope of the linear regression, estimator of the growth rate.
#' @seealso FeatFWHM, FeatLagGrow
#' @export
#'
FeatHalfMaxGrow <- function(y, ...){
  left <- FeatFWHM(y = y, ...)$left
  if(length(left) == 0) return(NA)
  if(is.na(left)) return(NA)
  tmaxi <- which.max(y)
  # Use the point closest to half-max as additional point for regression (if it is interpolated then it is not an integer)
  if(left %% 1 !=0){
    splinetraj <- mySpline(x = seq_along(y), y, 30*length(y))
    newpoint <- splinetraj$y[splinetraj$x == left]
    ytrim <- c(newpoint, y[(left+1):(tmaxi+1)]) - max(y)
    # Set interpolation time to 0, and shift all other times
    time.reg <- c(0,   seq(1, (length(ytrim)-1))-(1-left%%1) )
  } else {
    ytrim <- y[left:tmaxi] - max(y)
    time.reg <- seq_along(ytrim) - 1
  }
  fit <- lm(ytrim ~ time.reg)
  return(coef(fit)[2])
}

