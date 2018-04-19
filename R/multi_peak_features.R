################
# All features #
################

#' MPFeatAllFeat
#'
#' Compute a set of features aimed at characterizing multipeak/oscillating
#' signals.
#' @param x a numerical vector.
#' @param window.extrema integer, window width for detecting local extrema.
#' Used for detecting peaks maxima or minima. See ?MPFeatTrend_extrema.
#' @param window.rollmean integer, window width for classical decomposition.
#' In case of multi-peak signal, should correspond to the number of
#' points in one period. See ?classical.decomposition.
#' @param extrema one of c("mini", "maxi"). Which extrema should be
#' looked for and regressed on for trend estimation?
#' See ?MPFeatTrend_extrema.
#' @param trim.remain numeric, what is the minimum absolute remainder to
#' be considered for computing amplitude variation? See ?MPFeatAmp_variation.
#' @param robust.linreg logical, whether linear regression for trend estimation
#' should be performed robustly. See ?mblm::mblm
#' @param robust.decomp logical, whether classical decomposition for amplitude
#' features should be performed robustly. See ?classical.decomposition.
#'
#' @return Alist of 5 features of x:
#' \itemize{
#'   \item $trend.extrema: slope of linear regression performed on
#'   local extrema. Estimator for long-term trend.
#'   \item $trend.rollmean: slope of linear regression performed on
#'   rolling mean. Estimator for long-term trend.
#'   \item $amplt.dist.euclid: Euc lidean distance between x and
#'   its rolling mean. Not scaled! Estimator for oscillation amplitude.
#'   \item $amplt.season: amplitude (max-min) of the seasonal
#'   component. Estimator for oscillation amplitude.
#'   \item $amplt.variation: mean remainder of classical decomposition,
#'   after optionally trimming of small remainders.
#'   Estimator for oscillation amplitude variation.
#' }
#' @export
#'
#' @examples
#' # One sinusoid with linear trend
#' x <- sim_phase_shifted_with_fixed_trend(n = 1, noise = 0, slope = 0.1)
#' # How many points in one oscillation? Check maximum in power spectrum
#' x.spec <- spectrum(x$value)
#' x.period <- round(1/x.spec$freq[which.max(x.spec$spec)])
#' # Check by plotting extended rolling mean
#' plot(x$value, type = "b")
#' lines(seq_along(x$value), rollex(x$value, 31), col = "blue")
#' legend("topleft", legend="rolling mean", fill="blue")
#' # Extract features
#' x.features <- MPFeatAllFeat(x$value, x.period, 5)
#' # Visualize trend
#' plot(x$value, type = "b")
#' abline(a=0, b=x.features$trend.extrema, col = "blue")
#' abline(a=0, b=x.features$trend.rollmean, col = "red")
#' legend("topleft", legend=c("extrema", "rollmean"), fill=c("blue","red"),
#' title = "Trend estimation method")
#' # Visualize amplitude season
#' plot_decomposition(classical.decomposition(x$value, x.period))
#'
MPFeatAllFeat <- function(x, window.rollmean, window.extrema,
                          extrema = "maxi", trim.remain = NULL,
                          robust.linreg = FALSE, robust.decomp = TRUE){
  trend_extr <- MPFeatTrend_extrema(x, window.extrema, extrema, robust.linreg)
  trend_roll <- MPFeatTrend_rollmean(x, window.rollmean, robust.linreg)
  amp_euclid <- MPFeatAmp_euclidmean(x, window.rollmean)
  amp_season <- MPFeatAmp_seasonal(x, window.rollmean, robust.decomp)
  amp_vartin <- MPFeatAmp_variation(x, window.rollmean, trim.remain, robust.decomp)
  return(list(trend.extrema = trend_extr$trend,
              trend.rollmean = trend_roll$trend,
              amplt.dist.euclid = amp_euclid,
              amplt.season = amp_season,
              amplt.variation = amp_vartin))
}



#####################
# Analysis of trend #
#####################

#' MPFeatTrend_extrema
#'
#' Isolate local extrema and perform a linear regression.
#' Aims at giving an estimation of long-term trend in multi-peak signals.
#'
#' @param x a numerical vector.
#' @param window.size integer, width of window for maxima detection.
#'  See ?detect.peak.
#' @param what character, one of "minima" or "maxima". Define whether regression
#' should be performed on maxima or minima of peaks.
#' @param robust logical. If TRUE perform robusts linear regression.
#' Median-Based Linear Models are used, see ?mblm::mblm
#'
#' @return A list of 4.
#' \itemize{
#'  \item $trend, slope of the linear model
#'  \item $model, the full linear model
#'  \item $extremes, x extrema values and their index
#'  \item $type, type of extrema, "mini" or "maxi"
#'  }
#' @export
#' @seealso MPFeatTrend_rollmean
#'
#' @examples
#' x <- sim_phase_shifted_with_fixed_trend(n = 1, noise = 0, slope = 0.1)
#' x.trend <- MPFeatTrend_extrema(x = x$value, window.size = 5, what = "maxi")
#' plot(x$value, type = "b")
#' abline(x.trend$model, col = "blue", lwd = 2, lty = "dashed")
#'
MPFeatTrend_extrema <- function(x, window.size, what, robust = FALSE){
  require(mblm)
  # Extract local extrema
  extr.pos <- which(detect.peak(x, window.size, what))
  extr.x <- x[extr.pos]
  names(extr.x) <- extr.pos

  # Perform fitting
  if(robust){
    fit <- mblm(extr.x ~ extr.pos)
    return(list(trend=fit$coefficients[2], model=fit, extremes=extr.x, type=what))
  } else {
    fit <- lm(extr.x ~ extr.pos)
    return(list(trend=fit$coefficients[2], model=fit, extremes=extr.x, type=what))
  }
}


#' MPFeatTrend_rollmean
#'
#' Smooth by rolling mean and perform linear regression.
#' Rolling mean is extended at extremeties by linear extrapolation, see ?rollex.
#' Aims at giving an estimation of long-term trend in multi-peak signals.
#'
#' @param x a numerical vector.
#' @param window.size integer, width of window for rolling mean.
#' @param robust logical. If TRUE perform robusts linear regression.
#' Median-Based Linear Models are used, see ?mblm::mblm
#'
#' @return A list of 2: $model, the full linear model;
#' $trend, slope of the linear model.
#' @export
#' @seealso MPFeatTrend_extrema
#' @examples
#' width <- 30
#' x <- sim_phase_shifted_with_fixed_trend(n = 1, noise = 0, slope = 0.1)
#' # Rolling mean extended by linear extrapolation.
#' # This is also performed in MPFeatTrend_rollmean.
#' x.roll <- rollex(x$value, width)
#' plot(x$value, type = "b")
#' lines(x.roll, col ="red", lwd = 2, lty = "dashed")
#' x.trend <- MPFeatTrend_rollmean(x = x$value, window.size = width)
#' abline(x.trend$model, col = "blue", lwd = 2, lty = "dashed")
#'
MPFeatTrend_rollmean <- function(x, window.size, robust = FALSE){
  require(mblm)
  x.roll <- rollex(x, window.size)
  if(robust){
    fit <- mblm(x ~ seq_along(x))
    return(list(trend=fit$coefficients[2], model = fit))
  } else {
    fit <- lm(x ~ seq_along(x))
    return(list(trend=fit$coefficients[2], model = fit))
  }
}


############################
# Analysis of oscillations #
############################

#' MPFeatAmp_euclidmean
#'
#' Compute Euclidean distance between a signal and its rolling mean.
#' Rolling mean is extended at extremeties by linear extrapolation, see ?rollex.
#' Aims at giving an estimation of peak amplitude in oscillating signal.
#' @param x a numerical vector.
#' @param window.size integer, width of window for rolling mean.
#'
#' @return The Euclidean distance between x and its rolling mean.
#' @export
#' @seealso MPFeatAmp_seasonal, MPFeatAmp_variation
#'
MPFeatAmp_euclidmean <- function(x, window.size){
  return(sqrt(sum( (x - rollex(x, window.size) )^2 )))
}


#' MPFeatAmp_seasonal
#'
#' Decompose signal with classical decomposition and returns amplitude of the
#' seasonal component.
#' Aims at giving an estimation of peak amplitude in oscillating signal.
#' @param x a numerical vector.
#' @param window.size integer, width of window for rolling mean used to extract
#' the trend component. See ?classical.decomposition.
#' @param robust logical. If TRUE, seasonal component is estimated by taking
#' the median of values at each stage of a season. If FALSE, uses mean.
#' See ?classical.decomposition.
#'
#' @return The difference between maximum and minimum of the seasonal component.
#' @export
#' @seealso MPFeatAmp_euclidmean, MPFeatAmp_variation
#'
MPFeatAmp_seasonal <- function(x, window.size, robust = TRUE){
  x.decomp <- classical.decomposition(x, window.size, robust = robust)
  return(max(x.decomp[,"seasonal"]) - min(x.decomp[,"seasonal"]))
}


#' MPFeatAmp_variation
#'
#' Decompose signal with classical decomposition and returns mean of the
#' (trimmed) absolute remainder component.
#' Aims at at giving an estimation of the variation of peak amplitudes
#' in oscillating signal.
#' @param x a numerical vector.
#' @param window.size integer, width of window for rolling mean used to extract
#' the trend component. See ?classical.decomposition.
#' @param trim optional. Gives the minimum absolute value of remainder to be
#' considered.
#' @param robust logical. If TRUE, seasonal component is estimated by taking
#' the median of values at each stage of a season. If FALSE, uses mean.
#' See ?classical.decomposition.
#'
#' @return The mean (trimmed) absolute remainder.
#' @export
#' @seealso MPFeatAmp_euclidmean, MPFeatAmp_seasonal
#'
MPFeatAmp_variation <- function(x, window.size, trim = NULL, robust = TRUE){
  x.decomp <- classical.decomposition(x, window.size, robust = robust)
  x.remain <- x.decomp[,"remainder"]
  x.remain <- abs(x.remain)
  if(!is.null(trim)) x.remain <- x.remain[x.remain >= trim]
  return(mean(x.remain))
}
