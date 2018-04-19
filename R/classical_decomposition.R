#' Decompose a time series into trend, seasonal and remainder component
#'
#' The model is assumed to be additive, i.e. the time series is of the form:
#' y(t) = trend(t) + seasonal(t) + remainder(t)
#'
#' Trends is determined using a rolling mean with window width corresponding to
#' the signal period, extremities are padded with linear interpolation from the
#' first and last 2 measures.
#'
#' Season is determined by summarizing points "period-wise". Specifically,
#' let's consider a signal, which has 5 periods, each composed of 15 points.
#' The 1st point of a season is determined by taking the mean (or median if
#' robust=TRUE) of the 1st point of each of the 5 periods. Once all 15 points
#' have been summarized, this 'unit season' is repeated 5 times to form the
#' seasonal component.
#'
#' @param ts: a time series or numerical vector
#' @param period: number of time points in one cycle
#' @param robust: if T, the median of all points of a cycle are used to
#'   determined seasonal component; if F, the mean is used instead.
#'
#' @return a 3xn numeric matrix, with column trend, season and remainder
#' @note Good explanation of the decomposition at https://www.otexts.org/fpp/6/3
#' @seealso decompose
#' @examples
#' # Regular motif of length 15 repeated 5 times + Linear Trend
#' x <- rep(1:15, 5)
#' x <- x + seq(0, 10, length.out = length(x))
#' x_decomp <- classical.decomposition(ts = x, period = 15)
#' plot_decomposition(x_decomp)
#' @export
#'
classical.decomposition <- function(ts, period, robust = T){
  # 1) Get trend with extended rolling mean
  trend <- rollex(ts, period)
  # 2) Detrend
  detrend <- ts - trend
  # 3) Average every time point for each oscillation
  seasonal <- rep(NA, length(ts))
  if(robust){
    for(i in 1:period){
      index <- seq(i, length(ts), by = period)
      seasonal[index] <- median(detrend[index])
    }
  } else {
    for(i in 1:period){
      index <- seq(i, length(ts), by = period)
      seasonal[index] <- mean(detrend[index])
    }
  }
  # 4) Adjust to zero
  seasonal <- scale(seasonal, center = T, scale = F)
  # 5) Remainder
  remainder <- ts - (trend + seasonal)
  out <- cbind(ts, trend, seasonal, remainder)
  colnames(out) <- c("ts","trend", "seasonal", "remainder")
  return(out)
}

#' @describeIn classical.decomposition utility plot function
#' @export
plot_decomposition <- function(decomposition, main = "Time series decomposition"){
  xax <- seq_len(nrow(decomposition))
  layout(matrix(1:4, ncol = 1), widths = 1, heights = c(1.6,1,1,1.5), respect = F)
  par(mar = c(0, 5.1, 4.1, 2.1))
  plot(xax, decomposition[,"ts"], main = main, type = "l", xaxt = 'n', xlab = "", ylab = "Data", cex.lab = 1.4, cex.main = 1.5)
  par(mar = c(0, 5.1, 0, 2.1))
  plot(xax, decomposition[,"seasonal"], type = "l", xaxt = 'n', xlab = "", ylab = "Seasonal", cex.lab = 1.4)
  par(mar = c(0, 5.1, 0, 2.1))
  plot(xax, decomposition[,"trend"], type = "l", xaxt = 'n', xlab = "", ylab = "Trend", cex.lab = 1.4)
  par(mar = c(4.1, 5.1, 0, 2.1))
  plot(xax, decomposition[,"remainder"], type = "h", ylab = "Remainder", xlab = "Time", cex.lab = 1.4)
  abline(h=0, lty="dashed")
}

