#' synchrony.measures
#'
#' A set of measures for synchrony of (oscillating) trajectories.
#' @param x,y a numerical vector
#' @param window.size integer, size of the window used to compute
#' rolling mean. The latter is also used to get the trend component
#' in the classical decomposition.
#' @param method.clip one of c("rollmean","decomposition"). If "rollmean"
#' clipping is performed on x directly by comparing it to its
#' rolling mean. If "decomposition", x is decomposed via classical
#' decomposition (see ?classical.decomposition) and its seasonal
#' component is clipped according to its sign.
#' @param robust.decomp logical, whether classical decomposition should
#' be performed in a robust way. Default to TRUE.
#'
#' @details
#' Clipping method is robust to trends in the data, while correlations
#' are not. See example.
#' The clipping method aims at identifying ascending and
#' descending phase in signals and to see if these phases align
#' between two signals. Clipping transforms the signal in a binary
#' signal (0 or 1). An overlap of these clipped signals is then
#' returned. Rescaled between -1 (total asynchrony) and 1 (total
#' synchrony) for consistence with correlations.
#' @return A list of 4:
#' \itemize{
#'   \item $pearson: Pearson's correlation coefficient
#'   \item $spearman: Spearman's rank-based correlation coefficient
#'   \item $kendall: Kendall's rank-based correlation coefficient
#'   \item $overlap.clip: Overlap of clipped trajectories. Clipping
#'   performed according to method.clip.
#' }
#' @seealso cor, classical.decomposition
#' @export
#'
#' @examples
#' # Two phase-shifted sinusoids
#' x <- sin(seq(0, 20, 0.5))
#' y <- sin(seq(2, 22, 0.5))
#' # Number of points per oscillations
#' period <- get.period(x)
#' # Rollmean method
#' x.roll <- rollex(x, period)
#' y.roll <- rollex(y, period)
#' x.clip <- ifelse(x >= x.roll, 1, -1)
#' y.clip <- ifelse(y >= y.roll, 1, -1)
#' # Visualize
#' par(mfrow=c(2,1))
#' plot(x, type = "b", main = "Raw trajectory and rolling mean")
#' lines(x.roll, lty = "dashed", col = "darkgreen", lwd = 2)
#' lines(x.clip, type = "s", col = "blue")
#' plot(y, type = "b")
#' lines(y.roll, lty = "dashed", col = "darkgreen", lwd = 2)
#' lines(y.clip, type = "s", col = "red")
#' synchrony.measures(x, y, period, "rollmean")
#'
#' # With linear trend
#' x <- sin(seq(0, 20, 0.5))
#' y <- sin(seq(2, 22, 0.5))
#' trend <- seq(0,3,length.out=length(x))
#' x <- x + trend
#' y <- y + trend
#' # Rollmean method
#' x.roll <- rollex(x, period)
#' y.roll <- rollex(y, period)
#' x.clip <- ifelse(x >= x.roll, 3, 0)
#' y.clip <- ifelse(y >= y.roll, 3, 0)
#' # Visualize
#' par(mfrow=c(2,1))
#' plot(x, type = "b", main = "Raw trajectory and rolling mean,
#' with linear trend")
#' lines(x.roll, lty = "dashed", col = "darkgreen", lwd = 2)
#' lines(x.clip, type = "s", col = "blue")
#' plot(y, type = "b")
#' lines(y.roll, lty = "dashed", col = "darkgreen", lwd = 2)
#' lines(y.clip, type = "s", col = "red")
#' synchrony.measures(x, y, period, "rollmean")
#'
synchrony.measures <- function(x, y, window.size, method.clip, robust.decomp = TRUE){
  if(!method.clip %in% c("rollmean","decomposition")){
    stop("method.clip must be one of c('rollmean','decomposition')")
  }
  # Correlations
  cor.prs <- cor(x, y, method = "p")
  cor.spr <- cor(x, y, method = "s")
  cor.ken <- cor(x, y, method = "k")
  # Clip by comparing x to rolling mean
  if(method.clip=="rollmean"){
    x.roll <- rollex(x, window.size)
    y.roll <- rollex(y, window.size)
    x.clip <- ifelse(x >= x.roll, 1, 0)
    y.clip <- ifelse(y >= y.roll, 1, 0)
  # Clip by comparing seasonal component to 0
  } else if(method.clip=="decomposition"){
    x.decomp <- classical.decomposition(x, window.size, robust.decomp)
    y.decomp <- classical.decomposition(y, window.size, robust.decomp)
    x.clip <- ifelse(x.decomp[,"seasonal"] >= 0, 1, 0)
    y.clip <- ifelse(y.decomp[,"seasonal"] >= 0, 1, 0)
  }
  # Compute overlap between clipped trajectories. Scale on [-1,1]
  overlap <- sum(x.clip==y.clip)/length(x.clip)
  overlap <- 2*(overlap-1)+1
  return(list(pearson = cor.prs,
              spearman = cor.spr,
              kendall = cor.ken,
              overlap.clip = overlap))
}



#' max.cc
#'
#' Get maximal cross-correlation between 2 numerical vectors. Only 1st position
#' is returned.
#' @param x, y numerical vector
#' @param use.circular logical, whether to use circular cross-correlation.
#' @param plot logical. If TRUE, plots the cross-correlation histogram.
#' @param ... additional parameters for the cross-correlation if use.circular is FALSE.
#'  Of particular interest is lag.max. See ?ccf.
#'
#' @return A list of 2. $corr contains the maximum of correlation and
#' $lag the lag at which it was reached.
#' @export
#'
#' @examples
#' # 2 identical signals, shifted by a lag 3
#' pattern <- c(0.25, 0.5, 0.75, 1)
#' x <- rep(pattern, 10)
#' y <- c(x[4:length(x)],x[1:3])
#' par(mfrow=c(3,1))
#' plot(x, type = "h", main = "x", ylim = c(-0.2,1.2))
#' plot(y, type = "h", main = "y", ylim = c(-0.2,1.2))
#' # With circular cc, correlation is strictly 1
#' max_cc(x, y, use.circular = T, plot = T)
#'
max_cc <- function(x, y, use.circular = TRUE, plot = FALSE, ...){
  if(use.circular){
    cc <- circular.cc(x, y)
    if(plot){
      plot(names(cc), cc, type="h",
           xlab = "Lag", ylab = "Correlation")
      abline(h = 0, lty = "dashed")
    }
    max.indx <- which.max(cc)
    return(list(corr=cc[max.indx], lag=as.numeric(names(cc[max.indx]))))
  }else{
    cc <- ccf(x, y, plot = plot, ...)
    max.indx <- which.max(cc$acf)
    return(list(corr=cc$acf[max.indx], lag=cc$lag[max.indx]))
  }
}










all_pairwise_stats <- function(data, condition, label, measure, k_roll_mean = 5){
  # Compute all pairwise measures: DTW, overlap_clipping and correlations (Pearson, Spearman, Kendall)
  # The clipping of trajectories is performed as first step

  # data: a data table containing CLIPPED trajectories
  # condition: column name that fully define an experimental condition; MUST BE INTEGERS, NUMERIC, FACTOR OR CHARACTERS
  # measure: column name with clipped trajectories
  # label: column name with label of individual objects in each condition (cell label); LABELS MUST BE INTEGERS, NUMERIC, FACTOR OR CHARACTERS

  if(!(class(data[,get(condition)]) %in% c("integer", "numeric", "character", "factor"))){
    stop("Column 'condition' must be integer, numeric, factor or character.")
  }
  if(!(class(data[,get(label)]) %in% c("integer", "numeric", "character", "factor"))){
    stop("Column 'label' must be integer, numeric, factor or character.")
  }

  require(data.table)
  require(dtw)
  setkeyv(data, c(condition, label))

  # Perform clipping
  ClipData <- data[, .(clip_measure = wrap_clip(get(measure), k = k_roll_mean)), by = c(condition, label)]
  setkeyv(ClipData, c(condition, label))

  # Number of rows for data table initialization, sum of number of pairs in each condition
  nber_row <- 0
  for(i in unique(data[, get(condition)])){
    nber_row <- nber_row + choose(length(unique(data[.(i), get(label)])), 2)
  }


  # One row = one pair in one condition; set column types according to input types
  out <- data.table(matrix(ncol = 8, nrow = nber_row))
  colnames(out) <- c(condition, "Label1", "Label2", "Overlap", "Pearson", "Spearman", "Kendall", "DTW")
  if(class(data[,get(condition)]) == "integer"){out[[condition]] <- as.integer(out[[condition]])}
  else if(class(data[,get(condition)]) == "numeric"){out[[condition]] <- as.numeric(out[[condition]])}
  else if(class(data[,get(condition)]) == "character"){out[[condition]] <- as.character(out[[condition]])}
  else if(class(data[,get(condition)]) == "factor"){out[[condition]] <- as.factor(out[[condition]])}


  if(class(data[,get(label)]) == "integer"){
    out[, Label1 := as.integer(Label1)]
    out[, Label2 := as.integer(Label2)]
  }
  else if(class(data[,get(label)]) == "numeric"){
    out[, Label1 := as.numeric(Label1)]
    out[, Label2 := as.numeric(Label2)]
  }
  else if(class(data[,get(label)]) == "character"){
    out[, Label1 := as.character(Label1)]
    out[, Label2 := as.character(Label2)]
  }
  else if(class(data[,get(label)]) == "factor"){
    out[, Label1 := as.factor(Label1)]
    out[, Label2 := as.factor(Label2)]
  }

  out[, Overlap := as.numeric(Overlap)]
  out[, Pearson := as.numeric(Pearson)]
  out[, Spearman := as.numeric(Spearman)]
  out[, Kendall := as.numeric(Kendall)]
  out[, DTW := as.numeric(DTW)]

  curr_row <- 1L
  # Loop condition
  for(i in unique(data[, get(condition)])){
    labels <- unique(data[.(i), get(label)])
    # Loop 1st label
    for(j in 1:(length(labels)-1)){
      # Loop 2nd label
      for(k in (j+1):length(labels)){
        set(out, curr_row, 1L, i)
        set(out, curr_row, 2L, labels[j])
        set(out, curr_row, 3L, labels[k])
        # /!\ use clip data for overlap
        set(out, curr_row, 4L, overlap(ClipData[.(i, labels[j]), clip_measure], ClipData[.(i, labels[k]), clip_measure]))
        set(out, curr_row, 5L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)], method ="pearson"))
        set(out, curr_row, 6L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)], method ="spearman"))
        set(out, curr_row, 7L, cor(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)], method ="kendall"))
        set(out, curr_row, 8L, dtw(data[.(i, labels[j]), get(measure)], data[.(i, labels[k]), get(measure)])$distance)
        curr_row <- curr_row + 1L
      }
    }
  }
  return(out)
}


# # Solution to try
# library(data.table)
# dat <- data.table(Region = rep(LETTERS[24:26], each=3),
#                   Name = rep(LETTERS[1:3], 3),
#                   Val1 = rep(rnorm(3), 3),
#                   Val2 = rep(rnorm(3), 3))
# dat2 <- merge(dat, dat, by="Region", allow.cartesian = T)[Name.x < Name.y]
# dat2[, Val1Ratio := Val1.x / Val1.y]
# dat2[, Val2Ratio := Val2.x / Val2.y]
# dat2[, Diff := Val1Ratio - Val2Ratio]
