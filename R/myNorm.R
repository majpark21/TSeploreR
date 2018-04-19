#' myNorm
#'
#' Returns original dt with an additional column with normalized quantity.
#' Normalisation is based on part of the trajectory.
#' @param in.dt data.table, contain a column for: time, measurement(, grouping
#'   columns like object label)
#' @param in.meas.col character, name of column to be normalised.
#' @param in.rt.col character, name of the column with time (or index)
#' @param in.rt.min numeric, start time of normalization period.
#' @param in.rt.max numeric, end time of normalization period.
#' @param in.by.cols character vector with 'by' columns to calculate
#'   normalisation per group if NULL, no grouping is done
#' @param in.robust logical, whether robust measures should be used (median
#'   instead of mean, mad instead of sd).
#' @param in.type type of normalization: 'z.score' or 'mean' (fi.e. old change
#'   w.r.t. mean)
#'
#' @return A data.table with an additional column with normalized quantity. The
#'   name of additional column is the same as in.meas.col but with ".norm"
#'   suffix added.
#' @export
#'
myNorm <- function(in.dt,
                  in.meas.col,
                  in.rt.col = 'RealTime',
                  in.rt.min = 10,
                  in.rt.max = 20,
                  in.by.cols = NULL,
                  in.robust = TRUE,
                  in.type = 'z.score') {
  loc.dt <-
    copy(in.dt) # copy so as not to alter original dt object w intermediate assignments

  if (is.null(in.by.cols)) {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = mad(get(in.meas.col), na.rm = TRUE))]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = sd(get(in.meas.col), na.rm = TRUE))]

    loc.dt = cbind(loc.dt, loc.dt.pre.aggr)
  }  else {
    if (in.robust)
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = median(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = mad(get(in.meas.col), na.rm = TRUE)), by = in.by.cols]
    else
      loc.dt.pre.aggr = loc.dt[get(in.rt.col) > in.rt.min &
                                 get(in.rt.col) < in.rt.max, .(meas.md = mean(get(in.meas.col), na.rm = TRUE),
                                                               meas.mad = sd(get(in.meas.col), na.rm = TRUE)), by = in.by.cols]

    loc.dt = merge(loc.dt, loc.dt.pre.aggr, by = in.by.cols)
  }


  if (in.type == 'z.score') {
    loc.dt[, meas.norm := (get(in.meas.col) - meas.md) / meas.mad]
  } else {
    loc.dt[, meas.norm := (get(in.meas.col) / meas.md)]
  }

  setnames(loc.dt, 'meas.norm', paste0(in.meas.col, '.norm'))

  loc.dt[, c('meas.md', 'meas.mad') := NULL]
  return(loc.dt)
}
