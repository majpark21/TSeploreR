#' separability.measures
#'
#' Compute a set of measures for separability of 2 distributions. from:
#' https://stats.stackexchange.com/a/78855/43610
#' @param Vector.1,Vector.2  numerical vectors
#'
#' @return A vector of 5 measures: \itemize{ \item "jm": Jeffries-Matusita index.
#'   Bounded version of Bhattacharrya on (0,2). 0: identical dsitributions; 2:
#'   perfectly separable distributions. \item "bh": Bhattacharrya index Vary on
#'   (0, +inf). 0: identical distributions. \item "div": Kullback-Leibler
#'   Divergence. \item "tdiv": Transformed Divergence }
#' @export
#'
separability.measures <- function ( Vector.1 , Vector.2 ) {
  # convert vectors to matrices in case they are not
  Matrix.1 <- as.matrix (Vector.1)
  Matrix.2 <- as.matrix (Vector.2)
  # define means
  mean.Matrix.1 <- mean ( Matrix.1 )
  mean.Matrix.2 <- mean ( Matrix.2 )
  # define difference of means
  mean.difference <- mean.Matrix.1 - mean.Matrix.2
  # define covariances for supplied matrices
  cv.Matrix.1 <- cov ( Matrix.1 )
  cv.Matrix.2 <- cov ( Matrix.2 )
  # define the halfsum of cv's as "p"
  p <- ( cv.Matrix.1 + cv.Matrix.2 ) / 2
  # --%<------------------------------------------------------------------------
  # calculate the Bhattacharryya index
  bh.distance <- 0.125 *t ( mean.difference ) * p^ ( -1 ) * mean.difference +
    0.5 * log (det ( p ) / sqrt (det ( cv.Matrix.1 ) * det ( cv.Matrix.2 )
    )
    )
  # --%<------------------------------------------------------------------------
  # calculate Jeffries-Matusita
  # following formula is bound between 0 and 2.0
  jm.distance <- 2 * ( 1 - exp ( -bh.distance ) )
  # also found in the bibliography:
  # jm.distance <- 1000 * sqrt (   2 * ( 1 - exp ( -bh.distance ) )   )
  # the latter formula is bound between 0 and 1414.0
  # --%<------------------------------------------------------------------------
  # calculate the divergence
  # trace (is the sum of the diagonal elements) of a square matrix
  trace.of.matrix <- function ( SquareMatrix ) {
    sum ( diag ( SquareMatrix ) ) }
  # term 1
  divergence.term.1 <- 1/2 * trace.of.matrix (( cv.Matrix.1 - cv.Matrix.2 ) *
                                                ( cv.Matrix.2^ (-1) - cv.Matrix.1^ (-1) )
  )
  # term 2
  divergence.term.2 <- 1/2 * trace.of.matrix (( cv.Matrix.1^ (-1) + cv.Matrix.2^ (-1) ) *
                                                ( mean.Matrix.1 - mean.Matrix.2 ) *
                                                t ( mean.Matrix.1 - mean.Matrix.2 )
  )
  # divergence
  divergence <- divergence.term.1 + divergence.term.2
  # --%<------------------------------------------------------------------------
  # and the transformed divergence
  transformed.divergence  <- 2 * ( 1 - exp ( - ( divergence / 8 ) ) )

  # KS stat
  ks <- ks.test(Vector.1 , Vector.2)$statistic

  indices <- c(jm=jm.distance, bh=bh.distance, div=divergence,
               tdiv=transformed.divergence, ks=ks)
  return(indices)
}
