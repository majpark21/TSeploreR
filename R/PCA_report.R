#' Report PCA
#'
#' A function for running PCA and optionally get plots out of it and
#' Davies-Bouldin index.
#'
#' @param data A data table which contains time series in long format or a
#'   numeric matrix with time series per row.
#' @param what A character vector describing how the data should be handled
#'   before running pca, as well as which representations should be plotted.
#'   Valid isntructions are: c("cast.and.pca", "nocast.and.pca", "pca.only",
#'   "plot.pca", "plot.variance", "plot.extremes", "plot.loading", "DBindex") If
#'   data are a data.table in long format, use "cast.and.pca", if data are in a
#'   matrix use "nocast.and.pca". The "plot.xxx" settings call biplot.PCA and
#'   visualize.extremes.PCA.
#' @param na.fill Value to replace NA after casting data table from wide to
#'   long.
#' @param measure.var Character. Column name of the measurement that defines the
#'   time series in time.
#' @param condition.var Character vector. Column names used for casting from
#'   long to wide. Should also contain the name of the column used for coloring
#'   the PCA plot if requested. A combination of these variables must be
#'   sufficient to identify unambiguously a single trajectroy in long data
#'   table.
#' @param time.var Character. Column name of the time measure.
#' @param label.color.pca Character or Vector used for coloring PCA. If data is
#'   long data.table (i.e. 'what' is set to "cast.and.pca") should contain the
#'   name of the column used for coloring; note that this column should also be
#'   provided in 'condition.var'. If data is a matrix (i.e. 'what' is set to
#'   "nocast.and.pca') vector of length equal to number of rows in data.
#' @param center.pca Should variables be centered before running PCA. Default is
#'   TRUE.
#' @param scale.pca Should variables be scaled before running PCA. Default is
#'   TRUE, but is susceptible to be changed.
#' @param n.extremes Numeric. How many extremes trajectories to plot.
#' @param PC.biplot Numeric vector of length 2. PCs to use for biplot.
#' @param PC.extremes Numeric, PC from which to plot the extreme individuals.
#' @param PC.db Numeric vector. PC from which computing DBindex
#' @param PC.loading Numeric, plot loading ('composition') of these PCs.
#' @param log.loading Logical, should loading be log?
#' @param group.db A character. Variable to be use as grouping factor when
#'   computing Davies-Bouldin
#' @param var.db A numeric between 0 and 1. If provided, Davies-Bouldin will be
#'   computed on as many PCs as necessary to reach the value.
#' @param ... additional parameters for biplot.PCA, visualize.extremes and. For
#'   example var.axes=F to remove variable arrows or tails = "positive" to plot
#'   only extremes trajectories on positive tail of the PCs.
#'
#' @return If 'what' is set to "pca.only", returns PCA object. Otherwise plot
#'   PCA result.
#' @export
#'
#' @examples
#' library(data.table)
#' library(ggplot2)
#' # Create some dummy data, imagine 20 time series under four conditions A, B, C or D in a long data.table
#' number.measure <- 101
#' mydata <- data.table(Condition = rep(LETTERS[1:4], each = 20*number.measure),
#'  Label = rep(1:20, each = number.measure),
#'  Time=rep(seq(0,100), 80))
#'
#' # A: oscillate around 1 for 10 time units, then shift to oscillation around 1.3
#' # B: oscillate around 1 for 10 time units, then shift to oscillation around 1.25
#' # C: oscillate around 1 for 10 time units, then peak to 1.3 and gets back to 1
#' # D: oscillate around 1 all along trajectory
#'
#' mydata[Condition=="A", Measure := c(rnorm(10, 1, 0.05), rnorm(91, 1.3, 0.05)), by = "Label"]
#' mydata[Condition=="B", Measure := c(rnorm(10, 1, 0.05), rnorm(91, 1.25, 0.05)), by = "Label"]
#' mydata[Condition=="C", Measure := c(rnorm(10, 1, 0.05), rnorm(91, 1.3, 0.05) - seq(0, 0.3, length.out = 91)), by = "Label"]
#' mydata[Condition=="D", Measure := rnorm(101, 1, 0.05), by = "Label"]
#' ggplot(mydata, aes(x=Time, y=Measure)) + geom_line(aes(group=Label), alpha = 0.3) + facet_wrap("Condition") +
#'  stat_summary(fun.y = mean, geom = "line", col = "red", size = 1.25)
#'
#' report.PCA(mydata, what = c("cast.and.pca", "plot.variance", "plot.pca", "plot.extremes", "plot.loading", "DBindex"),
#' measure.var="Measure", condition.var=c("Condition", "Label"), time.var="Time",
#' center.pca = T, scale.pca = F,
#' PC.biplot=c(1,2), label.color.pca = "Condition", var.axes = F, n.extremes = 2,
#' PC.db = NULL, var.db = 0.8)
#'
report.PCA <- function(data, what = c("cast.and.pca", "plot.pca", "plot.variance", "plot.extremes", "plot.loading", "DBindex"),
                       measure.var = NULL, condition.var = NULL, time.var = NULL,
                       center.pca = T, scale.pca = F,
                       PC.db = c(1,2), var.db = NULL,
                       PC.biplot = c(1,2), PC.extremes = 1,
                       PC.loading = PC.biplot, log.loading = F,
                       na.fill = NULL, label.color.pca = NULL, group.db = label.color.pca, n.extremes = NULL,...){
  require(ggbiplot)
  require(heatmap3)
  # Argument check
  arg.what <- c("cast.and.pca", "nocast.and.pca", "pca.only", "plot.pca", "plot.variance", "plot.extremes", "plot.loading", "DBindex")
  if(!all(what %in% arg.what)) stop(paste("'what' arguments must be one of:", paste(arg.what, collapse = ", ")))
  if(!any(c("cast.and.pca", "nocast.and.pca") %in% what)) stop("One of c('cast.and.pca', 'nocast.and.pca') must be provided. Use cast.and.pca for data.table in long format. Use nocast.and.pca for numeric matrix ready for stats::prcomp")
  if(all(c("cast.and.pca", "nocast.and.pca") %in% what)) stop("Only one of c('cast.and.pca', 'nocast.and.pca') must be provided. Use cast.and.pca for data.table in long format. Use nocast.and.pca for numeric matrix ready for stats::prcomp")
  if(("data.frame" %in% class(data) & "nocast.and.pca" %in% what) | ("numeric" %in% class(data) & "cast.and.pca" %in% what)) warning("Use cast.and.pca for data.table in long format. Use nocast.and.pca for numeric matrix ready for stats::prcomp")

  # If data is a long data.table
  if("cast.and.pca" %in% what){
    cast <- cast.and.fill(data = data, condition.var = condition.var, time.var = time.var, na.fill = na.fill, measure.var = measure.var)
    pca <- prcomp(cast$mat2, center = center.pca, scale = scale.pca)
  }
  # If data is already a wide numeric matrix
  else if("nocast.and.pca" %in% what){
    pca <- prcomp(data, center = center.pca, scale. = scale.pca)
  }

  # Stop and return PCA object
  if("pca.only" %in% what) return(list(matrix=cast$mat, pca=pca))

  # Davies-Bouldin index
  if("DBindex" %in% what){
    explained.var <- cumsum(pca$sdev^2L)/sum(pca$sdev^2)
    # If threshold of explained variance is not provided
    if(is.null(var.db)){
      if("cast.and.pca" %in% what){
        db <- PCA.DB(pca = pca, PC = PC.db, labels.pca = cast$mat[, condition.var], label.cluster = group.db, PCvar = var.db)
      } else if ("nocast.and.pca" %in% what){
        db <- PCA.DB(pca = pca, PC = PC.db, labels.pca = label.color.pca, label.cluster = label.color.pca, PCvar = var.db)
      }
      explained.var <- explained.var[max(PC.db)]
    } else {
      if(!is.null(PC.db)) warning("Both 'PC.db' (with default = c(1,2)) and 'var.db' were provided. The former will be ignored.")
      # If threshold of explained variance is provided
      if("cast.and.pca" %in% what){
        db <- PCA.DB(pca = pca, PCvar = var.db, labels.pca = cast$mat[, condition.var], label.cluster = group.db)
      } else if ("nocast.and.pca" %in% what){
        db <- PCA.DB(pca = pca, PCvar = var.db, labels.pca = label.color.pca, label.cluster = label.color.pca)
      }
      PC.db <- seq(1, (max(which(explained.var < var.db)) + 1))
      explained.var <- explained.var[max(which(explained.var < var.db)) + 1]
    }
    cat(paste0("Davies-Bouldin Index: ", db),
        paste0("With grouping variable: ", group.db),
        paste0("Based on PCs: ", paste(PC.db, collapse = ", ")),
        paste0("Variance explained with these components: ", explained.var),
        sep = "\n")
  }

  # Plot explained variance
  if("plot.variance" %in% what) plot(pca)

  # Biplot PCA
  if("plot.pca" %in% what){
    if(is.null(label.color.pca)){
      stop("'label.color.pca' must be a name of a variable in 'data' and provided in 'condition.var' if data is a data.table. It must be a vector witht the groups if data is a matrix.")
    }

    if("nocast.and.pca" %in% what){
      plot(biplot.PCA(pca, label.color.pca, PC = PC.biplot, ...))
    }

    else if("cast.and.pca" %in% what & label.color.pca %in% colnames(cast$mat)){
      label.color.pca <- unlist(cast$mat[, label.color.pca])
      plot(biplot.PCA(pca, label.color.pca, PC = PC.biplot, ...))
    } else {
      stop("'label.color.pca' must be a name of a variable in 'data' and provided in 'condition.var' if data is a data.table. It must be a vector witht the groups if data is a matrix.")
    }
  }

  # Extremes plotting
  if("plot.extremes" %in% what){
    if(is.null(n.extremes)){
      n.extremes <- 3
      warning("n.extremes was not provided and automatically set to 3.")
    }
    if("cast.and.pca" %in% what) visualize.extremes.PCA(pca, cast$mat2, PC = PC.extremes, n = n.extremes, tails = "both", interact=F)
    else if("nocast.and.pca" %in% what) visualize.extremes.PCA(pca, data, PC = PC.extremes, n = n.extremes, tails = "both", interact=F)
  }

  # PC loadings plotting
  if("plot.loading" %in% what){
    PC.loading <- paste0("PC", PC.loading)
    loading.mat <- t(pca$rotation[, PC.loading])
    title <- "PC loadings - Contribution variables to PC"
    if(log.loading){loading.mat <- log(loading.mat); title <- paste("LOG", title)}
    heatmap3(loading.mat, Rowv = NA, Colv = NA, scale = "none", xlab = "Time variables", ylab = "Principal Components", main = title)
  }
}


#' cast.and.fill
#'
#' Cast a data table from long to wide and fill missing values. Output is
#' suitable for PCA.
#' @param data A data table in long format.
#' @param condition.var Character vector, names of variables used for rows in
#'   casting, define the conditions of the experiment.
#' @param time.var Character, name of Time variable used for columns in casting
#' @param measure.var Character, name of measurement variable, on which PCA is
#'   to be performed.
#' @param na.fill Numeric, values to replace NAs
#'
#' @return Two matrices. "mat" is a wide character matrix which contains the
#'   casted "condition.var" as first columns followed by casted measurements.
#'   "mat2" is a numeric wide matrix, which is essentially the same with
#'   conditions trimmed.
#' @export
#'
cast.and.fill <- function(data, condition.var, time.var, measure.var, na.fill){
  formula <- as.formula(paste0(paste(condition.var, collapse = " + "), " ~ ", time.var))
  mat <- as.matrix(dcast(data, formula, value.var = measure.var))
  mat2 <- mat[,-(1:length(condition.var))]
  class(mat2) <- "numeric"
  mat2[which(is.na(mat2))] <- na.fill
  return(list(mat = mat, mat2 = mat2))
}


biplot.PCA <- function(pca.obj, group.vec, ..., PC = c(1,2), obs.scale = 1, var.scale = 1, ellipse = T, circle = T, var.axes = T){
  require(ggbiplot)
  g <- ggbiplot(pca.obj, ..., obs.scale = obs.scale, var.scale = var.scale,
                groups = group.vec, ellipse = ellipse,
                circle = circle, var.axes = var.axes, choices = PC)
  g <- g + scale_color_discrete(name = '')
  g <- g + theme(legend.direction = 'horizontal',
                 legend.position = 'top')
  return(g)
}


visualize.extremes.PCA <- function(pca.object, matrix.data, PC, n, tails = "both", interact=T){
  # Visualize the n extreme individuals of the component PC. Tails indicates if
  # the extremes are to be picked from left, right or both tails of PC.
  # Matrix.data is the matrix in wide format that was used to run the pca
  # analysis with stats::prcomp.
  ord <- order(pca.object$x[,PC])
  ordd <- rev(ord)
  xaxis <- as.numeric(rownames(pca.object$rotation))
  if(tails=="both"){
    par(mfrow=c(1,2))
    for(i in 1:n){
      mini <-  min(matrix.data[ord[i],], matrix.data[ordd[i],])
      maxi <- max(matrix.data[ord[i],], matrix.data[ordd[i],])
      plot(xaxis, matrix.data[ord[i],], type = "l", ylim = c(mini, maxi), xlab = "Time", ylab = "Trajectory", main = paste0("Negative Tail - Coord PC", PC, ": ", round(pca.object$x[ord[i], PC], 4)))
      plot(xaxis, matrix.data[ordd[i],], type = "l", ylim = c(mini, maxi), xlab = "Time", ylab = "Trajectory", main = paste0("Positive Tail - Coord PC", PC, ": ", round(pca.object$x[ordd[i], PC], 4)))
      if(interact) readline(prompt="Press [enter] to see next plots")
    }
  }

  else if(tails=="positive"){
    for(i in 1:n){
      plot(xaxis, matrix.data[ordd[i],], type = "l", xlab = "Time", ylab = "Trajectory", main = paste0("Positive Tail - Coord PC", PC, ": ", round(pca.object$x[ordd[i], PC], 4)))
      if(interact) readline(prompt="Press [enter] to see next plots")
    }
  }

  else if(tails=="negative"){
    for(i in 1:n){
      plot(xaxis, matrix.data[ord[i],], type = "l", xlab = "Time", ylab = "Trajectory", main = paste0("Negative Tail - Coord PC", PC, ": ", round(pca.object$x[ord[i], PC], 4)))
      if(interact) readline(prompt="Press [enter] to see next plots")
    }
  }
}


#' PCA.DB
#'
#' Compute the Davies-Bouldin index for a pca clustering. This is used to
#' quantify how much groups are taken away in PCA space.
#'
#' @param pca A pca object.
#' @param PC Which PCs should be used to compute the DB index. Consider only a
#'   few components to use only leading trends in data.
#' @param PCvar If provided, replace PC. Instead of using a given set of PC,
#' will use as many as necessary to reach PCvar% of explained variance.
#' @param labels.pca A data.table or character matrix with a number of rows
#'   equal to the number of time series used in pca. Each column contains a
#'   different label for the associated trajectory. Note that column names must
#' @param label.cluster
#'   be provided.
#' @return A numeric, the Davies-Bouldin index of the clustering.
#' @export
#'
PCA.DB <- function(pca, PC = NULL, PCvar = NULL, labels.pca, label.cluster){
  require(RDRToolbox)
  if(!label.cluster %in% colnames(labels.pca)) stop("'cluster.label' must a single character contained in the column names of 'labels.pca'")
  if((is.null(PC) & is.null(PCvar)) | (!is.null(PC) & !is.null(PCvar))) stop("One of PC or PCvar (and one only!) must be provided")

  # Use provided PC
  if(is.null(PCvar)){
    PC <- paste0("PC", PC)
  } else {
    # Use as many PC as necessary to reach PCvar explained variance
    explained.var <- cumsum(pca$sdev^2L)/sum(pca$sdev^2)
    if(explained.var[1] > PCvar) stop("'PCvar' is inferior to variance explained by first component.")
    PC <- max(which(explained.var < PCvar)) + 1
    PC <- paste0("PC", seq(1,PC))
  }
  # Turn pca coordinates into a data table
  pca_coord <- as.data.table(pca$x)
  for(label in colnames(labels.pca)){
    pca_coord[, (label) := labels.pca[,label]]
  }
  db <- DBIndex(data = as.matrix(pca_coord[, ..PC]), labels = unlist(pca_coord[,..label.cluster]))
  return(db)
}
