#' @title bfastOpt
#'
#' @description Some test description
#'
#' @param x Numeric. Time series vector
#' @param dates Time (date) index of \code{x}
#' @param start Numeric. Start of monitoring period. See \code{\link{bfastmonitor}}
#' @param monend Numeric. Optional: end of the monitoring period. See \code{\link{bfmSpatial}}
#' @param mT Numeric. Threshold for breakpoints. All breakpoints with magnitude < \code{magnThresh} are ignored, all data \code{adjust} years preceding this breakpoint are deleted, and monitoring continues with a revised history period.
#' @param mP Numeric. Period of time (in years) from which to compute magnitude (as median of residuals)
#' @param adjust Numeric. Delete all observations within this many years preceding a breakpoint, if that breakpoint has a magnitude less than \code{magnThresh}
#' @param verbose Logical. Display error messages (if any)?
#' @param formula See \code{\link{bfastmonitor}}
#' @param order See \code{\link{bfastmonitor}}
#' @param returnBFM Logical. Return a regular \code{bfastmonitor} object?
#' @param ... Arguments to pass to \code{\link{bfastmonitor}}
#'
#' @return ... Test
#'
#' @author E. van Schaik
#'
#' @export


bfastOpt<- function(x, cells, start, shp, nIter = 10,
                    lmin = 0.001, lmax = 0.05,
                    hmin = 0.25, hmax = 1.0,
                    mtmin = -1000, mtmax = 1000,
                    mpmin = 0.5, mpmax = 2.0,
                    monend = NULL, adjust = 1, verbose = FALSE,
                    formula = response ~ harmon, mc.cores = 1,
                    order = 1, returnBFM = FALSE, ...) {

  # Check / correct parameter values
  if (lmin > lmax | mtmin > mtmax | mpmin > mpmax) {
    stop("Min value larger than max value")
  }

  lmin <- max(min(lmin, 0.05), 0.001)
  lmax <- min(max(lmax, 0.001), 0.05)
  mpmin <- max(mpmin, 0)

  # For all parameters create a random subsample (uniform)
  levelRand       <- runif(nIter, lmin, lmax)
  hRand           <- sample(c(0.25, 0.5, 1.0), nIter, replace=TRUE)
  magnThreshRand  <- runif(nIter, mtmin, mtmax)
  magnPeriodRand  <- runif(nIter, mpmin, mpmax)
  allRand         <- cbind(levelRand, hRand, magnThreshRand, magnPeriodRand)

  # Obtain dates information
  tsZoo <- zoo(x[cells[1]][1,], getSceneinfo(names(x))$date)
  dates <- time(tsZoo)

  # Create function that allows parallel processing
  mcFun <- function(i) {
    var <- bfmThresh(x = as.numeric(x[cells[i]]), start=start, dates = dates, level=level, h=h)
    return(var)
  }

  # Set
  nBreakpoint   <- c(1:nIter)
  locBreakpoint <- shp@data

  # Run for all iterations
  for (j in 1:nIter) {

    # Set parameters
    level      <- allRand[j,1]
    h          <- allRand[j,2]
    magnThresh <- allRand[j,3]
    magnPeriod <- allRand[j,4]

    # Run for one set of parameters
    data <- mclapply(1:length(cells), mcFun, mc.cores=mc.cores)

    # Put results in a matrix
    data <- matrix(unlist(data), ncol=9, byrow = TRUE)

    # Add to breakpoint
    nBreakpoint[j] <- sum(!is.na(data[,2]))
    locBreakpoint  <- cbind(locBreakpoint, data[,2])

    # Print iteration
    print(j)

  }

  # Validate
  omList  <- NULL
  comList <- NULL

  # Loop
  for (i in 1:length(nBreakpoint)) {

    bestFit <- i + 1

    # Create vectors for (no) breakpoints
    deforest <- valShp@data == 1 | valShp@data == 3 | valShp@data == 5
    forest   <- valShp@data == 0 | valShp@data == 2 | valShp@data == 4

    # Get accuracy
    omission   <- sum(!is.na(locBreakpoint[deforest,bestFit])) / length(locBreakpoint[deforest,bestFit]) * 100
    commission <- sum(is.na(locBreakpoint[forest,bestFit])) / length(locBreakpoint[forest,bestFit]) * 100

    # Store
    omList  <- c(omList, omission)
    comList <- c(comList, commission)
  }

  returnList <- list(allRand, omList, comList)

  return(returnList)


}
