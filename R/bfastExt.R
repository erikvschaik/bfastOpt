
#' @title bfastExt
#'
#' @description Multicore version of extended BFAST function. Should be replaced with B. DeVries's bfmThresh in due time.
#' @param TEST
#'
#' @return ... Test
#'
#' @author E. van Schaik
#'
#' @export


bfastExt <- function(x, start, monend = NULL, cells = NULL,
                     f = 1, min.thresh = NULL, sensor = NULL,
                     interactive = FALSE, plot = FALSE, order = 1,
                     h = 0.25, level = 0.05, mc.cores = 1,
                     magnThresh = 0, magnPeriod = 1,
                     verbose = FALSE,
                     ...) {

  # Create BFASTExt single-pixel function
  bfastExtSP <- function(x, start, cell, ...) {

    # Created timseries for cell
    tsZoo <- zoo(x[cell][1,], getSceneinfo(names(x))$date)
    bts <- bfastts(tsZoo, time(tsZoo), type='irregular')

    # Initialize iteration
    iterate <- TRUE

    # Start loop
    while(iterate) {

      bfm <- try(bfastmonitor(bts, start = start, level = level, h = h,
                              order = order, formula = response ~ harmon, ...), silent = !verbose)

      Coor <- xyFromCell(x, cell)
      xCoor <- Coor[1]
      yCoor <- Coor[2]

      if(class(bfm) == 'try-error') {
        bkp <- magn <- hist <- rsq <- adjrsq <- NA
        # coefs <- rep(NA, coef_len)
        err <- 1
        break
      } else {
        bkp <- bfm$breakpoint
        hist <- bfm$history[2] - bfm$history[1]
        rsq <- summary(bfm$model)$r.squared
        adjrsq <- summary(bfm$model)$adj.r.squared
        # coefs <- coef(bfm$model)
      }

      # breakpoint magnitude
      if(!is.na(bkp)) {
        postpp <- subset(bfm$tspp, time >= bkp & time <= (bkp + magnPeriod))
        magn <- median(postpp$response - postpp$prediction)
      } else {
        magn <- NA
      }

      if(!is.na(magn) & magn >= magnThresh) {
        bts[(time(bts) <= bkp) & (time(bts) > bkp - magnPeriod)] <- NA
        start <- c(floor(bkp), 1)
      } else {
        err <- 0
        iterate <- FALSE
      }
    }

    res <- c(cell, bkp, magn, hist, rsq, adjrsq, err, xCoor, yCoor)
    names(res) <- c("cell", "breakpoint", "magnitude", "history", "r.squared", "adj.r.squared",  "error", "x", "y")

    return(res)
  }


  # Create function that allows parallel processing
  mcFun <- function(i) {
    var <- bfastExtSP(x, start=start, cell = cells[i])
    return(var)
  }

  # Make sure the proper information is returned from this function
  return(mclapply(1:length(cells), mcFun, mc.cores=mc.cores))

}

