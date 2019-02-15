## make garchFit shut up
f_garch <- function(...) {
  capture.output({
    fit <- tryCatch(garchFit(...), error = function(e) e)
  })
  fit
}

## get AIC for successfully fitted GARCH models
get_aic <- function(x) {
  if (any(class(x) == "error")) {
    return(NA_real_)
  } else {
    x@fit$ics[["AIC"]]
  }
}

## simple RMSE function
rmse <- function(pred, obs) {
  sqrt(mean((obs - pred)^2))
}

powt1 <- function(x) {
  ## transform a ts object back into some 1-d dendro data.frame...
  .x <- data.frame(x)
  rownames(.x) <- attributes(x)$tsp[1]:attributes(x)$tsp[2]
  class(.x) <- c("rwl", "data.frame")
  ..x <- dplR::powt(.x)
  if (identical(..x, .x)) {
    return(x)
  } else {
    ## transform back into ts object
    ts(..x[, 1], start = attributes(x)$tsp[1])
  }
}

#' View individual series treatments of the prewhitening scheme
#'
#' @param x a data.frame as returned from \code{pwgarch}
#'
#' @return nothing, invoked for side effects (printing)
#' @export
treatments <- function(x) {
  if (any(class(x) == "pwgarch")) {
    print(attr(x, "treatments"))
  } else {
    cat("Nothing to do.\n")
  }
}

#' Prewhitening with GARCH model
#'
#' This is an experimental alternative prewhitening scheme for data with heavy
#' heteroscedasticity. As a first step, an ARMA model is fitted to the each
#' series, with the optimal order of MA and AR determined automatically. This
#' model is then used to prewhiten the series. If there is no significant
#' heteroscedasticity in the residuals, the series is returned. Otherwise, a new
#' ARMA model is fitted to the series after applying power transformation.
#' Again, if the residuals do not show significant heteroscedasticity, the
#' series is returned. Otherwise, a ARMA(ar, ma) + GARCH(alpha, beta) model is
#' fit to the series, where ar and ma are taken from the initial ARMA model, and
#' alpha and beta are variied according to user input. The best model is
#' selected according to AIC and the residuals are returned as prewhitended
#' series.
#' @param x a data.frame of tree-ring indices, dplR-style
#' @param alpha alpha parameter for arma(ar, ma) + garch(alpha, beta) model, can
#'   be given as integer skalar or integer vector
#' @param beta beta parameter for arma(ar, ma) + garch(alpha, beta) model, can
#'   be given as integer skalar or integer vector
#' @param lm should the Lagrange-Modifier test be used to check for
#'   heteroscedasticity?
#' @param verbose do you want to be informed what's going on?
#'
#' @return a dplR-style data.frame of prewhitened tree-ring series; in addition
#'   to a regular dplR data.frame, the return value has the attribute
#'   "treatments", which can easily be viewed using the function
#'   \code{treatments} on the object.
#' @import fGarch
#' @import forecast
#' @import aTSA
#' @import dplR
#' @export
#'
#' @examples
#' library(dplR)
#' data("cana533")
#' cana_detr <- detrend(cana533, method = "Spline", nyrs = 32)
#' cana_pwg <- pwgarch(cana_detr)
pwgarch <- function(x, alpha = 1, beta = 1:3, lm = FALSE, verbose = TRUE) {

  ## create container for results; make it compatible with dplR; we
  ## might, however, store further bits of information in the
  ## attributes
  n <- ncol(x)
  yrs <- as.numeric(rownames(x))
  .x <- as.data.frame(matrix(NA_real_, ncol = ncol(x), nrow = nrow(x)))
  rownames(.x) <- rownames(x)
  colnames(.x) <- colnames(x)
  class(.x) <- c("pwgarch", "rwl", "data.frame")
  treatments <- data.frame(series = character(n),
                           treatment = character(n),
                           stringsAsFactors = FALSE)

  ## loop through individual series (= trees/cores)
  for (i in 1:n) {

    if (verbose) {
      cat("Processing series", i, "of", n, ":\n")
      cat("...fitting ARMA model...\n")
    }

    .name <- names(x)[i]
    treatments[i, 1]  <- .name

    ## reduce to non-NA observations and convert into time series
    x_i <- x[, i]
    x_nna <- which(!is.na(x_i))
    start_year <- yrs[x_nna][1]
    x_i <- ts(x_i[x_nna], start = start_year)

    #############
    ## 1. ARMA ##
    #############

    ## first, we fit an ARMA model to the data; we expect that in the
    ## default case, an ARMA in combination with detrended data should
    ## be sufficient to describe the memory structure in the data

    ## here we use an automated way for getting the best ARIMA model
    aa <- forecast::auto.arima(x_i, stationary = FALSE)
    arma_ar <- aa$arma[1]
    arma_ma <- aa$arma[2]
    arma_cal <- arima(x_i, c(arma_ar, 0, arma_ma))
    mlt <- aTSA::arch.test(arma_cal, output = FALSE)
    ## we decide if there is a significant heteroscedasticity in the
    ## data, based on the Portmanteau Q test by default; this seems to
    ## be relatively insensitive; we could the more aggressive
    ## Lagrange Multiplier test instead by setting arg lm = TRUE
    if (lm) {
      need_powt <- any(mlt[, 5] < 0.05)
    } else {
      ## the default!
      need_powt <- any(mlt[, 3] < 0.05)
    }

    ## TODO sneak in power transformation

    if (!need_powt) {
      ## we are finished with this series and just use ARMA for
      ## prewhitening
      res <- residuals(arma_cal)
      cut_start <- arma_ar
      ## write into result data.frame
      if (cut_start > 0) {
        res[1:cut_start] <- NA_real_
      }
      .x[x_nna, i] <- res
      treatments[i, 2]  <- "ARMA"
    } else {

      #############################
      ## 2. Power transformation ##
      #############################

      if (verbose) {
        cat("...fitting ARMA on power-transformed series...\n")
      }

      ## in case of heteroscedasticity, we assume that a simple
      ## power-transformation will work in many cases, so we try that
      ## first and check again for heteroscedasticity after the
      ## power-transformation.

      x_i_powt <- powt1(x_i)

      ## check if power transformation has been performed
      has_powt <- !identical(x_i, x_i_powt)
      ## go to GARCH immediately if data was not transformed
      if (!has_powt) need_garch <- TRUE

      if (has_powt) {
        ## do arima again for power-transformed data, explanations see
        ## above
        aa2 <- forecast::auto.arima(x_i_powt, stationary = FALSE)
        arma_ar2 <- aa2$arma[1]
        arma_ma2 <- aa2$arma[2]
        arma_cal2 <- arima(x_i_powt, c(arma_ar2, 0, arma_ma2))
        mlt <- aTSA::arch.test(arma_cal2, output = FALSE)
        ## we decide if there is a significant heteroscedasticity in the
        ## data, based on the Portmanteau Q test by default; this seems to
        ## be relatively insensitive; we could the more aggressive
        ## Lagrange Multiplier test instead by setting arg lm = TRUE
        if (lm) {
          need_garch <- any(mlt[, 5] < 0.05)
        } else {
          ## the default!
          need_garch <- any(mlt[, 3] < 0.05)
        }
      }

      if (!need_garch) {
        ## we are finished with this series and use
        ## power-transformation and ARMA prewhitening
        res <- residuals(arma_cal2)
        cut_start <- arma_ar2
        ## write into result data.frame
        if (cut_start > 0) {
          res[1:cut_start] <- NA_real_
        }
        .x[x_nna, i] <- res
        treatments[i, 2]  <- "ARMA + Power-Transform"
      } else {

        ##############
        ## 3. GARCH ##
        ##############

        if (verbose) {
          cat("...fitting GARCH on ARMA-modelled series...\n")
        }

        ## we need a full GARCH to get rid of the heteroscedasticity
        ## (at least, we try to!); we model the raw data (w/o
        ## power-transformation) as an arma(ar, ma) + garch(alpha,
        ## beta) process, where ar and ma come from the previously
        ## fitted model (before power-transformation)
        ## the alpha and beta are varied according to function
        ## parameters alpha and beta

        mparams <- expand.grid(arma_ar, arma_ma, alpha, beta)
        formulas <- apply(mparams, 1, function(x) {
          f <- paste("~ arma(", x[1], ",", x[2], ") + garch(",
                     x[3], ",", x[4], ")")
          as.formula(f)
        })
        fits <- lapply(formulas, function(x) f_garch(x, data = x_i))
        aics <- sapply(fits, get_aic)
        best_index <- which.min(aics)[1]

        if (is.na(best_index)) {
          # FIXME GARCH could not be fitted, we skip the series!
          if (verbose) {
            cat("...GARCH could not be fitted, skipping series...\n")
          }
          .x[, i] <- NA_real_
          treatments[i, 2]  <- "Skipped"
        } else {
          best_model <- fits[[best_index]]
          res <- residuals(best_model)
          cut_start <- mparams[best_index, 1]
          ## write into result data.frame
          if (cut_start > 0) {
            res[1:cut_start] <- NA_real_
          }
          .x[x_nna, i] <- res
          treatments[i, 2]  <- "ARMA + GARCH"
        }
      }
    }
  }
  attr(.x, "treatments") <- treatments
  .x
}
