printCoefmat3 <- function(x, digits = max(3L, getOption("digits") - 2L),
                         signif.stars = getOption("show.signif.stars"),
                         signif.legend = signif.stars,
                         dig.tst = max(1L, min(5L, digits - 1L)),
                         cs.ind = 1:k,
                         tst.ind = k + 1,
                         zap.ind = integer(),
                         P.values = NULL,
                         has.Pvalue = nc >= 4L && length(cn <- colnames(x)) &&
                           substr(cn[nc], 1L, 3L) %in% c("Pr(", "p-v"),
                         eps.Pvalue = .Machine$double.eps,
                         na.print = "NA",
                         quote = FALSE,
                         right = TRUE,
                         ...) {
  if (is.null(d <- dim(x)) || length(d) != 2L) {
    stop("'x' must be coefficient matrix/data frame")
  }
  nc <- d[2L]
  if (is.null(P.values)) {
    scp <- getOption("show.coef.Pvalues")
    if (!is.logical(scp) || is.na(scp)) {
      warning("option \"show.coef.Pvalues\" is invalid: assuming TRUE")
      scp <- TRUE
    }
    P.values <- has.Pvalue && scp
  } else if (P.values && !has.Pvalue) {
    stop("'P.values' is TRUE, but 'has.Pvalue' is not")
  }
  if (has.Pvalue && !P.values) {
    d <- dim(xm <- data.matrix(x[, -nc, drop = FALSE]))
    nc <- nc - 1
    has.Pvalue <- FALSE
  } else {
    xm <- data.matrix(x)
  }
  k <- nc - has.Pvalue - (if (missing(tst.ind)) {
    1
  } else {
    length(tst.ind)
  })
  if (!missing(cs.ind) && length(cs.ind) > k) {
    stop("wrong k / cs.ind")
  }
  Cf <- array("", dim = d, dimnames = dimnames(xm))
  ok <- !(ina <- is.na(xm))
  digmin <- 1
  if (length(cs.ind)) {
    nozap.ind <- setdiff(cs.ind, zap.ind)
    acs <- abs(xm[, cs.ind, drop = FALSE])
    larger_acs <- abs(acs) >= 10^(-digits+1)
    larger_acs[, nozap.ind] <- TRUE
    if (any(ia <- is.finite(acs)) && length(acs <- acs[ia & larger_acs])) {
      digmin <- 1 + floor(log10(range(acs, finite = TRUE)))
    }
  }
  for (i in zap.ind) xm[, i] <- zapsmall(xm[, i], max(1L, digits - digmin))
  if (length(cs.ind)) {
    coef.se <- xm[, cs.ind, drop = FALSE]
    Cf[, cs.ind] <- format(round(coef.se, max(1L, digits -
                                                digmin)), digits = digits)
  }

  if (length(tst.ind)) {
    Cf[, tst.ind] <- format(round(xm[, tst.ind], digits = dig.tst),
      digits = digits
    )
  }
  if (any(r.ind <- !((1L:nc) %in% c(cs.ind, tst.ind, if (has.Pvalue) nc)))) {
    for (i in which(r.ind)) Cf[, i] <- format(xm[, i], digits = digits)
  }
  ok[, tst.ind] <- FALSE
  okP <- if (has.Pvalue) {
    ok[, -nc]
  } else {
    ok
  }
  x1 <- Cf[okP]
  dec <- getOption("OutDec")
  if (dec != ".") {
    x1 <- chartr(dec, ".", x1)
  }
  x0 <- (xm[okP] == 0) != (as.numeric(x1) == 0)
  if (length(not.both.0 <- which(x0 & !is.na(x0)))) {
    Cf[okP][not.both.0] <- format(xm[okP][not.both.0], digits = max(
      1L,
      digits - 1L
    ))
  }
  if (any(ina)) {
    Cf[ina] <- na.print
  }
  if (any(inan <- is.nan(xm))) {
    Cf[inan] <- "NaN"
  }
  if (P.values) {
    if (!is.logical(signif.stars) || is.na(signif.stars)) {
      warning("option \"show.signif.stars\" is invalid: assuming TRUE")
      signif.stars <- TRUE
    }
    if (any(okP <- ok[, nc])) {
      pv <- as.vector(xm[, nc])
      Cf[okP, nc] <- format.pval(pv[okP],
        digits = dig.tst,
        eps = eps.Pvalue
      )
      signif.stars <- signif.stars && any(pv[okP] < 0.1)
      if (signif.stars) {
        Signif <- symnum(pv,
          corr = FALSE, na = FALSE,
          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
          symbols = c("***", "**", "*", ".", " ")
        )
        Cf <- cbind(Cf, format(Signif))
      }
    } else {
      signif.stars <- FALSE
    }
  } else {
    signif.stars <- FALSE
  }
  print.default(Cf,
    quote = quote, right = right, na.print = na.print,
    ...
  )
  if (signif.stars && signif.legend) {
    if ((w <- getOption("width")) < nchar(sleg <- attr(
      Signif,
      "legend"
    ))) {
      sleg <- strwrap(sleg, width = w - 2, prefix = "  ")
    }
    cat("---\nSignif. codes:  ", sleg, sep = "", fill = w +
      4 + max(nchar(sleg, "bytes") - nchar(sleg)))
  }
  invisible(x)
}

zapsmall <- function(x, digits = getOption("digits")) {
  if (length(digits) == 0L) {
    stop("invalid 'digits'")
  }
  if (all(ina <- is.na(x))) {
    return(x)
  }
  mx <- max(abs(x[!ina]))
  round(x, digits = if (mx > 0) max(0L, digits - as.numeric(log10(mx))) else digits)
}

x = c(3.760667e-07, 4.760667e-06, 5.760667e-05,
      6.760667e-04,
      7.760667e-03, 8.760667e-02,
      0.4992263, 0.6009902, - 0.0001)
zapsmall(x, digits = 4)
zapsmall(4.760667e-06, digits = 4
         )
