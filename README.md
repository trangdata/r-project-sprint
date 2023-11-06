
<!-- README.md is generated from README.Rmd. Please edit that file -->

# r-project-sprint

<!-- badges: start -->
<!-- badges: end -->

<https://contributor.r-project.org/r-project-sprint-2023/projects/tweak-printCoefmat/>

<https://github.com/r-devel/r-svn/blob/236f2202c3ffb739fbd2b99bdb6c2670c00ce53b/src/library/stats/R/anova.R#L112-L114>

`printCoefmat3`:

``` r
  if (length(cs.ind)) {
    nozap.ind <- setdiff(cs.ind, zap.ind)
    acs <- abs(xm[, cs.ind, drop = FALSE])
    larger_acs <- abs(acs) >= 10^(-digits+1)
    larger_acs[, nozap.ind] <- TRUE
    if (any(ia <- is.finite(acs)) && length(acs <- acs[ia & larger_acs])) {
      digmin <- 1 + floor(log10(range(acs, finite = TRUE)))
    }
  }
```
