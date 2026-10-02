## Test provided by Arkadiusz Gladki 2026-10-02
## drm() must leave the na.action option as the caller set it.

library(drc)

dat <- data.frame(conc = rep(c(0.01, 0.1, 1, 10, 100), 3))
set.seed(1)
dat$response <- 1 / (1 + (dat$conc / 1)^1) + rnorm(nrow(dat), sd = 0.02)

op <- options(na.action = "na.exclude")

drm(response ~ conc, data = dat, fct = LL.4(), na.action = stats::na.omit)

stopifnot(identical(getOption("na.action"), "na.exclude"))

## A namespace-qualified na.action used to be left in the option, where get()
## cannot resolve it, so every later model.frame() in the session failed.
lm(response ~ conc, data = dat)

options(op)
