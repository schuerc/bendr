# bendr 0.4.1 (2026-09-29)

* `findNTC()`: new arguments `measured.criterion` and `narrow.ci.correction` (both default TRUE,
  as before) and `details`: with `details = TRUE` a list with NtC_lower, NtC_upper, their effects,
  the rule that chose the NtC, NtC_measured and the confidence band used is returned.
* `findNTC()`: always evaluates the criteria on the 0.001-log-unit grid of `extend_ci()`, so the
  NtC no longer depends on the resolution of the plotting grid (delta and bootstrap fits now give
  the same NtC).
* `extend_ci()`: for very flat curves the default lower end (99.9% effect) lay extremely far below
  the data and the grid could not be allocated; the default extension is now limited to
  `max.default = 3` log units, and ranges needing more than 1e6 grid points give an error.

# bendr 0.4.0 (2026-09-29)

Fork of https://github.com/alxbetz/bendr (0.3.0, commit 0f2d407), maintained at
https://github.com/schuerc/bendr.

Fixes, tested with R 4.6.0 (tests/testthat/test-fixes.R).

* `bootstrap_ci()`: removed the line that overwrote row 1 of the parameter table on every iteration.
* `bootstrap_ci()`: parameter columns were filled positionally from `coef()`, whose order follows the
  start vector (logEC50, slope), so slope and logEC50 were swapped; now assigned by name.
* `bootstrap_ci()`: resamples are fitted with `nlsLM` from the original estimates instead of
  `nls2(grid-search)` (which returned grid points, not optimised fits). Failed fits are caught,
  stored as NA, excluded from the quantiles and counted (`n.failed`, with a warning).
* `bootstrap_ci()`: resampling is stratified by concentration, so no level is dropped.
  The original fit is no longer included among the bootstrap replicates.
* `fitdr()`: `level` is now passed to `bootstrap_ci()` (was ignored). The result also holds
  `ci.method` and `boot.params`.
* `fitdr()`: EC10 CI now comes from a fit of the reparametrised model (logEC10 as parameter)
  with `nlsLM`. Previously the defined EC10 formula was unused and `nls2(brute-force)` refitted
  the original model from a single start, so the "EC10 CI" was meaningless.
* `extend_ci()`: uses the fit's CI method (bootstrap bands from stored bootstrap parameters)
  and confidence level; the extended range never shrinks the original range.
* `findNTC()`: clear error if no CI bands exist; `NA` with a warning instead of `-Inf` indexing
  when the lower band never exceeds 90%, or when both criterion effects are 0 (0/0).
* `findNTC()`: the concentration grid is extended (1 log unit at a time, up to `max.extension = 10`)
  until both CI criteria hold at its lowest concentration. Before, the upper-CI crossing could lie
  below the grid (the +/-20% range; `extend_ci` only ran for the lower criterion); it was then
  treated as effect 0, the ratio became Inf and the narrow-CI correction fired, giving a too low NtC.
* `findNTC()`: always uses delta-method bands. With `ci_method = "bootstrap"` the bootstrap is
  used for the EC50/EC10/slope CIs only: bootstrap bands stay below 100%, so the upper-CI criterion
  could never be met (previously the NtC then came from the /10 correction by accident).
* `extend_ci()`: new argument `ci.method` to override the method of the fit.
* `findNTC()`: the measured-data criterion always uses the individual replicates.
* `extend_ci()`: new arguments `lox` and `step` (grid spacing 0.001 log units instead of 1000 points).
* `fitdr_replicates()`: new argument `replicates = c("independent", "same_plate")`. "same_plate"
  fits the concentration means with df = n_conc - 2; bootstrap is refused for it. Default stays
  "independent" (previous behaviour).
