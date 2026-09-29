drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))

# Synthetic data (6 concentrations x 3 replicates): survival never exceeds 90%, and the
# upper-CI crossing lies below the default +-20% plotting range, so findNTC has to extend it
syn = data.frame(
  conc = c(0.5, 1, 2, 4, 8, 16),
  rep1 = c(85.8, 89.9, 37.3, 29.7, 23.4, 2.9),
  rep2 = c(87.8, 70.4, 49.9, 30.5, 20.3, 2.3),
  rep3 = c(84.6, 75.0, 39.7, 27.5, 8.8, 0.5)
)

test_that("EC10 CI comes from the reparametrised fit and brackets EC10", {
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc)
  expect_equal(unname(fo$ec10), unname(calc_ecx(90, fo$slope, fo$ec50)), tolerance = 1e-4)
  expect_true(fo$ec10 < fo$ec50)
  expect_true(fo$ec10.ci["lower"] < fo$ec10 && fo$ec10 < fo$ec10.ci["upper"])
})

test_that("bootstrap keeps parameter columns aligned and reports failed fits", {
  set.seed(1)
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc, ci_method = "bootstrap")
  bp = fo$boot.params
  expect_equal(nrow(bp), 1000)
  ok = !is.na(bp$slope)
  expect_true(all(bp$slope[ok] < 0))
  expect_equal(bp$EC50[ok], 10^bp$logEC50[ok])
  # estimates are optimised, not quantised to a grid
  expect_gt(length(unique(round(bp$logEC50[ok], 6))), 100)
  expect_true(fo$slope.ci["lower"] < fo$slope && fo$slope < fo$slope.ci["upper"])
})

test_that("extend_ci uses the method and level of the fit", {
  set.seed(1)
  fo.boot = fitdr_replicates(drc.formula, syn, 2:4, conc, ci_method = "bootstrap")
  ext = extend_ci(fo.boot)
  expect_true(all(is.finite(ext$ci.values.lower)))
  expect_true(min(ext$log.concentration) <= min(fo.boot$plot.data$log.concentration))

  fo90 = fitdr_replicates(drc.formula, syn, 2:4, conc, level = 0.90)
  fo95 = fitdr_replicates(drc.formula, syn, 2:4, conc, level = 0.95)
  w90 = with(extend_ci(fo90), ci.values.upper - ci.values.lower)
  w95 = with(extend_ci(fo95), ci.values.upper - ci.values.lower)
  expect_true(all(w90 < w95))
})

test_that("findNTC works when survival never exceeds 90%", {
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc)
  ntc = findNTC(fo)
  expect_true(is.finite(ntc) && ntc > 0 && ntc < 16)
})

test_that("findNTC fails clearly when no CI is available", {
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc)
  fo$plot.data$ci.values.lower = NA
  fo$plot.data$ci.values.upper = NA
  expect_error(findNTC(fo), "no confidence bands")
})

test_that("findNTC extends the grid until the upper-CI crossing is found", {
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc)
  # the crossing is not inside the default range ...
  i = which.min(fo$plot.data$log.concentration)
  expect_lte(fo$plot.data$ci.values.upper[i], 100)
  # ... but is found after extension (bendr 0.3.0 missed it and returned 0.0531)
  expect_equal(findNTC(fo), 0.0648, tolerance = 0.01)
})

test_that("same_plate fits concentration means with df = n_conc - 2", {
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc, replicates = "same_plate")
  expect_equal(nrow(fo$data), 6)
  expect_equal(nrow(fo$data.long), 18)
  expect_equal(fo$replicates, "same_plate")
  # df 4 instead of 16: wider band. The NtC goes UP (0.227 vs 0.0648): the upper criterion
  # (highest conc with upper CI > 100) moves to higher concentrations with a wider band.
  expect_equal(findNTC(fo), 0.227, tolerance = 0.01)
})

test_that("bootstrap with same_plate replicates is refused", {
  expect_error(fitdr_replicates(drc.formula, syn, 2:4, conc, ci_method = "bootstrap", replicates = "same_plate"),
               "independent replicates")
})

test_that("findNTC uses delta-method bands also for a bootstrap fit", {
  set.seed(1)
  fb = suppressWarnings(fitdr_replicates(drc.formula, syn, 2:4, conc, ci_method = "bootstrap"))
  fd = fitdr_replicates(drc.formula, syn, 2:4, conc)
  expect_equal(findNTC(fb), findNTC(fd))
})

test_that("findNTC switches and details", {
  fo = fitdr_replicates(drc.formula, syn, 2:4, conc)
  d = findNTC(fo, details = TRUE)
  expect_equal(d$NtC, findNTC(fo))
  expect_true(d$rule %in% c("lower", "upper", "narrow_ci_correction"))
  expect_equal(d$NtC_measured, 0.5)   # all replicates at conc 0.5 are below 90%
  expect_true(is.data.frame(d$plot.data))
  # without the correction the NtC can only come from one of the two criteria
  d0 = findNTC(fo, narrow.ci.correction = FALSE, measured.criterion = FALSE, details = TRUE)
  expect_true(d0$NtC %in% c(d0$NtC_lower, d0$NtC_upper))
  expect_false(d0$measured.applied)
})
