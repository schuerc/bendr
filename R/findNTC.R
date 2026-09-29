#' Find the non-toxic concentration.
#'
#' If you use the result of this calculation in your research, please cite Stadnicka-Michalak et al.
#'
#' @param bendrObj
#' A bendr result object
#' @param max.extension
#' maximum number of log10 units by which the concentration grid is extended below its start
#' @param measured.criterion
#' apply the measured value criterion (default TRUE)
#' @param narrow.ci.correction
#' apply the correction for narrow confidence intervals (default TRUE)
#' @param details
#' if TRUE, return a list with the intermediate values (NtC_lower, NtC_upper, their effects,
#' the rule that chose the NtC, NtC_measured, the confidence band used) instead of the NtC only
#'
#' @return
#' The absolute (non-logged) non-toxic concentration value, or a list if details = TRUE
#' @export
#'
#' @examples
#' require(tibble)
#'drc.formula = effect ~ 100 / (1 + 10^((logEC50-logconc) * slope))
#'rdata = as_tibble(
#'data.frame(
#'    conc=c(0.1,1,10,100,1000,5000),
#'    replicateA=c(1,12,27,65,88,100),
#'    replicateB=c(2,16,21,62,81,100),
#'    replicateC=c(4,15,21,61,85,96)
#'  ))
#'
#'fitObjectRep = fitdr_replicates(drc.formula,rdata,2:4,conc)
#'ntc = findNTC(fitObjectRep)
#'
findNTC = function(bendrObj, max.extension = 10, measured.criterion = TRUE,
                   narrow.ci.correction = TRUE, details = FALSE) {
  # individual replicates, also when the curve was fitted to concentration means
  measured = if(!is.null(bendrObj$data.long)) bendrObj$data.long else bendrObj$data
  # lowest tested concentration with >= 10% effect in any replicate
  toxic = measured$conc[measured$effect <= 90]
  NtC_measured = if(length(toxic) > 0) min(toxic) else NA_real_

  out = list(NtC = NA_real_, NtC_model = NA_real_, Effect_NtC = NA_real_,
             NtC_lower = NA_real_, NtC_upper = NA_real_,
             Effect_NtC_lower = NA_real_, Effect_NtC_upper = NA_real_,
             rule = NA_character_, measured.applied = FALSE, NtC_measured = NtC_measured,
             n.extensions = 0, plot.data = NULL, message = NA_character_)
  fail = function(msg) {
    warning(paste("findNTC:", msg))
    out$message = msg
    if(details) out else NA_real_
  }

  if(all(is.na(bendrObj$plot.data$ci.values.lower))) {
    stop("findNTC: no confidence bands available (the CI calculation failed), so the NtC cannot be determined")
  }
  # Both criteria look for the highest concentration of a region that starts at the low end
  # (upper CI > 100, lower CI > 90). That region is only found correctly if the grid starts
  # inside it; otherwise a criterion is "not found" and would be treated as effect 0, which
  # wrongly triggers the narrow-CI correction. Extend the grid 1 log unit at a time until
  # both criteria hold at its lowest concentration.
  # The NtC criteria need the delta-method bands: bootstrap bands are quantiles of fitted
  # curves, which stay below 100%, so the upper criterion could never be met. For a bootstrap
  # fit the bootstrap is used only for the EC50/EC10/slope CIs.
  covers_both = function(pd) {
    i = which.min(pd$log.concentration)
    isTRUE(pd$ci.values.upper[i] > 100) && isTRUE(pd$ci.values.lower[i] > 90)
  }
  # always evaluate on the fine grid of extend_ci (0.001 log units), so the result does not
  # depend on the resolution of the plotting grid
  plot.data = extend_ci(bendrObj, ci.method = "delta")
  lox = min(plot.data$log.concentration)
  n.ext = 0
  while(!covers_both(plot.data) && n.ext < max.extension) {
    lox = lox - 1
    n.ext = n.ext + 1
    plot.data = extend_ci(bendrObj, lox = lox, ci.method = "delta")
  }
  out$n.extensions = n.ext
  out$plot.data = plot.data
  if(!covers_both(plot.data)) {
    return(fail(sprintf("the CI criteria are not reached within %d log units below the tested range; NtC is undetermined", max.extension)))
  }

  upperCI = plot.data$ci.values.upper
  lowerCI = plot.data$ci.values.lower
  curve.predict = plot.data$curve.predict
  log.concentration = plot.data$log.concentration
  logEC50 = coef(bendrObj$curve.fit)['logEC50']
  slope = coef(bendrObj$curve.fit)['slope']

  ## NtC upper
  # highest concentration at which the upper confidence interval is above 100
  uci100idx = which(upperCI > 100)
  if(length(uci100idx) > 0) {
    inds = max(uci100idx)
    Effect_NtC_upper = 100 - curve.predict[inds]
    if (Effect_NtC_upper > 0){
      NtC_upper = 10^(log.concentration[inds])
    } else {
      NtC_upper = 0
    }
  } else {
    Effect_NtC_upper = 0
    NtC_upper = 0
  }

  ## NtC lower
  # highest concentration at which the lower confidence interval is above 90
  lci90idx = which(lowerCI > 90)
  if(length(lci90idx) == 0) {
    return(fail("the lower confidence band never exceeds 90% survival, even on the extended range; NtC is undetermined"))
  }
  inds = max(lci90idx)
  Effect_NtC_lower = 100 - curve.predict[inds]
  if (Effect_NtC_lower > 0){
    NtC_lower = 10^(log.concentration[inds])
  } else {
    NtC_lower = 0
  }
  out[c("NtC_lower", "NtC_upper", "Effect_NtC_lower", "Effect_NtC_upper")] =
    list(unname(NtC_lower), unname(NtC_upper), unname(Effect_NtC_lower), unname(Effect_NtC_upper))

  ## Choose NtC based on upper and lower CI criteria
  # In case of very narrow confidence intervals, caused e.g. by the non-optimal concentration range:
  #   If Effect_NtC_lowerCI > 10 * Effect_NtC_upperCI -> Effect_NtC = Effect_NtC_lowerCI / 10
  # The NtC is then calculated from Effect_NtC with the dose-response curve.
  if (NtC_upper > NtC_lower) {
    NtC = NtC_lower
    Effect_NtC = Effect_NtC_lower
    rule = "lower"
  } else if (Effect_NtC_lower == 0 && Effect_NtC_upper == 0) {
    # both effects 0 would give 0/0 in the ratio below
    return(fail("predicted effect is 0 at both CI criteria; NtC is undetermined"))
  } else if (narrow.ci.correction && Effect_NtC_lower/Effect_NtC_upper > 10) {
    Effect_NtC = Effect_NtC_lower/10
    NtC = 10^(logEC50-log10(100/(100-Effect_NtC)-1)/slope)
    rule = "narrow_ci_correction"
  } else {
    NtC = NtC_upper
    Effect_NtC = Effect_NtC_upper
    rule = "upper"
  }
  NtC = unname(NtC)
  out[c("NtC_model", "Effect_NtC", "rule")] = list(NtC, unname(Effect_NtC), rule)

  ## Measured value criterion: if any replicate at a concentration below the NtC shows
  # >= 10% effect, the NtC is the lowest such concentration
  if(measured.criterion) {
    belowNtC = measured %>%
      dplyr::filter(conc < NtC,effect <= 90)
    if(nrow(belowNtC) > 0) {
      NtC = belowNtC %>% pull(conc) %>% min()
      out$measured.applied = TRUE
    }
  }
  out$NtC = NtC

  if(details) out else NtC
}
