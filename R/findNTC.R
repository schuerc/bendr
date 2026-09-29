#' Find the non-toxic concentration.
#'
#' If you use the result of this calculation in your research, please cite Stadnicka-Michalak et al.
#'
#' @param bendrObj
#' A bendr result object
#' @param max.extension
#' maximum number of log10 units by which the concentration grid is extended below its start
#'
#' @return
#' The absolute (non-logged) non-toxic concentration value
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
findNTC = function(bendrObj, max.extension = 10) {
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
  delta.bands = !identical(bendrObj$ci.method, "bootstrap")
  covers_both = function(pd) {
    i = which.min(pd$log.concentration)
    isTRUE(pd$ci.values.upper[i] > 100) && isTRUE(pd$ci.values.lower[i] > 90)
  }
  plot.data = if(delta.bands) bendrObj$plot.data else extend_ci(bendrObj, ci.method = "delta")
  if(!covers_both(plot.data)) {
    plot.data = extend_ci(bendrObj, ci.method = "delta")
  }
  lox = min(plot.data$log.concentration)
  n.ext = 0
  while(!covers_both(plot.data) && n.ext < max.extension) {
    lox = lox - 1
    n.ext = n.ext + 1
    plot.data = extend_ci(bendrObj, lox = lox, ci.method = "delta")
  }
  if(!covers_both(plot.data)) {
    warning(sprintf("findNTC: the CI criteria are not reached within %d log units below the tested range; NtC is undetermined", max.extension))
    return(NA_real_)
  }

  upperCI = plot.data$ci.values.upper
  lowerCI = plot.data$ci.values.lower
  curve.predict = plot.data$curve.predict
  log.concentration = plot.data$log.concentration
  logEC50 = coef(bendrObj$curve.fit)['logEC50']
  slope = coef(bendrObj$curve.fit)['slope']

  # In case of very narrow confidence intervals, caused e.g. by the non-optimal concentration range:
  #   If Effect_NtClower CI> 10 * Effect_NtCupperCI -> Effect_NtC = Effect_NtClowerCI / 10 where: Effect_NtCi (%) is the effect caused by the NtC
  # NtC can be calculated then from Effect_NtC (and sigmoidal dose-response curve).
  #


  ## NtC upper
  # find concentration at which the upper confidence interval falls below
  # 100
  uci100 = upperCI > 100
  uci100idx = which(uci100)
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
  # find concentration at which lower confidence interval is 90
  lci90idx = which(lowerCI > 90)
  if(length(lci90idx) == 0) {
    warning("findNTC: the lower confidence band never exceeds 90% survival, even on the extended range; NtC is undetermined")
    return(NA_real_)
  }
  inds = max(lci90idx)
  Effect_NtC_lower = 100 - curve.predict[inds]
  if (Effect_NtC_lower > 0){
    NtC_lower = 10^(log.concentration[inds])
  } else {
    NtC_lower = 0
  }




  ## Choose NtC based on upper and lower CI criteria
  if (NtC_upper > NtC_lower) {
    NtC = NtC_lower
    Effect_NtC = Effect_NtC_lower
  } else {


    # both effects 0 would give 0/0 in the ratio below
    if (Effect_NtC_lower == 0 && Effect_NtC_upper == 0) {
      warning("findNTC: predicted effect is 0 at both CI criteria; NtC is undetermined")
      return(NA_real_)
    } else if (Effect_NtC_lower/Effect_NtC_upper > 10) {
      Effect_NtC = Effect_NtC_lower/10
      NtC = 10^(logEC50-log10(100/(100-Effect_NtC)-1)/slope)

    } else {
      NtC = NtC_upper
      Effect_NtC = Effect_NtC_upper
    }
  }

  ## Include measured data if the effect of any replicate (at concentration below
  # the concentration, which has been thus far chosen as the NtC) is > 10%
  # 1. which measured concentrations are below NtC
  # 2. does any of the effects measured at these concentrations have an effect > 10%
  # 3. if yes, select the lowest concentration that is below the NtC


  # individual replicates, also when the curve was fitted to concentration means
  measured = if(!is.null(bendrObj$data.long)) bendrObj$data.long else bendrObj$data
  belowNtC = measured %>%
    dplyr::filter(conc < NtC,effect <= 90)
  if(nrow(belowNtC) > 0) {
    NtC = belowNtC %>% pull(conc) %>% min()
  }

  NtC
}


