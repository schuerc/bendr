#' Calculate a vector of x-values for plotting
#'
#' @param m.formula
#' model formula
#' @param fdata
#' dataframe with concentration values
#'
#' @return vector of concentration values for plotting
calc_xrange = function(m.formula,fdata) {
  concName = intersect(all.vars(m.formula),colnames(fdata))[2]
  concV = dplyr::pull(fdata,concName)
  loc = min(concV)
  hic = max(concV)
  fext = 0.2
  aext = (hic-loc) * fext
  x.start = loc - aext
  x.end = hic + aext

  x.values = seq(x.start,x.end,length.out=1E3)
}

#' Inverse hill function
#'
#' @param effect effect value
#' @param logEC50 parameter; log10(EC50)
#' @param slope parameter; slope
#'
#' @return concentration for given effect
inverse.hill = function(effect,logEC50,slope) {
  return(logEC50 - (log10((100/effect)-1)/slope))
}


#' Extends the plot data down to a lower log concentration
#'
#' The prediction bands are computed with the same method (delta or bootstrap)
#' and confidence level as the original fit. For the bootstrap, the stored
#' bootstrap parameter estimates are re-evaluated on the extended grid.
#'
#' @param fo bendr fit object
#' @param lox lowest log10 concentration of the extended grid; by default the
#' concentration at which the fitted effect is 99.9 percent. The existing range is never shrunk.
#' @param step grid spacing in log10 units
#' @param ci.method "delta" or "bootstrap"; by default the method of the fit
#'
#' @return dataframe
#' @export
extend_ci = function(fo, lox = NULL, step = 0.001, ci.method = NULL) {
  m.formula = fo$curve.fit$m$formula()
  level = if(is.null(fo$confidence.level)) 0.95 else fo$confidence.level
  if(is.null(ci.method)) ci.method = if(is.null(fo$ci.method)) "delta" else fo$ci.method
  if(is.null(lox)) {
    #get concentration at which effect is 99.9%
    lox = inverse.hill(99.9,log10(fo$ec50),fo$slope)
  }
  lox = min(lox, min(fo$plot.data$log.concentration))
  #generate sequence from new minimal log concentration value to previous max
  xtend_x = seq(from = lox,
                to = max(fo$plot.data$log.concentration),
                by = step)

  #predict effect values along new log concentration value range
  curve.predict.xtend = predict(fo$curve.fit,newdata = data.frame(logconc=xtend_x))

  if(ci.method == "bootstrap") {
    bp = fo$boot.params
    bp = bp[!is.na(bp$slope),]
    if(is.null(bp) || nrow(bp) == 0) {
      stop("extend_ci: no bootstrap parameter estimates stored in the fit object")
    }
    preds = t(sapply(seq_len(nrow(bp)), function(i) {
      eval(m.formula[[3]], list(logEC50 = bp$logEC50[i], slope = bp$slope[i], logconc = xtend_x))
    }))
    level_lower = (1-level)/2
    ci.lower = apply(preds,2,quantile,probs=level_lower)
    ci.upper = apply(preds,2,quantile,probs=1-level_lower)
  } else {
    #calculate confidence interval with the delta method
    ci.values.xtend = prediction_ci(m.formula,
                                    x.values = xtend_x,
                                    curvefit = fo$curve.fit,
                                    curve.predict = curve.predict.xtend,
                                    dof = nrow(fo$data)-2,
                                    level = level)
    ci.lower = curve.predict.xtend - ci.values.xtend
    ci.upper = curve.predict.xtend + ci.values.xtend
  }

  return(data.frame(
    curve.predict = curve.predict.xtend,
    log.concentration = xtend_x,
    ci.values.upper = ci.upper,
    ci.values.lower = ci.lower
    ))

}


#' calculate ECx
#'
#' @param x desired effective concentration level (ECx)
#' @param slope slope of the hill curve
#' @param ec50 Effective concentration at 50 percent
#'
#' @return ECx value
#' @export
#'
#' @examples
#' x = 10
#' slope = 2
#' ec50 = 12
#' calc_ecx(x,slope,ec50)
calc_ecx = function(x,slope,ec50) {
  ecx = ((x/(100-x))^(1/slope)) * ec50
  names(ecx) = paste0("ec",as.character(x))
  ecx
}

#' Remove NA columns from dataframe and convert to long format
#'
#' @param df dataframe with dose response data. The first column needs to be the concentration.
#'
#' @return long dataframe
#' @export
cleanAndPivot = function(df) {
  if(is.na(df) || nrow(df) == 0 || ncol(df) <2) {
    return(NA)
  }

  df = df[colSums(!is.na(df)) > 0]
  nrep = ncol(df)-1



  if(nrep > 1) {
    repnames= sapply(as.character(seq(nrep)),function(x) paste0("replicate",x))
    colnames(df) = c("conc",repnames)
    df = tidyr::pivot_longer(df,2:ncol(df),names_to = "replicateID",values_to = "effect")
  } else {
    colnames(df) = c("conc","effect")
  }

  df = df  %>% dplyr::filter(conc > 0) %>% mutate(logconc = log10(conc)) %>% drop_na()
  if(ncol(df) <3) {
    return(NA)
  } else {
    return(df)
  }

}


#' Calculate quantile difference of a vector
#'
#' @param x vector of numeric values
#' @param lo lower quantile [0,1]
#' @param hi upper quantile [0,1]
#'
#' @return difference between lo and hi
#' @export
#'
#' @examples
#' qd = qdiff(seq(200),0.05,0.95)
#'
qdiff = function(x,lo,hi) {
  if(length(x) < 2 && is.na(x)) {
    NA
  }
  qlo = quantile(x,lo)
  qhi = quantile(x,hi)
  qdiff = qhi - qlo
  qdiff
}

#' Guess the logEC50 by fitting a linear model through the datapoints closest to 50 percent effect concentration
#'
#' @param df dataframe containing effect and log concentration columns
#'
#' @return guess for EC50 value
#' @export
guessEC50_lm = function(df) {
  df = df %>% arrange(logconc)
  if("replicateID" %in% colnames(df)) {
    df = df %>% group_by(logconc) %>% summarise(effect = mean(effect),.groups = 'drop')
  }
  eff = df %>% pull(effect)
  maxOver = max(df[df$effect > 50,]$logconc)
  #minOver = min(eff[eff > 50])
  maxOverIdx = which(df$logconc == maxOver)

  minUnder = min(df[df$effect < 50,]$logconc)
  minUnderIdx = which(df$logconc == minUnder)
  dfs = df[c(maxOverIdx,minUnderIdx),]
  coe = lm(effect ~ logconc,dfs)$coefficients
  guess = (50 -coe[1])/coe[2]
  as.double(guess)
}







