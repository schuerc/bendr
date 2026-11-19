#' two stage nonlinear least squares fit
#'
#' @param m.formula
#' model formula
#' @param start
#' vector of start parameters
#' @param fdata
#' data frame of data to fit
#'
#' @return
#' nls2 fit object
#' @export
nlsfit = function(m.formula,start,fdata) {
  #fit1 = nls2::nls2(m.formula, data = fdata,start = start, algorithm = "plinear-brute")
  #curvefit <- nls2::nls2(m.formula, data = fdata, start = coef(fit1)[1:2],
  #                 algorithm = "brute-force")
  #fit1 = nlmrt::nlxb(m.formula,
  #                   start = start,
  #                   trace = F,
  #                   data = fdata)
  #curvefit <- nls2::nls2(m.formula, data = fdata, start = fit1$coefficients,
  #                      algorithm = "brute-force")
  #curvefit <- minpack.lm::nlsLM(m.formula, data = fdata, start = start)
  curvefit = tryCatch(
    {
      minpack.lm::nlsLM(m.formula,fdata,start=start,trace=F)

    }, error = function(cond) {
      warning("Switching to grid-search")
      inits <- data.frame(slope = sort(c(0.5, 100)*start["slope"]),
                          logEC50 = sort(c(0.5, 3)*start["logEC50"]))
      nls2::nls2(m.formula, data = fdata, start = inits,
                 algorithm = "grid-search",
                 trace=F,
                 control=list(maxiter=100))
    }
  )

  return(curvefit)
}


#' Calculate the confidence intervals using the delta method
#'
#'
#'
#'
#' @param m.formula
#' model formula
#' @param x.values
#' concentration values at which the curve will be plotted
#' @param curvefit
#' model fit object from nonlinear least squares
#' @param curve.predict
#' predicted effect values at the concentrations specified in x.values
#' @param dof
#' degrees of freedom for the inverse student's t distribution ; default is (# measurements)-(number of parameters)
#' @param level
#' Level of confidence interval to use (0.95 by default). [0,1]
#' @param verbose
#' Print and plot intermediate results for verboseging
#'
#' @return
#' Confidence interval values for all concentrations in x.values
#' @export
prediction_ci = function(m.formula,x.values,curvefit,curve.predict,dof,level = 0.95,verbose=FALSE) {

  logEC50.value = coef(curvefit)['logEC50']
  slope.value = coef(curvefit)['slope']

  # get the covariance matrix for the parameters
  covmatrix= vcov(curvefit)

  # get the jacobian for all x.values
  residual.jacfun = nlmrt::model2jacfun(m.formula,coef(curvefit))
  #nameY = all.vars(m.formula)['effect']
  #nameX = all.vars(m.formula)['logconc']

  paramY = list(curve.predict)
  paramX = list(x.values)
  #names(paramY) = nameY
  #names(paramX) = nameX
  names(paramY) = "effect"
  names(paramX) = "logconc"
  residual.jacobian = residual.jacfun(coef(curvefit), paramX, paramY)
  if(verbose) {
    jac_plot = ggplot(as_tibble(residual.jacobian)) + geom_point(aes(x=logEC50,y=slope))
    ggsave("./residual_jacobian.png")
  }

  # get the confidence interval for all x_values
  # only the diagonal of J V J' is needed: rowSums((J V) * J) gives it without building the
  # n x n matrix (n = grid points, up to ~1e4 in extend_ci/findNTC -> ~800 MB before)
  error = rowSums((residual.jacobian %*% covmatrix) * residual.jacobian)
  error = sqrt(error)

  # t for the inverse studen't distribution depends on the number of measurements
  p = 1- ((1-level)/2)
  t_student = qt(p, df =dof)
  ci.values = error*t_student

  ci.values
}



#' Calculate parameter and prediction confidence intervals using the bootstrap
#'
#' @param m.formula model formula
#' @param data data for fitting
#' @param logC.values log concentration values for plotting
#' @param curvefit initial model fit
#' @param curve.predict predicted effect values from model fit
#' @param n.boot number of bootstrap repetitions
#' @param level confidence level
#' @param verbose print intermediate results
#'
#' Replicates are resampled within each concentration level, so every bootstrap
#' sample contains all tested concentrations. Each resample is fitted with
#' Levenberg-Marquardt starting from the original estimates; failed fits are
#' recorded as NA and excluded from the quantiles.
#'
#' @return list: bootstrapped confidence intervals for prediction parameters and prediction bands,
#' the raw bootstrap parameter estimates (boot.params) and the number of failed fits (n.failed)
bootstrap_ci = function(m.formula,data,logC.values,curvefit,curve.predict,n.boot=1000,level = 0.95,verbose=FALSE) {
  start = coef(curvefit)[c('logEC50','slope')]

  ## arrays to hold predictions and parameters (bootstrap resamples only)
  predictions <- matrix(NA, ncol=length(logC.values), nrow=n.boot)
  parameters <- data.frame(slope=rep(NA_real_, n.boot), logEC50=NA_real_, EC50=NA_real_, EC10=NA_real_)

  ## row indices grouped by concentration, for stratified resampling
  conc.groups = split(seq_len(nrow(data)), data$logconc)

  ## --- fit to bootstrap samples
  for(i in seq_len(n.boot)){
    ## resample replicates within each concentration level
    idx = unlist(lapply(conc.groups, function(g) g[sample.int(length(g), replace = TRUE)]),
                 use.names = FALSE)
    data2 <- data[idx,]

    ## fit model to resampled data; a failed fit leaves the row NA
    fit <- tryCatch(minpack.lm::nlsLM(m.formula, data = data2, start = start),
                    error = function(cond) NULL)
    if(is.null(fit)) next

    ## store predictions and parameter values (explicit names: coef() order depends on start)
    predictions[i,] <- predict(fit, newdata=list(logconc=logC.values))
    slope = unname(coef(fit)['slope'])
    logEC50 = unname(coef(fit)['logEC50'])
    EC50 = 10^logEC50
    # EC10 = 10% effect = 90% survival
    EC10 = unname(calc_ecx(90, slope, EC50))
    parameters[i,] <- c(slope, logEC50, EC50, EC10)
  }

  n.failed = sum(is.na(parameters$slope))
  if(n.failed > 0) {
    warning(sprintf("%d of %d bootstrap fits failed and were excluded", n.failed, n.boot))
  }

  level_lower = (1-level)/2
  level_upper = 1- level_lower

  param_ci = apply(parameters,2,quantile,probs=c(level_lower,level_upper),na.rm=TRUE)
  rownames(param_ci) = c("lower","upper")
  pred_ci =  apply(predictions,2,quantile,probs=c(level_lower,level_upper),na.rm=TRUE)
  rownames(pred_ci) = c("lower","upper")
  return(list(params=param_ci,preds=pred_ci,boot.params=parameters,n.failed=n.failed))
}




#' Dose Response curve fit
#'
#' @param m.formula
#' Model formula.
#' @param fdata
#' Data to fit the model to.
#' @param level
#' Level of confidence interval to use (0.95 by default). [0,1]
#' @param ci_method
#' Method for computing confidence intervals ["delta","bootstrap"] default: delta
#' @param start
#' initial parameter values for the dose response fitting. c(logEC50,slope)
#' @param verbose
#' logical, prints intermediate results from fitting and print and plot intermediate results from confidence interval calculation
#'
#'
#' @return list: fit object
#' @importFrom nlstools confint2
#' @export
fitdr = function(m.formula,fdata,level=0.95,ci_method="delta",start=vector(),verbose=FALSE) {
  #responseName = all.vars(m.formula)['effect']
  #concName = all.vars(m.formula)['logconc']
  responseName = "effect"
  concName = "logconc"
  if(length(start) == 0) {
    #start = c(logEC50 = median(fdata$logconc), slope = diff(range(fdata$logconc))/ diff(range(fdata$effect)))
    responseV = fdata %>% dplyr::pull(responseName)
    dependentV = fdata %>% dplyr::pull(concName)

    #set slopeGuess to +1 if the response values are ascending and -1 if they are descending
    #slopeGuess = sign(mean(diff(responseV)))
    slopeGuess = as.double(sign(responseV[length(responseV)]-responseV[1]))
    #if the concentration values are descending in the dataframe, negate the slope guess
    #if(mean(diff(dependentV)) < 0)
    if(dependentV[length(dependentV)] - dependentV[1] < 0 ) {
      slopeGuess = - slopeGuess
    }
    # if the majority of the tested concentration range is close to 100 or 0 % effect,
    # use a linear fit through the measured points that are closest to 50% effect, instead of the median
    meanEffect = mean(dplyr::pull(fdata,responseName))
    repMeans = fdata %>% group_by(conc) %>%
      summarise(effect = mean(effect),.groups="drop") %>% pull(effect)
    if((any(repMeans > 50) && any(repMeans < 50 )) &&
       (meanEffect > 80 || meanEffect < 20)) {
      if(verbose) {
        "Using linear model for EC50 guess"
      }
      logEC50Guess = guessEC50_lm(fdata)
    } else {
      logEC50Guess = median(dplyr::pull(fdata,concName))
    }

    start = c(logEC50 = logEC50Guess,slope=slopeGuess)


    if(verbose) {
      print("Initial values are")
      print(start)
    }

  }
  curvefit = nlsfit(m.formula,start=start,fdata = fdata)
  x.values = calc_xrange(m.formula,fdata)
  xdf = data.frame(x.values)
  names(xdf) = concName
  curve.predict = predict(curvefit, newdata = xdf)
  logec50 = coefficients(curvefit)['logEC50']
  slope = coefficients(curvefit)['slope']
  ec50 = 10^logec50
  ec10 = calc_ecx(90,slope,ec50)

  #define default values for confidence intervals
  slope.ci = NA
  ec50.ci = NA
  ec10.ci = NA
  ci.values.lower = NA
  ci.values.upper = NA
  boot.params = NULL


  if (ci_method == "delta") {
    tryCatch({
      ci.values = prediction_ci(m.formula,x.values,curvefit,curve.predict,dof=nrow(fdata)-2,level=level,verbose=verbose)
      ci.values.lower = curve.predict - ci.values
      ci.values.upper = curve.predict + ci.values
      },
    error= function(cond) {
      message("Failed to compute prediction bands")
      message("Original error message:")
      message(cond)

    }
      )

    #plot.data = data.frame(log.concentration = x.values,curve.predict = curve.predict,ci.values = ci.values)

    #ec10 = calc_ecx(10,slope,ec50)
    #logec10 = log10(ec10)
    tryCatch({
      ci.par = confint2(curvefit)
      colnames(ci.par) = c('lower','upper')
      ec50.ci = 10^ci.par['logEC50',]
      slope.ci = ci.par['slope',]

      # same curve, reparametrised so that logEC10 (90% survival) is a model parameter:
      # at logconc = logEC10 the effect is 100/(1 + 1/9) = 90
      m.formula.ec10 = effect ~ 100 / (1 + 10^(((logEC10-logconc) * slope) - log10(0.9/0.1)))

      start.ec10 = c(logEC10 = as.double(log10(ec10)),slope=as.double(slope))
      curvefit.ec10 = minpack.lm::nlsLM(m.formula.ec10, data = fdata, start = start.ec10)
      logec10 = coefficients(curvefit.ec10)['logEC10']
      ec10 = 10^logec10
      ci.par.ec10 = confint2(curvefit.ec10)
      colnames(ci.par.ec10) = c('lower','upper')
      ec10.ci = 10^ci.par.ec10['logEC10',]



    },
    error= function(cond) {
      message("Failed to compute parameter confidence intervals")
      message("Original error message:")
      message(cond)
    }
    )
    names(ec50) = c('')
    names(ec10) = c('')
    names(slope) = c('')

    plot.data = data.frame(log.concentration = x.values,curve.predict = curve.predict,ci.values.upper = ci.values.upper,
                           ci.values.lower = ci.values.lower)

  } else if (ci_method == "bootstrap") {

    boot = bootstrap_ci(m.formula,data=fdata,x.values,curvefit,curve.predict,n.boot = 1000,level = level)
    boot.params = boot$boot.params
    slope.ci = boot$params[,"slope"]
    ec50.ci = boot$params[,"EC50"]
    ec10.ci = boot$params[,"EC10"]
    ci.values.lower = boot$preds[1,]
    ci.values.upper = boot$preds[2,]
    plot.data = data.frame(log.concentration = x.values,curve.predict = curve.predict,ci.values.upper = ci.values.upper,
                           ci.values.lower = ci.values.lower)
  } else {
    print("Invalid confidence interval calculation method")
    return(NA)
  }




  list(plot.data = plot.data ,
       curve.fit=curvefit,
       confidence.level=level,
       ci.method=ci_method,
       boot.params=boot.params,
       data=fdata,
       ec50=ec50,
       ec50.ci = ec50.ci,
       ec10 = ec10,
       ec10.ci = ec10.ci,
       slope = slope,
       slope.ci = slope.ci,
       aic = AIC(curvefit),
       xname = concName,
       yname = responseName
  )
}



#' Dose response curve fit with replicates
#'
#' @param m.formula model formula
#' @param fdata dataframe with data to fit
#' @param effectColumns effect column names or indices in data frame
#' @param concColumn concentration column name in dataframe
#' @param level confidence level
#' @param start starting values
#' @param verbose print fitting intermediate results and confidence interval calculation intermediate  results
#' @param ci_method method for confidence interval calculation, one of ["delta","bootstrap"]
#' @param replicates what the replicate columns are:
#' "independent" (separate experiments): the curve is fitted to all replicate values, df = n_values - 2;
#' "same_plate" (wells on one plate, pseudo-replicates): the curve is fitted to the concentration
#' means, df = n_concentrations - 2. Bootstrap CIs need independent replicates.
#'
#' @return list: fit object
#' @importFrom tidyr drop_na
#' @export
fitdr_replicates= function(m.formula,fdata,effectColumns,concColumn,level=0.95,start=vector(),ci_method="delta",verbose=F,
                           replicates=c("independent","same_plate")) {
  replicates = match.arg(replicates)
  if(replicates == "same_plate" && ci_method == "bootstrap") {
    stop("Bootstrap confidence intervals need independent replicates; use ci_method = \"delta\" for replicates from the same plate")
  }
  quo_concColumn = enquo(concColumn)
  data.logged = fdata %>% dplyr::mutate(logconc := log10(!! quo_concColumn))
  data.long = data.logged %>% tidyr::gather(replicateID,effect,effectColumns) %>%
    drop_na() %>%
    arrange(desc(logconc))
  data.for.fit = data.long
  if(replicates == "same_plate") {
    data.for.fit = data.long %>%
      dplyr::group_by(!! quo_concColumn, logconc) %>%
      dplyr::summarise(effect = mean(effect), .groups = "drop") %>%
      arrange(desc(logconc))
  }
  data.fit = fitdr(m.formula,data.for.fit,level=level,start=start,ci_method=ci_method,verbose=verbose)
  data.fit$replicates = replicates
  data.fit$data.long = data.long
  data.fit$nreplicates = length(effectColumns)
  data.fit
}



#' Fit multiple dose response curves
#'
#' @param m.formula
#' model formula
#' @param fdata
#' data frame with one column corresponding to the (non-logged) concentration values
#' and an arbitrary number of columns corresponding to effects
#' @param effectColumns
#' vector of effect columns; can be either numeric indices or column names
#' @param concColumn
#' name of the concentration column; this parameter is lazily evaluated, so do not use quotes
#'
#' @return data frame: fit objects and corresponding summaries
#' @import rlang
#' @import dplyr
#' @export
fitdr_multi= function(m.formula,fdata,effectColumns,concColumn) {
  quo_concColumn = enquo(concColumn)
  data.logged = fdata %>% dplyr::mutate(logconc := log10(!! quo_concColumn))
  data.long = data.logged %>% tidyr::gather(sampleID,effect,effectColumns)
  data.fit = data.long %>% dplyr::group_by(sampleID) %>% tidyr::nest() %>%
    dplyr::mutate(m.fit = purrr::map(data,fitdr,m.formula = m.formula)) %>%
    dplyr::mutate(
      ec50=purrr::map_dbl(m.fit,function(mod) 10^(coefficients(mod$curve.fit)['logEC50'])),
      slope=purrr::map_dbl(m.fit,function(mod) coefficients(mod$curve.fit)['slope']),
      aic=purrr::map_dbl(m.fit,function(mod) mod$aic))
  data.fit
}
