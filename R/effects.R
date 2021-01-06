
#' Calculate effects for plotting
#'
#' @title Adapts package \code{effects}
#'
#' @method Effect fitTMB
#' @export
Effect.fitTMB <- function (focal.predictors, mod, which_formula="taxis", ...) {

  #assign("which_formula", 100, envir=.GlobalEnv)

  # Basic
  if( which_formula=="taxis" ){
    formula_orig = mod$formula_taxis
      formula_full = update.formula(formula_orig,y~.+0)
    whichnum = which(names(mod$parameter_estimates$par)=="alpha_logit_ratio_k")
    mod$parhat = mod$parameter_estimates$par[whichnum]
      mod$covhat = mod$parameter_estimates$SD$cov.fixed[whichnum,whichnum]
    mod$parhat = mod$parhat[ mod$Obj$env$map$alpha_logit_ratio_k ]
      mod$covhat = mod$covhat[ mod$Obj$env$map$alpha_logit_ratio_k, mod$Obj$env$map$alpha_logit_ratio_k ]
    mod$parhat = ifelse(is.na(mod$parhat),0,mod$parhat)
      mod$covhat = ifelse(is.na(mod$covhat),0,mod$covhat)
    names(mod$parhat)[] = "alpha_logit_ratio_k"
      rownames(mod$covhat) = colnames(mod$covhat) = names(mod$parhat)
    mod$call = mod$effects$call_taxis
  }else{
    stop("Other options are not implemented")
  }

  # Augment stuff
  mod$coefficients = mod$parhat
  mod$vcov = mod$covhat
  mod$formula = formula_full
  mod$family = gaussian(link = "identity")

  # Functions for package
  family.fitTMB = function(x,...) x$family
  vcov.fitTMB = function(x,...) x$vcov

  # Do stuff
  fam <- mod$family
  ## dummy functions to make Effect.default work
  dummyfuns <- list(variance=function(mu) mu,
                    initialize=expression(mustart <- y + 0.1),
                    dev.resids=function(...) poisson()$dev.res(...) )
  for (i in names(dummyfuns)) {
    if (is.null(fam[[i]])) fam[[i]] <- dummyfuns[[i]]
  }
  ## allow calculation of effects ...
  if (length(formals(fam$variance))>1) {
    warning("overriding variance function for effects: computed variances may be incorrect")
    fam$variance <- dummyfuns$variance
  }
  args <- list(call = mod$call,
               coefficients = mod$coefficients,
               vcov = mod$vcov,
               family = fam,
               formula = formula_full)
  if(!requireNamespace("effects")) stop("please install the effects package")
  effects::Effect.default(focal.predictors,
    mod,
    ...,
    sources = args)
}
