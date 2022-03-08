
#' @title
#' Simulate new dynamics and sampling data
#'
#' @description
#' \code{simulate_data} conducts a parametric bootstrap to simulate new data and potentially simulate new population dynamics and associated variables
#'
#' Simulate new data given various potential procedures to propagate uncertainty about parameters.
#'
#' @param fit output form \code{fit_model(.)}
#' @param type integer stating what type of simulation to use from the following options:
#' \itemize{
#' \item \code{type=1} is a "measurement error" or "conditional" simulator that simulates new data conditional upon estimated fixed and random effects.
#' \item \code{type=2} is an "unconditional" simulator that simulates new random effects conditional upon fixed effects
#' (but not otherwise conditioning upon original data), and new data conditional upon both.
#' \item \code{type=3} simulates new fixed and random effects from the joint precision matrix (i.e., conditioning upon the original data), and new data conditional upon these values.
#' \item \code{type=4} simulates new random effects from the internal Hessian matrix evaluated at the MLE (i.e., conditional on fixed effects estimates and the original data),
#' and new data conditional upon these values.
#' }
#' @param random_seed integer passed to \code{set.seed}, where the default value \code{random_seed=NULL} resets the random-number seed.
#' @param parvec vector of parameters to use as fixed values (or use as mean when simulating new parameters),
#'        to be used when simulating new data.
#'        default \code{par=fit$Obj$env$last.par.best} uses the best previous value,
#'        i.e., MLE if the model is optimized.
#'

#' @return Report object containing new data and population variables including
#' \describe{
#'   \item{b_i}{New simulated data}
#'   \item{D_gcy}{Density for each grid cell g, category c, and year y}
#'   \item{Index_cyl}{Index of abundance for each category c, year y, and stratum l}
#' }

#' @export
simulate_data <-
function( fit,
          type = 1,
          random_seed = NULL,
          parvec = fit$Obj$env$last.par.best ){

  # Sample from GMRF using sparse precision
  rmvnorm_prec <- function(mu, prec, n.sims, random_seed ) {
    set.seed( random_seed )
    z = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
    L = Matrix::Cholesky(prec, super=TRUE)
    z = Matrix::solve(L, z, system = "Lt") ## z = Lt^-1 %*% z
    z = Matrix::solve(L, z, system = "Pt") ## z = Pt    %*% z
    z = as.matrix(z)
    return(mu + z)
  }

  # Check for loaded VAST
  # Modified from TMB:::getUserDLL
  dlls <- getLoadedDLLs()
  isTMBdll <- function(dll) !is(try(getNativeSymbolInfo("MakeADFunObject",dll), TRUE), "try-error")
  TMBdll <- sapply(dlls, isTMBdll)
  if( sum(TMBdll)==0 ){
    stop("VAST is not linked as a DLL, so `simulate_data` will not work.
    Please re-run model (potentially from informative starting values to save time) to use `simulate_data`")
  }else if(sum(TMBdll)>=2){
    warning("VAST is linked to multiple DLLs. Please consider using dyn.unload() to unload
    earlier VAST runs to avoid potentially ambiguous behavior when running `simulate_data`")
  }

  # Extract stuff
  Obj = fit$Obj
  simulate_random = Obj$env$data$simulate_random
  if(length(parvec)!=length(fit$Obj$env$last.par.best)) stop("Check supplied `par`")

  # Revert settings when done
  revert_settings = function(simulate_random){Obj$env$data$simulate_random = simulate_random}
  on.exit( revert_settings(simulate_random) )

  # Simulate conditional upon fixed and random effect estimates
  if( type==1 ){
    # Change and revert settings
    set.seed(random_seed)
    Obj$env$data$simulate_random = 0
    Return = Obj$simulate( par=parvec, complete=TRUE )
  }

  # Simulate new random effects and data
  if( type==2 ){
    set.seed(random_seed)
    Obj$env$data$simulate_random = 1
    Return = Obj$simulate( par=parvec, complete=TRUE )
  }

  # Simulate from predictive distribution of fixed AND random effects, and then new data
  # Could instead explore fit$Obj$env$MC(.) for sampling-importance-resampling approach
  if( type==3 ){
    # Informative error messages
    if( !("jointPrecision" %in% names(fit$parameter_estimates$SD)) ){
      stop("jointPrecision not present in fit$parameter_estimates$SD; please re-run with `getJointPrecision=TRUE`")
    }

    # Sample from joint distribution
    newpar = rmvnorm_prec( mu = parvec,
                           prec = fit$parameter_estimates$SD$jointPrecision,
                           n.sims = 1,
                           random_seed = random_seed )[,1]

    # Simulate
    Obj$env$data$simulate_random = 0
    Return = Obj$simulate( par = newpar, complete=TRUE )
  }

  # Simulate from predictive distribution of random effects and NOT fixed effects, and then new data
  if( type==4 ){
    set.seed( random_seed )
    warning( "Type-4 residuals are still under development, please use with care and note that they may change at any point.")
    newpar = parvec
    MC = Obj$env$MC( keep = TRUE,
                     n = 1,
                     antithetic = FALSE )
    newpar[Obj$env$random] = attr(MC, "samples")

    # Simulate
    Obj$env$data$simulate_random = 0
    Return = Obj$simulate( par=newpar, complete=TRUE )
  }

  # Fix tag
  #  These aren't simulated in the CPP due to not having a rmultinom function in TMB
  for( i in 1:nrow(Return$conventional_hz) ){
    Return$conventional_hz[i,2] = sample( x=1:ncol(Return$like_conventional_hg), size=1, prob=Return$like_conventional_hg[i,] )
  }
  Return$like_conventional_hg[i,2] = Return$like_conventional_hg[i,2] - 1 # Must subtract 1 to do CPP indexing!
  if( nrow(Return$satellite_iz)>0 ) stop("Not implemented")

  # Return
  return( Return )
}


