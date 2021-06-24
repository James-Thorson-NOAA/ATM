
#' Fit ATM model in TMB
#'
#' \code{fitTMB} fits the Advection-taxis model in TMB
#'
#' @return Object of class \code{fitTMB}, containing parameter estimates and predicted time-integrated movement
#'
#' @author James Thorson
#' @export
fitTMB <-
function( Cov_stars,
      formula_taxis,
      formula_diffusion,
      coords_gz,
      #t_uy,
      satellite_iz = NULL,
      conventional_hz = NULL,
      survey_jz = NULL,
      fishery_fz = NULL,
      E_guy = NULL,
      uy_tz = NULL,
      duration_u = NULL,
      cpp_version = FishStatsUtils::get_latest_version(package="ATM"),
      tmb_dir = system.file("executables",package="ATM"),
      run_dir = getwd(),
      compile_dir = tmb_dir,
      log2steps = 20,
      sigma2 = 4^2,
      spde_aniso = NULL,
      use_REML = TRUE,
      build_model = TRUE,
      run_model = TRUE,
      start_param_list = NULL,
      alpha_ratio_bounds = 1,
      diffusion_bounds = 0,
      movement_penalty = 0,
      constant_tail_probability = 1e-8,
      ... ){

  data_list = make_data( Cov_stars = Cov_stars,
    formula_taxis = formula_taxis,
    formula_diffusion = formula_diffusion,
    coords_gz = coords_gz,
    uy_tz = uy_tz,
    satellite_iz = satellite_iz,
    conventional_hz = conventional_hz,
    survey_jz = survey_jz,
    fishery_fz = fishery_fz,
    E_guy = E_guy,
    duration_u = duration_u,
    cpp_version = cpp_version,
    log2steps = log2steps,
    spde_aniso = spde_aniso,
    constant_tail_probability = constant_tail_probability,
    alpha_ratio_bounds = alpha_ratio_bounds,
    diffusion_bounds = diffusion_bounds,
    movement_penalty = movement_penalty
  )

  # Make parameters
  rnorm_array = function( dim, mean=0, sd=0.1 ){
    array( rnorm(prod(dim),mean=mean,sd=sd), dim=dim )
  }
  if( cpp_version %in% c("ATM_v1_0_0","ATM_v0_9_0") ){
    param_list = list(
      "ln_sigma" = log(sqrt(sigma2)),
      "alpha_logit_ratio_k" = 0.01 * rnorm(dim(data_list$X_guyk)[4])
    )
  }
  if( cpp_version %in% c("ATM_v2_0_0") ){
    param_list = list(
      "ln_sigma" = log(sqrt(sigma2)),
      "alpha_logit_ratio_k" = 0.01 * rnorm(dim(data_list$X_guyk)[4]),
      "ln_H_input" = c(0,0),
      "ln_kappa" = -3,
      "ln_sigma_omega" = log(1),
      "ln_sigma_epsilon" = log(1),
      "ln_phi" = log(1),
      "power_prime" = qlogis( 1.5 - 1),
      "Beta_t" = rep(0, nrow(data_list$uy_tz)),
      "Omegainput_s" = rnorm_array(spatial_list$n_s),
      "Epsiloninput_st" = rnorm_array( c(spatial_list$n_s,nrow(data_list$uy_tz)) )
    )
  }
  if( cpp_version %in% c("ATM_v3_0_0") ){
    param_list = list(
      "ln_sigma" = log(sqrt(sigma2)),
      "alpha_logit_ratio_k" = 0.01 * rnorm(dim(data_list$X_guyk)[4]),
      "ln_H_input" = c(0,0),
      "ln_kappa" = -3,
      "ln_sigma_epsilon0" = log(3),
      "ln_sigma_epsilon" = log(1),
      "ln_phi" = log(1),
      "power_prime" = qlogis( 1.5 - 1 ),
      "Beta_t" = rep(0, nrow(data_list$uy_tz)),
      "ln_d_st" = rnorm_array( c(spatial_list$n_s,nrow(data_list$uy_tz)) )
    )
  }
  if( cpp_version %in% c("ATM_v4_0_0") ){
    param_list = list(
      "ln_sigma_l" = c( log(sqrt(sigma2)), rep(0,dim(data_list$Z_guyl)[4]-1) ),
      "alpha_logit_ratio_k" = 0.01 * rnorm(dim(data_list$X_guyk)[4]),
      "ln_H_input" = c(0,0),
      "ln_kappa" = -3,
      "ln_sigma_epsilon0" = log(3),
      "ln_sigma_epsilon" = log(1),
      "ln_phi" = log(1),
      "power_prime" = qlogis( 1.5 - 1 ),
      "Beta_t" = rep(0, nrow(data_list$uy_tz)),
      "ln_d_st" = rnorm_array( c(spatial_list$n_s,nrow(data_list$uy_tz)) )
    )
  }
  if( cpp_version %in% c("ATM_v6_0_0","ATM_v5_0_0") ){
    param_list = list(
      "ln_sigma_l" = c( log(sqrt(sigma2)), rep(0,dim(data_list$Z_guyl)[4]-1) ),
      "alpha_logit_ratio_k" = 0.01 * rnorm(dim(data_list$X_guyk)[4]),
      "ln_H_input" = c(0,0),
      "ln_kappa" = -3,
      "ln_sigma_epsilon0" = log(3),
      "ln_sigma_epsilon" = log(1),
      "ln_phi" = log(1),
      "power_prime" = qlogis( 1.5 - 1 ),
      "ln_CV" = log(1),
      "lambda" = 0,
      "Beta_t" = rep(0, nrow(data_list$uy_tz)),
      "ln_d_st" = rnorm_array( c(spatial_list$n_s,nrow(data_list$uy_tz)) )
    )
  }

  # User supplied start values
  if( !is.null(start_param_list) ){
    if( all(sapply(start_param_list,length) == sapply(param_list,length)) ){
      param_list = start_param_list
    }else{
      stop("Check `user_param_list` for issues")
    }
  }

  # Which map
  map = NULL
  # Map off betas for years without survey data
  if( "Beta_t" %in% names(param_list) ){
    Beta_t = 1:length(param_list$Beta_t) - 1
    Beta_t = ifelse( Beta_t %in% c(data_list$t_j,data_list$t_f), Beta_t, NA )
    map$Beta_t = factor(Beta_t)
  }
  # Map off density parameters if no survey data AND no fishery
  if( length(data_list$b_j)==0 & length(data_list$b_f)==0 ){
    if("ln_H_input"%in%names(param_list)) map$ln_H_input = factor(c(NA,NA))
    if("ln_kappa"%in%names(param_list)) map$ln_kappa = factor(NA)
    if("ln_sigma_epsilon0"%in%names(param_list)) map$ln_sigma_epsilon0 = factor(NA)
    if("ln_sigma_epsilon"%in%names(param_list)) map$ln_sigma_epsilon = factor(NA)
    if("ln_d_st"%in%names(param_list)) map$ln_d_st = factor(array(NA,dim=dim(param_list$ln_d_st)))
  }
  # Map off survey measurement-error if no survey data
  if( length(data_list$b_j) == 0 ){
    if("ln_phi"%in%names(param_list)) map$ln_phi = factor(NA)
    if("power_prime"%in%names(param_list)) map$power_prime = factor(NA)
  }
  # Map off fishery measurement-error if no fishery data
  if( length(data_list$b_f) == 0 ){
    if("ln_CV"%in%names(param_list)) map$ln_CV = factor(NA)
  }
  # Map of catchability ratio if no fishery or no survey
  if( length(data_list$b_f)==0 | length(data_list$b_j)==0){
    if("lambda"%in%names(param_list)) map$lambda = factor(NA)
  }
  # Map off density effects for the intercept
  if( any(dimnames(data_list$X_guyk)[[4]]=="(Intercept)") ){
    map$alpha_logit_ratio_k[which(dimnames(data_list$X_guyk)[[4]]=="(Intercept)")] = NA
  }
  # Map off constants in X_guyk
  turnoff_k = apply( data_list$X_guyk, MARGIN=4, FUN=function(x){var(as.vector(x))==0} )
  if( any(turnoff_k) ){
    map$alpha_logit_ratio_k = 1:length(param_list$alpha_logit_ratio_k)
    map$alpha_logit_ratio_k[which(turnoff_k)] = NA
    map$alpha_logit_ratio_k = factor(map$alpha_logit_ratio_k)
    param_list$alpha_logit_ratio_k[which(dimnames(data_list$X_guyk)[[4]]=="(Intercept)")] = 0
  }
  # Map off constants in Z_guyl except first term
  turnoff_l = apply( data_list$Z_guyl, MARGIN=4, FUN=function(x){var(as.vector(x))==0} )
  turnoff_l[1] = FALSE
  if( any(turnoff_l) ){
    map$ln_sigma_l = 1:length(param_list$ln_sigma_l)
    map$ln_sigma_l[which(turnoff_l)] = NA
    map$ln_sigma_l = factor(map$ln_sigma_l)
    param_list$ln_sigma_l[which(turnoff_l)] = 0
  }

  #
  random = c("Omegainput_s", "Epsiloninput_st", "ln_d_st")
  if( use_REML==TRUE ){
    random = union( random, c("Beta_t","ln_phi","lambda") ) # ,"alpha_logit_ratio_k","ln_sigma","ln_sigma_l"
  }
  random = random[which(random %in% names(param_list))]
  if( length(random)==0) random = NULL
  #random = NULL

  # Return inputs
  Return = list( "data_list"=data_list, "param_list"=param_list, "random"=random )

  # Optionally build and run
  if( build_model == TRUE ){
    # Compile
    # dyn.unload( paste0(compile_dir,"/",TMB::dynlib(TMB:::getUserDLL())) )
    file.copy( from=paste0(tmb_dir,"/",cpp_version,".cpp"), to=paste0(compile_dir,"/",cpp_version,".cpp"), overwrite=FALSE)
    origwd = getwd()
    on.exit(setwd(origwd),add=TRUE)
    setwd( compile_dir )
    TMB::compile( paste0(cpp_version,".cpp"), CPPFLAGS="-Wno-ignored-attributes" )

    # Build object
    dyn.load( paste0(compile_dir,"/",TMB::dynlib(cpp_version)) ) # random=random,
    Return$Obj = TMB::MakeADFun( data=data_list, parameters=param_list, hessian=FALSE, map=map, random=random, DLL=cpp_version )  #
    # Report = Obj$report()

    # Print number of parameters
    ThorsonUtilities::list_parameters( Return$Obj )

    # Optionally run model
    if( run_model == TRUE ){
      Return$Obj$env$beSilent()
      Return$parameter_estimates = TMBhelper::fit_tmb( Return$Obj, control=list(trace=1),
        savedir=run_dir, getReportCovariance=TRUE, ... ) #
      Return$parhat = Return$Obj$env$parList( Return$parameter_estimates$par )
    }

    # Extract stuff
    Return$Report = Return$Obj$report( )
  }

  # Add stuff for effects package integration
  Data_zp = cbind( format_covariates( Cov_stars, formula=formula_taxis )$Data_zp, "dummy_zeroes"=0)
  formula_taxis_full = update.formula(formula_taxis, dummy_zeroes~.+0)
  formula_diffusion_full = update.formula(formula_diffusion, dummy_zeroes~.+0)
  call_taxis = lm( formula_taxis_full, data=Data_zp)$call
  call_diffusion = lm( formula_diffusion_full, data=Data_zp)$call
  Return$effects = list( Data_zp=Data_zp, call_taxis=call_taxis, call_diffusion=call_diffusion )
  Return$formula_taxis = formula_taxis
  Return$formula_diffusion = formula_diffusion

  # return
  class(Return) = "fitTMB"
  return(Return)
}

#' Print parameter estimates and standard errors.
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{fitTMB}}
#' @param ... Not used
#' @return NULL
#' @method print fitTMB
#' @export
print.fitTMB <- function(x, ...)
{
  cat("fitTMB(.) result\n")
  if( "parameter_estimates" %in% names(x) ){
    print( x$parameter_estimates )
  }else{
    cat("`parameter_estimates` not available in `fitTMB`\n")
  }
  invisible(x$parameter_estimates)
}

#' Predict habitat preference
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{fitTMB}}
#' @param ... Not used
#' @return NULL
#' @method predict fitTMB
#' @export
predict.fitTMB <- function(x,
               newdata,
               formula,
               varname = "alpha_k",
               origdata = NULL,
               prediction_type = 1,
               seed = NULL )
{
  #message("Running `predict.fitTMB`")
  if( !is.null(origdata) ){
    for(cI in 1:ncol(newdata)){
      if(is.factor(origdata[,cI])){
        newdata[,cI] = factor(newdata[,cI], levels=levels(origdata[,cI]))
      }
    }
  }
  fulldata = rbind( newdata, origdata )

  # MLE for covariance predictions
  if( prediction_type==1 ){
    if( "parameter_estimates" %in% names(x) ){
      if(varname=="alpha_k") var_vec = c( 0, x$Report$alpha_k )
      if(varname=="ln_sigma_l") var_vec = x$Report$ln_sigma_l
    }else{
      stop("`parameter_estimates` not available in `fitTMB`\n")
    }
  }
  if( prediction_type==2 ){
    if( "SD" %in% names(x$parameter_estimates) ){
      # Simulate new fixed-effect values
      set.seed(seed)
      fixed_sim = mvtnorm::rmvnorm( n=1, mean=x$parameter_estimates$par, sigma=x$parameter_estimates$SD$cov.fixed )
      # Add simulated fixed to joint parameter vector
      joint_par = x$Obj$env$last.par
      if(length(x$Obj$env$random)>0){
        joint_par[-x$Obj$env$random] = fixed_sim
      }else{
        joint_par = fixed_sim
      }
      # Use cheap report
      x$Obj$env$data$report_early = 1
      Report_sim = x$Obj$report(joint_par)
      if(varname=="alpha_k") var_vec = c( 0, Report_sim$alpha_k )  # Add zero because formula used here has an intercept, but fitted one doesn't
      if(varname=="ln_sigma_l") var_vec = Report_sim$ln_sigma_l
      x$Obj$env$data$report_early = 0
    }else{
      stop("`parameter_estimates` doesn't include covariance in `fit_TMB`")
    }
  }

  #print( x$parameter_estimates )
  new_matrix = model.matrix( update.formula(formula, ~.+1), data=fulldata )
  pred = ( new_matrix %*% var_vec )[,1]
  newpred = pred[1:nrow(newdata)]
  return(newpred)
}

#' Predict habitat preference
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{fitTMB}}
#' @param ... Not used
#' @return NULL
#' @method simulate fitTMB
#' @export
simulate.fitTMB <- function(x, random_seed=NULL, ...)
{

  # Check for loaded VAST
  # Modified from TMB:::getUserDLL
  dlls <- getLoadedDLLs()
  isTMBdll <- function(dll) !is(try(getNativeSymbolInfo("MakeADFunObject",dll), TRUE), "try-error")
  TMBdll <- sapply(dlls, isTMBdll)
  if( sum(TMBdll)==0 ){
    stop("ADM is not linked as a DLL, so `simulate_data` will not work.
    Please re-run model (potentially from informative starting values to save time) to use `simulate_data`")
  }else if(sum(TMBdll)>=2){
    warning("ADM is linked to multiple DLLs. Please consider using dyn.unload() to unload
    earlier ADM runs to avoid potentially ambiguous behavior when running `simulate_data`")
  }

  # Extract stuff
  Obj = x$Obj

  # Simulate
  set.seed(random_seed)
  Return = Obj$simulate( complete=TRUE )

  # return stuff
  return( Return )
}

#' Check residuals etc.
#'
#' \code{summary.fitTMB} extracts commonly used quantities derived from a fitted VAST model
#'
#' \code{summary.fitTMB} faciliates common queries for model output including:
#'
#' @return NULL
#' @method summary fitTMB
#' @export
summary.fitTMB <-
function( x,
          what="survey_residuals",
          n_samples=250,
          working_dir=NULL, ...)
{
  ans = NULL

  # Residuals
  if( tolower(what) %in% c("survey_residuals","fishery_residuals") ){
    # extract objects
    Obj = x$Obj

    message( "Sampling from the distribution of data conditional on estimated fixed and random effects" )
    if( what == "survey_residuals" ){
      x_z = x$data_list$b_j
      xhat_z = fit$Report$bhat_j
      x_zr = matrix(NA, nrow=length(x$data_list$b_j), ncol=n_samples)
      varname = "b_j"
    }
    if( what == "fishery_residuals" ){
      x_z = x$data_list$b_f
      xhat_z = fit$Report$bhat_f
      x_zr = matrix(NA, nrow=length(x$data_list$b_f), ncol=n_samples)
      varname = "b_f"
    }
    for( rI in 1:n_samples ){
      if( rI%%max(1,floor(n_samples/10)) == 0 ){
        message( "  Finished sample ", rI, " of ",n_samples )
      }
      x_zr[,rI] = Obj$simulate()[[varname]]
    }
    if( any(is.na(x_zr)) ){
      stop("Check simulated residuals for NA values")
    }

    # Run DHARMa
    dharmaRes = DHARMa::createDHARMa(simulatedResponse = x_zr, # + 1e-10*array(rnorm(prod(dim(b_iz))),dim=dim(b_iz)),
      observedResponse = x_z,
      fittedPredictedResponse = xhat_z,
      integer = FALSE)

    # Calculate probability-integral-transform (PIT) residuals
    message( "Substituting probability-integral-transform (PIT) residuals for DHARMa-calculated residuals" )
    prop_lessthan_j = apply( x_zr<outer(x_z,rep(1,n_samples)),
      MARGIN=1,
      FUN=mean )
    prop_lessthanorequalto_j = apply( x_zr<=outer(x_z,rep(1,n_samples)),
      MARGIN=1,
      FUN=mean )
    PIT_j = runif(min=prop_lessthan_j, max=prop_lessthanorequalto_j, n=length(prop_lessthan_j) )
    dharmaRes$scaledResiduals = PIT_j

    # do plot
    if( is.null(working_dir) ){
      plot(dharmaRes, ...)
    }else if(!is.na(working_dir) ){
      png(file=paste0(working_dir,"quantile_",what,".png"), width=8, height=4, res=200, units='in')
        plot(dharmaRes, ...)
      dev.off()
    }

    # Return stuff
    ans = dharmaRes
    message( "Invisibly returning output from `DHARMa::createDHARMa`, e.g., to apply `plot.DHARMa` to this output")
  }

  # diagnostic plots
  return(invisible(ans))
}

#' Get formula
#'
#' \code{formula.fitTMB} is necessary for \code{effects::predictorEffects}
#'
#' @method formula fitTMB
#' @export
formula.fitTMB = function(x){
  #if(which_formula=="taxis" | !exists(which_formula)){
    return(x$formula_taxis)
  #}else{
  #  stop("No `which_formula` in global environment")
  #}
}

