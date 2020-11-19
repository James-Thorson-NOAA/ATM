
#' Fit ATM model in TMB
#'
#' \code{fitTMB} fits the Advection-taxis model in TMB
#'
#' @return Object of class \code{fitTMB}, containing parameter estimates and predicted time-integrated movement
#'
#' @author James Thorson
#' @export
fitTMB <-
function( X_guyk,
  coords_gz,
  #t_uy,
  uy_tz = NULL,
  satellite_iz = NULL,
  survey_jz = NULL,
  duration_u = NULL,
  cpp_version = FishStatsUtils::get_latest_version(package="ATM"),
  tmb_dir = system.file("executables",package="ATM"),
  run_dir = getwd(),
  compile_dir = run_dir,
  log2steps = 20,
  sigma2 = 4^2,
  spde_aniso = NULL,
  use_REML = TRUE,
  ... ){

  # Build data
  if( FALSE ){
    coords_gz = loc_gz@coords
    duration_u = NULL
    uy_tz = NULL
    log2steps = 0
    cpp_version = "ATM_v3_0_0"
    sigma2 = 0.1
    tmb_dir = "C:/Users/James.Thorson/Desktop/Git/ATM/inst/executables/"
    run_dir = "C:/Users/James.Thorson/Desktop/Work files/Collaborations/2020 -- Advection-taxis movement/"
    compile_dir = run_dir
    use_REML = TRUE
  }
  data_list = make_data( X_guyk = X_guyk,
    coords_gz = coords_gz,
    uy_tz = uy_tz,
    satellite_iz = satellite_iz,
    survey_jz = survey_jz,
    duration_u = duration_u,
    cpp_version = cpp_version,
    log2steps = log2steps,
    spde_aniso = spde_aniso
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

  # Which map
  map = NULL
  if( "Beta_t" %in% names(param_list) ){
    Beta_t = 1:length(param_list$Beta_t) - 1
    Beta_t = ifelse( Beta_t %in% data_list$t_j, Beta_t, NA )
    map$Beta_t = factor(Beta_t)
  }

  #
  Random = c("Omegainput_s", "Epsiloninput_st", "ln_d_st")
  if( use_REML==TRUE ){
    Random = union( Random, c("Beta_t","ln_phi","alpha_logit_ratio_k") )
  }
  Random = Random[which(Random %in% names(param_list))]
  if( length(Random)==0) Random = NULL
  #Random = NULL

  # Compile
  # dyn.unload( paste0(compile_dir,"/",TMB::dynlib(TMB:::getUserDLL())) )
  file.copy( from=paste0(tmb_dir,"/",cpp_version,".cpp"), to=paste0(compile_dir,"/",cpp_version,".cpp"), overwrite=FALSE)
  origwd = getwd()
  on.exit(setwd(origwd),add=TRUE)
  setwd( compile_dir )
  TMB::compile( paste0(cpp_version,".cpp"), CPPFLAGS="-Wno-ignored-attributes" )

  # Build object
  dyn.load( paste0(compile_dir,"/",TMB::dynlib(cpp_version)) ) # random=Random,
  Obj = TMB::MakeADFun( data=data_list, parameters=param_list, hessian=FALSE, map=map, random=Random, DLL=cpp_version )  #
  # Report = Obj$report()

  # Print number of parameters
  ThorsonUtilities::list_parameters( Obj )

  # Optimize
  Obj$env$beSilent()
  parameter_estimates = TMBhelper::fit_tmb( Obj, control=list(trace=1), ... )

  # Extract stuff
  Return = list("parameter_estimates"=parameter_estimates, "data_list"=data_list)
  class(Return) = "fitTMB"
  Return$Report = Obj$report( )
  Return$parhat = Obj$env$parList( parameter_estimates$par )

  # return
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

