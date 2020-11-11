
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
  loc_gz,
  #t_uy,
  uy_tz,
  satellite_iz = NULL,
  survey_jz = NULL,
  duration_u = NULL,
  cpp_version = FishStatsUtils::get_latest_version(package="ATM"),
  tmb_dir = system.file("executables",package="ATM"),
  run_dir = getwd(),
  compile_dir = run_dir,
  log2steps = 20,
  ... ){

  # Build data
  data_list = make_data( X_guyk = X_guyk,
    loc_gz = loc_gz,
    uy_tz = uy_tz,
    satellite_iz = satellite_iz,
    survey_jz = survey_jz,
    duration_u = duration_u,
    cpp_version = cpp_version,
    log2steps = log2steps
  )
  #data_list$A_gg = as(data_list$A_gg, "dgTMatrix")

  # Convert indexing from R to CPP convention
  data_list$satellite_iz = data_list$satellite_iz - 1
  data_list$uy_tz = data_list$uy_tz - 1

  # Make parameters
  param_list = list(
    "ln_sigma" = log(100),
    "alpha_logit_ratio_k" = 0.1 * rnorm(dim(data_list$X_guyk)[4])
  )

  # Which map
  map = NULL

  # Compile
  # dyn.unload( paste0(compile_dir,"/",TMB::dynlib(TMB:::getUserDLL())) )
  file.copy( from=paste0(tmb_dir,"/",cpp_version,".cpp"), to=paste0(compile_dir,"/",cpp_version,".cpp"), overwrite=FALSE)
  origwd = getwd()
  on.exit(setwd(origwd),add=TRUE)
  setwd( compile_dir )
  TMB::compile( paste0(cpp_version,".cpp"), CPPFLAGS="-Wno-ignored-attributes" )

  # Build object
  dyn.load( paste0(compile_dir,"/",TMB::dynlib(cpp_version)) ) # random=Random,
  Obj = TMB::MakeADFun( data=data_list, parameters=param_list, hessian=FALSE, map=map, DLL=cpp_version )  #
  Report = Obj$report()

  # Optimize
  parameter_estimates = TMBhelper::fit_tmb( Obj, control=list(trace=1), ... )

  # Extract stuff
  Return = list("parameter_estimates"=parameter_estimates, "data_list"=data_list)
  class(Return) = "fitTMB"
  Return$Report = Obj$report( parameter_estimates$par )
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

