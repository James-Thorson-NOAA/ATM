
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
  ... ){

  # Build data
  data_list = make_data( X_guyk = X_guyk,
    loc_gz = loc_gz,
    uy_tz = uy_tz,
    satellite_iz = satellite_iz,
    survey_jz = survey_jz,
    duration_u = duration_u,
    cpp_version = cpp_version
  )
  #data_list$A_gg = as(data_list$A_gg, "dgTMatrix")

  # Convert indexing from R to CPP convention
  data_list$satellite_iz = data_list$satellite_iz - 1

  # Make parameters
  param_list = list(
    "ln_sigma" = log(1),
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

  if( FALSE ){
    Report = Obj$report()
    image(Report$Taxis_gg)
    image(Report$Diffusion_gg)
    summary(as.vector(Report$Diffusion_gg))
  }

  # Optimize
  parameter_estimates = TMBhelper::fit_tmb( Obj, ..., control=list(trace=1) )

  # Extract stuff
  Report = Obj$report()
  parhat = Obj$env$parList()

  # return
  Return = list("parameter_estimates"=parameter_estimates, "parhat"=parhat, "data_list"=data_list)
  class(Return) = "fitTMB"
  return(Return)
}
