
#' Fit ATM model in R
#'
#' \code{fitR} fits the Advection-taxis model in base-R
#'
#' @return Object of class \code{fitR}, containing parameter estimates and predicted time-integrated movement
#'
#' @author James Thorson
#' @export
fitR <-
function( X_guyk,
  loc_gz,
  #t_uy,
  uy_tz,
  satellite_iz = NULL,
  survey_jz = NULL,
  #cpp_version = FishStatsUtils::get_latest_version(package="ATM"),
  duration_u = NULL ){

  # Build data
  data_list = make_data( X_guyk = X_guyk,
    loc_gz = loc_gz,
    uy_tz = uy_tz,
    satellite_iz = satellite_iz,
    survey_jz = survey_jz,
    duration_u = duration_u #,
    #cpp_version = cpp_version
  )

  # Objective function
  Obj = function( param_vec, skeleton, what="NLL", data_list ){
    # Extract parameters
    param_list = relist(param_vec,skeleton)
    sigma2 = exp(2 * param_list$ln_sigma)
    alpha = sigma2 * (-1+2*plogis(param_list$alpha_logit_ratio))
    if(!is.null(dimnames(data_list$X_guyk)[[4]])){
      names(alpha) = dimnames(data_list$X_guyk)[[4]]
    }
    if(what=="params") return(list("alpha"=alpha, "sigma2"=sigma2))

    # Initiatize cumulator
    Mprimesum_gg = data_list$A_gg
    Mprimesum_gg[]= 0

    #
    satellite_igt = array(0, dim=c(nrow(data_list$satellite_iz),dim(data_list$X_guyk)[1],nrow(data_list$uy_tz)))

    # Check for issues
    Problem_tz = matrix( FALSE, nrow=nrow(data_list$uy_tz), ncol=3)

    # Loop through times
    for( tI in 1:nrow(data_list$uy_tz) ){
      # Diffusion-rate matrix
      Diffusion_gg = sigma2 * data_list$A_gg
      diag(Diffusion_gg) = -colSums(Diffusion_gg) ## Mass conservation

      # Advection-rate matrix
      Preference_g = as.vector( X_guyk[,data_list$uy_tz[tI,'u'],data_list$uy_tz[tI,'y'],] %*% alpha )
      Taxis_gg = data_list$A_gg * outer( Preference_g, Preference_g, FUN=function(a,b)a-b )
      diag(Taxis_gg) = -colSums(Taxis_gg) ## Mass conservation

      # Movement-rate matrix
      Mprime_gg = Diffusion_gg + Taxis_gg
      Mprimesum_gg = Mprimesum_gg + Mprime_gg

      # Movement probability matrix
      Movement_gg = expm(Mprime_gg)

      # Check
      #assign(x="par", value=par, envir=.GlobalEnv)
      if( any(is.na(Movement_gg)) ){
        Problem_tz[tI,1] = TRUE
      }else if( any(Movement_gg<0) ){
        Problem_tz[tI,2] = TRUE
      }#else if( abs(matlib::powerMethod(as.matrix(M))$value-1)>0.01 ){
      #  Problem[3] = TRUE
      #}

      # Apply to satellite tags
      for( iI in 1:nrow(satellite_iz) ){
        if( data_list$satellite_iz[iI,'t_release']==tI ){
          init = rep(0,dim(data_list$X_guyk)[1])
          init[satellite_iz[iI,'g_release']] = 1
          satellite_igt[iI,,tI] = as.vector(Movement_gg %*% init)
        }else if(satellite_iz[iI,'t_release']<tI & satellite_iz[iI,'t_recovery']>=tI){
          satellite_igt[iI,,tI] = as.vector(Movement_gg %*% satellite_igt[iI,,tI-1])
        }
      }
    }

    #
    NLL_satellite = 0
    for( iI in 1:nrow(satellite_iz) ){
      NLL_satellite = NLL_satellite - log(satellite_igt[iI,satellite_iz[iI,'g_recovery'],satellite_iz[iI,'t_recovery']])
    }

    #
    NLL = NLL_satellite
    if(any(Problem_tz)) NLL = 1000000
    if(what=="Msum_gg") return( expm(Mprimesum_gg) )

    return(NLL)
  }

  # Make parameter list
  sigma2 = 4^2
  param_list = list( "ln_sigma"=log(4),
    "alpha_logit_ratio"=rep(0,dim(data_list$X_guyk)[4]) )

  # Compress to vector
  param_vec = unlist(param_list)
  #return(Obj(param_vec, skeleton=param_list, data_list=data_list))

  # Optimize
  parameter_estimates = nlminb( objective=Obj, start=param_vec, skeleton=param_list, data_list=data_list, control=list(trace=1) )
  parhat = Obj(parameter_estimates$par, what="params", skeleton=param_list, data_list=data_list)
  Msum_gg = Obj(parameter_estimates$par, what="Msum_gg", skeleton=param_list, data_list=data_list)

  # return
  Return = list("parameter_estimates"=parameter_estimates, "parhat"=parhat, "Msum_gg"=Msum_gg, "data_list"=data_list)
  class(Return) = "fitR"
  return(Return)
}
