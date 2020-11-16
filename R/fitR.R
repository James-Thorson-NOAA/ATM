
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
  uy_tz = NULL,
  satellite_iz = NULL,
  survey_jz = NULL,
  #cpp_version = FishStatsUtils::get_latest_version(package="ATM"),
  duration_u = NULL,
  log2steps = 20,
  sigma2 = 4^2 ){

  #
  start_time = Sys.time()

  # Build data
  data_list = make_data( X_guyk = X_guyk,
    loc_gz = loc_gz,
    uy_tz = uy_tz,
    satellite_iz = satellite_iz,
    survey_jz = survey_jz,
    duration_u = duration_u,
    #cpp_version = cpp_version,
    log2steps = log2steps
  )

  # Objective function
  Obj = function( param_vec, skeleton, what="NLL", data_list ){
    # Extract parameters
    param_list = relist(param_vec,skeleton)
    sigma2 = exp(2 * param_list$ln_sigma)
    alpha = sigma2 * (2 * plogis(param_list$alpha_logit_ratio) - 1)
    if(!is.null(dimnames(data_list$X_guyk)[[4]])){
      names(alpha) = dimnames(data_list$X_guyk)[[4]]
    }
    if(what=="params") return(list("alpha"=alpha, "sigma2"=sigma2))

    # Initiatize cumulator
    Mprime_ggt = array( NA, dim=c(dim(data_list$A_gg),nrow(data_list$uy_tz)) )
    Preference_gt = array( NA, dim=c(dim(data_list$A_gg)[1],nrow(data_list$uy_tz)) )

    #
    prob_satellite_igt = array(0, dim=c(nrow(data_list$satellite_iz),dim(data_list$X_guyk)[1],nrow(data_list$uy_tz)))

    # Check for issues
    Problem_tz = matrix( FALSE, nrow=nrow(data_list$uy_tz), ncol=3)

    # Loop through times
    for( tI in 1:nrow(data_list$uy_tz) ){
      # Diffusion-rate matrix
      Diffusion_gg = sigma2 * data_list$A_gg
      diag(Diffusion_gg) = -colSums(Diffusion_gg) ## Mass conservation

      # Advection-rate matrix
      Preference_gt[,tI] = as.vector( X_guyk[,data_list$uy_tz[tI,'u'],data_list$uy_tz[tI,'y'],] %*% alpha )
      Taxis_gg = data_list$A_gg * outer( Preference_gt[,tI], Preference_gt[,tI], FUN=function(a,b)a-b )
      diag(Taxis_gg) = -colSums(Taxis_gg) ## Mass conservation

      # Movement-rate matrix
      Mprime_ggt[,,tI] = Diffusion_gg + Taxis_gg

      # Movement probability matrix
      Movement_gg = matexp(Mprime_ggt[,,tI], log2steps=data_list$log2steps)

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
          prob_satellite_igt[iI,,tI] = as.vector(Movement_gg %*% init)
        }else if(satellite_iz[iI,'t_release']<tI & satellite_iz[iI,'t_recovery']>=tI){
          prob_satellite_igt[iI,,tI] = as.vector(Movement_gg %*% prob_satellite_igt[iI,,tI-1])
        }
      }
    }

    #
    NLL_i = rep( NA, nrow(satellite_iz) )
    for( iI in 1:nrow(satellite_iz) ){
      NLL_i[iI] = -1 * log(prob_satellite_igt[iI,satellite_iz[iI,'g_recovery'],satellite_iz[iI,'t_recovery']])
    }

    #
    if(what=="Msum_gg") return( expm(apply(Mprime_ggt,MARGIN=1:2,FUN=sum)) )
    if(what=="Mprimesum_gg") return( apply(Mprime_ggt,MARGIN=1:2,FUN=sum) )
    if(what=="Mprime_ggt") return( Mprime_ggt )
    if(what=="final_D") return( Diffusion_gg )
    if(what=="final_T") return( Taxis_gg )
    if(what=="Preference_gt") return( Preference_gt )
    if(what=="NLL_i") return( NLL_i )
    if(what=="Mannual_ggt"){
      Mannual_ggt = array(NA, dim=dim(Mprime_ggt) )
      for( tI in 1:nrow(data_list$uy_tz) ){
        t_set = tI + 1:dim(X_guyk)[2] - 1
        if( all(t_set <= nrow(data_list$uy_tz)) ){
          Mannual_ggt[,,tI] = matrix(expm(apply(Mprime_ggt[,,t_set,drop=FALSE],MARGIN=1:2,FUN=sum)))
        }
      }
      return( Mannual_ggt )
    }

    NLL = sum(NLL_i)
    if(any(Problem_tz)) NLL = 1000000
    return(NLL)
  }

  # Make parameter list
  if(missing(sigma2) || is.na(sigma2)) sigma2 = 4^2
  param_list = list( "ln_sigma"=log(sqrt(sigma2)),
    "alpha_logit_ratio"=0.1*rnorm(dim(data_list$X_guyk)[4]) )

  # Compress to vector
  param_vec = unlist(param_list)
  #return(Obj(param_vec, skeleton=param_list, data_list=data_list))

  # Optimize
  parameter_estimates = nlminb( objective=Obj, start=param_vec, skeleton=param_list, data_list=data_list, control=list(trace=1) )
  parameter_estimates[["time_for_MLE"]] = Sys.time() - start_time

  # Add slots
  Return = list("parameter_estimates"=parameter_estimates, "data_list"=data_list)
  class(Return) = "fitR"
  Return$parhat = Obj(parameter_estimates$par, what="params", skeleton=param_list, data_list=data_list)
  Return$Msum_gg = Obj(parameter_estimates$par, what="Msum_gg", skeleton=param_list, data_list=data_list)
  Return$Mannual_ggt = Obj(parameter_estimates$par, what="Mannual_ggt", skeleton=param_list, data_list=data_list)
  Return$Preference_gt = Obj(parameter_estimates$par, what="Preference_gt", skeleton=param_list, data_list=data_list)
  #Return$Mprime_ggt = Obj(parameter_estimates$par, what="Mprime_ggt", skeleton=param_list, data_list=data_list)
  #Return$Mprimesum_gg = Obj(parameter_estimates$par, what="Mprimesum_gg", skeleton=param_list, data_list=data_list)
  #Return$Diffusion_gg = Obj(parameter_estimates$par, what="final_D", skeleton=param_list, data_list=data_list)
  #Return$Taxis_gg = Obj(parameter_estimates$par, what="final_T", skeleton=param_list, data_list=data_list)
  #Return$NLL_i = Obj(parameter_estimates$par, what="NLL_i", skeleton=param_list, data_list=data_list)

  # return
  return(Return)
}

#' Print parameter estimates and standard errors.
#'
#' @title Print parameter estimates
#' @param x Output from \code{\link{fitR}}
#' @param ... Not used
#' @return NULL
#' @method print fitR
#' @export
print.fitR <- function(x, ...)
{
  cat("fitR(.) result\n")
  if( "parameter_estimates" %in% names(x) ){
    print( x$parameter_estimates )
  }else{
    cat("`parameter_estimates` not available in `fitR`\n")
  }
  invisible(x$parameter_estimates)
}

