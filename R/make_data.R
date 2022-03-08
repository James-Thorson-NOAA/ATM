
#' Build data input for ATM model
#'
#' \code{make_data} builds a tagged list of data inputs used by TMB for running the model
#'
#' @return Object of class \code{make_data}, containing inputs to function \code{\link{make_model}}
#'
#' @author James Thorson
#' @export
make_data <-
function( Cov_stars,
      formula_taxis,
      formula_diffusion,
      coords_gz,
      #t_uy,
      uy_tz = NULL,
      conventional_hz = NULL,
      satellite_iz = NULL,
      survey_jz = NULL,
      fishery_fz = NULL,
      E_guy = NULL,
      duration_u = NULL,
      cpp_version = "R",
      log2steps = 20,
      spde_aniso = NULL,
      alpha_ratio_bounds = 1,
      diffusion_bounds = 0,
      constant_tail_probability = 1e-8,
      movement_penalty = 0 ){

  #
  X_guyk = format_covariates( Cov_stars, formula=formula_taxis, remove_intercept=TRUE )$Data_guyk
  Z_guyl = format_covariates( Cov_stars, formula=formula_diffusion, remove_intercept=FALSE )$Data_guyk

  # Check for issues
  if( !is.na(log2steps) && abs(log2steps)==Inf ) stop("`log2steps` cannot be Inf")
  if( any(is.na(X_guyk)) ) stop("`X_guyk` includes NA values; please fix")
  if( any(is.na(Z_guyl)) ) stop("`Z_guyl` includes NA values; please fix")

  # Get dimensions
  n_g = dim(X_guyk)[1]
  n_u = dim(X_guyk)[2]
  n_y = dim(X_guyk)[3]

  # Fill in missing
  #if(missing(t_uy)){
  #  t_uy = matrix( 1:(n_u*n_t), nrow=n_u, ncol=n_y) # Keep with R indexing until bundling for TMB list
  #}
  if(missing(uy_tz) | is.null(uy_tz)){
    uy_tz = cbind( "u"=rep(1:n_u,n_y), "y"=rep(1:n_y,each=n_u) )
  }
  if(missing(satellite_iz) | is.null(satellite_iz)){
    satellite_iz = matrix(NA, nrow=0, ncol=4, dimnames=list(NULL,c("g_release","g_recovery","t_release","t_recovery")) )
  }else{
    if( any(satellite_iz[,'t_recovery'] < satellite_iz[,'t_release']) ){
      stop("Check for problem in satellite tags")
    }
  }
  if(missing(conventional_hz) | is.null(conventional_hz)){
    conventional_hz = matrix(NA, nrow=0, ncol=4, dimnames=list(NULL,c("g_release","g_recovery","t_release","t_recovery")) )
  }else{
    if( any(conventional_hz[,'t_recovery'] < conventional_hz[,'t_release']) ){
      stop("Check for problem in conventional tags")
    }
  }
  if(missing(duration_u) | is.null(duration_u)){
    duration_u = 1 / n_u
  }
  if(missing(survey_jz) | is.null(survey_jz)){
    survey_jz = matrix(NA, nrow=0, ncol=3, dimnames=list(NULL,c("b_j","t_j","g_j")) )
  }
  if(missing(fishery_fz) | is.null(fishery_fz)){
    fishery_fz = matrix(NA, nrow=0, ncol=3, dimnames=list(NULL,c("b_f","t_f","g_f")) )
  }
  if(missing(E_guy) | is.null(E_guy)){
    E_guy = array(1, dim=c(n_g,n_u,n_y))
  }

  # Number of time intervals
  n_t = nrow(uy_tz)

  # Rescale E_guy to be a proportion by location in each season-year
  Eprime_guy = E_guy / outer( rep(1,n_g), apply(E_guy,MARGIN=2:3,FUN=sum) )

  # Calculate adjacency matrix
  distance_gg = as.matrix(dist(coords_gz[,c('x','y')]))
  min_distance = min( ifelse(distance_gg==0,NA,distance_gg), na.rm=TRUE )
  Adense_gg = ifelse( abs(as.matrix(distance_gg)-min_distance)/min_distance < 1e-5, 1, 0 )
  Asparse_gg = as(Adense_gg, "dgTMatrix")

  # Bundle while switching indexing from R to CPP
  if( cpp_version %in% c("R") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz, "satellite_iz"=satellite_iz,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Adense_gg, "log2steps"=log2steps )
  }
  if( cpp_version %in% c("ATM_v0_9_0") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz-1, "satellite_iz"=satellite_iz-1,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Asparse_gg, "log2steps"=log2steps,
      "A_ij"=cbind(Asparse_gg@i,Asparse_gg@j) )
  }
  if( cpp_version %in% c("ATM_v1_0_0") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz-1, "satellite_iz"=satellite_iz-1,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Adense_gg, "log2steps"=log2steps )
  }
  if( cpp_version %in% c("ATM_v3_0_0","ATM_v2_0_0") ){
    data_list = list( "log2steps"=log2steps,"alpha_ratio_bounds"=alpha_ratio_bounds,
      "constant_tail_probability"=constant_tail_probability, "report_early"=FALSE,
      "X_guyk"=X_guyk, "uy_tz"=uy_tz-1, "satellite_iz"=satellite_iz-1,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Adense_gg,
      "spde_aniso"=spde_aniso, "b_j"=survey_jz[,'b_j'], "t_j"=survey_jz[,'t_j']-1, "g_j"=survey_jz[,'g_j']-1 )
  }
  if( cpp_version %in% c("ATM_v4_0_0") ){
    data_list = list( "log2steps"=log2steps,"alpha_ratio_bounds"=alpha_ratio_bounds, "diffusion_bounds"=diffusion_bounds,
      "constant_tail_probability"=constant_tail_probability, "report_early"=FALSE,
      "X_guyk"=X_guyk, "Z_guyl"=Z_guyl, "uy_tz"=uy_tz-1,
      "satellite_iz"=satellite_iz-1, "conventional_hz"=conventional_hz-1,
      "survey_jz"=survey_jz, "E_guy"=Eprime_guy, "duration_u"=duration_u, "A_gg"=Adense_gg,
      "spde_aniso"=spde_aniso, "b_j"=survey_jz[,'b_j'], "t_j"=survey_jz[,'t_j']-1, "g_j"=survey_jz[,'g_j']-1 )
  }
  if( cpp_version %in% c("ATM_v5_0_0") ){
    data_list = list( "log2steps"=log2steps,"alpha_ratio_bounds"=alpha_ratio_bounds,
      "diffusion_bounds"=diffusion_bounds, "movement_penalty"=movement_penalty,
      "constant_tail_probability"=constant_tail_probability, "report_early"=FALSE,
      "X_guyk"=X_guyk, "Z_guyl"=Z_guyl, "uy_tz"=uy_tz-1,
      "satellite_iz"=satellite_iz-1, "conventional_hz"=conventional_hz-1,
      "survey_jz"=survey_jz, "fishery_fz"=fishery_fz,
      "E_guy"=Eprime_guy, "duration_u"=duration_u, "A_gg"=Adense_gg, "spde_aniso"=spde_aniso,
      "b_j"=survey_jz[,'b_j'], "t_j"=survey_jz[,'t_j']-1, "g_j"=survey_jz[,'g_j']-1,
      "b_f"=fishery_fz[,'b_f'], "t_f"=fishery_fz[,'t_f']-1, "g_f"=fishery_fz[,'g_f']-1 )
  }
  if( cpp_version %in% c("ATM_v6_0_0") ){
    data_list = list( "log2steps"=log2steps,"alpha_ratio_bounds"=alpha_ratio_bounds,
      "diffusion_bounds"=diffusion_bounds, "movement_penalty"=movement_penalty,
      "constant_tail_probability"=constant_tail_probability, "report_early"=FALSE,
      "X_guyk"=X_guyk, "Z_guyl"=Z_guyl, "uy_tz"=uy_tz-1,
      "satellite_iz"=satellite_iz-1, "conventional_hz"=conventional_hz-1,
      "survey_jz"=survey_jz, "fishery_fz"=fishery_fz,
      "Eprime_guy"=Eprime_guy, "duration_u"=duration_u, "A_gg"=Adense_gg, "spde_aniso"=spde_aniso,
      "b_j"=survey_jz[,'b_j'], "t_j"=survey_jz[,'t_j']-1, "g_j"=survey_jz[,'g_j']-1,
      "b_f"=fishery_fz[,'b_f'], "t_f"=fishery_fz[,'t_f']-1, "g_f"=fishery_fz[,'g_f']-1 )
  }
  if( cpp_version %in% c("ATM_v7_0_0") ){
    data_list = list( "log2steps"=log2steps,"alpha_ratio_bounds"=alpha_ratio_bounds,
      "diffusion_bounds"=diffusion_bounds, "movement_penalty"=movement_penalty,
      "constant_tail_probability"=constant_tail_probability, "report_early"=FALSE,
      "simulate_random" = 0,
      "X_guyk"=X_guyk, "Z_guyl"=Z_guyl, "uy_tz"=uy_tz-1,
      "satellite_iz"=satellite_iz-1, "conventional_hz"=conventional_hz-1,
      "survey_jz"=survey_jz, "fishery_fz"=fishery_fz,
      "Eprime_guy"=Eprime_guy, "duration_u"=duration_u, "A_gg"=Adense_gg, "spde_aniso"=spde_aniso,
      "b_j"=survey_jz[,'b_j'], "t_j"=survey_jz[,'t_j']-1, "g_j"=survey_jz[,'g_j']-1,
      "b_f"=fishery_fz[,'b_f'], "t_f"=fishery_fz[,'t_f']-1, "g_f"=fishery_fz[,'g_f']-1 )
  }

  # return
  class(data_list) = "make_data"
  return(data_list)
}
