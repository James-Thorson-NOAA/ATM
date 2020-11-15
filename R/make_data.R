
#' Build data input for ATM model
#'
#' \code{make_data} builds a tagged list of data inputs used by TMB for running the model
#'
#' @return Object of class \code{make_data}, containing inputs to function \code{\link{make_model}}
#'
#' @author James Thorson
#' @export
make_data <-
function( X_guyk,
  loc_gz,
  #t_uy,
  uy_tz = NULL,
  satellite_iz = NULL,
  survey_jz = NULL,
  duration_u = NULL,
  cpp_version = "R",
  log2steps = 20 ){

  # Check for issues
  if( !is.na(log2steps) && abs(log2steps)==Inf ) stop("`log2steps` cannot be Inf")
  if( any(is.na(X_guyk)) ) stop("`X_guyk` includes NA values; please fix")

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
  }
  if(missing(duration_u) | is.null(duration_u)){
    duration_u = 1 / n_u
  }
  if(missing(survey_jz) | is.null(survey_jz)){
    survey_jz = matrix(NA, nrow=0, ncol=3, dimnames=list(NULL,c("t_i","g_i","b_i")) )
  }

  # Number of time intervals
  n_t = nrow(uy_tz)

  # Calculate adjacency matrix
  distance_gg = as.matrix(dist(loc_gz[,c('x','y')]))
  min_distance = min( ifelse(distance_gg==0,NA,distance_gg), na.rm=TRUE )
  Adense_gg = ifelse( abs(as.matrix(distance_gg)-min_distance)/min_distance < 1e-5, 1, 0 )
  Asparse_gg = as(Adense_gg, "dgTMatrix")

  # bundle
  if( cpp_version %in% c("R") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz, "satellite_iz"=satellite_iz,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Adense_gg, "log2steps"=log2steps )
  }
  if( cpp_version %in% c("ATM_v1_0_0") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz-1, "satellite_iz"=satellite_iz-1,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Adense_gg, "log2steps"=log2steps )
  }
  if( cpp_version %in% c("ATM_v2_0_0") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz-1, "satellite_iz"=satellite_iz-1,
      "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=Asparse_gg, "log2steps"=log2steps,
      "A_ij"=cbind(Asparse_gg@i,Asparse_gg@j) )
  }

  # return
  class(data_list) = "make_data"
  return(data_list)
}
