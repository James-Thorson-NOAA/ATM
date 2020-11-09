
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
  uy_tz,
  satellite_iz,
  survey_jz,
  cpp_version = "R",
  duration_u ){

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
  A_gg = ifelse( abs(as.matrix(distance_gg)-min_distance)/min_distance < 1e-5, 1, 0 )
  A_gg = as(A_gg, "dgTMatrix")

  # bundle
  if( cpp_version %in% c("ATM_v1_0_0.cpp", "R") ){
    data_list = list( "X_guyk"=X_guyk, "uy_tz"=uy_tz, "satellite_iz"=satellite_iz, "survey_jz"=survey_jz, "duration_u"=duration_u, "A_gg"=A_gg )
  }

  # return
  class(data_list) = "make_data"
  return(data_list)
}
