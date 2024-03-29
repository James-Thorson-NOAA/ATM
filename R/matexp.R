
#' Matrix exponential and Euler approximation
#'
#' \code{matexp} implements the matrix exponential and/or sparse Euler approximation
#'
#' @author James Thorson
#' @export
matexp <-
function( mat_gg, log2steps=0, zap_small=FALSE ){

  # Make sparse
  if( (log2steps <=0 ) || (log2steps > 100) ){
    return( Matrix::expm(mat_gg) )
  }else{
    mat_gg = Matrix::Diagonal(nrow(mat_gg)) + mat_gg / (2^log2steps)
    for(stepI in seq(1,log2steps,length=log2steps)){
      mat_gg = mat_gg %*% mat_gg
    }
    if(zap_small) mat_gg = zapsmall(mat_gg)
    return( mat_gg )
  }
}
