
#' @export
format_covariates = function( Cov_stars, formula=~1, remove_intercept=TRUE ){

  ## S3 method for class 'stars'
  #Data_zk = as.data.frame(Cov_stars)
  Data_zp = na.omit(as.data.frame(Cov_stars))
  if("Bathy" %in% colnames(Data_zp)) Data_zp$Bathy = Data_zp$Bathy / 1000
  if("x" %in% colnames(Data_zp)) Data_zp$x = Data_zp$x / 1000000
  if("y" %in% colnames(Data_zp)) Data_zp$y = Data_zp$y / 1000000

  if( FALSE ){
    # See number of rows
    nrow(unique(Data_zp[,c('x','y')]))

    # Check season averages
    tapply( Data_zp[,'BT'], INDEX=Data_zp[,'season'], FUN=mean )
    tapply( Data_zp[,'BT'], INDEX=list(Data_zp[,'year'],Data_zp[,'season']), FUN=mean )
  }

  # Expand with model.matrix(.)
  formula_to_use = update.formula(formula, ~.+1)
  Model_matrix = model.matrix( formula_to_use, data=Data_zp )
  if( remove_intercept==TRUE ){
    Columns_to_keep = which( attr(Model_matrix,"assign") != 0 )
  }else{
    Columns_to_keep = 1:ncol(Model_matrix)
  }
  coefficient_names = attr(Model_matrix,"dimnames")[[2]][Columns_to_keep]
  Data_zk = Model_matrix[,Columns_to_keep,drop=FALSE]
  dimnames(Data_zk) = list(NULL, coefficient_names)

  # Remove all-0 columns (which occur when using two basis-splines for reasons I don't understand)
  # Avoid doing this to simplify predict.fitTMB, where model.matrix(formula) is also used
  # Instead mapping off parameters for constant covariates
  #all_zero_k = apply(Data_zk, MARGIN=2, FUN=function(vec){all(vec==0)})
  #if(any(all_zero_k)){
  #  Data_zk = Data_zk[,setdiff(1:ncol(Data_zk),which(all_zero_k))]
  #}

  # Reformat:   Season-Year must be in same order as along in st_as_stars(.)
  Data_uygk = array(NA, dim=c(length(unique(Data_zp[,'season'])), length(unique(Data_zp[,'year'])), nrow(unique(Data_zp[,c('x','y')])), ncol(Data_zk)) )
  Xvec = as.vector(Data_zk)
  #Xvec = c( Xvec, rep(NA,prod(dim(X_yugk))-length(Xvec)) )
  # Xvec = 1:length(Xvec)
  Data_uygk[] = Xvec
  Data_guyk = aperm( Data_uygk, c(3,1,2,4) )
  dimnames(Data_guyk) = list( NULL,
    "season"=st_dimensions(Cov_stars)$season$values,
    "year"=st_dimensions(Cov_stars)$year$values,
    "covariates"=colnames(Data_zk) )

  # If none present, then add back an intercept
  if( dim(Data_guyk)[4] == 0 ){
    Data_guyk = array( 1, dim=c(dim(Data_guyk)[1:3],1), dimnames=c(dimnames(Data_guyk)[1:3],"(Intercept)") )
  }

  Return = list("Data_guyk"=Data_guyk, "Data_zp"=Data_zp, "Data_zk"=Data_zk)
  return(Return)
}
