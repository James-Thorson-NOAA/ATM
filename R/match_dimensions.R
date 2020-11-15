
#' @export
match_dimensions = function( X_guyk,
  loc_gz,
  lat_i,
  lon_i,
  season_i,
  year_i,
  CRSobj = sp::CRS('+proj=longlat') ){

  # Coordinates
  lonlat_i = SpatialPoints( coords=cbind(lon_i,lat_i), proj4string=sp::CRS('+proj=longlat') )
  locproj_i = sp::spTransform( x=lonlat_i, CRSobj=CRSobj )
  locproj_g = sp::spTransform( loc_gz, CRSobj=CRSobj )

  # Match
  g_i = RANN::nn2( data=locproj_g@coords, query=locproj_i@coords, k=1 )$nn.idx[,1]
  u_i = match( season_i, dimnames(X_guyk)[[2]] )
  y_i = match( year_i, dimnames(X_guyk)[[3]] )
  t_i = u_i + (y_i-1)*dim(X_guyk)[2]

  #
  indices_i = data.frame( "g_i"=g_i, "u_i"=u_i, "y_i"=y_i, "t_i"=t_i )
  return(indices_i)
}

