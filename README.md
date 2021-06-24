# ATM
Advection-taxis movement model

### Demo usage
```R
# Load package and case-study data
library(ATM)
library(splines)
data( Cod_data )

# Generate spatial domain
extrapolation_list = FishStatsUtils::make_extrapolation_info(
  Region = c("eastern_bering_sea","northern_bering_sea"),
  projargs = NA,
  #projargs = CRSobj@projargs,
  create_strata_per_region = TRUE,
  max_cells = Inf )

# Generate spatial objects
lonlat_gz = sp::spTransform(
  Cod_data$loc_gz,
  CRSobj = CRS("+proj=longlat") )
spatial_list = FishStatsUtils::make_spatial_info(
  fine_scale = FALSE,
  n_x = nrow(lonlat_gz@coords),
  Extrapolation_List = extrapolation_list,
  Lon_i = lonlat_gz@coords[,'x'],
  Lat_i = lonlat_gz@coords[,'y'],
  knot_method = "samples",
  grid_size_km = 5000 )

# Format SPDE object representing spatial smoother
spde_aniso = list(
  "n_s" = spatial_list$MeshList$anisotropic_spde$n.spde,
  "n_tri" = nrow(spatial_list$MeshList$anisotropic_mesh$graph$tv),
  "Tri_Area" = spatial_list$MeshList$Tri_Area,
  "E0" = spatial_list$MeshList$E0,
  "E1" = spatial_list$MeshList$E1,
  "E2" = spatial_list$MeshList$E2,
  "TV" = spatial_list$MeshList$TV-1,
  "G0" = spatial_list$MeshList$anisotropic_spde$param.inla$M0,
  "G0_inv" = INLA::inla.as.dgTMatrix(solve(spatial_list$MeshList$anisotropic_spde$param.inla$M0)) )

# Fit model
# NOTE:  turning off fishery data, which after truncating to avoid dislosing confidential data was resulting
# in poor fit (not PD-Hessian) with the covariates used in the paper.
fit = fitTMB(
  Cov_stars = Cod_data$Cov_stars,
  formula_diffusion = ~ 1,
  formula_taxis = ~ bs(sqrt(Bathy), knots=3, intercept=FALSE):season + bs(BT, knots=3, intercept=FALSE),
  coords_gz = Cod_data$loc_gz@coords,
  satellite_iz = NULL,
  survey_jz = Cod_data$survey_jz,
  fishery_fz = NULL,
  conventional_hz = Cod_data$conventional_hz,
  E_guy = Cod_data$E_guy,
  log2steps = 0,
  sigma2 = 0.1,
  spde_aniso = spde_aniso,
  alpha_ratio_bounds = 0,
  movement_penalty = 10000,
  constant_tail_probability = 1e-7 )
```
