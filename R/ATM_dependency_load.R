################################################################################
##### Rtools4.2 install
# Install wizard
# https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html

################################################################################
##### TMB install
install.packages("TMB")

##### TMB example
library(TMB)
## Optionally:
## precompile()
runExample(all=TRUE)

################################################################################
##### VAST install
# Install and load devtools package
install.packages("devtools")
library("devtools")

# Install INLA using currently recommended method
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

devtools::install_local(path="C:/Users/Matt Siskey/desktop/TMB_contrib_R-master/TMBhelper", dependencies=FALSE)
devtools::install_local(path="C:/Users/Matt Siskey/desktop/utilities-master", dependencies=FALSE)
devtools::install_local(path="C:/Users/Matt Siskey/desktop/FishStatsUtils-main", dependencies=FALSE)
devtools::install_local(path="C:/Users/Matt Siskey/desktop/geostatistical_delta-GLMM-master", dependencies=FALSE)
devtools::install_local(path="C:/Users/Matt Siskey/desktop/VAST-main", dependencies=FALSE)

##### VAST example
# Load package
library(VAST)

# load data set
# see `?load_example` for list of stocks with example data 
# that are installed automatically with `FishStatsUtils`. 
example = load_example( data_set="EBS_pollock" )

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 100, 
                          Region = example$Region, 
                          purpose = "index2", 
                          bias.correct = FALSE )

# Run model
fit = fit_model( settings = settings, 
                 Lat_i = example$sampling_data[,'Lat'], 
                 Lon_i = example$sampling_data[,'Lon'], 
                 t_i = example$sampling_data[,'Year'], 
                 b_i = example$sampling_data[,'Catch_KG'], 
                 a_i = example$sampling_data[,'AreaSwept_km2'] )

# Plot results
plot( fit )

################################################################################
##### ATM install
devtools::load_all("C:/Users/Matt Siskey/Desktop/ATM-main")

# Next, it is currently required that the user go to ATM_v6_0_0.cpp
# and comment-out lines 51-59 for fitTMB to run in the demo below

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
fit = ATM::fitTMB(
  Cov_stars = Cod_data$Cov_stars,
  formula_diffusion = ~ 1,
  formula_taxis = ~ bs(sqrt(Bathy), knots=3, intercept=FALSE):season + bs(BT, knots=3, intercept=FALSE),
  coords_gz = Cod_data$loc_gz@coords,
  satellite_iz = NULL,
  survey_jz = Cod_data$survey_jz,
  fishery_fz = NULL,
  conventional_hz = Cod_data$conventional_hz,
  E_guy = Cod_data$E_guy,
  tmb_dir                   = "C:/Users/Matt Siskey/Desktop/ATM-main/inst/executables/",
  cpp_version               = "ATM_v6_0_0",
  log2steps = 0,
  sigma2 = 0.1,
  spde_aniso = spde_aniso,
  alpha_ratio_bounds = 0,
  movement_penalty = 10000,
  constant_tail_probability = 1e-7 )
