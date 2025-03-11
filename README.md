# ATM
Advection-taxis movement model

`ATM` can be used to fit realistic movement models to tags, survey, and fishery data.

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
  Lon_i = lonlat_gz@coords[,1],
  Lat_i = lonlat_gz@coords[,2],
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

### Preference landscapes
Users can then visualize the landscape of habitat preference:

```R
X = fit$Report$Preference_gt - (rep(1,length(Cod_data$loc_gz)) %o% colMeans(fit$Report$Preference_gt))
DF = SpatialPointsDataFrame( coords=Cod_data$loc_gz@coords, data=data.frame(X[1:length(Cod_data$loc_gz),]) )
spplot(DF)
```

### Covariate responses
Similarly, users can visualize covariate-response curves using a partial-dependence-plot with or without uncertainty intervals:

```R
library(pdp)
library(ggplot2)
Data_zp = format_covariates( Cod_data$Cov_stars,
                             formula=fit$formula_taxis )$Data_zp
# MLE of effects, using predict.fit_model
PDP = partial( fit,
              train = Data_zp,
              origdata = Data_zp,
              pred.var = c("Bathy","season"),
              type = "regression",
              varname = "alpha_k",
              formula = fit$formula_taxis )
autoplot(PDP)

# Sample distribution for effects
Interval = plot_partial_dependence( x = fit,
                                    n_samples = 100,
                                    train = Data_zp,
                                    pred.var = "BT",
                                    type = "regression",
                                    filename = paste0(run_dir,"taxis_BT_CI"),
                                    varname = "alpha_k",
                                    formula = fit$formula_taxis,
                                    bounds_type = "whiskers" )
Interval = plot_partial_dependence( x = fit,
                                    n_samples = 100,
                                    train = Data_zp,
                                    pred.var = c("Bathy","season"),
                                    type = "regression",
                                    filename = paste0(run_dir,"taxis_Bathy_CI"),
                                    varname = "alpha_k",
                                    formula = fit$formula_taxis,
                                    bounds_type = "whiskers" )
```

### Simulating new data
Finally, `ATM` can be used as an operating model, either conditional upon estimated values (a "conditional simulator") or simulating new random effects (an "unconditional simulator").  The latter is helpful to generate a new history of the population, e.g., when manually changing movement or other parameters.

```R
# Modify newpar as needed
# e.g., newpar[1:3] = [CHANGES]
# Here I'm just randomizing to make a point
newpar = fit$Obj$env$last.par.best
newpar[grep("alpha_logit_ratio_k",names(newpar))] =
  newpar[grep("alpha_logit_ratio_k",names(newpar))] + rnorm(length(grep("alpha_logit_ratio_k",names(newpar))))

# Plot changes in preference function
PDP = partial( fit,
              train = Data_zp,
              origdata = Data_zp,
              pred.var = c("Bathy","season"),
              type = "regression",
              varname = "alpha_k",
              formula = fit$formula_taxis,
              parvec = newpar )
autoplot(PDP)

# "unconditional" simulation of random effects
# using the new preference function
out = simulate_data( fit,
                     type = 2,
                     parvec = newpar )
out$tmp_s
colMeans(exp(out$ln_d_st))
```

# NOAA Enterprise GitHub disclaimer
This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
