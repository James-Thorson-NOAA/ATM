
#' Data to demonstrate advection-diffusion-taxis model
#'
#' Data sufficient to replicate some features of the application of a
#'      Advection-Diffusion-Taxis model to Bering Sea Pacific cod
#'
#' \itemize{
#'   \item \code{Cov_stars} Stars object containing bathymetry and bottom temperature
#'   \item \code{loc_gz} SpatialPoints object containing projected location of grid cells
#'   \item \code{conventional_hz} Release and recapture location (as cell numbers) and times (as year-season index)
#'   \item \code{survey_jz} Survey samples containing cod numbers per hectare, grid cell, and year-season index
#'   \item \code{fishery_fz} Fishery samples containing cod numbers-per-hook, grid cell, and year-season index
#'   \item \code{E_guy} Fishing effort for every grid-cell, season, and year (included as dummy array containing 1's to avoid disclosing business information)
#' }
#'
#' @name Cod_data
#' @docType data
#' @usage data(Cod_data)
#' @keywords data
NULL

