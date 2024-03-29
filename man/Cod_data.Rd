% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{Cod_data}
\alias{Cod_data}
\title{Data to demonstrate advection-diffusion-taxis model}
\usage{
data(Cod_data)
}
\description{
Data sufficient to replicate some features of the application of a
     Advection-Diffusion-Taxis model to Bering Sea Pacific cod
}
\details{
\itemize{
  \item \code{Cov_stars} Stars object containing bathymetry and bottom temperature
  \item \code{loc_gz} SpatialPoints object containing projected location of grid cells
  \item \code{conventional_hz} Release and recapture location (as cell numbers) and times (as year-season index)
  \item \code{survey_jz} Survey samples containing cod numbers per hectare, grid cell, and year-season index
  \item \code{fishery_fz} Fishery samples containing cod numbers-per-hook, grid cell, and year-season index
  \item \code{E_guy} Fishing effort for every grid-cell, season, and year (included as dummy array containing 1's to avoid disclosing business information)
}
}
\keyword{data}
