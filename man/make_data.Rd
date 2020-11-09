\name{make_data}
\alias{make_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_data(X_guyk, loc_gz, uy_tz, satellite_iz, survey_jz, duration_u, cpp_version = FishStatsUtils::get_latest_version(package = "ATM"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X_guyk}{
%%     ~~Describe \code{X_guyk} here~~
}
  \item{loc_gz}{
%%     ~~Describe \code{loc_gz} here~~
}
  \item{uy_tz}{
%%     ~~Describe \code{uy_tz} here~~
}
  \item{satellite_iz}{
%%     ~~Describe \code{satellite_iz} here~~
}
  \item{survey_jz}{
%%     ~~Describe \code{survey_jz} here~~
}
  \item{duration_u}{
%%     ~~Describe \code{duration_u} here~~
}
  \item{cpp_version}{
%%     ~~Describe \code{cpp_version} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X_guyk, loc_gz, uy_tz, satellite_iz, survey_jz, duration_u, 
    cpp_version = FishStatsUtils::get_latest_version(package = "ATM")) 
{
    n_g = dim(X_guyk)[1]
    n_u = dim(X_guyk)[2]
    n_y = dim(X_guyk)[3]
    if (missing(uy_tz)) {
        uy_tz = cbind(u = rep(1:n_u, n_y), y = rep(1:n_y, each = n_u))
    }
    if (missing(satellite_iz)) {
        satellite_iz = matrix(NA, nrow = 0, ncol = 4, dimnames = list(NULL, 
            c("g_release", "g_recovery", "t_release", "t_recovery")))
    }
    if (missing(duration_u)) {
        duration_u = 1/n_u
    }
    if (missing(survey_jz)) {
        survey_jz = matrix(NA, nrow = 0, ncol = 3, dimnames = list(NULL, 
            c("t_i", "g_i", "b_i")))
    }
    n_t = nrow(uy_tz)
    distance_gg = as.matrix(dist(loc_gz[, c("x", "y")]))
    min_distance = min(ifelse(distance_gg == 0, NA, distance_gg), 
        na.rm = TRUE)
    A_gg = ifelse(abs(as.matrix(distance_gg) - min_distance)/min_distance < 
        1e-05, 1, 0)
    A_gg = as(A_gg, "dgTMatrix")
    if (cpp_version \%in\% c("ATM_v1_0_0.cpp")) {
        data_list = list(X_guyk = X_guyk, uy_tz = uy_tz, sattelite_iz = sattelite_iz, 
            survey_jz = survey_iz, duration_u = duration_u, A_gg = A_gg)
    }
    return(data_list)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
