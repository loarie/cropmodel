% Generated by roxygen2 (4.1.1.9000): do not edit by hand
% Please edit documentation in R/crop_function.R
\name{crop_function}
\alias{crop_function}
\title{Crop Function}
\usage{
crop_function(b0, Vbcoef, s1, s2, y, x1, x2)
}
\arguments{
\item{b0,}{Vbcoef, s1, s2, y, x1, x2}
}
\description{
This function takes responses (y), and 2 environmental variables (x1 & x2) along with prior parameters for beta (b0, Vbcoef) and sigma2 (s1,s2) and produces estimates for beta and sigma using a Gibb's sampler and quadradic Bayesian regression
}
\examples{
crop_function()
}
\keyword{Bayesian}
\keyword{crop}
\keyword{model}

