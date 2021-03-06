% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simjm_benchmark.R
\name{simjm_benchmark}
\alias{simjm_benchmark}
\title{Benchmark "True" analysis of data simulated by \code{simjm} function}
\usage{
simjm_benchmark(data, surv_model = "Cox", marker_model = "RE",
  corr = "Low")
}
\arguments{
\item{data}{A data.frame produced by \code{\link[survtd]{simjm}}.}

\item{surv_model}{Time-to-event model that was used to generate \code{data}. Options are "Cox" for Cox proportional
hazard model and "Add" for semi-parametric additive hazards model.}

\item{marker_model}{Multiple marker model that was used to generate \code{data}. Options are "RE" for the correlated random effects model and "PN" for the
product normal model.}

\item{corr}{Degree of marginal pairwise correlations between the markers that was used to generate \code{data}. Options are "Low", "Mod" and "High".}
}
\value{
Returns regression coefficient estimates for each covariate based on the "True" analysis, along with 95%
  confidence intervals and p-values.
}
\description{
\code{simjm_benchmark} fits Cox proportional hazards and semi-parametric additive hazards models
to data generated by \code{\link[survtd]{simjm}} using the perfect data, i.e. the true values of the time-dependent markers.
}
\details{
The function is used to perform what is referred to as the "True" analysis in Moreno-Betancur et al. (2017) which is based on perfect data,
i.e. the true values of the markers at each of the event times as drawn from the simulation model. This provides a benchmark for analyses performed
with error-polluted and incomplete data using \code{\link[survtd]{survtd}}.
}
\examples{

  dat<-simjm(n=200,surv_model="Cox",marker_model="PN",
             MErr="High",Miss="None",effects="Weak",corr="Low")

  simjm_benchmark(dat,surv_model="Cox",marker_model="PN",corr="Low")

}
\references{
Moreno-Betancur M, Carlin JB, Brilleman SL, Tanamas S, Peeters A, Wolfe R (2017). Survival analysis
with time-dependent covariates subject to missing data or measurement error: Multiple Imputation for Joint Modeling (MIJM).
\emph{Biostatistics} [Epub ahead of print 12 Oct 2017].
}
