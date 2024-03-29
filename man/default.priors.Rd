% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write.jags.R
\name{default.priors}
\alias{default.priors}
\title{Sets default priors for JAGS model code}
\usage{
default.priors(
  fun = dloglin(),
  UME = FALSE,
  regress.mat = NULL,
  regress.effect = "common",
  om = list(rel = 5, abs = 10)
)
}
\arguments{
\item{fun}{An object of \code{class("dosefun")} that specifies a functional form to be assigned to the
dose-response. See Details.}

\item{UME}{A boolean object to indicate whether to fit an Unrelated Mean Effects model
that does not assume consistency and so can be used to test if the consistency
assumption is valid.}

\item{regress.mat}{A Nstudy x Ncovariate design matrix of meta-regression covariates}

\item{regress.effect}{Indicates whether effect modification should be assumed to be
\code{"common"} (assumed to be equal versus Placebo throughout the network),
\code{"random"} (assumed to be exchangeable versus Placebo throughout the network),
\code{"agent"} (assumed to be equal versus Placebo within each agent), or
\code{"class"} (assumed to be equal versus Placebo within each class).}

\item{om}{a list with two elements that report the maximum relative (\code{"rel"}) and maximum absolute (\code{"abs"}) efficacies
on the link scale.}
}
\value{
A list, each element of which is a named JAGS snippet
corresponding to a prior in the MBNMA JAGS code.
}
\description{
This function creates JAGS code snippets for default MBNMA model priors.
}
\examples{
\donttest{
default.priors(fun=demax())
}

}
