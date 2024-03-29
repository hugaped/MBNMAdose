% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.functions.R
\name{ref.synth}
\alias{ref.synth}
\title{Synthesise single arm dose = 0 / placebo studies to estimate E0}
\usage{
ref.synth(
  data.ab,
  mbnma,
  synth = "fixed",
  n.iter = mbnma$BUGSoutput$n.iter,
  n.burnin = mbnma$BUGSoutput$n.burnin,
  n.thin = mbnma$BUGSoutput$n.thin,
  n.chains = mbnma$BUGSoutput$n.chains,
  ...
)
}
\arguments{
\item{data.ab}{A data frame of arm-level data in "long" format containing the
columns:
\itemize{
\item \code{studyID} Study identifiers
\item \code{y} Numeric data indicating the aggregate response for a continuous outcome. Required for
continuous data.
\item \code{se} Numeric data indicating the standard error for a given observation. Required for
continuous data.
\item \code{r} Numeric data indicating the number of responders within a study arm. Required for
binomial or poisson data.
\item \code{n} Numeric data indicating the total number of participants within a study arm. Required for
binomial data
\item \code{E} Numeric data indicating the total exposure time for participants within a study arm. Required
for poisson data.
}}

\item{mbnma}{An S3 object of class \code{"mbnma"} generated by running
a dose-response MBNMA model}

\item{synth}{A character object that can take the value \code{"fixed"} or \code{"random"} to
specify the the type of pooling to use for synthesis of \code{E0} if a data frame
has been provided for it. Using \code{"random"} rather
than \code{"fixed"} for \code{synth} will result in wider 95\\\% CrI for predictions.}

\item{n.iter}{number of total iterations per chain (including burn in;
    default: 2000)}

\item{n.burnin}{length of burn in, i.e. number of iterations to
    discard at the beginning. Default is \code{n.iter/2}, that is,
    discarding the first half of the simulations. If n.burnin is 0,
    \code{jags()} will run 100 iterations for adaption.}

\item{n.thin}{thinning rate. Must be a positive integer.  Set
    \code{n.thin} > 1 to save memory and computation time if
    \code{n.iter} is large.  Default is \code{max(1, floor(n.chains *
    (n.iter-n.burnin) / 1000))} which will only thin if there are at
    least 2000 simulations.}

\item{n.chains}{number of Markov chains (default: 3)}

\item{...}{Arguments to be sent to \code{\link[R2jags:jags]{R2jags::jags()}} for synthesis of the network
reference treatment effect (using \code{\link[=ref.synth]{ref.synth()}})}
}
\value{
A list of named elements corresponding to E0 and the between-study standard deviation for
E0 if \code{synth="random"}. Each element contains the full MCMC results from the synthesis.
}
\description{
Synthesises single arm studies to estimate E0. Used in predicting responses from a
dose-response MBNMA.
}
\details{
\code{data.ab} can be a collection of studies that closely resemble the
population of interest intended for the prediction, which could be
different to those used to estimate the MBNMA model, and could include
single arms of RCTs or observational studies. If other data is not
available, the data used to estimate the MBNMA model can be used by
selecting only the studies and arms that investigate dose = 0 (placebo).

Defaults for \code{n.iter}, \code{n.burnin}, \code{n.thin} and \code{n.chains} are those used to estimate
\code{mbnma}.
}
\examples{
\donttest{
# Using the triptans data
network <- mbnma.network(triptans)

# Run an Emax dose-response MBNMA
emax <- mbnma.run(network, fun=demax(), method="random")

# Data frame for synthesis can be taken from placebo arms
ref.df <- triptans[triptans$agent=="placebo",]

# Meta-analyse placebo studies using fixed treatment effects
E0 <- ref.synth(ref.df, emax, synth="fixed")
names(E0)

# Meta-analyse placebo studies using random treatment effects
E0 <- ref.synth(ref.df, emax, synth="random")
names(E0)
}

}
