% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MBNMAdose.R
\docType{package}
\name{MBNMAdose-package}
\alias{MBNMAdose}
\alias{MBNMAdose-package}
\title{MBNMAdose for dose-response Model-Based Network Meta-Analysis}
\description{
\code{MBNMAdose} provides a collection of useful commands that allow users to run dose-response
Model-Based Network Meta-Analyses (MBNMA).
}
\section{Introduction}{

\code{MBNMAdose} allows meta-analysis of studies that compare multiple doses of different agents in a way that can
account for the dose-response relationship.

Whilst making use of all the available evidence in a statistically robust and biologically plausible framework,
this also can help connect networks at the agent level that may otherwise be disconnected at the dose/treatment
level, and help improve precision of estimates. It avoids "lumping" of doses that is often done in standard
Network Meta-Analysis (NMA). All models and analyses are implemented
in a Bayesian framework, following an extension of the standard NMA methodology presented by
\insertCite{lu2004}{MBNMAdose} and are run in \insertCite{jags;textual}{MBNMAdose}. For full details of
dose-response MBNMA methodology see \insertCite{mawdsley2016;textual}{MBNMAdose}. Within this package we
refer to a \strong{treatment} as a specific \strong{dose} or a specific \strong{agent}.
}

\section{Workflow}{

Functions within \code{MBNMAdose} follow a clear pattern of use:
\enumerate{
\item Load your data into the correct format using \code{mbnma.network()}
\item Analyse your data using \code{mbnma.run()}, or any of the available wrapper dose-response functions
\item Test for consistency at the treatment-level using functions like \code{nma.nodesplit()} and \code{nma.run()}
\item Examine model results using forest plots and treatment rankings
\item Use your model to predict responses using \code{predict()}
}

At each of these stages there are a number of informative plots that can be generated to help understand
the data and to make decisions regarding model fitting.
}

\examples{
# Generate an "mbnma.network" object that stores data in the correct format
network <- mbnma.network(HF2PPITT)

# Generate a network plot at the dose/treatment level
plot(network, level="treatment")

# Generate a network plot at the agent level
plot(network, level="agent", remove.loops=TRUE)

\donttest{
# Perform "split" NMA to examine dose-response relationship
nma <- nma.run(network)
plot(nma)

# Analyse data using mbnma.run()
result <- mbnma.run(network, fun="emax",
  beta.1="rel", beta.2="rel",
  method="common")

# ...or achieve the same result by using a wrapper function for mbnma.run()
result <- mbnma.emax(network,
  emax="rel", ed50="rel",
  method="common")

# Generate forest plots for model results
plot(result)

# Rank results and plot rankograms
ranks <- rank(result)
plot(ranks, params="d.emax")

# Predict responses
pred <- predict(result, E0=0.2)

# Plot predicted response with "split" NMA results displayed
plot(pred, disp.obs=TRUE, network=network)
}

}
\references{
\insertAllCited{}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4999602/}
}

}
\author{
\strong{Maintainer}: Hugo Pedder \email{hugopedder@gmail.com}

Other contributors:
\itemize{
  \item David Mawdsley [reviewer]
  \item Nicky Welton [reviewer]
  \item Sofia Dias [reviewer]
  \item Meg Bennetts [reviewer]
  \item Martin Boucher [reviewer]
  \item Adil Karim [contributor]
}

}
\keyword{internal}
