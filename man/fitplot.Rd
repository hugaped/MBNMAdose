% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.functions.R
\name{fitplot}
\alias{fitplot}
\title{Plot fitted values from MBNMA model}
\usage{
fitplot(
  mbnma,
  disp.obs = TRUE,
  n.iter = mbnma$BUGSoutput$n.iter,
  n.thin = mbnma$BUGSoutput$n.thin,
  ...
)
}
\arguments{
\item{mbnma}{An S3 object of class \code{"mbnma"} generated by running
a dose-response MBNMA model}

\item{disp.obs}{A boolean object to indicate whether raw data responses should be
plotted as points on the graph}

\item{n.iter}{number of total iterations per chain (including burn in;
    default: 2000)}

\item{n.thin}{thinning rate. Must be a positive integer.  Set
    \code{n.thin} > 1 to save memory and computation time if
    \code{n.iter} is large.  Default is \code{max(1, floor(n.chains *
    (n.iter-n.burnin) / 1000))} which will only thin if there are at
    least 2000 simulations.}

\item{...}{Arguments to be sent to \code{ggplot2::geom_point()} or \code{ggplot2::geom_line()}}
}
\value{
Generates a plot of fitted values from the MBNMA model and returns a list containing
the plot (as an object of \code{class(c("gg", "ggplot"))}), and a data.frame of posterior mean
fitted values for each observation.
}
\description{
Plot fitted values from MBNMA model
}
\details{
Fitted values should only be plotted for models that have converged successfully.
If fitted values (\code{theta}) have not been monitored in \code{mbnma$parameters.to.save}
then additional iterations will have to be run to get results for these.
}
\examples{
\donttest{
# Using the triptans data
network <- mbnma.network(triptans)

# Run an Emax dose-response MBNMA and predict responses
emax <- mbnma.run(network, fun=demax(), method="random")

# Plot fitted values and observed values
fitplot(emax)

# Plot fitted values only
fitplot(emax, disp.obs=FALSE)

# A data frame of fitted values can be obtained from the object
#returned by `fitplot`
fits <- fitplot(emax)
head(fits$fv)
}

}
