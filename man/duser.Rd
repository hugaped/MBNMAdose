% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dose.functions.R
\name{duser}
\alias{duser}
\title{User-defined dose-response function}
\usage{
duser(fun, beta.1 = "rel", beta.2 = "rel", beta.3 = "rel", beta.4 = "rel")
}
\arguments{
\item{fun}{A formula specifying any relationship including \code{dose} and
one/several of: \code{beta.1}, \code{beta.2}, \code{beta.3}, \code{beta.4}.}

\item{beta.1}{Pooling for the 1st coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}

\item{beta.2}{Pooling for the 2nd coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}

\item{beta.3}{Pooling for the 3rd coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}

\item{beta.4}{Pooling for the 4th coefficient. Can take \code{"rel"}, \code{"common"}, \code{"random"} or be
assigned a numeric value (see details).}
}
\value{
An object of \code{class("dosefun")}
}
\description{
User-defined dose-response function
}
\section{Dose-response parameters}{
\tabular{ll}{
   \strong{Argument} \tab \strong{Model specification} \cr
   \code{"rel"} \tab Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. \cr
   \code{"common"} \tab Implies that all agents share the same common effect for this dose-response parameter. \cr
   \code{"random"} \tab Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. \cr
   \code{numeric()} \tab Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. \cr
}


When relative effects are modelled on more than one dose-response parameter,
correlation between them is automatically estimated using a vague inverse-Wishart prior.
This prior can be made slightly more informative by specifying the scale matrix \code{omega}
and by changing the degrees of freedom of the inverse-Wishart prior
using the \code{priors} argument in \code{mbnma.run()}.
}

\examples{

dr <- ~ beta.1 * (1/(dose+1)) + beta.2 * dose^2
duser(fun=dr,
  beta.1="common", beta.2="rel")

}
\references{
\insertAllCited
}
