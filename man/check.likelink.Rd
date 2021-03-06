% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run.functions.R
\name{check.likelink}
\alias{check.likelink}
\title{Check likelihood and link function}
\usage{
check.likelink(data.ab, likelihood = NULL, link = NULL)
}
\arguments{
\item{data.ab}{A data frame of arm-level data in "long" format containing the columns:
\itemize{
\item \code{studyID} Study identifiers
\item \code{dose} Numeric data indicating the dose (must take positive values)
\item \code{agent} Agent identifiers (can be numeric, factor or character)
\item \code{y} Numeric data indicating the aggregate response for a continuous outcome. Required for
continuous data.
\item \code{se} Numeric data indicating the standard error for a given observation. Required for
continuous data.
\item \code{r} Numeric data indicating the number of responders within a study arm. Required for
binomial or poisson data.
\item \code{N} Numeric data indicating the total number of participants within a study arm. Required for
binomial data
\item \code{E} Numeric data indicating the total exposure time for participants within a study arm. Required
for poisson data.
\item \code{class} An optional column indicating a particular class code. Agents with the same identifier
must also have the same class code.
}}

\item{likelihood}{A string indicating the likelihood to use in the model. Can take either \code{"binomial"},
\code{"normal"} or \code{"poisson"}. If left as \code{NULL} the likelihood will be inferred from the data.}

\item{link}{A string indicating the link function to use in the model. Can take any link function
defined within JAGS (e.g. \code{"logit"}, \code{"log"}, \code{"probit"}, \code{"cloglog"}) or be assigned the value \code{"identity"} for
and identity link function. If left as \code{NULL} the link function will be automatically assigned based
on the likelihood.}
}
\description{
Checks that likelihood and link function is provided and confirm that the correct
form of data is provided.
}
