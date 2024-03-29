% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inconsistency.functions.R
\name{inconsistency.loops}
\alias{inconsistency.loops}
\title{Identify comparisons in loops that fulfill criteria for node-splitting}
\usage{
inconsistency.loops(df, checkindirect = TRUE, incldr = FALSE)
}
\arguments{
\item{df}{A data frame containing variables \code{studyID} and \code{treatment} (as
numeric codes) that indicate which treatments are used in which studies. If \code{checkindirect = TRUE}
then variables \code{agent} and \code{dose} are also required.}

\item{checkindirect}{A boolean object to indicate whether or not to perform an additional
check to ensure network remains connected even after dropping direct evidence on a comparison.
Default is \code{TRUE} and should be kept as \code{TRUE} if working with dose-response data, though this requires
further computational iterations to confirm. If set to \code{FALSE}, additional comparisons may be identified, though computation will be much more
rapid.}

\item{incldr}{A boolean object indicating whether or not to allow for indirect evidence contributions via
the dose-response relationship. This can be used when node-splitting in dose-response MBNMA to allow
for a greater number of potential loops in which to check for consistency.}
}
\value{
A data frame of comparisons that are informed by direct and indirect
evidence from independent sources. Each row of the data frame is a
different treatment comparison. Numerical codes in \code{t1} and \code{t2} correspond
to treatment codes. \code{path} indicates the treatment codes that connect the
shortest path of indirect evidence.

If \code{incldr=TRUE} then \code{path} may indicate \code{doseresp} for some comparisons.
These are comparisons for which indirect evidence is only available via the
dose-response relationship. The two numbers given after (e.g. \verb{3 2}) indicate the
number of doses available in the indirect evidence with which to estimate the
dose-response function for the treatments in \code{t1} and \code{t2} respectively/
}
\description{
Identify comparisons informed by both direct and indirect evidence from
independent sources, which therefore fulfill the criteria for testing for
inconsistency via node-splitting.
}
\details{
Similar to \code{gemtc::mtc.nodesplit.comparisons()} but uses a fixed
reference treatment and therefore identifies fewer loops in which to test for
inconsistency. Heterogeneity can also be parameterised as inconsistency and
so testing for inconsistency in additional loops whilst changing the
reference treatment would also be identifying heterogeneity. Depends on
\code{\link[igraph]{igraph}}.
}
\examples{
\donttest{
# Identify comparisons informed by direct and indirect evidence
#in triptans dataset
network <- mbnma.network(triptans)
inconsistency.loops(network$data.ab)

# Include indirect evidence via dose-response relationship
inconsistency.loops(network$data.ab, incldr=TRUE)
}


# Do not perform additional connectivity check on data
data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
            treatment=c(1,2,1,3,2,3,3,4,1,2,4)
            )
inconsistency.loops(data, checkindirect=FALSE)
}
\references{
\insertAllCited{}
}
