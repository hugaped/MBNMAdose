# Functions for assessing inconsistency in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-30


#' Unrelated Mean Effects model for testing consistency
MBNMA.ume <- function(level="treatment", ...) {

  # Fit UME split NMA and compare to split NMA

  splitNMA <- NMA.run(network=network, method=method,
                      likelihood=likelihood, link=link, ...)
}
