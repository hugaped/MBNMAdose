# Functions for assessing inconsistency in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-30


#' Node-splitting model for testing consistency
#'
#' Uses GeMTC
MBNMA.nodesplit <- function(network, level="treatment", ...) {

  # Change names to allow for GeMTC
  data.ab <- network$data.ab
  data.ab$treatment <- paste(factor(data.ab$agent, labels=network$agents),
                             data.ab$dose,
                             sep="_")

  varlist <- list(
    c("studyID", "study"),
    c("y", "mean"),
    c("se", "std.err"),
    c("r", "responders"),
    c("N", "sampleSize"),
    c("E", "exposure")
    )

  for (i in seq_along(varlist)) {
    if (varlist[[i]][1] %in% names(data.ab)) {
      names(data.ab)[names(data.ab)==varlist[[i]][1]] <- varlist[[i]][2]
    }
  }

  # Convert to GeMTC
  mtc <- mtc.network(data.ab)


}
