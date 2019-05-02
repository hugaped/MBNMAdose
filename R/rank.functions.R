# Functions for ranking in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-26


rank.MBNMA.predict <- function() {

}




#' Rank parameter estimates
#'
#' Only parameters that vary by agent can be ranked
rank.MBNMA <- function(mbnma, params=NULL, direction=1, agents=NULL) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, classes="MBNMA", add=argcheck)
  checkmate::assertChoice(params, choices = , null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(direction, choices = c(-1,1), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # If treats have not been specified then select all of them - WONT WORK IF PLACEBO NOT INCLUDED
  if (is.null(agents)) {
    agents <- c(2:mbnma[["model"]][["data"]]()[["agent"]])
  }

  if (direction==-1) {
    decreasing <- FALSE
  } else if (direction==1) {
    decreasing <- TRUE
  } else {stop("`direction` must be either -1 or 1 for ranking")}


  rank.result <- list()
  for (i in seq_along(params)) {
    if (params[i] %in% mbnma[["parameters.to.save"]]) {
      param.mod <- mbnma[["BUGSoutput"]][["sims.list"]][[params[i]]]

      # Check that selected parameter is different over multiple treatments
      if (!is.matrix(param.mod) | ncol(param.mod)<=1) {
        msg <- paste0(params[i], " does not vary by treatment and therefore cannot be ranked by treatment")
        stop(msg)
      }

      param.mod <- param.mod[,treats]
      rank.mat <- t(apply(param.mod, MARGIN=1, FUN=function(x) {
        order(order(x, decreasing = decreasing), decreasing=FALSE)
      }))
      colnames(rank.mat) <- treats

      rank.result[[params[i]]] <-
        list("summary"=sumrank(rank.mat),
             "prob.matrix"=calcprob(rank.mat, treats=treats),
             "rank.matrix"=rank.mat)

    }
}
