# Functions for ranking in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-26

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "do", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))

#' Set rank as a method
#'
#' @param x An object on which to apply the rank method
#' @param ... Arguments to be passed to methods
#'
#' @export
rank <- function (x, ...) {
  UseMethod("rank", x)
}



#' Calculates a matrix of ranking probabilities from a matrix of treatment/agent/class
#' rankings
#'
#' @param rank.mat Numeric matrix of treatment/agent/class rankings
#' @param treats A numeric vector of treatment codes for which to calculate ranking probabilities
#' @noRd
calcprob <- function(rank.mat, treats=NULL) {
  NT <- ncol(rank.mat)
  rank.prob <- vector(length=NT)

  for (c in 1:NT) {
    pos.vec <- vector()
    for (r in 1:NT) {
      pos.vec <- append(pos.vec,
                        length(rank.mat[rank.mat[,c]==r,c])/nrow(rank.mat))
    }
    rank.prob <- cbind(rank.prob, pos.vec)
  }
  rank.prob <- rank.prob[,-1]

  if (!is.null(treats)) {
    colnames(rank.prob) <- treats
  }

  return("rank.prob"=rank.prob)
}





#' Generates a summary data frame from a matrix of treatment/agent/class rankings
#'
#' @inheritParams calcprob
#'
#' @noRd
sumrank <- function(rank.mat) {
  if (is.null(colnames(rank.mat))) {
    colnames(rank.mat) <- c(1:ncol(rank.mat))
  }

  quantiles.rank <- apply(X=rank.mat, MARGIN = 2,
                          function(x) stats::quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  summary.rank <- data.frame(
    "rank.param"=colnames(rank.mat),
    "mean"= apply(X=rank.mat, MARGIN = 2, mean),
    "sd"= apply(X=rank.mat, MARGIN = 2, stats::sd)
  )
  summary.rank <- cbind(summary.rank, t(quantiles.rank))
  rownames(summary.rank) <- NULL

  return(summary.rank)
}






calcauc <- function(df) {
  str <- paste(paste(c(df$Var1[1], df$Var1, df$Var1[nrow(df)], df$Var1[1]),
                     c(0, df$value, 0, 0),
                     sep=" "),
               collapse=",")
  str <- paste0("POLYGON((", str, "))")

  polygon <- rgeos::readWKT(str)
  #temp$auc <- rgeos::gArea(polygon)
  auc <- rgeos::gArea(polygon)
  #auc <- rep(auc, nrow(temp))
  return(auc)
}
