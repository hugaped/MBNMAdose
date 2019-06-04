# Functions for ranking in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-26

#' Rank predicted doses of different agents
#'
#' @inheritParams rank.MBNMA
#' @param rank.doses A list of numeric vectors. Each named element corresponds to an
#' agent (as named/coded in `predict`), and each number within the vector for that element corresponds to the dose
#' for that agent. Doses of agents specified in `rank.doses` *must* be a subset of those
#' for which responses have been predicted in `predict`. If left as `NULL` (the default)
#' then all doses of all agents in `predict` will be ranked.
#'
#' @details
#' If `predict` contains multiple predictions at dose=0, then only the first of these
#' will be included, to avoid duplicating rankings.
rank.MBNMA.predict <- function(predict, direction=1, rank.doses=NULL) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(predict, classes="MBNMA.predict", add=argcheck)
  checkmate::assertChoice(direction, choices = c(-1,1), add=argcheck)
  #checkmate::assertList(rank.doses, types="numeric", null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (direction==-1) {
    decreasing <- FALSE
  } else if (direction==1) {
    decreasing <- TRUE
  } else {stop("`direction` must be either -1 or 1 for ranking")}

  # If rank doses is null, get values from predict
  if (is.null(rank.doses)) {
    rank.doses <- list()
    incl.zero <- FALSE
    for (i in seq_along(predict$predicts)) {
      doses <- as.numeric(names(predict$predicts[[i]]))

      # Drop zero doses from all but one agent
      if (incl.zero==TRUE) {
        doses <- doses[doses!=0]
      } else if (0 %in% doses) {
        incl.zero <- TRUE
      }

      rank.doses[[names(predict$predicts)[i]]] <- doses
    }

  } else {
    # Check that rank.doses is a subset of predict
    agents.missing <- vector()
    doses.msg <- vector()
    new.ranks <- rank.doses
    for (i in seq_along(rank.doses)) {
      doses.missing <- vector()
      if (!(names(rank.doses)[i] %in% names(predict[["predicts"]]))) {
        warning(paste0("Agent ", names(rank.doses)[i], " not in `predicts` so will not be included in ranking"))
        new.ranks[[i]] <- NULL
      } else {
        doses <- as.numeric(names(predict[["predicts"]][[names(rank.doses)[i]]]))
        for (k in seq_along(rank.doses[[i]])) {
          if (!(rank.doses[[i]][k] %in% doses)) {
            doses.missing <- append(doses.missing, rank.doses[[i]][k])
            new.ranks[[i]] <- new.ranks[[i]][new.ranks[[i]]!=rank.doses[[i]][k]]
          }
        }
      }

      if (length(doses.missing)>0) {
        warning(paste0("For ", names(rank.doses)[i], " in `rank.doses`, the following doses are missing from `predict` and will not be included in ranking: ", paste(doses.missing, collapse=", ")))
      }
    }
    rank.doses <- new.ranks
  }


  treats <- vector()
  rank.mat <- NULL
  for (i in seq_along(rank.doses)) {
    for (k in seq_along(rank.doses[[i]])) {
      treats <- append(treats, paste(names(rank.doses)[i], rank.doses[[i]][k], sep="_"))
      temp <- predict$predicts[[
        names(rank.doses)[i]
        ]][[
          as.character(rank.doses[[i]][k])
          ]]
      if (is.null(rank.mat)) {
        rank.mat <- temp
      } else {
        rank.mat <- cbind(rank.mat, temp)
      }
    }
  }

  # Assign ranks
  rank.mat <- t(apply(rank.mat, MARGIN=1, FUN=function(x) {
    order(order(x, decreasing = decreasing), decreasing=FALSE)
  }))
  colnames(rank.mat) <- treats

  result <- list("summary"=sumrank(rank.mat),
                 "prob.matrix"=calcprob(rank.mat, treats=treats),
                 "rank.matrix"=rank.mat)
  result <- list("Predictions"=result)

  class(result) <- "MBNMA.rank"

  return(result)

}




#' Rank parameter estimates
#'
#' Only parameters that vary by agent/class can be ranked.
#'
#' @param to.rank A numeric vector containing the codes for the agents/classes you wish to rank.
#' If left `NULL` then all agents/classes (depending on the value assigned to `level`) in
#' the model will be ranked. Numbers must be greater than
#' 2 if placebo has been modelled, since placebo will not be included in the ranking
#' @param level Can be set to `"agent"` to rank across different agents or `"class"` to rank
#' across different classes.
#' @param params A character vector of named parameters in the model that vary by either agent
#' or class (depending on the value assigned to `level`). If left as `NULL` (the default), then
#' ranking will be calculated for all available parameters that vary by agent/class.
#'
#' @export
rank.MBNMA <- function(mbnma, params=NULL, direction=1, to.rank=NULL, level="agent") {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, classes="MBNMA", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(direction, choices = c(-1,1), add=argcheck)
  checkmate::assertNumeric(to.rank, lower = 2, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(level, choices = c("agent","class"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (mbnma$model.arg$fun %in% c("nonparam.up", "nonparam.down")) {
    stop("Ranking cannot currently be performed for non-parametric models")
  }

  if (level=="class") {
    if (is.null(mbnma[["model"]][["data"]]()[["class"]])) {
      stop("`level` has been set to `class` but classes have not been used in the model")
    }
    if (!is.null(to.rank)) {
      warning("Codes in `to.rank` correspond to class codes rather than treatment codes")
    }
  }

  # If treats have not been specified then select all of them - WONT WORK IF PLACEBO NOT INCLUDED
  starttrt <- ifelse(mbnma$agents[1]=="Placebo", 2, 1)
  codes.mod <- c(starttrt:max(mbnma[["model"]][["data"]]()[[level]], na.rm=TRUE))
  if (is.null(to.rank)) {
    to.rank <- codes.mod
  }
  if (mbnma$agents[1]=="Placebo") {
    to.rank <- to.rank-1
  }

  if (level=="agent") {
    agents <- mbnma$agents[to.rank]
  } else if (level=="class") {
    agents <- mbnma$classes[to.rank]
  }


  if (direction==-1) {
    decreasing <- FALSE
  } else if (direction==1) {
    decreasing <- TRUE
  } else {stop("`direction` must be either -1 or 1 for ranking")}

  if (is.null(params)) {
    for (i in seq_along(mbnma[["parameters.to.save"]])) {
      if (length(mbnma$BUGSoutput$long.short[i][[1]])==length(codes.mod)) {
        params <- append(params, mbnma[["parameters.to.save"]][i])
      }
    }
  }

  rank.result <- list()
  for (i in seq_along(params)) {
    if (params[i] %in% mbnma[["parameters.to.save"]]) {
      param.mod <- mbnma[["BUGSoutput"]][["sims.list"]][[params[i]]]

      # Check that selected parameter is different over multiple treatments
      if (!is.matrix(param.mod) | ncol(param.mod)!=length(codes.mod)) {
        msg <- paste0(params[i], " does not vary by ", level, " and therefore cannot be ranked")
        stop(msg)
      }

      param.mod <- param.mod[,to.rank]
      rank.mat <- t(apply(param.mod, MARGIN=1, FUN=function(x) {
        order(order(x, decreasing = decreasing), decreasing=FALSE)
      }))
      #colnames(rank.mat) <- to.rank
      colnames(rank.mat) <- agents

      rank.result[[params[i]]] <-
        list("summary"=sumrank(rank.mat),
             #"prob.matrix"=calcprob(rank.mat, treats=to.rank),
             "prob.matrix"=calcprob(rank.mat, treats=agents),
             "rank.matrix"=rank.mat)

    }
  }
  class(rank.result) <- "MBNMA.rank"

  if (length(rank.result)==0) {
    stop(paste0("There are no parameters saved in the model that vary by ", level))
  }

  return(rank.result)
}





#' Calculates a matrix of ranking probabilities from a matrix of treatment/agent/class
#' rankings
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
sumrank <- function(rank.mat) {
  if (is.null(colnames(rank.mat))) {
    colnames(rank.mat) <- c(1:ncol(rank.mat))
  }

  quantiles.rank <- apply(X=rank.mat, MARGIN = 2,
                          function(x) quantile(x, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
  summary.rank <- data.frame(
    "rank.param"=colnames(rank.mat),
    "mean"= apply(X=rank.mat, MARGIN = 2, mean),
    "sd"= apply(X=rank.mat, MARGIN = 2, sd)
  )
  summary.rank <- cbind(summary.rank, t(quantiles.rank))
  rownames(summary.rank) <- NULL

  return(summary.rank)
}
