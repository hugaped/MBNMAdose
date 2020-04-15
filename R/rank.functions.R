# Functions for ranking in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-26

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Set rank as a method
#'
#' @param x An object on which to apply the rank method
#' @param ... Arguments to be passed to methods
#'
#' @export
rank <- function (x, ...) {
  UseMethod("rank", x)
}

#' Rank predicted doses of different agents
#'
#' @inheritParams rank.mbnma
#' @inheritParams predict.mbnma
#' @param rank.doses A list of numeric vectors. Each named element corresponds to an
#' agent (as named/coded in `predict`), and each number within the vector for that element corresponds to the dose
#' for that agent. Doses of agents specified in `rank.doses` **must** be a subset of those
#' for which responses have been predicted in `predict`. If left as `NULL` (the default)
#' then all doses of all agents in `predict` will be ranked.
#' @inheritParams plot.mbnma.predict
#'
#' @details
#' If `predict` contains multiple predictions at dose=0, then only the first of these
#' will be included, to avoid duplicating rankings.
#'
#' @return An object of `class("mbnma.rank")` which is a list containing a summary data
#' frame, a matrix of rankings for each MCMC iteration, and a matrix of probabilities
#' that each agent has a particular rank, for each parameter that has been ranked.
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Rank predictions from a linear dose-response MBNMA
#' linear <- mbnma.run(network, fun="linear")
#' pred <- predict(linear, E0 = 0.5)
#' rank <- rank(pred)
#' summary(rank)
#'
#' # Rank selected predictions from an Emax dose-response MBNMA
#' emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#' doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
#' pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
#'             exact.doses=doses)
#' rank <- rank(pred,
#'             rank.doses=list("eletriptan"=c(0,2), "rizatriptan"=2))
#'
#' # Print and generate summary data frame for `mbnma.rank` object
#' summary(rank)
#' print(rank)
#'
#' # Plot `mbnma.rank` object
#' plot(rank)
#' }
#'
#' @export
rank.mbnma.predict <- function(x, direction=1, rank.doses=NULL, ...) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, classes="mbnma.predict", add=argcheck)
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
    for (i in seq_along(x$predicts)) {
      doses <- as.numeric(names(x$predicts[[i]]))

      # Drop zero doses from all but one agent
      if (incl.zero==TRUE) {
        doses <- doses[doses!=0]
      } else if (0 %in% doses) {
        incl.zero <- TRUE
      }

      rank.doses[[names(x$predicts)[i]]] <- doses
    }

  } else if (!is.null(rank.doses)) {
    # Check that rank.doses is a subset of predict
    agents.missing <- vector()
    doses.msg <- vector()
    new.ranks <- rank.doses
    for (i in seq_along(rank.doses)) {
      doses.missing <- vector()
      if (!(names(rank.doses)[i] %in% names(x[["predicts"]]))) {
        stop(paste0("Agent ", names(rank.doses)[i], " not in `predicts` so cannot be included in ranking"))
        new.ranks[[i]] <- NULL
      } else {
        doses <- as.numeric(names(x[["predicts"]][[names(rank.doses)[i]]]))
        for (k in seq_along(rank.doses[[i]])) {
          if (!(rank.doses[[i]][k] %in% doses)) {
            doses.missing <- append(doses.missing, rank.doses[[i]][k])
            new.ranks[[i]] <- new.ranks[[i]][new.ranks[[i]]!=rank.doses[[i]][k]]
          }
        }
      }

      if (length(doses.missing)>0) {
        stop(paste0("For ", names(rank.doses)[i], " in `rank.doses`, the following doses are missing from `x` and cannot be included in ranking: ", paste(doses.missing, collapse=", ")))
      }
    }
    rank.doses <- new.ranks
  }


  treats <- vector()
  rank.mat <- NULL
  for (i in seq_along(rank.doses)) {
    for (k in seq_along(rank.doses[[i]])) {
      treats <- append(treats, paste(names(rank.doses)[i], rank.doses[[i]][k], sep="_"))
      temp <- x$predicts[[
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
                 "rank.matrix"=rank.mat,
                 "direction"=direction)
  result <- list("Predictions"=result)

  class(result) <- "mbnma.rank"

  return(result)

}




#' Rank parameter estimates
#'
#' Only parameters that vary by agent/class can be ranked.
#'
#' @param direction Indicates whether negative responses are better (taking the
#'   value `-1`) or positive responses are better (taking the value `1`)
#' @param to.rank A numeric vector containing the codes for the agents/classes you wish to rank.
#' If left `NULL` then all agents/classes (depending on the value assigned to `level`) in
#' the model will be ranked. Included codes must be greater than
#' `2` if placebo has been modelled, since placebo cannot be included in the ranking
#' @param level Can be set to `"agent"` to rank across different agents or `"class"` to rank
#' across different classes.
#' @param params A character vector of named parameters in the model that vary by either agent
#' or class (depending on the value assigned to `level`). If left as `NULL` (the default), then
#' ranking will be calculated for all available parameters that vary by agent/class.
#' @param ... Arguments to be passed to methods
#' @inheritParams predict.mbnma
#' @inheritParams rank
#'
#' @details Ranking cannot currently be performed on non-parametric dose-response MBNMA
#'
#' @return An object of `class("mbnma.rank")` which is a list containing a summary data
#' frame, a matrix of rankings for each MCMC iteration, a matrix of probabilities
#' that each agent has a particular rank, and a matrix of cumulative ranking probabilities
#' for each agent, for each parameter that has been ranked.
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Rank selected agents from a linear dose-response MBNMA
#' linear <- mbnma.run(network, fun="linear")
#' ranks <- rank(linear, to.rank=c("zolmitriptan", "eletriptan", "sumatriptan"))
#' summary(ranks)
#'
#' # Rank only ED50 parameters from an Emax dose-response MBNMA
#' emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#' ranks <- rank(emax, params="d.ed50")
#' summary(ranks)
#'
#'
#' #### Ranking by class ####
#' # Generate some classes for the data
#' class.df <- HF2PPITT
#' class.df$class <- ifelse(class.df$agent=="placebo", "placebo", "active1")
#' class.df$class <- ifelse(class.df$agent=="eletriptan", "active2", class.df$class)
#' netclass <- mbnma.network(class.df)
#' emax <- mbnma.emax(netclass, emax="rel", ed50="rel", method="random",
#'             class.effect=list("ed50"="common"))
#'
#' # Rank by class, with negative responses being "better"
#' ranks <- rank(emax, level="class", direction=-1)
#' print(ranks)
#'
#' # Print and generate summary data frame for `mbnma.rank` object
#' summary(ranks)
#' print(ranks)
#'
#' # Plot `mbnma.rank` object
#' plot(ranks)
#' }
#'
#' @export
rank.mbnma <- function(x, params=NULL, direction=1, level="agent", to.rank=NULL, ...) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, classes="mbnma", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(direction, choices = c(-1,1), add=argcheck)
  #checkmate::assertNumeric(to.rank, lower = 2, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(level, choices = c("agent","class"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (any(x$model.arg$fun %in% c("nonparam.up", "nonparam.down"))) {
    stop("Ranking cannot currently be performed for non-parametric models")
  }
  if (length(x$model.arg$fun)>1) {
    stop("Ranking cannot currently be performed for models with multiple dose-response functions")
  }

  # Change agent/class to agents/classes
  levels <- ifelse(level=="agent", "agents", "classes")

  if (level=="class") {
    if (is.null(x[["model"]][["data"]]()[["class"]])) {
      stop("`level` has been set to `class` but classes have not been used in the model")
    }
    if (!is.null(to.rank)) {
      warning("Codes in `to.rank` correspond to class codes rather than treatment codes")
    }
  }

  # If treats have not been specified then select all of them - WONT WORK IF PLACEBO NOT INCLUDED
  starttrt <- ifelse(x$network$agents[1]=="Placebo", 2, 1)
  codes.mod <- c(starttrt:max(x[["model"]][["data"]]()[[level]], na.rm=TRUE))
  if (is.null(to.rank)) {
    to.rank <- codes.mod
  } else if (is.numeric(to.rank)) {
    if (!all(to.rank %in% seq(1:max(x[["model"]][["data"]]()[[level]], na.rm=TRUE)))) {
      stop("`to.rank` codes must match those in the dataset for either `agent` or `class`")
    }
  } else if (is.character(to.rank)) {
    if (!all(to.rank %in% x$network[[levels]])) {
      stop("`to.rank` agent/class names must match those in the network for either `agent` or `class`")
    }
    to.rank <- as.numeric(factor(to.rank, levels=x$network[[levels]]))
  }

  if (x$network$agents[1]=="Placebo") {
    to.rank <- to.rank-1
    if (any(to.rank==0)) {
      warning("Placebo (d[1] or D[1]) cannot be included in the ranking for relative effects and will therefore be excluded")
      to.rank <- to.rank[to.rank!=0]
    }
    agents <- x$network[[levels]][to.rank+1]
  } else {
    agents <- x$network[[levels]][to.rank]
  }


  if (direction==-1) {
    decreasing <- FALSE
  } else if (direction==1) {
    decreasing <- TRUE
  } else {stop("`direction` must be either -1 or 1 for ranking")}

  if (is.null(params)) {
    for (i in seq_along(x$BUGSoutput$root.short)) {
      if (length(x$BUGSoutput$long.short[i][[1]])==length(codes.mod)) {
        params <- append(params, x$BUGSoutput$root.short[i])
      }
    }
  } else {
    for (i in seq_along(params)) {
      if (!(params[i] %in% x[["parameters.to.save"]])) {
        stop(paste0(params[i], " has not been monitored by the model. `params` can only include model parameters that have been monitored and vary by agent/class"))
      }
    }
  }

  rank.result <- list()
  for (i in seq_along(params)) {
    if (params[i] %in% x[["parameters.to.save"]]) {
      param.mod <- x[["BUGSoutput"]][["sims.list"]][[params[i]]]

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

      # Calculate ranking probabilities
      prob.mat <- calcprob(rank.mat, treats=agents)

      # Calculate cumulative ranking probabilities
      cum.mat <- apply(prob.mat, MARGIN=2,
                       FUN=function(col) {cumsum(col)})

      rank.result[[params[i]]] <-
        list("summary"=sumrank(rank.mat),
             #"prob.matrix"=calcprob(rank.mat, treats=to.rank),
             "prob.matrix"=prob.mat,
             "rank.matrix"=rank.mat,
             "cum.matrix"=cum.mat,
             "direction"=direction)

    }
  }
  class(rank.result) <- "mbnma.rank"

  if (length(rank.result)==0) {
    stop(paste0("There are no parameters saved in the model that vary by ", level))
  }

  return(rank.result)
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
