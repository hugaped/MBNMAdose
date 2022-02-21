##############################################
#### Functions for class("mbnma.predict") ####
##############################################

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))

#' Rank predicted doses of different agents
#'
#' Ranks predictions at different doses from best to worst.
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
#' network <- mbnma.network(triptans)
#'
#' # Rank all predictions from a log-linear dose-response MBNMA
#' loglin <- mbnma.run(network, fun=dloglin())
#' pred <- predict(loglin, E0 = 0.5)
#' rank <- rank(pred)
#' summary(rank)
#'
#' # Rank selected predictions from an Emax dose-response MBNMA
#' emax <- mbnma.run(network, fun=demax(), method="random")
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
rank.mbnma.predict <- function(x, lower_better=TRUE, rank.doses=NULL, ...) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, classes="mbnma.predict", add=argcheck)
  checkmate::assertLogical(lower_better, add=argcheck)
  checkmate::reportAssertions(argcheck)

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


  # Generate matrix of rankings
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
    order(order(x, decreasing = lower_better), decreasing=FALSE)
  }))
  colnames(rank.mat) <- treats

  result <- list("summary"=sumrank(rank.mat),
                 "prob.matrix"=calcprob(rank.mat, treats=treats),
                 "rank.matrix"=rank.mat,
                 "lower_better"=lower_better)
  result <- list("Predictions"=result)

  class(result) <- "mbnma.rank"

  return(result)

}




#' Plots predicted responses from a dose-response MBNMA model
#'
#' Plots predicted responses on the natural scale from a dose-response MBNMA model.
#'
#' @param x An object of class `"mbnma.predict"` generated by
#'   `predict("mbnma")`
#' @param disp.obs A boolean object to indicate whether to show the location of observed doses
#'   in the data on the 95\\% credible intervals of the predicted dose-response curves as shaded regions (`TRUE`)
#'   or not (`FALSE`). If set to `TRUE` the original network object used for the model
#'   **must** be specified in `network`.
#' @param overlay.split A boolean object indicating whether to overlay a line
#'   showing the split (treatment-level) NMA results on the plot (`TRUE`) or not (`FALSE`). This will
#'   require automatic running of a split NMA model.
#'   For `overlay.split=TRUE` the original network object used for the model
#'   **must** be specified in `network`.
#' @param method Indicates the type of split (treatment-level) NMA to perform when `overlay.split=TRUE`. Can
#'   take either `"common"` or `"random"`.
#' @param agent.labs A character vector of agent labels to display on plots. If
#'   left as `NULL` (the default) the names of agents will be taken from `predict`. The position of
#'   each label corresponds to each element of `predict`. The number of labels must equal
#'   the number of active agents in `predict`. If placebo / dose=0 data is included in the predictions
#'   then a label for placebo **should not** be included in `agent.labs`. It will not be shown
#'   in the final plot since placebo is the point within each plot at which dose = 0 (rather
#'   than a separate agent).
#' @param ... Arguments for `ggplot2`
#' @inheritParams mbnma.run
#' @inheritParams plot.mbnma.rank
#' @inheritParams ggplot2::facet_wrap
#'
#' @details For the S3 method `plot()`, it is advisable to ensure predictions in
#'   `predict` are estimated using a sufficient number of doses to ensure a smooth
#'   predicted dose-response curve. If `disp.obs = TRUE` it is
#'   advisable to ensure predictions in `predict` are estimated using an even
#'   sequence of time points to avoid misrepresentation of shaded densities.
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(triptans)
#'
#' # Run an Emax dose-response MBNMA and predict responses
#' emax <- mbnma.run(network, fun=demax(), method="random")
#' pred <- predict(emax, E0 = 0.5)
#' plot(pred)
#'
#' # Display observed doses on the plot
#' plot(pred, disp.obs=TRUE)
#'
#' # Display split NMA results on the plot
#' plot(pred, overlay.split=TRUE)
#'
#' # Split NMA results estimated using random treatment effects model
#' plot(pred, overlay.split=TRUE, method="random")
#'
#' # Add agent labels
#' plot(pred, agent.labs=c("Elet", "Suma", "Frov", "Almo", "Zolmi",
#'       "Nara", "Riza"))
#'
#' # These labels will throw an error because "Placebo" is included in agent.labs when
#' #it will not be plotted as a separate panel
#' #### ERROR ####
#' #plot(pred, agent.labs=c("Placebo", "Elet", "Suma", "Frov", "Almo", "Zolmi",
#' #      "Nara", "Riza"))
#'
#'
#' # If insufficient predictions are made across dose-response function
#' # then the plotted responses are less smooth and can be misleading
#' pred <- predict(emax, E0 = 0.5, n.doses=3)
#' plot(pred)
#' }
#'
#' @export
plot.mbnma.predict <- function(x, disp.obs=FALSE,
                               overlay.split=FALSE, method="common",
                               agent.labs=NULL, scales="free_x", ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma.predict", add=argcheck)
  checkmate::assertLogical(disp.obs, len=1, add=argcheck)
  checkmate::assertLogical(overlay.split, len=1, add=argcheck)
  checkmate::assertCharacter(agent.labs, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(scales, c("free_x", "fixed"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  sum.pred <- summary(x)
  sum.df <- sum.pred

  # Check agent.labs and that the number of labels there are is correct
  if (!is.null(agent.labs)) {
    if ("Placebo" %in% names(x[["predicts"]])) {
      agent.labs <- c("Placebo", agent.labs)
    }
    if (length(agent.labs)!=
        length(x[["predicts"]])) {
      stop("The length of `agent.labs` does not equal the number of active agents that responses have been predicted for in `x`. `Placebo` will be added automatically and should not be given a label")
    }
    sum.df$agent <- factor(sum.df$agent, labels=agent.labs)
  }

  sum.df <- sum.df[!(sum.df$agent %in% c("1", "Placebo")),]


  # Plot predictions
  g <- ggplot2::ggplot(sum.df, ggplot2::aes(x=as.numeric(as.character(dose)),
                                            y=`50%`, ymin=`2.5%`, ymax=`97.5%`), ...)


  # Plot observed data as shaded regions
  if (disp.obs==TRUE) {
    network <- x$network
    checkmate::assertClass(network, "mbnma.network", null.ok=TRUE)

    # Check that predict labels and agent labels in network are consistent
    if (!all(sum.pred$agent %in% network$agents)) {
      stop("Agent labels in `network` differ from those in `pred`")
    }

    gobs <- disp.obs(g=g, network=network, predict=x,
                     col="green", max.col.scale=NULL)

    g <- gobs[[1]]
    cols <- gobs[[2]]

  }
  if (overlay.split==TRUE) {
    network <- x$network
    checkmate::assertClass(network, "mbnma.network", null.ok=TRUE)

    # Check that placebo is included (or dose=0 in networks without placebo)
    if (network$agents[1]!="Placebo") {
      stop("Placebo required in `network` for calculation of relative effects in split NMA")
    }

    # Check that at least one prediction is at a dose=0
    if (!(0 %in% sum.pred$dose)) {
      stop("`x` must include a predicted response at dose = 0 for at least one agent")
    }

    g <- overlay.split(g=g, network=network, E0=x$E0, method=method,
                       likelihood = x[["likelihood"]],
                       link = x[["link"]])

  }

  # Add overlayed lines and legends
  g <- g +
    ggplot2::geom_line(ggplot2::aes(y=`2.5%`, linetype="95% Interval")) +
    ggplot2::geom_line(ggplot2::aes(y=`97.5%`, linetype="95% Interval")) +
    ggplot2::geom_line(ggplot2::aes(linetype="Posterior Median"))

  g <- g + ggplot2::facet_wrap(~agent, scales=scales) +
    ggplot2::labs(y="Predicted response", x="Dose")

  g <- g + ggplot2::scale_linetype_manual(name="",
                                          values=c("Posterior Median"="solid",
                                                   "95% Interval"="dashed")) +
    theme_mbnma()

  return(g)
}





#' Produces a summary data frame from an mbnma.predict object
#'
#' @param object An object of `class("mbnma.predict)"` generated by
#'   `predict("mbnma")`
#' @param ... additional arguments affecting the summary produced.
#'
#' @return A data frame containing posterior summary statistics from predicted responses
#'   from a dose-response MBNMA model
#' @export
summary.mbnma.predict <- function(object, ...) {

  checkmate::assertClass(object, "mbnma.predict")

  predict <- object[["predicts"]]

  output <- data.frame()
  for (i in seq_along(predict)) {
    for (k in seq_along(predict[[i]])) {
      quant <- stats::quantile(predict[[i]][[k]], probs=c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)
      df <- data.frame("agent"=names(predict)[i],
                       "dose"=as.numeric(as.character(names(predict[[i]])[k])),
                       "mean"=mean(predict[[i]][[k]]),
                       "sd"=stats::sd(predict[[i]][[k]]), stringsAsFactors = TRUE,
                       ...
      )
      output <- rbind(output, cbind(df, t(quant)))
    }
  }

  return(output)
}


#' Print summary information from an mbnma.predict object
#'
#' @param x An object of `class("mbnma.predict")` generated by `predict.mbnma()`
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.mbnma.predict <- function(x, ...) {

  checkmate::assertClass(x, "mbnma.predict")

  #x <- x[["predicts"]]

  sum.df <- summary(x)

  agents <- unique(sum.df$agent)

  # Check if agent labels are numeric or not
  named.a <- TRUE
  if (any(is.na(suppressWarnings(as.numeric(agents))))) {
    named.a <- FALSE
  }

  head <- crayon::bold("=====================\nPredicted doses\n=====================\n")
  info <- vector()
  for (i in seq_along(agents)) {
    doses <- sum.df$dose[sum.df$agent==agents[i]]
    if (named.a==FALSE) {
      info <- append(info, paste0(crayon::bold(paste0("Agent ", i, ": ")), paste(doses, collapse=", ")))
    } else {
      info <- append(info, paste0(crayon::bold(paste0(agents[i], ": ")), paste(doses, collapse=", ")))
    }
  }
  out <- c(head, info)
  cat(paste(out, collapse="\n"), ...)
}
