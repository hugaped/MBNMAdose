# Functions for printing in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-25



#' Produces a summary data frame from an MBNMA.predict object
#'
#' @return A data frame containing posterior summary statistics from predicted responses
#'   from a dose-response MBNMA model
#' @export
summary.MBNMA.predict <- function(predict) {

  checkmate::assertClass(predict, "MBNMA.predict")

  predict <- predict[["predicts"]]

  output <- data.frame()
  for (i in seq_along(predict)) {
    for (k in seq_along(predict[[i]])) {
      quant <- quantile(predict[[i]][[k]], probs=c(0.025,0.25,0.5,0.75,0.975))
      df <- data.frame("agent"=names(predict)[i],
                       "dose"=as.numeric(as.character(names(predict[[i]])[k])),
                       "mean"=mean(predict[[i]][[k]]),
                       "sd"=sd(predict[[i]][[k]])
      )
      output <- rbind(output, cbind(df, t(quant)))
    }
  }

  return(output)
}


#' Print summary information from an MBNMA.predict object
#' @export
print.MBNMA.predict <- function(predict) {

  checkmate::assertClass(predict, "MBNMA.predict")

  #predict <- predict[["predicts"]]

  sum.df <- summary(predict)

  agents <- unique(sum.df$agent)

  # Check if agent labels are numeric or not
  named.a <- TRUE
  if (any(is.na(suppressWarnings(as.numeric(agents))))) {
    named.a <- FALSE
  }

  head <- "#### Predicted doses ####"
  info <- vector()
  for (i in seq_along(agents)) {
    doses <- sum.df$dose[sum.df$agent==agents[i]]
    if (named.a==FALSE) {
      info <- append(info, paste0("Agent ", i, ": ", paste(doses, collapse=", ")))
    } else {
      info <- append(info, paste0(agents[i], ": ", paste(doses, collapse=", ")))
    }
  }
  out <- c(head, info)
  cat(paste(out, collapse="\n"))
}




summary.MBNMA.rank <- function(mbnma.rank) {
  checkmate::assertClass(mbnma.rank, "MBNMA.rank")

  output <- list()
  for (i in seq_along(mbnma.rank)) {
    output[[names(mbnma.rank)[i]]] <- mbnma.rank[[i]]$summary
  }
  return(output)
}

print.MBNMA.rank <- function(mbnma.rank) {
  head <- "\n###### Ranking of dose-response MBNMA ######"

  intro <- c()
  if ("Predictions" %in% names(mbnma.rank)) {
    intro <- c(intro, "Includes ranking of predictions from dose-response MBNMA")
  }
  if (any(grepl("^d\\.", names(mbnma.rank)))) {
    add <- "Includes ranking of relative effects from dose-response MBNMA:"
    add <- paste(add,
                 paste(names(mbnma.rank)[grepl("^d\\.", names(mbnma.rank))], collapse="\t"),
                 sep="\n")
    intro <- c(intro, add)
  }
  if (any(grepl("^D\\.", names(mbnma.rank)))) {
    add <- "Includes ranking of class effects from dose-response MBNMA:"
    add <- paste(add,
                 paste(names(mbnma.rank)[grepl("^D\\.", names(mbnma.rank))], collapse="\t"),
                 sep="\n")
    intro <- c(intro, add)
  }

  rankinfo <- paste(nrow(mbnma.rank[[1]]$summary), "parameters ranked", sep=" ")
  if (mbnma.rank[[1]]$direction==1) {
    rankinfo <- paste(rankinfo, "with positive responses ranked as `better`", sep=" ")
  } else if (mbnma.rank[[1]]$direction==-1) {
    rankinfo <- paste(rankinfo, "with negative responses ranked as `better`", sep=" ")
  }

  intro <- paste(intro, collapse="\n")
  out <- paste(head, intro, rankinfo, sep="\n\n")

  return(cat(out))
}



print.MBNMA.nodesplit <- function(nodesplit) {

}

summary.MBNMA.nodesplit <- function(nodesplit) {

}



print.MBNMA <- function(mbnma) {

}

summary.MBNMA <- function(mbnma) {

}






rhat.warning <- function(mbnma, cutoff=1.02) {
  rhats <- mbnma$BUGSoutput$summary[,colnames(mbnma$BUGSoutput$summary)=="Rhat"]
  rhats <- names(rhats)[rhats>cutoff]
  if (length(rhats)>0) {
    msg <- paste0("The following parameters have Rhat values > ",
                  cutoff,
                  " which could be due to convergence issues:\n")
    warning(paste0(msg, paste(rhats, collapse="\n")))
  }
}
