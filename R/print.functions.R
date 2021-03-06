# Functions for printing in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-25


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

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
      quant <- stats::quantile(predict[[i]][[k]], probs=c(0.025,0.25,0.5,0.75,0.975))
      df <- data.frame("agent"=names(predict)[i],
                       "dose"=as.numeric(as.character(names(predict[[i]])[k])),
                       "mean"=mean(predict[[i]][[k]]),
                       "sd"=stats::sd(predict[[i]][[k]]),
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
#' @inheritParams summary.mbnma.predict
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
  cat(paste(out, collapse="\n"), ...)
}



#' Generates summary data frames for an mbnma.rank object
#'
#' @param object An object of `class("mbnma.rank")` generated by `rank.mbnma()`
#' @param ... additional arguments affecting the summary produced
#'
#' @return A list in which each element represents a parameter that has been ranked
#' in `mbnma.rank` and contains a data frame of summary ranking results.
#'
#' @export
summary.mbnma.rank <- function(object, ...) {
  checkmate::assertClass(object, "mbnma.rank")

  output <- list(...)
  for (i in seq_along(object)) {
    output[[names(object)[i]]] <- object[[i]]$summary
  }
  return(output)
}




#' Prints summary information about an mbnma.rank object
#'
#' @param x An object of class `"mbnma.rank"` generated by `rank.mbnma()`
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.mbnma.rank <- function(x, ...) {
  checkmate::assertClass(x, "mbnma.rank")

  head <- "\n###### Ranking of dose-response MBNMA ######"

  intro <- c()
  if ("Predictions" %in% names(x)) {
    intro <- c(intro, "Includes ranking of predictions from dose-response MBNMA")
  }
  if (any(grepl("^d\\.", names(x)))) {
    add <- "Includes ranking of relative effects from dose-response MBNMA:"
    add <- paste(add,
                 paste(names(x)[grepl("^d\\.", names(x))], collapse="\t"),
                 sep="\n")
    intro <- c(intro, add)
  }
  if (any(grepl("^D\\.", names(x)))) {
    add <- "Includes ranking of class effects from dose-response MBNMA:"
    add <- paste(add,
                 paste(names(x)[grepl("^D\\.", names(x))], collapse="\t"),
                 sep="\n")
    intro <- c(intro, add)
  }

  rankinfo <- paste(nrow(x[[1]]$summary), "parameters ranked", sep=" ")
  if (x[[1]]$direction==1) {
    rankinfo <- paste(rankinfo, "with positive responses ranked as `better`", sep=" ")
  } else if (x[[1]]$direction==-1) {
    rankinfo <- paste(rankinfo, "with negative responses ranked as `better`", sep=" ")
  }

  intro <- paste(intro, collapse="\n")
  out <- paste(head, intro, rankinfo, sep="\n\n")

  return(cat(out, ...))
}


neatCrI <- function(vals, digits=3) {
  vals <- signif(vals, digits = digits)
  neat <- paste0(vals[2], " (", vals[1], ", ", vals[3], ")")
  return(neat)
}


#' Prints summary results from an nma.nodesplit object
#'
#' @param x An object of `class("nma.nodesplit")`
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.nma.nodesplit <- function(x, ...) {
  checkmate::assertClass(x, "nma.nodesplit")

  width <- "\t\t"
  output <- "========================================\nNode-splitting analysis of inconsistency\n========================================\n"

  comparisons <- names(x)
  colnam <- "comparison\tp.value\t\t\tMedian (95% CrI)"
  paramsect <- colnam
  for (i in seq_along(comparisons)) {
    pval <- signif(x[[i]]$p.values,
                   max(3L, getOption("digits") - 3L))
    tab <- x[[i]]$quantiles

    heading <- paste(names(x)[i], pval, sep=width)
    direct <- paste("-> direct", "", neatCrI(tab$direct), sep=width)
    indirect <- paste("-> indirect", "", neatCrI(tab$indirect), sep=width)
    nma <- paste("-> NMA", "\t", neatCrI(tab$nma), sep=width)

    out <- paste(heading, direct, indirect, nma, "", sep="\n")

    paramsect <- paste(paramsect, out, sep="\n")
  }
  output <- append(output, paramsect)
  cat(output, ...)
}



#' Generates a summary data frame for nma.nodesplit objects
#'
#' @param object An object of `class("nma.nodesplit")`
#' @param ... further arguments passed to or from other methods
#'
#' @export
summary.nma.nodesplit <- function(object, ...) {
  checkmate::assertClass(object, "nma.nodesplit")

  if ("quantiles" %in% names(object[[1]])) {
    type <- "dose"
  } else if ("quantiles" %in% names(object[[1]][[1]])) {
    type <- "time"
  } else {stop("Type of MBNMA cannot be identified (dose or time)")}

  sum.mat <- matrix(ncol=3)
  comp <- vector()
  time.param <- vector()
  evidence <- vector()
  pvals <- vector()

  for (i in seq_along(object)) {

    if (type=="dose") {
      post <- object[[i]]$quantiles

      sum.mat <- rbind(sum.mat, post$direct)
      evidence <- c(evidence, "Direct")

      sum.mat <- rbind(sum.mat, post$indirect)
      evidence <- c(evidence, "Indirect")

      sum.mat <- rbind(sum.mat, post$nma)
      evidence <- c(evidence, "NMA")

      pvals <- c(pvals, rep(object[[i]]$p.values, 3))
      comp <- c(comp, rep(names(object)[i], 3))
    } else if (type=="time") {
      for (k in seq_along(object[[i]])) {
        post <- object[[i]][[k]]$quantiles

        sum.mat <- rbind(sum.mat, post$direct)
        evidence <- c(evidence, "Direct")

        sum.mat <- rbind(sum.mat, post$indirect)
        evidence <- c(evidence, "Indirect")

        pvals <- c(pvals, rep(object[[i]][[k]]$p.values, 2))
        time.param <- c(time.param, rep(names(object[[i]])[k], 2))
        comp <- c(comp, rep(names(object)[i], 2))
      }
    }
  }
  sum.mat <- round(sum.mat[-1,], digits = max(3L, getOption("digits") - 5L))
  pvals <- round(pvals, max(3L, getOption("digits") - 5L))

  if (type=="time") {
    sum.df <- data.frame(comp, time.param,
                         evidence, sum.mat[,2],
                         sum.mat[,1], sum.mat[,3],
                         pvals
    )

    names(sum.df) <- c("Comparison", "Parameter", "Evidence", "Median",
                       "2.5%", "97.5%", "p.value")
  } else if (type=="dose") {
    sum.df <- data.frame(comp,
                         evidence, sum.mat[,2],
                         sum.mat[,1], sum.mat[,3],
                         pvals, ...
    )

    names(sum.df) <- c("Comparison", "Evidence", "Median",
                       "2.5%", "97.5%", "p.value")
  }
  return(sum.df)
}



get.beta.names <- function(mbnma) {
  betanames <- list()
  for (i in 1:4) {
    if (!is.null(mbnma$model.arg[[paste0("beta.",i)]])) {
      if (is.null(mbnma$model.arg$arg.params)) {
        betanames[[paste0("beta.", i)]] <- paste0("beta.", i)
      } else if (!is.null(mbnma$model.arg$arg.params)) {
        temp <- mbnma$model.arg$arg.params$wrap.params[
          mbnma$model.arg$arg.params$run.params==paste0("beta.", i)
          ]
        betanames[[paste0("beta.", i)]] <- temp
      }
    }
  }

  return(betanames)
}



get.timeparam.str <- function(mbnma, beta=NULL, param="d") {
  betanames <- get.beta.names(mbnma)

  if (grepl("beta", betanames[[beta]])) {
    temp <- strsplit(betanames[[beta]], split="\\.")[[1]][2]
  } else {
    temp <- betanames[[beta]]
  }

  match <- paste0("^", param, "\\.", temp, "(\\[[0-9]+\\])?")

  sum.mat <- mbnma$BUGSoutput$summary[grepl(match, rownames(mbnma$BUGSoutput$summary)),
                                      c(3,5,7)]

  # Check for UME
  if (any(grepl("\\[[0-9]+,[0-9]+\\]", rownames(sum.mat)))) {
    sum.mat <- NULL
  }

  if (length(sum.mat)>0) {
    tab.str <- c()

    if (is.matrix(sum.mat)) {
      if (any(grepl("^d\\..+\\[1\\]", rownames(sum.mat)[1]) |
              grepl("^D\\..+\\[1\\]", rownames(sum.mat)[1]))) {
        tab.str <- paste(tab.str,
                         paste(rownames(sum.mat)[1], "Network reference",
                               sep="\t"),
                         sep="\n"
        )
        count <- 2
      } else if (any(grepl("^beta\\..+", rownames(sum.mat)[1]) |
                     grepl("^BETA\\..+", rownames(sum.mat)[1]))) {
        tab.str <- paste(tab.str,
                         paste(rownames(sum.mat)[1],
                               neatCrI(sum.mat),
                               sep="\t"),
                         sep="\n"
        )
        count <- 1
      } else {count <-1}

      for (i in count:nrow(sum.mat)) {
        tab.str <- paste(tab.str,
                         paste(rownames(sum.mat)[i], neatCrI(sum.mat[i,]),
                               sep="\t"),
                         sep="\n"
        )
      }
    } else if (is.vector(sum.mat)) {
      tab.str <- paste(tab.str,
                       paste(rownames(mbnma$BUGSoutput$summary)[grepl(match, rownames(mbnma$BUGSoutput$summary))],
                             neatCrI(sum.mat),
                             sep="\t\t"),
                       sep="\n"
      )
    }

    return(tab.str)
  } else {
    return(NULL)
  }
}




#' Neatly prints a summary of results
#'
#' @inheritParams predict.mbnma
#' @noRd
print.treat.str <- function(mbnma) {
  betanames <- get.beta.names(mbnma)

  treat.sect <- c()
  # DR parameters for each agent (generate treat.str)
  for (i in seq_along(betanames)) {
    if (!(names(betanames)[i] %in% names(mbnma$model.arg$class.effect))) {

      sect.head <- paste("####", betanames[[i]], "dose-response parameter results ####", sep=" ")

      data.head <- paste("Parameter", "Median (95%CrI)", sep="\t")
      data.head <- paste(data.head, "---------------------------------", sep="\n")

      if (is.character(mbnma$model.arg[[names(betanames)[i]]])) {

        if (mbnma$model.arg[[names(betanames)[i]]]=="rel") {
          param <- "d"
        } else if (mbnma$model.arg[[names(betanames)[i]]] %in% c("common", "random")) {
          param <- "beta"
        }

        data.tab <- get.timeparam.str(mbnma, beta=paste0("beta.",i), param = param)
        if (!is.null(data.tab)) {
          data.str <- paste(data.head,
                            data.tab,
                            sep="")
        }

        if (mbnma$model.arg[[names(betanames)[i]]]=="random") {
          sd.name <- paste0("sd.", betanames[[i]])
          sd.vals <- neatCrI(mbnma$BUGSoutput$summary[
            rownames(mbnma$BUGSoutput$summary)==sd.name, c(3,5,7)])
          sd.str <- paste(sd.name, sd.vals, sep="\t")
          data.str <- paste(data.str, sd.str, sep="\n")
        }

        # Parameters on exponential scale
        if (mbnma$model.arg$fun=="emax" | mbnma$model.arg$fun=="emax.hill") {
          if (names(betanames)[i] %in% c("beta.2", "et50")) {
            sect.head <- paste(sect.head,
                               "Parameter modelled on exponential scale to ensure it takes positive values on the natural scale", sep="\n")
          }
        }

        # String for pooling
        if (mbnma$model.arg[[names(betanames)[i]]]=="rel") {
          pool <- "relative effects"
        } else if (mbnma$model.arg[[names(betanames)[i]]] %in% c("common", "random")) {
          pool <- "absolute single parameter"
        }
        pool.str <- paste("Pooling:", pool, "\n", sep=" ")

        treat.str <- paste(sect.head, pool.str, data.str, sep="\n")

      } else if (is.numeric(mbnma$model.arg[[names(betanames)[i]]])) {
        data.str <- paste("Assigned a numeric value:",
                          mbnma$model.arg[[names(betanames)[i]]],
                          sep=" ")

        treat.str <- paste(sect.head, data.str)
      }

      treat.sect <- paste(treat.sect, treat.str, "", sep="\n\n")
    }
  }
  return(treat.sect)
}



print.method.sect <- function(mbnma) {
  # String for method
  data.head <- paste("Parameter", "Median (95%CrI)", sep="\t\t\t\t\t")
  data.head <- paste(data.head, "-----------------------------------------------------------------------", sep="\n")

  if (mbnma$model.arg$method=="common") {
    method <- "Common (fixed) effects estimated for relative effects"
    method <- paste0(method, "\n\n")
  } else if (mbnma$model.arg$method=="random") {
    method <- "Random effects estimated for relative effects"

    # Check if >1 relative effect
    betas <- c(mbnma$model.arg$beta.1, mbnma$model.arg$beta.2, mbnma$model.arg$beta.3)
    if (table(betas)[names(table(betas))=="rel"]>=2) {
      method <- paste(method, "Correlation modelled between relative effect dose-response parameters", sep="\n")
    }

    temp <- mbnma$BUGSoutput$summary[grepl("^sd$", rownames(mbnma$BUGSoutput$summary)),
                                     c(3,5,7)]
    if (!is.vector(temp)) {stop("temp should only be length 1")}

    sd.str <- paste(data.head,
                    paste("Between-study SD for relative effects", neatCrI(temp), sep="\t\t"),
                    sep="\n")

    method <- paste0(method, "\n\n")
    method <- paste0(method, sd.str)
  }

  method.str <- paste("Method:", method, sep=" ")
  method.str <- paste("\n\n#### Pooling method ####", method.str, sep="\n\n")
  return(method.str)
}



print.class.str <- function(mbnma) {
  if (length(mbnma$model.arg$class.effect)>0) {
    classes <- mbnma$model.arg$class.effect

    head <- "\n#### Class effects ####\n"
    data.head <- paste("Parameter", "Median (95%CrI)", sep="\t")
    data.head <- paste(data.head, "---------------------------------", sep="\n")

    data.str <- c(data.head)
    class.str <- c()
    sd.str <- NULL
    for (i in seq_along(classes)) {
      betaparam <- names(classes)[i]

      if (!is.null(mbnma$model.arg$arg.params)) {
        wrapparam <- mbnma$model.arg$arg.params$wrap.params[
          mbnma$model.arg$arg.params$run.params==names(classes)[i]]

        class.str <- c(class.str, paste("Class effect on",
                                        paste0(wrapparam,":"), classes[[i]], "\n",
                                        sep=" "))
      } else {
        wrapparam <- ""
        class.str <- c(class.str, paste("Class effect on",
                                        paste0(betaparam,":"), classes[[i]], "\n",
                                        sep=" "))
      }

      if (mbnma$model.arg[[betaparam]]=="rel") {
        data.str <- paste0(data.str,
                           get.timeparam.str(mbnma, beta = names(classes)[i], param = "D"))
      } else {
        stop("`class.effects` seem to be specified on dose-response parameters not modelled by relative effects")
      }

      if (classes[[i]]=="random") {
        sd.str <- "\n# Within-class SD\n"
        match.1 <- paste0("(", strsplit(names(classes)[i], split="\\.")[[1]][2], ")?")
        match.2 <- paste0("(", wrapparam, ")?")
        match <- paste0("^sd\\.[A-Z]+\\.", match.1, match.2)
        temp <- mbnma$BUGSoutput$summary[grepl(match, rownames(mbnma$BUGSoutput$summary)),
                                         c(3,5,7)]
        if (!is.vector(temp)) {stop("temp should only be length 1")}

        sd.str <- paste(sd.str, data.head,
                        paste(rownames(mbnma$BUGSoutput$summary)[grepl(match, rownames(mbnma$BUGSoutput$summary))],
                              neatCrI(temp), sep="\t\t"),
                        sep="\n")
      }
    }
    class.sect <- c(head,
                    paste(class.str, collapse="\n"),
                    paste(data.str, collapse="\n"))

    if (!is.null(sd.str)) {
      class.sect <- c(class.sect, sd.str, "\n")
    }

    class.sect <- paste(class.sect, collapse="\n")

    return(class.sect)
  }
}



print.modfit.str <- function(mbnma) {
  totresdev.str <- c()

  head <- "#### Model Fit Statistics ####\n"

  # pD
  pd.str <- "Effective number of parameters:"
  if (mbnma$model.arg$pd=="pv") {
    pd <- "pD (pV) calculated using the rule, pD = var(deviance)/2 ="
  } else if (mbnma$model.arg$pd=="plugin") {
    pd <- "pD calculated using the plug-in method ="
  } else if (mbnma$model.arg$pd=="pd.kl") {
    pd <- "pD calculated using the Kullback-Leibler divergence ="
  } else if (mbnma$model.arg$pd=="popt") {
    pd <- "pD calculated using an optimism adjustment ="
  }
  pd.str <- paste(pd.str, paste(pd, round(mbnma$BUGSoutput$pD,0), sep=" "), sep="\n")

  # Deviance
  dev <- mbnma$BUGSoutput$summary[
    rownames(mbnma$BUGSoutput$summary)=="deviance", 5]
  dev.str <- paste("Deviance =", round(dev, 0), sep=" ")

  # Totresdev
  if ("totresdev" %in% mbnma$parameters.to.save) {
    totresdev <- round(
      mbnma$BUGSoutput$summary[
        rownames(mbnma$BUGSoutput$summary)=="totresdev", 5],
      0)
  } else {
    totresdev <- "NOT MONITORED IN MODEL"
  }
  totresdev.str <- paste("Residual deviance =", totresdev, sep=" ")

  dic <- mbnma$BUGSoutput$DIC
  dic.str <- paste("Deviance Information Criterion (DIC) =", round(dic, 0), "\n", sep=" ")

  modfit.sect <- paste(head, pd.str, dev.str, totresdev.str, dic.str, sep="\n")
  return(modfit.sect)
}



#' Print summary of MBNMA results to the console
#' @param object An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param ... additional arguments affecting the summary produced
#'
#' @export
summary.mbnma <- function(object, ...) {
  checkmate::assertClass(object, "mbnma")

  # State that function does not work if "parameters.to.save" has been specified
  if (!is.null(object$model.arg$parameters.to.save)) {
    stop("Cannot use `summary()` method if `parameters.to.save` have been assigned. Use `print()` instead.")
  }
  if (object$model.arg$fun %in% c("nonparam.up", "nonparam.down")) {
    stop("Cannot use `summary()` method for non-parametric dose-response functions. Use `print()` instead.")
  }

  # Check for rhat < 1.02
  rhat.warning(object)

  ##### Overall section #####
  # Print title
  title <- "========================================\nDose-response MBNMA\n========================================\n"

  # Print DR function
  overall.sect <- paste("Dose-response function:", object$model.arg$fun, sep=" ")
  overall.sect <- paste(title, overall.sect, sep="\n")

  # Print method section
  method.sect <- print.method.sect(object)

  # Print treatment-level section
  treat.sect <- print.treat.str(object)

  # Class-effect section
  class.sect <- print.class.str(object)

  # Model fit statistics section
  modfit.sect <- print.modfit.str(object)

  output <- paste(overall.sect, treat.sect, method.sect, "\n", class.sect, "\n\n", modfit.sect, sep="")
  cat(output, ...)
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





#' Print mbnma.network information to the console
#'
#' @param x An object of class `mbnma.network`.
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.mbnma.network <- function(x,...) {
  nn <- names(x)
  ll <- length(x)
  if (length(nn) != ll)
    nn <- paste("Component", seq.int(ll))
  for (i in seq_len(ll)) {
    cat(nn[i], ":\n")
    if (is.data.frame((x[[i]]))) {
      print(x[[i]], max=ncol(x[[i]])*6, ...)
    } else {
      print(x[[i]], ...)
    }
    cat("\n")
  }
  invisible(x)
}
