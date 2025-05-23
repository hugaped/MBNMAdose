# Functions for printing in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-25


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))



#' Print results as string with credible intervals in brackets
#' @noRd
neatCrI <- function(vals, digits=3) {
  vals <- signif(vals, digits = digits)
  neat <- paste0(vals[2], " (", vals[1], ", ", vals[3], ")")
  return(neat)
}



#' Print a warning if any monitored parameters have rhat above the cutoff in an MBNMA model
#' @noRd
rhat.warning <- function(mbnma, cutoff=1.2) {
  rhats <- mbnma$BUGSoutput$summary[,colnames(mbnma$BUGSoutput$summary)=="Rhat"]
  rhats <- names(rhats)[rhats>cutoff]
  if (length(rhats)>0) {
    msg <- paste0("The following parameters have Rhat values > ",
                  cutoff,
                  "\nwhich could be due to convergence issues:\n")
    warning(paste0(msg, paste(rhats, collapse="\n")))
  }
}




#' Neatly prints a summary of results
#'
#' @inheritParams predict.mbnma
#' @noRd
treat.str <- function(mbnma, digits=3, ...) {

  fun <- mbnma$model.arg$fun
  datasum <- as.data.frame(cbind(mbnma$BUGSoutput$summary[,5],
                                 mbnma$BUGSoutput$summary[,3],
                                 mbnma$BUGSoutput$summary[,7]))
  datasum$param <- rownames(datasum)

  treat.sect <- c()

  # For each dose-response parameter generate results
  for (i in seq_along(fun$params)) {
    headbeta <- fun$params[i]
    sect.head <- crayon::bold(paste(headbeta, "dose-response parameter results\n\n", sep=" "))
    cat(crayon::underline(sect.head))

    # if (fun$params[i] %in% c("ed50", "hill")) {
    #   cat("Parameter modelled on exponential scale to ensure it takes positive values\non the natural scale\n")
    # }

    agents <- mbnma$network$agents[mbnma$network$agents!="Placebo"]

    if (fun$apool[i] %in% "rel") {
      cat("Pooling: relative effects for each agent")

      if (TRUE %in% mbnma$model.arg$UME | fun$params[i] %in% mbnma$model.arg$UME) {
        cat("Unrelated Mean Effects modelled for this parameter")
        cat("\n\n---- RESULTS TOO LONG FOR SUMMARY ----\n\n\n")
      } else {
        trt.df <- datasum[grepl(paste0("^", fun$params[i], "\\["), datasum$param),]

        if (nrow(trt.df)==length(agents)) {
          trt.df$agents <- agents
        } else if (length(fun$name)>1) {

          listlen <- lengths(fun$paramlist)
          count <- 0
          k <- 0
          while (count<i) {
            count <- count + listlen[k+1]
            k <- k+1
          }
          if ("Placebo" %in% mbnma$network$agents) {
            temppos <- fun$posvec[-1]
          } else {
            temppos <- fun$posvec
          }
          subagents <- agents[which(temppos==k)]
          trt.df$agents <- subagents

        } else {
          stop("Agent named cannot be matched to parameters")
        }

        trt.df <- trt.df[,c(5,4,1,2,3)]
        rownames(trt.df) <- NULL
        print(knitr::kable(trt.df, col.names = c("Agent", "Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
        cat("\n\n")
      }

    } else if (fun$apool[i] %in% c("common", "random")) {
      cat("Pooling: single parameter across all agents in the network")
      trt.df <- datasum[grepl(paste0("^", fun$params[i], "$"), datasum$param),]

      trt.df <- trt.df[,c(4,1,2,3)]
      rownames(trt.df) <- NULL
      print(knitr::kable(trt.df, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
      cat("\n\n")

    } else if (suppressWarnings(!is.na(as.numeric(fun$apool[i])))) {
      data.str <- paste("Assigned a numeric value:",
                        fun$apool[i],
                        sep=" ")
      cat(data.str)
      cat("\n\n")
    }

    if (fun$apool[i] %in% "random") {

      trt.df <- datasum[grepl(paste0("sd\\.", fun$params[i]), datasum$param),]
      trt.df <- trt.df[,c(4,1,2,3)]
      rownames(trt.df) <- NULL
      cat("Between-study SD for random (exchangeable) dose-response parameter:")
      print(knitr::kable(trt.df, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
      cat("\n\n")

    }


  }

  invisible(mbnma)
}


#' Neatly prints a summary of the method used in the model
#'
#' @noRd
method.str <- function(mbnma) {
  # String for method
  data.head <- paste("Parameter", "Median (95%CrI)", sep="\t\t\t\t\t")
  data.head <- paste(crayon::bold(data.head, "-----------------------------------------------------------------------", sep="\n"))

  fun <- mbnma$model.arg$fun
  method <- vector()

  # Check if >1 relative effect
  if (length(fun$name)==1 & mbnma$model.arg$cor==TRUE) {
    if (sum(fun$apool %in% "rel")>=2) {
      method <- paste0(method, "Correlation modelled between relative effect dose-response parameters")
    }
  }

  if (mbnma$model.arg$method=="common") {
    method <- paste0(method, "Common (fixed) effects estimated for relative effects\n")
  } else if (mbnma$model.arg$method=="random") {
    method <- paste0(method, "Random effects estimated for relative effects")

    temp <- mbnma$BUGSoutput$summary[grepl("^sd$", rownames(mbnma$BUGSoutput$summary)),
                                     c(3,5,7)]
    if (!is.vector(temp)) {stop("temp should only be length 1")}

    sd.str <- paste(data.head,
                    paste("Between-study SD for relative effects", neatCrI(temp), sep="\t\t"),
                    sep="\n")

    method <- paste0(method, "\n\n")
    method <- paste0(method, sd.str)
  }

  method.sect <- paste("Method:", method, sep=" ")
  method.sect <- paste(crayon::bold(crayon::underline("\n\nPooling method")), method.sect, "", sep="\n\n")
  return(method.sect)
}







#' Neatly prints class results
#' @noRd
class.str <- function(mbnma, digits=4, ...) {

  if (length(mbnma$model.arg$class.effect)>0) {

    cat(crayon::bold(crayon::underline("\nClass effects\n")))

    classef <- mbnma$model.arg$class.effect

    classes <- mbnma$network$classes
    if ("Placebo" %in% mbnma$network$agents) {
      classes <- classes[-1]
    }

    fun <- mbnma$model.arg$fun
    datasum <- as.data.frame(cbind(mbnma$BUGSoutput$summary[,5],
                                   mbnma$BUGSoutput$summary[,3],
                                   mbnma$BUGSoutput$summary[,7]))
    datasum$param <- rownames(datasum)

    # For each dose-response parameter generate results
    for (i in seq_along(classef)) {
      sect.head <- paste("Class effect results for:", toupper(names(classef)[i]), "\n", sep=" ")
      cat(sect.head)
      cat(paste0("Model assumes ", classef[[i]], " class effects"))

      class.df <- datasum[grepl(paste0("^",toupper(names(classef)[i])), datasum$param),]

      if (nrow(class.df)==length(classes)) {
        class.df$classes <- classes
      } else if (length(fun$name)>1) {
        stop("Class names cannot be matched to parameters from agent-specific dose-response functions")
      } else {
        stop("Class names cannot be matched to parameters")
      }

      class.df <- class.df[,c(5,4,1,2,3)]
      rownames(class.df) <- NULL
      print(knitr::kable(class.df, col.names = c("Class", "Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
      cat("\n\n")

      if (classef[[i]] %in% "random") {
        class.df <- datasum[grepl(paste0("sd\\.", toupper(names(classes)[i])), datasum$param),]
        class.df <- class.df[,c(4,1,2,3)]
        rownames(class.df) <- NULL
        cat("Between-study SD for random (exchangeable) class effects:")
        print(knitr::kable(class.df, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
        cat("\n\n")
      }

    }

    invisible(mbnma)
  }
}


#' Neatly prints model fit details
#' @noRd
modfit.str <- function(mbnma) {
  totresdev.str <- c()

  cat(crayon::bold(crayon::underline("Model Fit Statistics\n")))

  # pD
  pd.str <- "Effective number of parameters:"
  if (mbnma$model.arg$pD==FALSE) {
    pd <- "pD (pV) calculated using the rule, pD = var(deviance)/2 ="
    pd.str <- paste(pd.str, paste(pd, round(mbnma$BUGSoutput$pV,0), sep=" "), sep="\n")
  } else if (mbnma$model.arg$pD==TRUE) {
    pd <- "pD calculated using the Kullback-Leibler divergence ="
    pd.str <- paste(pd.str, paste(pd, round(mbnma$BUGSoutput$pD,0), sep=" "), sep="\n")
  }

  # Deviance
  dev <- mbnma$BUGSoutput$summary[
    rownames(mbnma$BUGSoutput$summary)=="deviance", 5]
  cat(paste("Deviance =", round(dev, 1), sep=" "))
  cat("\n")

  # Totresdev
  if ("totresdev" %in% mbnma$parameters.to.save) {
    totresdev <- round(
      mbnma$BUGSoutput$summary[
        rownames(mbnma$BUGSoutput$summary)=="totresdev", 5],
      1)
  } else {
    totresdev <- crayon::red("NOT MONITORED IN MODEL")
  }
  cat(paste("Residual deviance =", totresdev, sep=" "))
  cat("\n")

  dic <- mbnma$BUGSoutput$DIC
  cat(paste("Deviance Information Criterion (DIC) =", round(dic, 1), "\n", sep=" "))

}




#' Neatly prints model regression details
#' @noRd
regress.str <- function(mbnma, digits=4, ...) {

  if (!is.null(mbnma$model.arg$regress)) {

    cat(crayon::bold(crayon::underline("Meta-regression\n\n")))

    cat(paste0("Covariates interacting with study-level relative effects: ", paste(crayon::bold(colnames(mbnma$model.arg$regress.mat)), collapse=", "), "\n"))

    # Create summary data frame
    datasum <- as.data.frame(cbind(mbnma$BUGSoutput$summary[,5],
                                   mbnma$BUGSoutput$summary[,3],
                                   mbnma$BUGSoutput$summary[,7]))
    datasum$param <- rownames(datasum)

    if ("common" %in% mbnma$model.arg$regress.effect) {
      str <- "Common (identical) covariate-by-treatment effects"
      labs <- "Common effect"
      labs.head <- "Regression effect"

    } else if ("random" %in% mbnma$model.arg$regress.effect) {
      str <- "Random (exchangeable) covariate-by-treatment effects"
      labs <- "Random effect"
      labs.head <- "Regression effect"

    } else if ("independent" %in% mbnma$model.arg$regress.effect) {
      str <- "Independent covariate-by-treatment effects"
      labs <- mbnma$network$treatments[mbnma$network$treatments!="Placebo_0"]
      labs.head <- "Treatment"

    } else if ("agent" %in% mbnma$model.arg$regress.effect) {
      str <- "Common (identical) covariate-by-agent effects"
      labs <- mbnma$network$agents[mbnma$network$agents!="Placebo"]
      labs.head <- "Agent"

    } else if ("class" %in% mbnma$model.arg$regress.effect) {
      str <- "Common (identical) covariate-by-class effects"
      labs <- mbnma$network$classes[mbnma$network$classes!="Placebo"]
      labs.head <- "Class"

    }
    cat(paste0(str, "\n"))

    trt.df <- datasum[grepl("^B\\.", datasum$param),]

    trt.df$labs <- rep(labs, ncol(mbnma$model.arg$regress.mat))

    trt.df <- trt.df[,c(5,4,1,2,3)]
    rownames(trt.df) <- NULL
    print(knitr::kable(trt.df, col.names = c(labs.head, "Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
    cat("\n\n")


    # Random effects
    if ("random" %in% mbnma$model.arg$regress.effect) {
      method <- cat("Standard deviation for random covariate-by-treatment effects")

      trt.df <- datasum[grepl("^sd\\.B\\.", datasum$param),]

      trt.df <- trt.df[,c(4,1,2,3)]
      rownames(trt.df) <- NULL
      print(knitr::kable(trt.df, col.names = c("Parameter", "Median", "2.5%", "97.5%"), digits=digits, ...))
      cat("\n\n")
    }
  }
}




#' Neatly prints heading section
#' @noRd
overall.str <- function(mbnma) {

  # Print title
  cat(crayon::bold("========================================\nDose-response MBNMA\n========================================\n\n"))

  # Print likelihood and link function
  cat(paste0("Likelihood: ", mbnma$model.arg$likelihood, "\n"))

  if ("smd" %in% mbnma$model.arg$link) {
    cat("Link function: identity\nTreatment effects synthesised using standardised mean differences\n")
  } else {
    cat(paste0("Link function: ", mbnma$model.arg$link, "\n"))
  }

  # Print DR function
  if (length(mbnma$model.arg$fun$name)==1) {
    cat(paste("Dose-response function:", mbnma$model.arg$fun$name, sep=" "))
  } else if (length(mbnma$model.arg$fun$name)>1) {
    drtab <- matrix(mbnma$model.arg$fun$name[mbnma$model.arg$fun$posvec],
                    ncol=1)

    paramvec <- sapply(mbnma$model.arg$fun$posvec, FUN = function(x) {
      paste(names(mbnma$model.arg$fun$paramlist[[x]]), collapse=", ")
    })

    drtab <- cbind(drtab, paramvec)

    fun.df <- data.frame(Agents=mbnma$network$agents,
                         Function=drtab[,1],
                         Parameters=drtab[,2])

    cat("Dose-response functions:")
    print(knitr::kable(fun.df))
  }

  if (any(mbnma$model.arg$fun$name == "user")) {
    cat("\nuser.fun:", unique(mbnma$model.arg$jags[mbnma$model.arg$fun$name %in% "user"]))
  }
}
