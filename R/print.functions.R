# Functions for printing in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-25


## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))



#' Print results as string with credible intervals in brackets
neatCrI <- function(vals, digits=3) {
  vals <- signif(vals, digits = digits)
  neat <- paste0(vals[2], " (", vals[1], ", ", vals[3], ")")
  return(neat)
}



#' Print a warning if any monitored parameters have rhat above the cutoff in an MBNMA model
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
print.treat.str <- function(mbnma) {

  if (!is.null(mbnma$model.arg$arg.params)) {
    wrapper <- TRUE
  } else {wrapper <- FALSE}

  betas <- assignfuns(fun=mbnma$model.arg$fun, agents=mbnma$network$agents, user.fun=mbnma$model.arg$user.fun,
                      wrapper=wrapper,
                      knots=ifelse(length(mbnma$model.arg$knots)==1, mbnma$model.arg$knots, length(mbnma$model.arg$knots)))

  datasum <- as.data.frame(cbind(mbnma$BUGSoutput$summary[,5],
                                 mbnma$BUGSoutput$summary[,3],
                                 mbnma$BUGSoutput$summary[,7]))

  treat.sect <- c()

  # For each dose-response parameter generate results
  for (i in seq_along(betas)) {
    if (!(betas[[i]]$betaname %in% names(mbnma$model.arg$class.effect))) {

      # Create section heading
      if (wrapper) {
        headbeta <- betas[[i]]$betaname
      } else {
        headbeta <- paste0(betas[[i]]$betaname, " (", betas[[i]]$fun, ", ", betas[[i]]$param, ")")
      }

      sect.head <- paste(crayon::bold(headbeta), "dose-response parameter results\n\n", sep=" ")
      cat(crayon::underline(sect.head))

      if (betas[[i]]$param %in% c("lambda", "ed50")) {
        cat("Parameter modelled on exponential scale to ensure it takes positive values\non the natural scale\n")
      }

      agents <- mbnma$network$agents[mbnma$network$agents!="Placebo"]

      if (is.character(mbnma$model.arg[[names(betas)[i]]])) {

        if (mbnma$model.arg[[names(betas)[i]]]=="rel") {
          param <- "d"
          cat("Pooling: relative effects\n\n")
        } else if (mbnma$model.arg[[names(betas)[i]]] %in% c("common", "random")) {
          param <- "beta"
          cat("Pooling: single parameter shared across the network\n\n")
        }

        # Create table of results
        datai <- vector()
        datai <- append(datai,
                        which(grepl(paste0("^", param, "\\.", betas[[i]]$param), rownames(mbnma$BUGSoutput$summary))))
        datai <- append(datai,
                        which(grepl(paste0("^", param, "\\.", i), rownames(mbnma$BUGSoutput$summary))))
        datai <- datai[datai!=0]

        datatab <- datasum
        names(datatab) <- c("Median", "2.5%", "97.5%")
        datatab$Parameter <- rownames(datatab)
        datatab <- datatab[datai,c(4,1,2,3)]


        if (param=="d") {
          # Drop rows that aren't relevant for multi-fun models
          datatab <- datatab[betas[[i]]$agents,]
          rownames(datatab) <- agents[betas[[i]]$agents]
        } else if (param=="beta") {
          rownames(datatab) <- ""
        }

        # Print table of results
        print(datatab)
        cat("\n\n")


        # Print between-study SD
        if (mbnma$model.arg[[names(betas)[i]]]=="random") {

          datai <- vector()
          datai <- append(datai, which(grepl(paste0("^sd\\.", i), rownames(mbnma$BUGSoutput$summary))))
          datai <- append(datai, which(grepl(paste0("^sd\\.", betas[[i]]$param), rownames(mbnma$BUGSoutput$summary))))
          datai <- datai[datai!=0]

          datatab <- as.data.frame(datasum)
          colnames(datatab) <- c("Median", "2.5%", "97.5%")
          datatab$Parameter <- rownames(datatab)
          datatab <- datatab[datai,c(4,1,2,3)]

          rownames(datatab) <- ""
          print(datatab)
          cat("\n\n")

        }

        # Print result if parameter assigned a numeric value
      } else if (is.numeric(mbnma$model.arg[[names(betas)[i]]])) {
        data.str <- paste("Assigned a numeric value:",
                          mbnma$model.arg[[names(betas)[i]]],
                          sep=" ")
        cat(data.str)
        cat("\n\n")

      }

    }
  }
  invisible(mbnma)
}


#' Neatly prints a summary of the method used in the model
#'
#' @noRd
print.method.sect <- function(mbnma) {
  # String for method
  data.head <- paste("Parameter", "Median (95%CrI)", sep="\t\t\t\t\t")
  data.head <- paste(crayon::bold(data.head, "-----------------------------------------------------------------------", sep="\n"))

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
  method.str <- paste(crayon::bold(crayon::underline("\n\nPooling method")), method.str, "\n", sep="\n\n")
  return(method.str)
}


#' Neatly prints class results
print.class.str <- function(mbnma) {

  if (length(mbnma$model.arg$class.effect)>0) {
    if (!is.null(mbnma$model.arg$arg.params)) {
      wrapper <- TRUE
    } else {wrapper <- FALSE}

    classes <- mbnma$model.arg$class.effect
    betas <- assignfuns(fun=mbnma$model.arg$fun, agents=mbnma$network$agents, user.fun=mbnma$model.arg$user.fun,
                        wrapper=wrapper)

    datasum <- as.data.frame(cbind(mbnma$BUGSoutput$summary[,5],
                                   mbnma$BUGSoutput$summary[,3],
                                   mbnma$BUGSoutput$summary[,7]))

    head <- crayon::bold(crayon::underline("\nClass effects\n"))

    # For each class
    for (i in seq_along(classes)) {
      if (wrapper) {
        for (k in seq_along(betas)) {
          if (names(classes)[i]==betas[[k]]$betaname) {

            sect.head <- paste("Class effect results for:", names(classes)[i], "\n\n", sep=" ")
            cat(sect.head)

            if (mbnma$model.arg[[names(betas)[k]]]=="rel") {
              param <- "D"
            } else if (mbnma$model.arg[[names(betas)[k]]] %in% c("common", "random")) {
              param <- "BETA"
            }

            # Generate table of class results
            datai <- vector()
            datai <- append(datai,
                            which(grepl(paste0("^", param, "\\.", betas[[i]]$param), rownames(mbnma$BUGSoutput$summary))))
            datai <- append(datai,
                            which(grepl(paste0("^", param, "\\.", i), rownames(mbnma$BUGSoutput$summary))))
            datai <- datai[datai!=0]

            datatab <- datasum
            names(datatab) <- c("Median", "2.5%", "97.5%")
            datatab$Parameter <- rownames(datatab)
            datatab <- datatab[datai,c(4,1,2,3)]

            if (param=="D") {
              if (nrow(datatab)==length(mbnma$network$classes)) {
                rownames(datatab) <- mbnma$network$classes
              } else {
                rownames(datatab) <- mbnma$network$classes[-1]
              }
            }

            # Print table of class results
            print(datatab)
            cat("\n\n")


            # Add between-agent SD
            if (classes[[i]]=="random") {

              cat(paste0("Within-class SD for ", names(classes)[i], "\n\n"))

              datai <- vector()
              datai <- append(datai, which(grepl(paste0("^sd\\.", param, "\\.", i), rownames(mbnma$BUGSoutput$summary))))
              datai <- append(datai, which(grepl(paste0("^sd\\.", param, "\\.", betas[[i]]$param), rownames(mbnma$BUGSoutput$summary))))
              datai <- datai[datai!=0]

              datatab <- as.matrix(datasum)
              colnames(datatab) <- c("Median", "2.5%", "97.5%")
              datatab <- datatab[datai,]
              print(datatab)
              cat("\n\n")

            }
          }
        }
      }
    }

    invisible(mbnma)
  }
}


#' Neatly prints model fit details
print.modfit.str <- function(mbnma) {
  totresdev.str <- c()

  head <- crayon::bold(crayon::underline("Model Fit Statistics\n"))

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
  pd.str <- paste(pd.str, paste(pd, round(mbnma$BUGSoutput$pD,1), sep=" "), sep="\n")

  # Deviance
  dev <- mbnma$BUGSoutput$summary[
    rownames(mbnma$BUGSoutput$summary)=="deviance", 5]
  dev.str <- paste("Deviance =", round(dev, 1), sep=" ")

  # Totresdev
  if ("totresdev" %in% mbnma$parameters.to.save) {
    totresdev <- round(
      mbnma$BUGSoutput$summary[
        rownames(mbnma$BUGSoutput$summary)=="totresdev", 5],
      1)
  } else {
    totresdev <- crayon::red("NOT MONITORED IN MODEL")
  }
  totresdev.str <- paste("Residual deviance =", totresdev, sep=" ")

  dic <- mbnma$BUGSoutput$DIC
  dic.str <- paste("Deviance Information Criterion (DIC) =", round(dic, 1), "\n", sep=" ")

  modfit.sect <- paste(head, pd.str, "", dev.str, totresdev.str, dic.str, sep="\n")
  return(modfit.sect)
}
