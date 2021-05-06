##############################################
#### Functions for class("nodesplit") ####
##############################################

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))


#' Prints summary results from a nodesplit object
#'
#' @param x An object of `class("nodesplit")`
#' @param ... further arguments passed to or from other methods
#'
#' @export
print.nodesplit <- function(x, ...) {
  checkmate::assertClass(x, "nodesplit")

  cat(crayon::bold("========================================\nNode-splitting analysis of inconsistency\n========================================\n"))
  sum.df <- summary.nodesplit(x)
  comps <- unique(sum.df$Comparison)

  out.df <- sum.df[1,]
  for (i in seq_along(comps)) {

    head <- sum.df[sum.df$Comparison==comps[i],][1,]
    head$Evidence <- NA
    head$Median <- NA
    head$`2.5%` <- NA
    head$`97.5%` <- NA

    tail <- sum.df[sum.df$Comparison==comps[i],]
    tail$Comparison <- c("-> direct", "-> indirect", "-> MBNMA")
    tail$p.value <- rep(NA,3)

    tail <- rbind(tail, rep(NA, ncol(tail)))

    out.df <- rbind(out.df, rbind(head,tail))
  }
  out.df <- out.df[-1,c(1,6,3:5)]
  out <- knitr::kable(out.df, row.names = FALSE,
                      col.names = c("Comparison", "p-value", "Median", "2.5%", "97.5%"), ...)
  cat(gsub('\\bNA\\b', '  ', out), sep='\n')
}





#' Generates a summary data frame for nodesplit objects
#'
#' @param object An object of `class("nodesplit")`
#' @param ... further arguments passed to or from other methods
#'
#' @export
summary.nodesplit <- function(object, ...) {
  checkmate::assertClass(object, "nodesplit")

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

      if ("mbnma" %in% names(object[[i]]$quantiles)) {
        sum.mat <- rbind(sum.mat, post$mbnma)
        evidence <- c(evidence, "MBNMA")
      } else if ("nma" %in% names(object[[i]]$quantiles)) {
        sum.mat <- rbind(sum.mat, post$nma)
        evidence <- c(evidence, "NMA")
      }

      pvals <- c(pvals, rep(object[[i]]$p.values, 3))
      comp <- c(comp, rep(paste(object[[i]]$comparison[1],
                                object[[i]]$comparison[2], sep=" vs "), 3))
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




#' @describeIn mbnma.nodesplit Plot outputs from treatment-level nodesplit MBNMA models
#'
#' @param x An object of `class("nodesplit")`
#' @param plot.type A character string that can take the value of `"forest"` to plot
#' forest plots or `"density"` to plot posterior density plots.
#' @param ... Arguments to be sent to [ggplot2::ggplot()] or [forestplot::forestplot()]
#'
#' @details The S3 method `plot()` on an `nodesplit` object generates either
#' forest plots of posterior medians and 95\\% credible intervals, or density plots
#' of posterior densities for direct and indirect evidence.
#'
#' @return Plots the desired graph if `plot.type="forest"` and plots and returns an object
#' of `class(c("gg", "ggplot"))` if `plot.type="density"`.
#'
#' @export
plot.nodesplit <- function(x, plot.type="forest", ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "nodesplit", add=argcheck)
  checkmate::assertChoice(plot.type, choices = c("density", "forest"), null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (plot.type == "forest") {

    forest <- forest.splits(x, ...)

  } else if (plot.type == "density") {
    densitydata <- x[[1]]$density.plot$data[0,]
    densityfacet <- vector()
    for (i in seq_along(x)) {
      comp <- paste(x[[i]]$comparison, collapse=" vs ")

      temp <- x[[i]]$density.plot$data
      densityfacet <- append(densityfacet, rep(comp, nrow(temp)))
      densitydata <- rbind(densitydata, temp)
    }
    densitydata$comp <- densityfacet

    # Density plot
    density <- ggplot2::ggplot(densitydata, ggplot2::aes(x=value,
                                                         linetype=Estimate, fill=Estimate), ...) +
      ggplot2::geom_density(alpha=0.2) +
      ggplot2::xlab("Treatment effect (95% CrI)") +
      ggplot2::ylab("Posterior density") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                     axis.title = ggplot2::element_text(size=14)) +
      ggplot2::facet_wrap(~factor(comp)) +
      theme_mbnma() +
      ggplot2::labs(linetype="Evidence", fill="Evidence")

    # Check for size of y-axis and allow scale=free if large difference (e.g. 50x)
    g <- ggplot2::ggplot_build(density)
    dens <- vector()
    for (i in seq_along(unique(g$data[[1]]$PANEL))) {
      dens <- append(dens, max(g$data[[1]]$density[g$data[[1]]$PANEL==unique(g$data[[1]]$PANEL)[i]]))
    }
    if (max(dens) > (min(dens) * 50)) {
      density <- density + ggplot2::facet_wrap(~factor(densitydata$comp), scales="free_y")
    }

    graphics::plot(density)
    return(invisible(density))
  }

}
