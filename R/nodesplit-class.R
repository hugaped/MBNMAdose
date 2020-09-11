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

  width <- "\t\t"
  output <- crayon::bold("========================================\nNode-splitting analysis of inconsistency\n========================================\n")

  comparisons <- names(x)
  colnam <- "comparison\tp.value\t\t\tMedian (95% CrI)\n"
  paramsect <- colnam
  for (i in seq_along(comparisons)) {
    pval <- signif(x[[i]]$p.values,
                   max(3L, getOption("digits") - 3L))
    tab <- x[[i]]$quantiles

    heading <- paste(crayon::bold(names(x)[i]), pval, sep=width)
    direct <- paste("-> direct", "", neatCrI(tab$direct), sep=width)
    indirect <- paste("-> indirect", "", neatCrI(tab$indirect), sep=width)

    if ("mbnma" %in% names(x[[i]]$quantiles)) {
      nma <- paste("-> MBNMA", "", neatCrI(tab$mbnma), sep=width)
    } else if ("nma" %in% names(x[[i]]$quantiles)) {
      nma <- paste("-> NMA", "\t", neatCrI(tab$nma), sep=width)
    }

    out <- paste(heading, direct, indirect, nma, "", sep="\n")

    paramsect <- paste(paramsect, out, sep="\n")
  }
  output <- append(output, paramsect)
  cat(output, ...)
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




#' @describeIn mbnma.nodesplit Plot outputs from treatment-level nodesplit MBNMA models
#'
#' @param x An object of `class("nodesplit")`
#' @param plot.type A character string that can take the value of `"forest"` to plot
#' only forest plots, `"density"` to plot only density plots, or left as `NULL` (the
#' default) to plot both types of plot.
#' @param ... Arguments to be sent to [ggplot2::ggplot()]
#'
#' @details The S3 method `plot()` on an `nodesplit` object generates either
#' forest plots of posterior medians and 95\\% credible intervals, or density plots
#' of posterior densities for direct and indirect evidence.
#'
#' @return Plots the desired graph(s) and returns an object (or list of object if
#' `plot.type=NULL`) of `class(c("gg", "ggplot"))`
#'
#' @export
plot.nodesplit <- function(x, plot.type=NULL, ...) {
  # ... are commands to be sent to geom_histogram

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "nodesplit", add=argcheck)
  checkmate::assertChoice(plot.type, choices = c("density", "forest"), null.ok=TRUE, add=argcheck)
  #checkmate::assertLogical(facet, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (is.null(plot.type)) {
    plot.type <- c("density", "forest")
  }

  forestdata <- x[[1]]$forest.plot$data[0,]
  densitydata <- x[[1]]$density.plot$data[0,]
  forestfacet <- vector()
  densityfacet <- vector()
  for (i in seq_along(x)) {
    comp <- paste(x[[i]]$comparison, collapse=" vs ")
    temp <- x[[i]]$forest.plot$data
    forestfacet <- append(forestfacet, rep(comp, nrow(temp)))
    forestdata <- rbind(forestdata, temp)

    temp <- x[[i]]$density.plot$data
    densityfacet <- append(densityfacet, rep(comp, nrow(temp)))
    densitydata <- rbind(densitydata, temp)
  }
  forestdata$comp <- forestfacet
  densitydata$comp <- densityfacet


  if ("forest" %in% plot.type) {
    forest <-
      ggplot2::ggplot(data=forestdata, ggplot2::aes(x=source, y=med,
                                                    ymin=l95, ymax=u95), ...) +
      ggplot2::geom_pointrange() +
      ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
      ggplot2::xlab("") + ggplot2::ylab("Treatment effect (95% CrI)") +
      ggplot2::theme(axis.text = ggplot2::element_text(size=10),
                     axis.title = ggplot2::element_text(size=12),
                     title=ggplot2::element_text(size=18)) +
      ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,1),"cm")) +
      ggplot2::facet_wrap(~factor(comp)) +
      ggplot2::theme_bw()
  }
  if ("density" %in% plot.type) {

    density <- ggplot2::ggplot(densitydata, ggplot2::aes(x=value,
                                                         linetype=Estimate, fill=Estimate), ...) +
      ggplot2::geom_density(alpha=0.2) +
      ggplot2::xlab("Treatment effect (95% CrI)") +
      ggplot2::ylab("Posterior density") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                     axis.title = ggplot2::element_text(size=14)) +
      ggplot2::facet_wrap(~factor(comp)) +
      ggplot2::theme_bw() +
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

  }

  if (identical(sort(plot.type), c("density", "forest"))) {
    graphics::plot(forest)
    graphics::plot(density)
    return(invisible(list(forest, density)))
  } else {
    if (plot.type=="forest") {
      graphics::plot(forest)
      return(invisible(forest))
    } else if (plot.type=="density") {
      graphics::plot(density)
      return(invisible(density))
    }
  }

}
