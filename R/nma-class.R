######################################
#### Functions for class("nma") ####
######################################

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))


#' @describeIn nma.run Plot outputs from treatment-level NMA models
#'
#' Results can be plotted either as a single forest plot, or facetted by agent
#' and plotted with increasing dose in order to identify potential dose-response
#' relationships.
#'
#' @param x An object of `class("nma")`
#' @param bydose A boolean object indicating whether to plot responses with dose
#' on the x-axis (`TRUE`) to be able to examine potential dose-response shapes, or
#' to plot a conventional forest plot with all treatments on the same plot (`FALSE`)
#' @param ... Arguments to be sent to [ggplot2::ggplot()]
#' @inheritParams plot.mbnma.predict
#'
#' @export
plot.nma <- function(x, bydose=TRUE, scales="free_x", ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "nma", add=argcheck)
  checkmate::assertLogical(bydose, len=1, add=argcheck)
  checkmate::assertChoice(scales, c("free_x", "fixed"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  intercept <- 0 # Leaving here in case want to allow user to change it at later date
  split.df <- x[["jagsresult"]]$BUGSoutput$summary

  # Check if NMA is from UME model
  if (any(grepl("^d\\[[0-9+],[0-9]+\\]", rownames(x$jagsresult$BUGSoutput$summary)))) {
    split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+,1\\]", rownames(split.df)), c(3,5,7)])
  } else {
    split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+\\]", rownames(split.df)), c(3,5,7)])
  }


  # Get doses, treatments and agent codes
  split.df$treatment <- x[["trt.labs"]]
  split.df$agent <- sapply(x[["trt.labs"]],
                           function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])
  split.df$dose <- as.numeric(sapply(x[["trt.labs"]],
                                     function(x) strsplit(x, split="_", fixed=TRUE)[[1]][2]))

  if (split.df$`50%`[1]!=0 & split.df$`2.5%`[1]!=0) {
    row <- split.df[0,]
    row[,1:3] <- 0
    row$treatment <- "Placebo_0"
    row$agent <- "Placebo"
    row$dose <- 1
    split.df <- rbind(row, split.df)
  }

  # Add intercept for all agents
  agents <- unique(split.df$agent)
  ref.agent <- "Placebo"
  ref.trt <- "Placebo_0"
  ylab.es <- "Effect size on link scale versus Placebo"

  if (!"Placebo" %in% agents) {
    ref.agent <- split.df$agent[1]
    ref.trt <- split.df$treatment[1]
    message(paste0("Placebo not included in dataset - reference treatment for plot is ", ref.trt))
    ylab.es <- paste0("Effect size on link scale versus ", ref.trt)
  }

  # Plot faceted by agent as dose-response splitplot
  if (bydose==TRUE) {

    # agents <- agents[agents!="Placebo"]
    # for (i in seq_along(agents)) {
    #   row <- split.df[split.df$agent=="Placebo",]
    #   row$agent <- agents[i]
    #   split.df <- rbind(row, split.df)
    # }
    # split.df <- split.df[split.df$agent!="Placebo",]

    agents <- agents[agents!="Placebo"]
    for (i in seq_along(agents)) {
      row <- split.df[split.df$treatment==ref.trt,]
      row$agent <- agents[i]
      split.df <- rbind(row, split.df)
    }

    if (ref.agent=="Placebo") {
      split.df <- split.df[split.df$agent!=ref.agent,]
    } else {
      split.df <- split.df[!(split.df$treatment==ref.trt & split.df$agent!=ref.agent),]
    }

    if (intercept!=0) {
      split.df[,1:3] <- split.df[,1:3] + intercept
    }

    # Generate forest plot
    g <- ggplot2::ggplot(split.df, ggplot2::aes(x=dose, y=`50%`), ...) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`, width=.05)) +
      ggplot2::facet_wrap(~factor(agent), scales = scales) +
      ggplot2::xlab("Dose") +
      ggplot2::ylab(ylab.es) +
      theme_mbnma()

  } else if (bydose==FALSE) {
    # Plot conventional forest plot
    split.df$treatment <- factor(split.df$treatment, levels=x[["trt.labs"]])

    g <- ggplot2::ggplot(split.df, ggplot2::aes(y=`50%`, x=treatment), ...) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`), width=.2) +
      ggplot2::coord_flip() +
      ggplot2::ylab(ylab.es) +
      ggplot2::xlab("Treatment") +
      theme_mbnma()
  }

  graphics::plot(g)
  return(invisible(g))
}
