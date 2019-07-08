# Functions for plotting in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-18

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' @describeIn mbnma.network Generate a network plot
#' @inheritParams mbnma.run
#'
#' @param x An object of class `mbnma.network`.
#' @param layout_in_circle A boolean value indicating whether the network plot
#'   should be shown in a circle or left as the igraph default layout.
#' @param edge.scale A number to scale the thickness of connecting lines
#'   (edges). Line thickness is proportional to the number of studies for a
#'   given comparison. Set to 0 to make thickness equal for all comparisons.
#' @param v.color Can take either `"connect"` (the default) to indicate that nodes should
#'   only be coloured if they are connected to the network reference treatment (indicates
#'   network connectivity) or `"agent"` to colour nodes by agent.
#' @param v.scale A number with which to scale the size of the nodes. If the variable `N`
#'   (to indicate the numbers of participants in each study arm) is included in the
#'   dataset then the size of the nodes will be proportional to the number of participants
#'   within a treatment/agent in the network.
#' @param label.distance A number scaling the distance of labels from the nodes
#'   to improve readability. The labels will be directly on top of the nodes if
#'   the defauls of 0 is used. Option only applicable if `layout_in_circle` is
#'   set to `TRUE`.
#' @param level A string indicating whether nodes/facets should represent `"treatment"`
#'   or `"agent"` in the plot. Can be used to examine the expected impact of modelling
#'   dose-response in terms of network connectivity.
#' @param remove.loops A boolean value indicating whether to include loops that
#'   indicate comparisons within a node.
#' @param doselink If given an integer value it indicates that connections via the dose-response
#' relationship with placebo should be plotted. The integer represents the minimum number of doses
#' from which a dose-response function could be estimated and is equivalent to the number of
#' parameters in the desired dose-response function plus one. If left as `NULL` (the default), connections
#' to placebo via dose-response relationships will not be included.
#' @param ... Options for plotting in `igraph`.
#'
#' @details The S3 method `plot()` on an `mbnma.network` object generates a
#'   network plot that shows how different treatments are connected within the
#'   network via study comparisons. This can be used to identify how direct and
#'   indirect evidence are informing different treatment comparisons. Depends on
#'   \code{\link[igraph]{igraph}}.
#'
#' @return `plot()`: An object of `class("igraph")` - any functions from the `igraph` package
#' can be applied to this object to change its characteristics.
#'
#' @examples
#' # Create an mbnma.network object from the data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Generate a network plot from the data
#' plot(network)
#'
#' # Generate a network plot at the agent level that removes loops indicating comparisons
#' #within a node
#' plot(network, level="agent", remove.loops=TRUE)
#'
#' # Generate a network plot at the treatment level that colours nodes by agent
#' plot(network, v.color="agent", remove.loops=TRUE)
#'
#' # Generate a network plot that includes connections via the dose-response function
#' # For a one parameter dose-response function (e.g. exponential)
#' plot(network, level="treatment", doselink=1, remove.loops=TRUE)
#'
#' # For a two parameter dose-response function (e.g. Emax)
#' plot(network, level="treatment", doselink=2, remove.loops=TRUE)
#'
#'
#' #### Plot a network with no placebo data included ####
#' # Make data with no placebo
#' noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
#' net.noplac <- mbnma.network(noplac.df)
#'
#' # Plotting network automatically plots connections to Placebo via dose-response
#' plot(net.noplac)
#' @export
plot.mbnma.network <- function(x, level="treatment", v.color="connect", doselink=NULL,
                               layout_in_circle = TRUE, remove.loops=FALSE,
                               edge.scale=1, v.scale=NULL, label.distance=0,
                               ...)
  # Requires igraph
  #S3method(plot, mbnma.network)

  # x is an object of class mbnma.network

{
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma.network", add=argcheck)
  checkmate::assertLogical(layout_in_circle, len=1, add=argcheck)
  checkmate::assertNumeric(edge.scale, finite=TRUE, len=1, add=argcheck)
  checkmate::assertNumeric(label.distance, finite=TRUE, len=1, add=argcheck)
  checkmate::assertNumeric(v.scale, lower = 0, finite=TRUE, null.ok=TRUE, len=1, add=argcheck)
  checkmate::assertChoice(level, choices = c("treatment", "agent"), add=argcheck)
  checkmate::assertChoice(v.color, choices = c("connect", "agent"), add=argcheck)
  checkmate::assertLogical(remove.loops, len=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Generate comparisons (using get.latest.time and mbnma.contrast?
  data.ab <- x$data.ab

  # Add "s" onto level to make consistent with network names
  levels <- paste0(level, "s")

  # Check if level="agent" that agents are present in dataset
  if (!(levels %in% names(x))) {
    stop(paste0("`level` has been set to `",
                level,
                "` but ", levels, " is not a variable within the dataset"))
  }

  #nodes <- x[[levels]]
  #data.ab$node <- as.character(factor(data.ab[[level]], labels=x[[levels]]))

  # if (!(nodes[1] %in% c("Placebo", "Placebo_0"))) {
  #   plac.incl <- FALSE
  #   net.lbls <- c("Placebo", x[[levels]])
  #   data.ab <- add.plac.row(data.ab)
  #
  # } else {
  #   plac.incl <- TRUE
  #   net.lbls <- x[[levels]]
  # }
  net.lbls <- x[[levels]]
  nodes <- net.lbls
  data.ab$node <- as.character(factor(data.ab[[level]], labels=net.lbls))

  # Calculate participant numbers (if v.scale not NULL)
  if (!is.null(v.scale)) {
    if (!("N" %in% names(data.ab))) {
      warning("`N` not included as a column in dataset. Vertices/nodes will all be scaled to be the same size.")
    }

    size.vec <- vector()
    for (i in seq_along(nodes)) {
      size.vec <- append(size.vec, sum(data.ab$N[data.ab$node==nodes[i]]))
    }
    # Scale size.vec by the max node.size
    size.vec <- size.vec/ (max(size.vec)/20)

    node.size <-
      stats::setNames(size.vec, nodes)
    node.size <- as.matrix(node.size*v.scale)
  } else {
    node.size <- NULL
  }

  # Change treatment column for agent if necessary
  if (level=="agent") {
    data.ab$treatment <- data.ab$agent
  }

  comparisons <- mbnma.comparisons(data.ab)

  # Add coloured vertices for plac if plac.incl!=TRUE
  if ((x$agents[1] != "Placebo" & x$treatments[1]!="Placebo_0")) {
    plac.incl <- FALSE
    # if (is.null(doselink)) {
    #   doselink <- 2
    # }
  } else {plac.incl <- TRUE}

  if (!is.null(doselink)) {
    dr.comp <- DR.comparisons(x$data.ab, level=level, doselink=doselink)
    if (plac.incl==TRUE) {
      dr.comp$t1 <- dr.comp$t1 + 1
    } else if (plac.incl==FALSE & nrow(dr.comp)>0) {
      nodes <- c("Placebo", nodes)
      if (!is.null(node.size)) {
        node.size <- c(1, node.size)
      }
    }
    #comparisons <- rbind(comparisons, dr.comp)
    comparisons <- rbind(dr.comp, comparisons)
  }


  # Code to make graph.create as an MBNMA command if needed
  g <- igraph::graph.empty()
  g <- g + igraph::vertex(nodes)
  ed <- t(matrix(c(comparisons[["t1"]], comparisons[["t2"]]), ncol = 2))
  ed <- factor(as.vector(ed), labels=nodes)
  edges <- igraph::edges(ed, weight = comparisons[["nr"]], arrow.mode=0)
  #edges <- igraph::edges(as.vector(ed), weight = comparisons[["nr"]], arrow.mode=0)
  g <- g + edges


  igraph::E(g)$curved <- FALSE # ensure edges are straight

  if (!is.null(doselink)) {
    igraph::E(g)$color <- c(rep("red", nrow(dr.comp)),
                            rep("black", nrow(comparisons)-nrow(dr.comp)))
    # igraph::E(g)$lty <- c(rep("dashed", nrow(dr.comp)),
    #                       rep("solid", nrow(comparisons)-nrow(dr.comp)))
  }


  if (remove.loops==TRUE) {
    g <- igraph::simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
  }

  # Check network is connected and produce warning message if not
  disconnects <- check.network(g)
  if (v.color=="connect") {
    igraph::V(g)$color <- "SkyBlue2"
    igraph::V(g)$color[which(names(igraph::V(g)) %in% disconnects)] <- "white"
  } else if (v.color=="agent") {
    if (!("agents" %in% names(x))) {
      stop("`level` has been set to `agent` but there are no agent codes given in the dataset")
    }

    # Get large vector of distinct colours using Rcolorbrewer
    cols <- genmaxcols()

    if (level=="treatment") {
      #igraph::V(g)$color <- cols[1:length(sapply(nodes, function(x) strsplit(x, "[_]")[[1]][1]))]
      temp <- as.character(sapply(nodes, function(x) strsplit(x, "[_]")[[1]][1]))
      igraph::V(g)$color <- cols[as.numeric(factor(temp))]
    } else if (level=="agent") {
      igraph::V(g)$color <- cols[1:length(nodes)]
    }
  }

  # Plot netgraph
  if (layout_in_circle==TRUE) {
    lab.locs <- radian.rescale(x=seq(1:length(nodes)), direction=-1, start=0)
    igraph::plot.igraph(g,
                        edge.width = edge.scale * comparisons[["nr"]],
                        layout = igraph::layout_in_circle(g),
                        vertex.label.dist=label.distance,
                        vertex.label.degree=lab.locs,
                        vertex.size = node.size,
                        ...
    )
  } else {
    igraph::plot.igraph(g, edge.width = edge.scale * comparisons[["nr"]], vertex.size = node.size,
                        ...)
  }

  if (!is.null(doselink)) {
    message(paste0("Dose-response connections to placebo plotted based on a dose-response
                   function with ", (doselink-1),
                   " degrees of freedom"))
  }

  return(invisible(g))
}






#' Check if all nodes in the network are connected (identical to function in `MBNMAtime`)
#'
#' @param g An network plot of `class("igraph")`
#' @param reference A numeric value indicating which treatment code to use as the reference treatment for
#' testing that all other treatments connect to it
#'
#' @noRd
check.network <- function(g, reference=1) {

  # Can add component to test for if placebo is missing:
  # First use check.network with reference=2 (so that next treatment is reference)

  # Then any loops which aren't connected, check if there are X doses of ANY agent within
  #that loop AND within the main loop. Report the max number of doses that are common to
  #all loops in the network (perhaps even colour vertices of loops different depending on
  #the number of doses they can connect via)

  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               to=reference))
  treats <- rownames(connects)[connects==FALSE]

  if (length(treats>0)) {
    warning(paste0("The following treatments/agents are not connected to the network reference:\n",
                   paste(treats, collapse = "\n")))
  }
  return(treats)
}




#' Calculate position of label with respect to vertex location within a circle
#'
#' Useful for graphs drawn using `igraph` to reposition labels relative to vertices when vertices
#' are laid out in a circle (as is common in network plots). `igraph` interprets position within
#' `vertex.label.degree` as radians, so it is necessary to convert locations into radian values. This
#' is the main role of this function.
#'
#' @param x A numeric vector of positions around a circle, typically sequentially numbered.
#' @param start A number giving the offset from 12 o'clock in radians for the label locations.
#' @param direction Either `1` for clockwise numbering (based on the order of `x`) or `-1` for
#' anti-clockwise.
#'
#' @examples
#' radian.rescale(c(1:10), start=0, direction=1)
#'
#' @noRd
#'
#' @references
#' https://gist.github.com/kjhealy/834774/a4e677401fd6e4c319135dabeaf9894393f9392c
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}




#' Get large vector of distinct colours using Rcolorbrewer
#' @noRd
genmaxcols <- function() {

  cols1 <- RColorBrewer::brewer.pal(9, "Set1")
  cols2 <- RColorBrewer::brewer.pal(9, "Pastel1")
  cols <- c(rbind(cols1, cols2))

  cols1 <- RColorBrewer::brewer.pal(8, "Set2")
  cols2 <- RColorBrewer::brewer.pal(8, "Pastel2")
  cols <- c(cols, c(rbind(cols1, cols2)))

  cols <- c(cols, RColorBrewer::brewer.pal(12, "Set3"))

  return(cols)
}





#' Forest plot for results from dose-response MBNMA models
#'
#' Generates a forest plot for dose-response parameters.
#'
#' @param x An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param params A character vector of dose-response parameters to plot.
#' Parameters must be given the same name as monitored nodes in `mbnma` and must be
#' modelled as relative effects (`"rel"`). Can be set to
#' `NULL` to include all available dose-response parameters estimated by `mbnma`.
#' @param agent.labs A character vector of agent labels (including `"Placebo"` if it
#' has been included in the original network). If left as `NULL` (the default) then
#' labels will be used as defined in the data.
#' @param class.labs A character vector of class labels if `mbnma` was modelled using class effects
#' (including `"Placebo"` if it
#' has been included in the original network). If left as `NULL`
#' (the default) then labels will be used as defined in the data.
#' @param ... Arguments to be passed to methods, such as graphical parameters
#'
#' @return A forest plot of class `c("gg", "ggplot")` that has separate panels for
#' different dose-response parameters. Results are plotted on the link scale.
#'
#' @examples
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an exponential dose-response MBNMA and generate the forest plot
#' exponential <- mbnma.run(network, fun="exponential")
#' plot(exponential)
#'
#' # Plot only Emax parameters from an Emax dose-response MBNMA
#' emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#' plot(emax, params=c("d.emax"))
#'
#'
#' #### Forest plots including class effects ####
#' # Generate some classes for the data
#' class.df <- HF2PPITT
#' class.df$class <- ifelse(df$agent=="placebo", "placebo", "active1")
#' class.df$class <- ifelse(df$agent=="eletriptan", "active2", df$class)
#' netclass <- mbnma.network(class.df)
#' emax <- mbnma.emax(netclass, emax="rel", ed50="rel", method="random",
#'             class.effect=list("ed50"="common"))
#'
#' # Plot forest plot with different labels for classes
#' plot(emax, class.labs=c("Placebo", "Other Active", "Eletriptan"))
#'
#' # Since "Placebo" is included in the network, it must be included in labels
#' # Failure to do so will cause an error
#' ## ERROR ## plot(emax, class.labs=c("Other Active", "Eletriptan"))
#'
#' @export
plot.mbnma <- function(x, params=NULL, agent.labs=NULL, class.labs=NULL, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma", add=argcheck)
  checkmate::assertChoice(params, choices=x[["parameters.to.save"]], null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check that specified params are modelled using relative effects
  for (i in seq_along(params)) {
    if (!(grepl("d\\.", params[i]) | grepl("D\\.", params[i]))) {
      stop(paste0(params[i], " has not been modelled using relative effects and does not vary by agent or class"))
    }
  }

  # Add all available params if is.null(params)
  if (is.null(params)) {
    params <- vector()

    # Add d
    params <- append(params,
                     x[["parameters.to.save"]][grep("^d\\.", x[["parameters.to.save"]])]
    )

    # Add D
    params <- append(params,
                     x[["parameters.to.save"]][grep("^D\\.", x[["parameters.to.save"]])]
    )

    if (length(params)==0) {
      stop("No dose-response relative effects can be identified from the model")
    }
  }


  # Compile parameter data into one data frame
  mbnma.sum <- as.data.frame(x[["BUGSoutput"]][["summary"]])
  plotdata <- mbnma.sum[0,]
  for (i in seq_along(params)) {
    paramdata <- mbnma.sum[grepl(paste0("^", params[i]),rownames(mbnma.sum)),]
    paramdata[["doseparam"]] <- rep(params[i], nrow(paramdata))
    plotdata <- rbind(plotdata, paramdata)
  }
  plotdata[["param"]] <- as.numeric(gsub("(.+\\[)([0-9]+)(\\])", "\\2", rownames(plotdata)))
  if (any(is.na(plotdata[["param"]]))) {
    plotdata[["param"]] <- c(1:nrow(plotdata))
  }

  # Change param labels for agents
  agentdat <- plotdata[grepl("^d\\.", rownames(plotdata)),]
  if (!is.null(agent.labs)) {
    agentcodes <- as.numeric(gsub("(^.+\\[)([0-9]+)(\\])", "\\2", rownames(agentdat)))
    if (length(agent.labs)!=max(agentcodes)) {
      stop("`agent.labs` length does not equal number of agents within the model")
    } else {
      a.labs <- agent.labs[sort(unique(agentcodes))]
    }
  } else if ("agents" %in% names(x)) {
    if (x$model.arg$fun %in% c("nonparam.up", "nonparam.down")) {
      a.labs <- x[["treatments"]]
    } else {
      a.labs <- x[["agents"]][x[["agents"]]!="Placebo"]
    }
  } else {
    a.labs <- sort(unique(agentdat$param))
  }

  # Change param labels for classes
  classdat <- plotdata[grepl("^D\\.", rownames(plotdata)),]
  c.labs <- vector()
  if (nrow(classdat)!=0) {
    if (!is.null(class.labs)) {
      classcodes <- as.numeric(gsub("(^.+\\[)([0-9]+)(\\])", "\\2", rownames(classdat)))
      c.labs <- class.labs[classcodes]
    } else if ("classes" %in% names(x)) {
      c.labs <- x[["classes"]][x[["classes"]]!="Placebo"]
    } else {
      c.labs <- sort(unique(classdat$param))
    }
  }

  # Increase param number for classes
  nagent <- ifelse(nrow(agentdat)>0, max(agentdat$param), 0)
  plotdata$param[grepl("^D\\.", rownames(plotdata))] <-
    plotdata$param[grepl("^D\\.", rownames(plotdata))] + nagent

  # Attach labels
  if (nrow(agentdat)>0) {
    all.labs <- c(a.labs, c.labs)
  } else {all.labs <- c.labs}
  plotdata$param <- factor(plotdata$param, labels=all.labs)

  if (any(is.na(levels(plotdata$param)))) {
    stop("`agent.labs` or `class.labs` have not been specified correctly. Perhaps include `Placebo` in labels")
  }

  g <- ggplot2::ggplot(plotdata, ggplot2::aes(y=plotdata$`50%`, x=plotdata$param)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=plotdata$`2.5%`, ymax=plotdata$`97.5%`)) +
    ggplot2::coord_flip()

  g <- g + ggplot2::facet_wrap(~plotdata$doseparam, scales="free")

  # Axis labels
  g <- g + ggplot2::xlab("Agent / Class") +
    ggplot2::ylab("Effect size")

  return(g)
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
#'   If `overlay.split = TRUE`, or `disp.obs = TRUE` then the original dataset must be specified
#'   by including the original `mbnma.network` object used to estimate the model as the `network`
#'   argument.
#'
#' @examples
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an Emax dose-response MBNMA and predict responses
#' emax <- mbnma.emax(network, method="random")
#' pred <- predict(emax, E0 = 0.5)
#' plot(pred)
#'
#' # Display observed doses on the plot (must include `network`)
#' plot(pred, disp.obs=TRUE, network=network)
#'
#' # Display split NMA results on the plot (must include `network`)
#' plot(pred, overlay.split=TRUE, network=network)
#'
#' # Split NMA results estimated using random treatment effects model
#' plot(pred, overlay.split=TRUE, network=network, method="random")
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
#'
#'
#' @export
plot.mbnma.predict <- function(x, network, disp.obs=FALSE,
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
    checkmate::assertClass(network, "mbnma.network", null.ok=TRUE)

    # Check that predict labels and agent labels in network are consistent
    if (!all(sum.pred$agent %in% network$agents)) {
      stop("Agent labels in `network` differ from those in `pred`")
    }

    g <- disp.obs(g=g, network=network, predict=x,
                  col="green", max.col.scale=NULL)

  }
  if (overlay.split==TRUE) {
    checkmate::assertClass(network, "mbnma.network", null.ok=TRUE)

    # Check that placebo is included (or dose=0 in networks without placebo)
    if (network$agents[1]!="Placebo") {
      stop("Placebo required in `network` for calculation of relative effects in split NMA")
    }

    # Check that at least one prediction is at a dose=0
    if (!(0 %in% sum.pred$dose)) {
      stop("`x` must include a predicted response at dose = 0 for at least one agent")
    }

    g <- overlay.split(g=g, network=network, method=method,
                       likelihood = x[["likelihood"]],
                       link = x[["link"]])

  }


  # Add overlayed lines and legends
  g <- g +
    ggplot2::geom_line(ggplot2::aes(y=`2.5%`, linetype="95% CrI")) +
    ggplot2::geom_line(ggplot2::aes(y=`97.5%`, linetype="95% CrI")) +
    ggplot2::geom_line(ggplot2::aes(linetype="Posterior Median"))

  g <- g + ggplot2::facet_wrap(~agent, scales=scales) +
    ggplot2::labs(y="Predicted response", x="Dose")

  g <- g + ggplot2::scale_linetype_manual(name="",
                                          values=c("Posterior Median"="solid",
                                                   "95% CrI"="dashed"))

  return(g)
}





#' Overlays observations as shaded regions on a time-course
#'
#' @inheritParams mbnma.run
#' @inheritParams plot.mbnma.predict
#' @param g An object of `class("ggplot")`
#' @param max.col.scale The maximum rgb numeric value to use for the colour scale
#' @noRd
disp.obs <- function(g, network, predict, col="red", max.col.scale=NULL) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(g, c("gg", "ggplot"), add=argcheck)
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertClass(predict, "mbnma.predict", add=argcheck)
  checkmate::assertInt(max.col.scale, lower=1, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  raw.data <- network[["data.ab"]]
  raw.data$agent <- factor(raw.data$agent, labels=network[["agents"]])
  predict.data <- summary(predict)
  predict.data[["count"]] <- NA
  predict.data[["cum.count"]] <- NA
  #predict.data <- predict.data[0,]

  # Change predict.data factors to numeric - REMOVE LATER
  #predict.data$agent <- as.numeric(as.character(predict.data$agent))
  #predict.data$dose <- as.numeric(as.character(predict.data$dose))


  # Identify counts of doses in raw data
  for (i in 1:nrow(predict.data)) {
    if (i==1) {
      dose.low <- -0.01
    } else if (predict.data$agent[i-1]!=predict.data$agent[i]) {
      dose.low <- -0.01
    } else {
      dose.low <- predict.data$dose[i-1]
    }

    dose.high <- predict.data$dose[i]
    p.agent <- predict.data$agent[i]

    count <- length(raw.data$studyID[as.character(raw.data$agent)==p.agent &
                                     raw.data$dose>dose.low & raw.data$dose<=dose.high])
    cum.count <- length(raw.data$studyID[as.character(raw.data$agent)==p.agent &
                                       raw.data$dose<=dose.high])

    predict.data$count[i] <- count
    predict.data$cum.count[i] <- cum.count
  }

  # Identify if placebo is in data
  if (any(unique(predict.data$agent) %in% c("1", "Placebo"))) {
    plac.count <- predict.data$count[predict.data$agent %in% c("1", "Placebo")]
    message(paste0(plac.count, " placebo arms in the dataset are not shown within the plots"))
    agents <- unique(predict.data$agent)[!(unique(predict.data$agent) %in% c("1", "Placebo"))]
  } else {
    agents <- unique(predict.data$agent)
  }


  # Check max.col.scale
  n.cut <- max(predict.data$count[predict.data$dose!=0])
  if (!is.null(max.col.scale)) {
    if (!is.numeric(max.col.scale)) {
      stop("max.col.scale must be a number greater than the maximum number of observed counts in the plotted treatments")
    }
    if (n.cut > max.col.scale) {
      stop("max.col.scale must be a number greater than the maximum number of observed counts in the plotted treatments")
    }
    n.cut <- max.col.scale
  }

  # Generate colours
  #cols <- col.scale(n.cut=max(predict.data$count), col=col)
  cols <- alpha.scale(n.cut=n.cut, col=col)

  for (agent in seq_along(agents)) {
    subset <- predict.data[predict.data$agent==agents[agent],]
    subset$agent <- factor(subset$agent, labels=levels(g$data$agent)[agents[agent]])

    # Start assuming lowest time = 0
    for (m in 2:nrow(subset)) {
      g <- g + ggplot2::geom_ribbon(data=subset[subset$dose<=subset$dose[m] &
                                                  subset$dose>=subset$dose[m-1]
                                                ,],
                                    ggplot2::aes(x=dose,
                                                 ymin=`2.5%`,
                                                 ymax=`97.5%`),
                                    fill=cols[subset$count[m]+1])
    }

  }

  return(g)
}




#' Generates colours with varying degrees of transparency
#'
#' Identical to function in `MBNMAtime` package
#'
#' @param ncut A number indicating the number of different counts in the dataset
#' @param col Colour to use for shading
#' @noRd
alpha.scale <- function(n.cut, col="blue") {
  # Run checks
  checkmate::assertIntegerish(n.cut, lower=1, len=1)

  # Set colour intensities
  if (is.character(col)) {
    if (col=="blue") {
      hue <- c(0,0,200)
    } else if (col=="red") {
      hue <- c(120,0,0)
    } else if (col=="green") {
      hue <- c(0,100,0)
    } else {
      stop("Permitted colours are `blue`, `red`, `green`, or an RGB code")
    }
  } else if (is.numeric(col)) {
    if (length(col)!=3) {
      stop("Specified RGB code must have length 3")
    }
    if (any(col>255) | any(col<0)) {
      stop("RGB values must be between 0 and 255")
    }
    hue <- col
  }

  #min.rgb <- rgb.end
  cut <- (255-0)/(n.cut)
  alpha <- 0
  alpha.vec <- alpha

  for (i in 1:n.cut) {
    alpha <- alpha+cut
    alpha.vec <- append(alpha.vec, alpha)
  }

  hexcol <- vector()
  for (i in 1:(n.cut+1)) {
    hexcol <- append(hexcol, grDevices::rgb(hue[1], hue[2],
                                 hue[3], alpha.vec[i], maxColorValue=255
    ))
  }

  return(hexcol)
}




#' Overlays results of split (treatment-level) NMA on a plot
#'
#' Requires automatically running an NMA
#'
#' @inheritParams plot.mbnma.predict
#' @inheritParams disp.obs
#' @noRd
overlay.split <- function(g, network, method="common",
                          likelihood="binomial", link="logit", ...) {

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  splitNMA <- nma.run(network=network, method=method,
                      likelihood=likelihood, link=link, ...)

  split.df <- splitNMA[["jagsresult"]]$BUGSoutput$summary
  split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+\\]", rownames(split.df)), c(3,5,7)])

  split.df$treatment <- splitNMA[["trt.labs"]]
  split.df$agent <- sapply(splitNMA[["trt.labs"]],
                       function(x) strsplit(x, split="_", fixed=TRUE)[[1]][1])
  split.df$dose <- as.numeric(sapply(splitNMA[["trt.labs"]],
                           function(x) strsplit(x, split="_", fixed=TRUE)[[1]][2]))

  split.df$agent <- factor(split.df$agent, levels=network$agents)
  #split.df$agent <- as.numeric(factor(split.df$agent, levels=network$agents))
  #split.df$agent <- factor(split.df$agent, labels=levels(g$data$agent))

  # Drop d[1] (reference)
  #split.df <- split.df[-1,]

  # Only include agents whose predictions are plotted
  split.df <- split.df[split.df$agent %in% g$data$agent,]

  # Add E0 (on correct scale)
  E0 <- unique(g$data$`50%`[g$data$dose==0])
  if (length(E0)>1) {
    warning("E0 in predictions appears to differ for different agents")
  }

  # Convert to scale using link function
  E0 <- rescale.link(E0, direction="link", link=link)

  split.df[,c(1:3)] <- split.df[,c(1:3)] + E0

  # Convert split.df to natural scale
  for (i in 1:3) {
    split.df[,i] <- rescale.link(split.df[,i], direction="natural", link=link)
  }



  # Add split NMAs to plot
  g <- g + ggplot2::geom_point(data=split.df, ggplot2::aes(x=split.df$dose, y=split.df$`50%`)) +
    ggplot2::geom_errorbar(data=split.df, ggplot2::aes(x=split.df$dose, ymin=split.df$`2.5%`, ymax=split.df$`97.5%`))


  #### Report split NMA model fit statistics ####
  nma.msg <- vector()
  nma.msg <- c(nma.msg, paste0("Split NMA residual Deviance: ",
                               round(splitNMA$jagsresult$BUGSoutput$median$totresdev, 1)))

  if (method=="random") {
    sd <- splitNMA$jagsresult$BUGSoutput$summary[rownames(splitNMA$jagsresult$BUGSoutput$summary)=="sd",]
    sd <- paste0(signif(sd[5], 3), " (", signif(sd[3], 3), ", ", signif(sd[7], 3), ")")

    nma.msg <- c(nma.msg, paste0("Split NMA between-study SD: ", sd
                                 ))
  }

  message(cat(paste(nma.msg, collapse="\n")))

  return(g)

}







#' Plot deviance contributions from an MBNMA model
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param dev.type *STILL IN DEVELOPMENT FOR MBNMAdose!* Deviances to plot - can be either residual
#' deviances (`"resdev"`, the default) or deviances (`"dev"`)
#' @param plot.type Deviances can be plotted either as scatter points (`"scatter"` - the default)
#' or as boxplots (`"box"`)
#' @param facet A boolean object that indicates whether or not to facet (by agent for `MBNMAdose`
#' and by treatment for `MBNMAtime`)
#' @param ... Arguments to be sent to `ggplot2::geom_point()` or `ggplot2::geom_boxplot`
#' @inheritParams predict.mbnma
#'
#' @return Generates a plot of deviance contributions and returns a list containing the
#' plot (as an object of `class(c("gg", "ggplot"))`), and a data.frame of posterior mean
#' deviance/residual deviance contributions for each observation.
#'
#' @details
#' Deviances should only be plotted for models that have converged successfully. If deviance
#' contributions have not been monitored in `mbnma$parameters.to.save` then additional
#' iterations will have to be run to get results for these.
#'
#' For `MBNMAtime`, deviance contributions cannot be calculated for models with a multivariate likelihood (i.e.
#' those that account for correlation between observations) because the covariance matrix in these
#' models is treated as unknown (if `rho = "estimate"`) and deviance contributions will be correlated.
#'
#' @examples
#' ###########################
#' ###### For MBNMAdose ######
#' ###########################
#'
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an Emax dose-response MBNMA and predict responses
#' emax <- mbnma.emax(network, method="random")
#'
#' # Plot deviances
#' devplot(emax)
#'
#' # Plot deviances using boxplots
#' devplot(emax, plot.type="box")
#'
#' # Plot deviances on a single plot (not facetted by agent)
#' devplot(emax, facet=FALSE, plot.type="box")
#'
#' # A data frame of deviance contributions can be obtained from the object
#' #returned by `devplot`
#' devs <- devplot(emax)
#' head(devs$dev.data)
#'
#' # Other deviance contributions not currently implemented but in future
#' #it will be possible to plot them like so
#' #devplot(emax, dev.type="dev")
#'
#'
#' ###########################
#' ###### For MBNMAtime ######
#' ###########################
#'
#' @export
devplot <- function(mbnma, plot.type="scatter", facet=TRUE, dev.type="resdev",
                    ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(dev.type, choices = c("dev", "resdev"), add=argcheck)
  checkmate::assertChoice(plot.type, choices = c("scatter", "box"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!is.null(mbnma$model.arg$rho)) {
    stop("Correlation between time points has been modelled using a multivariate likelihood.
         Deviances cannot be calculated for this model.")
  }

  if (!(dev.type %in% mbnma$parameters.to.save)) {
    msg <- paste0("`", dev.type, "` not monitored in mbnma$parameters.to.save.\n",
                  "additional iterations will be run in order to obtain results for `", dev.type, "`")
    message(msg)

    dev.df <- mbnma.update(mbnma, param=dev.type)

  } else {

    dev.df <- get.theta.dev(mbnma, param=dev.type)
  }

  # Plots the residual deviances
  if (mbnma$type=="time") {
    xlab <- "Follow-up count"
    facetscale <- "fixed"
  } else if (mbnma$type=="dose") {
    agents <- mbnma$agents

    # Remove placebo results if they are present
    if (mbnma$agents[1]=="Placebo") {
      dev.df <- dev.df[dev.df$facet!=1,]
      agents <- agents[-1]
    }

    xlab <- "Dose"
    facetscale <- "free_x"
  }

  if ("agents" %in% names(mbnma)) {
    dev.df$facet <- factor(dev.df$facet, labels=agents)
  } else if ("treatments" %in% names(mbnma)) {
    dev.df$facet <- factor(dev.df$facet, labels=mbnma$treatments)
  }

  if (plot.type=="scatter") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=dev.df$fupdose, y=mean), group=dev.df$groupvar) +
      ggplot2::geom_point(...)
  } else if (plot.type=="box") {
    g <- ggplot2::ggplot(dev.df, ggplot2::aes(x=factor(dev.df$fupdose), y=mean)) +
      ggplot2::geom_boxplot(...)
  }

  # Add axis labels
  g <- g + ggplot2::xlab(xlab) +
    ggplot2::ylab("Posterior mean")

  if (facet==TRUE) {
    g <- g + ggplot2::facet_wrap(~facet, scales = facetscale) +
      ggplot2::expand_limits(x=0)
  }

  graphics::plot(g)
  return(invisible(list("graph"=g, "dev.data"=dev.df)))
}





#' Extracts fitted values or deviance contributions into a data.frame with indices
#'
#' @inheritParams predict.mbnma
#' @param The parameter for which to extract values - can take either `"theta"`, `"dev"` or `"resdev"`
#' @noRd
get.theta.dev <- function(mbnma, param="theta") {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertChoice(param, choices = c("dev", "resdev", "theta"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  dev.df <- as.data.frame(
    mbnma$BUGSoutput$summary[grep(paste0("^", param, "\\["),
                                  rownames(mbnma$BUGSoutput$summary)),]
  )

  if (mbnma$type=="time") {
    # Takes the study, arm and follow-up measurement identifier for each residual deviance point
    id <- matrix(unlist(strsplit(
      gsub(
        "([a-z]+\\[)([0-9]+)(,)([0-9]+)(,)([0-9]+)(\\])",
        "\\2.\\4.\\6", rownames(dev.df)),
      split="\\.")),
      byrow=TRUE,
      ncol=3
    )

    dev.df$fupdose <- as.numeric(id[,3])
    dev.df$groupvar <- paste(as.numeric(id[,1]), as.numeric(id[,2]), sep="_")

    # Treatment as facet
    temp <- replicate(max(dev.df$fupdose), mbnma$model$data()$treatment)
    dev.df$facet <- as.vector(temp)[
      stats::complete.cases(as.vector(temp))
      ]

  } else if (mbnma$type=="dose") {
    # Takes the study, arm and follow-up measurement identifier for each residual deviance point
    id <- matrix(unlist(strsplit(
      gsub(
        "([a-z]+\\[)([0-9]+)(,)([0-9]+)(\\])",
        "\\2.\\4", rownames(dev.df)),
      split="\\.")),
      byrow=TRUE,
      ncol=2
    )

    # Agent as facet
    dev.df$facet <- as.vector(mbnma$model$data()$agent)[
      stats::complete.cases(as.vector(mbnma$model$data()$agent))
      ]

    dev.df$fupdose <- as.vector(mbnma$model$data()$dose)[
      stats::complete.cases(as.vector(mbnma$model$data()$dose))
      ]
    dev.df$groupvar <- as.numeric(id[,1])
  }

  dev.df$study <- as.numeric(id[,1])
  dev.df$arm <- as.numeric(id[,2])


  return(dev.df)
}





#' Plot fitted values from MBNMA model
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param disp.obs A boolean object to indicate whether raw data responses should be
#' plotted as points on the graph
#' @param ... Arguments to be sent to `ggplot2::geom_point()` or `ggplot2::geom_line()`
#' @inheritParams predict.mbnma
#'
#' @return Generates a plot of fitted values from the MBNMA model and returns a list containing
#' the plot (as an object of `class(c("gg", "ggplot"))`), and a data.frame of posterior mean
#' fitted values for each observation.
#'
#' @details
#' Fitted values should only be plotted for models that have converged successfully.
#' If fitted values (`theta`) have not been monitored in `mbnma$parameters.to.save`
#' then additional iterations will have to be run to get results for these.
#'
#' @examples
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an Emax dose-response MBNMA and predict responses
#' emax <- mbnma.emax(network, method="random")
#'
#' # Plot fitted values and observed values
#' fitplot(emax)
#'
#' # Plot fitted values only
#' fitplot(emax, disp.obs=FALSE)
#'
#' # A data frame of fitted values can be obtained from the object
#' #returned by `fitplot`
#' fits <- fitplot(emax)
#' head(fits$fv)
#'
#' @export
fitplot <- function(mbnma, disp.obs=TRUE, ...) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertLogical(disp.obs, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!("theta" %in% mbnma$parameters.to.save)) {
    msg <- paste0("`theta` not monitored in mbnma$parameters.to.save.\n",
                  "additional iterations will be run in order to obtain results")
    message(msg)

    theta.df <- mbnma.update(mbnma, param="theta")
  } else {
    theta.df <- get.theta.dev(mbnma, param="theta")
  }


  # Obtain raw responses to plot over fitted
  if (mbnma$model.arg$likelihood=="binomial") {
    theta <- mbnma$model$data()$r / mbnma$model$data()$N
  } else if (mbnma$model.arg$likelihood=="poisson") {
    theta <- mbnma$model$data()$r / mbnma$model$data()$E
  } else if (mbnma$model.arg$likelihood=="normal") {
    theta <- mbnma$model$data()$y
  }

  theta <- rescale.link(theta, direction="link", link=mbnma$model.arg$link)
  theta.df$y <- as.vector(theta)[stats::complete.cases(as.vector(theta))]

  # Add E0 for each agent and drop dose==0
  E0 <- mean(theta.df$mean[theta.df$fupdose==0])
  theta.df <- theta.df[theta.df$fupdose!=0,]

  for (i in seq_along(unique(theta.df$study))) {
    subset <- theta.df[theta.df$study==unique(theta.df$study)[i],]
    for (k in seq_along(unique(subset$facet))) {
      row <- subset[subset$facet==unique(subset$facet)[k],][1,]
      row$mean <- E0
      row$y <- NA
      row$fupdose <- 0
      theta.df <- rbind(theta.df, row)
    }
  }

  theta.df <- dplyr::arrange(theta.df, theta.df$study, theta.df$facet, theta.df$fupdose)

  # Axis labels
  if (mbnma$type=="time") {
    xlab <- "Follow-up count"
    facetscale <- "fixed"
  } else if (mbnma$type=="dose") {
    xlab <- "Dose"
    facetscale <- "free_x"
  }
  ylab <- "Response on link scale"

  # Add facet labels
  if ("agents" %in% names(mbnma)) {
    if (mbnma$agents[1]=="Placebo" & mbnma$treatments[1]=="Placebo_0") {
      labs <- mbnma$agents[-1]
    } else {labs <- mbnma$agents}
    theta.df$facet <- factor(theta.df$facet, labels=labs)
  } else if ("treatments" %in% names(mbnma)) {
    labs <- mbnma$treatments[-1]
    theta.df$facet <- factor(theta.df$facet, labels=labs)
  }

  # Generate plot
  g <- ggplot2::ggplot(theta.df,
                       ggplot2::aes(x=theta.df$fupdose, y=mean, group=theta.df$groupvar)) +
    ggplot2::geom_line()

  # Overlay observed responses
  if (disp.obs==TRUE) {
    g <- g + ggplot2::geom_point(ggplot2::aes(y=theta.df$y), size=1)
  }

  # Add facets
  g <- g + ggplot2::facet_wrap(~theta.df$facet, scales = facetscale)

  # Add axis labels
  g <- g + ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab)

  suppressWarnings(graphics::plot(g))

  return(invisible(list("graph"=g, "fv"=theta.df)))
}






#' Plot histograms of rankings from MBNMA models
#'
#' @param x An object of class `"mbnma.rank"` generated by `rank.mbnma()`
#' @param treat.labs A vector of treatment labels in the same order as treatment codes.
#' Easiest to use treatment labels stored by `mbnma.network()`
#' @param ... Arguments to be sent to `ggplot::geom_bar()`
#' @inheritParams rank.mbnma
#'
#' @return A series of histograms that show rankings for each treatment/agent/prediction, with a
#' separate panel for each parameter
#' The object returned is a list containing a separate element for each parameter in `params`
#' which is an object of `class(c("gg", "ggplot"))`.
#'
#' @examples
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Estimate rankings  from an Emax dose-response MBNMA
#' emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#' ranks <- rank(emax)
#'
#' # Plot rankings for both dose-response parameters (in two separate plots)
#' plot(ranks)
#'
#' # Plot rankings just for ED50
#' plot(ranks, params="d.ed50")
#'
#' # Plot rankings from prediction
#' doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
#' pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
#'             exact.doses=doses)
#' rank <- rank(pred)
#' plot(rank)
#' @export
plot.mbnma.rank <- function(x, params=NULL, treat.labs=NULL, ...) {
  # ... are commands to be sent to geom_histogram

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma.rank", add=argcheck)
  checkmate::assertCharacter(params, null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(treat.labs, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  output <- list()

  if (is.null(params)) {
    params <- names(x)
  }

  for (param in seq_along(params)) {

    rank.mat <- x[[params[param]]]$rank.matrix
    treats <- colnames(rank.mat)

    ranks <- vector()
    treat <- vector()
    for (i in seq_along(treats)) {
      treat <- append(treat, rep(treats[i], nrow(rank.mat)))
      ranks <- append(ranks, rank.mat[,i])
    }
    data <- data.frame("ranks"=ranks, "treat"=treat)

    if (!is.null(treat.labs)) {
      if (length(treat.labs)!=length(unique(data$treat))) {
        stop("`treat.labs` must be a character vector of the same length as the number of ranked tretments/agents")
      }
      data$treat <- factor(data$treat, labels=treat.labs)
    } else {
      #data$treat <- factor(as.numeric(as.character(data$treat)))
      data$treat <- factor(data$treat)
    }

    g <- ggplot2::ggplot(data, ggplot2::aes(x=ranks)) +
      ggplot2::geom_bar(...) +
      ggplot2::xlab("Rank (1 = best)") +
      ggplot2::ylab("MCMC iterations") +
      ggplot2::facet_wrap(~treat) +
      ggplot2::ggtitle(params[param])

    graphics::plot(g)

    output[[params[param]]] <- g
  }

  return(invisible(output))
}




#' @describeIn nma.nodesplit Plot outputs from treatment-level nodesplit models
#'
#' @param x An object of `class("nma.nodesplit")`
#' @param plot.type A character string that can take the value of `"forest"` to plot
#' only forest plots, `"density"` to plot only density plots, or left as `NULL` (the
#' default) to plot both types of plot.
#' @param ... Arguments to be sent to [ggplot2::ggplot()]
#'
#' @details The S3 method `plot()` on an `nma.nodesplit` object generates either
#' forest plots of posterior medians and 95\\% credible intervals, or density plots
#' of posterior densities for direct and indirect evidence.
#'
#' @return Plots the desired graph(s) and returns an object (or list of object if
#' `plot.type=NULL`) of `class(c("gg", "ggplot"))`
#'
#' @export
plot.nma.nodesplit <- function(x, plot.type=NULL, ...) {
  # ... are commands to be sent to geom_histogram

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "nma.nodesplit", add=argcheck)
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
      ggplot2::ggplot(data=forestdata, ggplot2::aes(x=forestdata$source, y=forestdata$med,
                                                    ymin=forestdata$l95, ymax=forestdata$u95), ...) +
      ggplot2::geom_pointrange() +
      ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
      ggplot2::xlab("") + ggplot2::ylab("Treatment effect (95% CrI)") +
      ggplot2::theme(axis.text = ggplot2::element_text(size=10),
                     axis.title = ggplot2::element_text(size=12),
                     title=ggplot2::element_text(size=18)) +
      ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,1),"cm")) +
      ggplot2::facet_wrap(~factor(forestdata$comp))
  }
  if ("density" %in% plot.type) {

    density <- ggplot2::ggplot(densitydata, ggplot2::aes(x=densitydata$value,
                                                         linetype=densitydata$Estimate, fill=densitydata$Estimate), ...) +
      ggplot2::geom_density(alpha=0.2) +
      ggplot2::xlab("Treatment effect (95% CrI)") +
      ggplot2::ylab("Posterior density") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                     axis.title = ggplot2::element_text(size=14)) +
      ggplot2::facet_wrap(~factor(densitydata$comp))
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
  #checkmate::assertNumeric(intercept, len=1, add=argcheck)
  checkmate::assertLogical(bydose, len=1, add=argcheck)
  checkmate::assertChoice(scales, c("free_x", "fixed"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  intercept <- 0 # Leaving here in case want to allow user to change it at later date
  split.df <- x[["jagsresult"]]$BUGSoutput$summary

  # Check if NMA is from UME model
  if (any(grepl("^d\\[[0-9+],[0-9]+\\]", rownames(x$jagsresult$BUGSoutput$summary)))) {
    split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+,1\\]", rownames(split.df)), c(3,5,7)])
    #split.df <- as.data.frame(split.df[grepl("^d\\[[0-9+],1\\]", rownames(split.df)), c(3,5,7)])
  } else {
    split.df <- as.data.frame(split.df[grepl("^d\\[[0-9]+\\]", rownames(split.df)), c(3,5,7)])
  }


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

  # Plot faceted by agent as dose-response splitplot
  if (bydose==TRUE) {

    # Add intercept for all agents
    agents <- unique(split.df$agent)
    agents <- agents[agents!="Placebo"]
    for (i in seq_along(agents)) {
      row <- split.df[split.df$agent=="Placebo",]
      row$agent <- agents[i]
      split.df <- rbind(row, split.df)
    }
    split.df <- split.df[split.df$agent!="Placebo",]

    if (intercept!=0) {
      split.df[,1:3] <- split.df[,1:3] + intercept
    }

    g <- ggplot2::ggplot(split.df, ggplot2::aes(x=dose, y=`50%`), ...) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`)) +
      ggplot2::facet_wrap(~factor(agent), scales = scales) +
      ggplot2::xlab("Dose") +
      ggplot2::ylab("Effect size on link scale")

  } else if (bydose==FALSE) {
    # Plot conventional forest plot
    split.df$treatment <- factor(split.df$treatment, levels=x[["trt.labs"]])

    g <- ggplot2::ggplot(split.df, ggplot2::aes(y=`50%`, x=treatment), ...) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`)) +
      ggplot2::coord_flip() +
      ggplot2::ylab("Effect size on link scale") +
      ggplot2::xlab("Treatment")
  }

  graphics::plot(g)
  return(invisible(g))
}
