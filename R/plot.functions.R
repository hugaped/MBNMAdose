# Functions for plotting in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-18



#' @describeIn MBNMA.network Generate a network plot
#' @inheritParams MBNMA.run
#'
#' @param layout_in_circle A boolean value indicating whether the network plot
#'   should be shown in a circle or left as the igraph default layout.
#' @param edge.scale A number to scale the thickness of connecting lines
#'   (edges). Line thickness is proportional to the number of studies for a
#'   given comparison. Set to `0` to make thickness equal for all comparisons.
#' @param v.color Can take either `"connect"` (the default) to indicate that nodes should
#'   only be coloured if they are connected to the network reference treatment (indicates
#'   network connectivity) or `"agent"` to colour nodes by agent.
#' @param v.scale A number with which to scale the size of the nodes. If the variable `N`
#'   (to indicate the numbers of participants in each study arm) is included in the
#'   dataset then the size of the nodes will be proportional to the number of participants
#'   within a treatment/agent in the network.
#' @param label.distance A number scaling the distance of labels from the nodes
#'   to improve readability. The labels will be directly on top of the nodes if
#'   the defauls of `0` is used. Option only applicable if `layout_in_circle` is
#'   set to `TRUE`.
#' @param level A string indicating whether nodes/facets should represent `treatment`
#'   or `agent` in the plot. Can be used to examine the expected impact of modelling
#'   dose-response in terms of network connectivity.
#' @param remove.loops A boolean value indicating whether to include loops that
#'   indicate comparisons within a node.
#' @param ... Options for plotting in `igraph`.
#'
#' @details The S3 method `plot()` on an `MBNMA.network` object generates a
#'   network plot that shows how different treatments are connected within the
#'   network via study comparisons. This can be used to identify how direct and
#'   indirect evidence are informing different treatment comparisons. Depends on
#'   \code{\link[igraph]{igraph}}.
#'
#' @examples
#' # Create an MBNMA.network object from the data
#' network <- MBNMA.network(HF2PPITT)
#'
#' # Generate a network plot from the data
#' plot(network, layout_in_circle=TRUE)
#'
#' # Generate a network plot at the agent level that removes loops indicating comparisons
#' #within a node
#' plot(network, level="agent", remove.loops=TRUE)
#'
#' # Generate a network plot at the treatment level that colours nodes by agent
#' plot(network, v.color="agent", remove.loops=TRUE)
#' @export
plot.MBNMA.network <- function(network, layout_in_circle = TRUE, edge.scale=1, label.distance=0,
                               level="treatment", remove.loops=FALSE, v.color="connect",
                               v.scale=NULL,
                               ...)
  # Requires igraph
  #S3method(plot, MBNMA.network)

  # network is an object of class MBNMA.network

{
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "MBNMA.network", add=argcheck)
  checkmate::assertLogical(layout_in_circle, len=1, add=argcheck)
  checkmate::assertNumeric(edge.scale, finite=TRUE, len=1, add=argcheck)
  checkmate::assertNumeric(label.distance, finite=TRUE, len=1, add=argcheck)
  checkmate::assertNumeric(v.scale, lower = 0, finite=TRUE, null.ok=TRUE, len=1, add=argcheck)
  checkmate::assertChoice(level, choices = c("treatment", "agent"), add=argcheck)
  checkmate::assertChoice(v.color, choices = c("connect", "agent"), add=argcheck)
  checkmate::assertLogical(remove.loops, len=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Generate comparisons (using get.latest.time and MBNMA.contrast?
  data.ab <- network$data.ab


  # Check if level="agent" that agents are present in dataset
  temp <- factor(data.ab$agent, labels=network$agents)
  if (level=="agent") {
    if (!("agents" %in% names(network))) {
      stop("`level` has been set to `agent` but there are no agent codes given in the dataset")
    }
    nodes <- sort(network[["agents"]])
    data.ab$treatment <- as.character(temp)

  } else if (level=="treatment") {
    data.ab$treatment <- paste(as.character(temp), data.ab$dose, sep="_")
    nodes <- sort(network[["treatments"]])
    #nodes <- sort(unique(data.ab$treatment))
  }

  # Calculate participant numbers (if v.scale not NULL)
  if (!is.null(v.scale)) {
    if (!("N" %in% names(data.ab))) {
      stop("`N` not included as a column in dataset. Vertices/nodes will all be scaled to be the same size.")
    }

    size.vec <- vector()
    for (i in seq_along(nodes)) {
      size.vec <- append(size.vec, sum(data.ab$N[data.ab$treatment==nodes[i]]))
    }
    # Scale size.vec by the max node.size
    size.vec <- size.vec/ (max(size.vec)/20)

    node.size <-
      setNames(size.vec, nodes)
    node.size <- as.matrix(node.size*v.scale)
  } else {
    node.size <- NULL
  }


  comparisons <- MBNMA.comparisons(data.ab)

  # Code to make graph.create as an MBNMA command if needed
  g <- igraph::graph.empty()
  g <- g + igraph::vertex(nodes)
  ed <- t(matrix(c(comparisons[["t1"]], comparisons[["t2"]]), ncol = 2))
  edges <- igraph::edges(as.vector(ed), weight = comparisons[["nr"]], arrow.mode=0)
  g <- g + edges

  if (remove.loops==TRUE) {
    g <- igraph::simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
  }

  # Check network is connected and produce warning message if not
  disconnects <- check.network(g)
  if (v.color=="connect") {
    igraph::V(g)$color <- "SkyBlue2"
    igraph::V(g)$color[which(names(igraph::V(g)) %in% disconnects)] <- "white"
  } else if (v.color=="agent") {
    if (!("agents" %in% names(network))) {
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

  return(invisible(g))
}






#' Check if all nodes in the network are connected (identical to MBNMAtime)
check.network <- function(g, reference=1) {
  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               reference))
  treats <- colnames(connects)[connects==FALSE]

  if (length(treats>0)) {
    warning(paste0("The following treatments are not connected to the network reference:\n",
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
#' @references
#' https://gist.github.com/kjhealy/834774/a4e677401fd6e4c319135dabeaf9894393f9392c
radian.rescale <- function(x, start=0, direction=1) {
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
}




#' Get large vector of distinct colours using Rcolorbrewer
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
#' @inheritParams predict.MBNMA
#' @param params A character vector of dose-response parameters to plot.
#' Parameters must be given the same name as monitored nodes in `mbnma` and must be
#' modelled as relative effects (`"rel"`). Can be set to
#' `NULL` to include all available dose-response parameters estimated by `mbnma`.
#' @param agent.labs A caracter vector of agent labels (including `"placebo"` only if it
#' has been included in the original network). This can be easily inputted by using the
#' object of `class("MBNMA.network")` on which the MBNMA model was run (i.e. `network$agents`).
#' @param class.labs A character vector of class labels if `mbnma` was modelled using class effects
#' (including `"placebo"` only if it has been included in the original network). This can be
#' easily inputted by using the object of `class("MBNMA.network")` on which the MBNMA model
#' was run (i.e. `network$classes`).
#'
#' @return A forest plot of class `c("gg", "ggplot")` that has separate panels for different dose-response parameters
#' @export
plot.MBNMA <- function(mbnma, params=NULL, agent.labs=NULL, class.labs=NULL) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "MBNMA", add=argcheck)
  checkmate::assertChoice(params, choices=mbnma[["parameters.to.save"]], null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check that specified params are modelled using relative effects
  for (i in seq_along(params)) {
    if (grepl("d\\.", params[i]) | grepl("D\\.", params[i])) {
      stop(paste0(params[i], " has not been modelled using relative effects and does not vary by agent or class"))
    }
  }

  # Add all available params if is.null(params)
  if (is.null(params)) {
    params <- vector()

    # Add d
    params <- append(params,
                     mbnma[["parameters.to.save"]][grep("^d\\.", mbnma[["parameters.to.save"]])]
    )

    # Add D
    params <- append(params,
                     mbnma[["parameters.to.save"]][grep("^D\\.", mbnma[["parameters.to.save"]])]
    )

    if (length(params)==0) {
      stop("No dose-response relative effects can be identified from the model")
    }
  }


  # Compile parameter data into one data frame
  mbnma.sum <- as.data.frame(mbnma[["BUGSoutput"]][["summary"]])
  plotdata <- mbnma.sum[0,]
  for (i in seq_along(params)) {
    paramdata <- mbnma.sum[grepl(paste0("^", params[i]),rownames(mbnma.sum)),]
    paramdata[["doseparam"]] <- rep(params[i], nrow(paramdata))
    plotdata <- rbind(plotdata, paramdata)
  }
  plotdata[["param"]] <- as.numeric(gsub("(.+\\[)([0-9]+)(\\])", "\\2", rownames(plotdata)))

  # Change param labels for agents
  agentdat <- plotdata[grepl("^d\\.", rownames(plotdata)),]
  if (!is.null(agent.labs)) {
    agentcodes <- as.numeric(gsub("(^.+\\[)([0-9]+)(\\])", "\\2", rownames(agentdat)))
    a.labs <- agent.labs[sort(unique(agentcodes))]
  } else {
    a.labs <- sort(unique(agentdat$param))
  }

  # Change param labels for classes
  classdat <- plotdata[grepl("^D\\.", rownames(plotdata)),]
  if (nrow(classdat)!=0) {
    if (!is.null(class.labs)) {
      classcodes <- as.numeric(gsub("(^.+\\[)([0-9]+)(\\])", "\\2", rownames(classdat)))
      c.labs <- class.labs[classcodes]
    } else {
      c.labs <- sort(unique(classdat$param))
    }
  }
  # Increase param number for classes
  plotdata$param[grepl("^D\\.", rownames(plotdata))] <-
    plotdata$param[grepl("^D\\.", rownames(plotdata))] + max(agentdat$param)

  # Attach labels
  all.labs <- c(a.labs, c.labs)
  plotdata$param <- factor(plotdata$param, labels=all.labs)

  g <- ggplot2::ggplot(plotdata, ggplot2::aes(y=`50%`, x=param)) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=`2.5%`, ymax=`97.5%`)) +
    ggplot2::coord_flip()

  g <- g + ggplot2::facet_wrap(~doseparam, scales="free")

  # Axis labels
  g <- g + ggplot2::xlab("Agent / Class") +
    ggplot2::ylab("Effect size")

  return(g)
}
