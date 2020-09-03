##############################################
#### Functions for class("mbnma.network") ####
##############################################

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "2.5%", "97.5%", "50%",
                                                        "treatment", "studyID", "agent", "dose"))

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




#' Print summary mbnma.network information to the console
#'
#' @param x An object of class `mbnma.network`.
#' @param ... further arguments passed to or from other methods
#'
#' @export
summary.mbnma.network <- function(x,...) {

  cat(crayon::underline(crayon::bold("Description:", x$description, "\n")))
  cat("Number of studies:", length(unique(x$data.ab$studyID)), "\n")
  cat("Number of treatments:", length(x$treatments), "\n")
  cat("Number of agents:", length(x$agents), "\n")

  if ("classes" %in% names(x)) {
    cat("Number of classes:", length(x$classes), "\n")
  }

  # Count doses per agent
  agentdf <- unique(x$data.ab[, c("treatment", "agent")])
  agentdf <- agentdf %>% dplyr::group_by(agent) %>%
    dplyr::mutate(ndose = dplyr::n())
  agentdf <- dplyr::arrange(agentdf, agent)[,c("agent", "ndose")]
  agentdf <- unique(agentdf)

  cat("Median (max, min) doses per agent: ", median(agentdf$ndose),
      " (", min(agentdf$ndose), ", ", max(agentdf$ndose), ")\n", sep="")

  # Check network is connected at agent-level
  g <- suppressWarnings(plot.invisible(x, level="agent"))
  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               to=1))
  if (any(connects==FALSE)) {
    cat("Agent-level network is", crayon::bold(crayon::red("DISCONNECTED"), "\n"))
  } else {
    cat("Agent-level network is", crayon::bold(crayon::green("CONNECTED"), "\n"))
  }

  # Check network is connected at treatment-level
  g <- suppressWarnings(plot.invisible(x, level="treatment"))
  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               to=1))
  if (any(connects==FALSE)) {
    cat("Treatment-level network is", crayon::bold(crayon::red("DISCONNECTED"), "\n"))
  } else {
    cat("Ttreatment-level network is", crayon::bold(crayon::green("CONNECTED"), "\n"))
  }
  invisible(x)
}





#' @describeIn mbnma.network Generate a network plot
#' @inheritParams mbnma.run
#'
#' @param x An object of class `mbnma.network`.
#' @param layout An igraph layout specification. This is a function specifying an igraph
#'   layout that determines the arrangement of the vertices (nodes). The default
#'   `igraph::as_circle()` arranged vertices in a circle. Two other useful layouts for
#'   network plots are: `igraph::as_star()`, `igraph::with_fr()`. Others can be found
#'   in \code{\link[igraph]{layout_}}
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
#'   the default of 0 is used. Option only applicable if `layout_in_circle` is
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
#' # Arrange network plot in a star with the reference treatment in the centre
#' plot(network, layout=igraph::as_star(), label.distance=3)
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
                               layout=igraph::in_circle(), remove.loops=FALSE,
                               edge.scale=1, v.scale=NULL, label.distance=0,
                               ...)
  # Requires igraph
  #S3method(plot, mbnma.network)

  # x is an object of class mbnma.network

{
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(x, "mbnma.network", add=argcheck)
  checkmate::assertClass(layout, "igraph_layout_spec", add=argcheck)
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
    comparisons <- comparisons[comparisons$t1!=comparisons$t2,]
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

  # Add attributes
  igraph::V(g)$label.dist <- label.distance
  if (!is.null(node.size)) {igraph::V(g)$size <- node.size}
  igraph::E(g)$width <- edge.scale * comparisons[["nr"]]

  # Change label locations if layout_in_circle
  laycheck <- as.character(layout)[2]
  if (any(
    grepl("layout_in_circle", laycheck) |
    grepl("layout_as_star", laycheck))) {
    lab.locs <- radian.rescale(x=seq(1:length(nodes)), direction=-1, start=0)
    igraph::V(g)$label.degree <- lab.locs
  }

  # Plot netgraph
  layout <- igraph::layout_(g, layout)
  igraph::plot.igraph(g,
                      layout = layout,
                      ...
  )

  if (v.color=="agent") {
    legend("bottomleft", x$agents, pt.bg=unique(igraph::V(g)$color), pch=21, pt.cex=1.5, cex=0.8)
  }

  # if (layout_in_circle==TRUE) {
  #   lab.locs <- radian.rescale(x=seq(1:length(nodes)), direction=-1, start=0)
  #   igraph::V(g)$label.degree <- lab.locs
  #   igraph::plot.igraph(g,
  #                       layout = igraph::layout_in_circle(g),
  #                       ...
  #   )
  # } else {
  #   igraph::plot.igraph(g, ...)
  # }

  if (!is.null(doselink)) {
    message(paste0("Dose-response connections to placebo plotted based on a dose-response
                   function with ", (doselink-1),
                   " degrees of freedom"))
  }

  return(invisible(g))
}
