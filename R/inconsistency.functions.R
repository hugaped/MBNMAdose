# Functions for assessing inconsistency in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-30


#' Node-splitting model for testing consistency at the treatment level
#'
#' Splits contributions for a given set of treatment comparisons into direct and
#' indirect evidence. A discrepancy between the two suggests that the consistency
#' assumption required for NMA and MBNMA may violated.
#'
#' @param drop.discon A boolean object that indicates whether to drop treatments
#' that are disconnected at the treatment level. Default is `TRUE`. If set to `FALSE` then
#' this could lead to indentification of nodesplit comparisons that are not connected
#' to the network reference treatment.
#' @param comparisons A matrix specifying the comparisons to be split (one row per comparison).
#' The matrix must have two columns indicating each treatment for each comparison. Values can
#' either be character (corresponding to the treatment names given in `network`) or
#' numerical (corresponding to treatment codes within the `network` - note that these
#' may change if `drop.discon=TRUE`).
#' @param ... Arguments to be sent to `R2jags`
#' @inheritParams MBNMA.run
#' @inheritParams MBNMA.network
#'
#' @examples
#' # Using the triptans data
#' network <- MBNMA.network(HF2PPITT)
#'
#' split <- MBNMA.nodesplit(network, likelihood = "binomial", link="logit",
#'   method="common")
#'
#'
#' #### To perform nodesplit on selected comparisons ####
#'
#' # Check for closed loops of treatments with independent evidence sources
#' loops <- inconsistency.loops(network$data.ab)
#'
#' split <- MBNMA.nodesplit(network, likelihood = "binomial", link="logit",
#'   method="random", comparisons=rbind(c(6,23), c(6,12)))
#'
#' # Drop treatments that are disconnected from the network in the analysis
#' split <- MBNMA.nodesplit(net.noplac, likelihood = "binomial", link="logit",
#'   method="random", drop.discon=TRUE)
#'
#' # Plot results
#' plot(split, plot.type="density") # Plot density plots of posterior densities
#' plot(split, plot.type="forest") # Plot forest plots of direct and indirect evidence
#'
#' # Print and summarise results
#' print(split)
#' summary(split) # Generate a data frame of summary results
#' @export
MBNMA.nodesplit <- function(network, likelihood="binomial", link="logit", method="common",
                            drop.discon=FALSE, comparisons=NULL,
                            ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "MBNMA.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::assertLogical(drop.discon, add=argcheck)
  #checkmate::assertDataFrame(comparisons, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  # Load data
    # Check treatments that are not connected and remove if not
  if (drop.discon==TRUE) {
    connect <- drop.disconnected(network)
    data.ab <- connect[["data.ab"]]
    trt.labs <- connect[["trt.labs"]]
  } else if (drop.discon==FALSE) {
    data.ab <- network$data.ab
    trt.labs <- network$treatments
  }

  # Identify closed loops of treatments
  if (is.null(comparisons)) {
    comparisons <- inconsistency.loops(data.ab)
  } else {
    if (class(comparisons)=="data.frame") {
      if (all(c("t1", "t2") %in% names(comparisons))) {
        comparisons <- data.frame(comparisons$t1, comparisons$t2)
      }
      comparisons <- as.matrix(comparisons)
    }
    if (ncol(comparisons)!=2) {
      stop("`comparisons` must be a matrix of comparisons on which to split containing exactly two columns")
    }
    if (is.character(comparisons)) {
      if (!all(comparisons %in% trt.labs)) {
        stop("Treatment names given in `comparisons` do not match those within `network` or match treatments that have been dropped from the network due to being disconnected")
      }
      comparisons <- matrix(as.numeric(factor(comparisons, levels=trt.labs)), ncol=2)
    }
    if (!all(comparisons %in% data.ab$treatment)) {
      stop("Treatment codes given in `comparisons` do not match those within `network` or match treatments that have been dropped from the network due to being disconnected")
    }

    # Sort comparisons so that lowest is t1
    comparisons <- t(apply(comparisons, MARGIN=1, FUN=function(x) {sort(x)}))

    # Ensure comparisons are nested within inconsistency.loops
    check <- paste(comparisons[,1], comparisons[,2], sep="_")
    fullcomp <- inconsistency.loops(data.ab)
    match <- match(check, paste(fullcomp[,1], fullcomp[,2], sep="_"))
    if (any(is.na(match))) {
      out <- comparisons[is.na(match),]
      out <- matrix(unlist(lapply(out, FUN=function(x) {trt.labs[x]})), ncol=2)
      printout <- c()
      for (i in 1:nrow(out)) {
        printout <- paste(printout, paste(out[i,], collapse=" "), sep="\n")
      }
      stop(cat(paste0("The following `comparisons` are not part of closed loops of treatments informed by direct and indirect evidence from independent sources:\n",
                      printout, "\n\n")))
    }
  }


  ##### Run NMA #####
  nma.net <- suppressMessages(MBNMA.network(data.ab))
  nma.jags <- NMA.run(nma.net, method=method,
                      likelihood=likelihood, link=link,
                      warn.rhat=FALSE, drop.discon = FALSE)#, ...)

  nodesplit.result <- list()
  for (split in seq_along(comparisons[,1])) {

    comp <- as.numeric(comparisons[split,1:2])
    print(paste0("Calculating nodesplit for: ",
                 paste0(trt.labs[comp[2]], " vs ", trt.labs[comp[1]])))

    ##### Estimate NMA #####
    nma.res <- nma.jags$jagsresult$BUGSoutput$sims.list$d[,comp[2]] -
      nma.jags$jagsresult$BUGSoutput$sims.list$d[,comp[1]]

    ##### Estimate Indirect #####

    ind.df <- data.ab

    # Drop studies/comparisons that compare comps
    dropID <- vector()
    dropcomp <- vector()
    studies <- unique(ind.df$studyID)
    for (study in seq_along(studies)) {
      subset <- ind.df[ind.df$studyID==studies[study],]
      if (all(comp %in% subset$treatment)) {
        if (subset$narm[1]<=2) {
          dropID <- append(dropID, subset$studyID[1])
        } else if (subset$narm[1]>2) {
          dropcomp <- append(dropcomp, subset$studyID[1])
        }
      }
    }

    # Drop studies
    ind.df <- ind.df[!(ind.df$studyID %in% dropID),]

    # Drop comparisons from studies
    ind.df <- drop.comp(ind.df, drops=dropcomp, comp=comp)
    # stoploop <- FALSE
    # while(stoploop==FALSE) {
    #   temp <- drop.comp(ind.df, drops=dropcomp, comp=comp)
    #   temp.net <- MBNMA.network(temp)
    #   nt <- length(temp.net$treatments)
    #   if (nt==length(nma.net$treatments)) {
    #     g <- plot(temp.net, doseparam=1000)
    #     connectcheck <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
    #                                                      to=1)[
    #                                                        c(comp[1], comp[2])
    #                                                        ])
    #     if (all(connectcheck==TRUE)) {
    #       ind.df <- temp
    #       stoploop <- TRUE
    #     }
    #   }
    #   #print("Restarting drop.comp")
    # }


    # Run NMA
    ind.net <- suppressMessages(MBNMA.network(ind.df))
    ind.jags <- NMA.run(ind.net, method=method, likelihood=likelihood, link=link,
                        warn.rhat=FALSE, drop.discon = FALSE)#, ...)
    ind.res <- ind.jags$jagsresult$BUGSoutput$sims.list$d[,comp[2]] -
      ind.jags$jagsresult$BUGSoutput$sims.list$d[,comp[1]]


    ##### Estimate Direct #####
    dir.net <- suppressMessages(change.netref(MBNMA.network(data.ab), ref=comp[1]))
    dir.jags <- NMA.run(dir.net, method=method, likelihood=likelihood, link=link,
                        warn.rhat=FALSE, drop.discon=FALSE, UME=TRUE)#, ...)
    dir.res <- dir.jags$jagsresult$BUGSoutput$sims.matrix[
      ,colnames(dir.jags$jagsresult$BUGSoutput$sims.matrix)==paste0("d[", comp[2],",1]")
    ]


    ##### Generate plots/results #####

    # Overlaps
    overlap.mat <- list("direct"=dir.res, "indirect"=ind.res)
    overlap <- overlapping::overlap(overlap.mat, plot=FALSE)
    p.values <- overlap$OV

    # Quantiles
    quantile_dif <- quantile(ind.res - dir.res, c(0.025, 0.5, 0.975))
    quantile_dir <- quantile(dir.res, c(0.025, 0.5, 0.975))
    quantile_ind <- quantile(ind.res, c(0.025, 0.5, 0.975))
    quantile_nma <- quantile(nma.res, c(0.025, 0.5, 0.975))
    quantiles <- list("difference" = quantile_dif, "nma"=quantile_nma,
                      "direct"=quantile_dir, "indirect"=quantile_ind)


    # GGplots
    source <- c("NMA", "Direct", "Indirect")
    l95 <- c(quantile_nma[1], quantile_dir[1], quantile_ind[1])
    med <- c(quantile_nma[2], quantile_dir[2], quantile_ind[2])
    u95 <- c(quantile_nma[3], quantile_dir[3], quantile_ind[3])
    plotdata <- data.frame(source, l95, med, u95)

    title <- paste0(trt.labs[comp[2]], " vs ", trt.labs[comp[1]])

    gg <-
      ggplot2::ggplot(data=plotdata, ggplot2::aes(x=source, y=med, ymin=l95, ymax=u95)) +
      ggplot2::geom_pointrange() +
      ggplot2::coord_flip() +  # flip coordinates (puts labels on y axis)
      ggplot2::xlab("") + ggplot2::ylab("Treatment effect (95% CrI)") + ggplot2::ggtitle(title) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                     axis.title = ggplot2::element_text(size=18),
                     title=ggplot2::element_text(size=18)) +
      ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,1),"cm"))

    # Density plots (with shaded area of overlap)
    molten <- data.frame(ind.res, dir.res)
    molten <- suppressMessages(reshape2::melt(molten))
    names(molten) <- c("Estimate", "value")
    linetypes <- c("solid", "dash")
    levels(molten$Estimate) <- c("Indirect", "Direct")

    dens <- ggplot2::ggplot(molten, ggplot2::aes(x=value, linetype=Estimate, fill=Estimate)) +
      ggplot2::geom_density(alpha=0.2) +
      ggplot2::xlab(title) +
      ggplot2::ylab("Posterior density") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                     axis.title = ggplot2::element_text(size=14))

    nodesplit <- list("comparison"= c(trt.labs[comp[2]], trt.labs[comp[1]]),
                      "direct"=dir.res, "indirect"=ind.res, "nma"=nma.res,
                      "overlap matrix"=overlap.mat,
                      "p.values"=p.values, "quantiles"=quantiles,
                      "forest.plot"=gg, "density.plot"=dens,
                      "direct.model"=dir.jags, "indirect.model"=ind.jags,
                      "nma.model"=nma.jags)

    nodesplit.result[[paste0(comp[2], "v", comp[1])]] <-
      nodesplit
  }

  class(nodesplit.result) <- "MBNMA.nodesplit"

  return(nodesplit.result)
}






#' Identify comparisons in loops that fulfill criteria for node-splitting
#'
#' Identify comparisons informed by both direct and indirect evidence from
#' independent sources, which therefore fulfill the criteria for testing for
#' inconsistency via node-splitting. Follows the method of van Valkenhoef \insertCite{RN35;textual}{MBNMAtime}.
#'
#' @param data A data frame containing variables `studyID` and `treatment` (as
#'   numeric codes) that indicate which treatments are used in which studies.
#'
#' @details Similar to \code{\link[gemtc]{mtc.nodesplit.comparisons}} but uses a fixed
#'   reference treatment and therefore suggests fewer loops in which to test for
#'   inconsistency. Heterogeneity can also be parameterised as inconsistency and
#'   so testing for inconsistency in additional loops whilst changing the
#'   reference treatment would also be identifying heterogeneity. Depends on
#'   \code{\link[igraph]{igraph}}.
#'
#' @return A data frame of comparisons that are informed by direct and indirect
#'   evidence from independent sources. Each row of the data frame is a
#'   different treatment comparison. Numerical codes in `t1` and `t2` correspond
#'   to treatment codes.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' inconsistency.loops(data)
#' @export
inconsistency.loops <- function(data)
{
  # Assert checks
  checkmate::assertDataFrame(data)

  treatments <- factor(unique(data$treatment))

  data <- data %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(design=list(as.numeric(treatment)))


  comparisons <- ref.comparisons(data)

  splits1 <- vector()
  splits2 <- vector()
  paths <- vector()
  loops <- vector()

  pb <- txtProgressBar(0, nrow(comparisons), style = 3)
  for (i in 1:nrow(comparisons)) {
    setTxtProgressBar(pb, i)

    drops <- comparisons[-i,]

    # Alternative graph create (non-gemtc)
    g <- igraph::graph.empty()
    g <- g + igraph::vertex(levels(treatments))
    #g <- g + igraph::edges.create(new.comparisons, arrow.mode=0)

    ed <- t(matrix(c(drops[["t1"]], drops[["t2"]]), ncol = 2))
    edges <- igraph::edges(as.vector(ed), arrow.mode=0)
    g <- g + edges

    # Check whether there is still an indirect connection once direct evidence studies are removed
    if (as.logical(is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                                    comparisons[i,1], comparisons[i,2]))) == TRUE) {

      # Check if dropping 2-arm studies with both treatments and then either arm from multi-arm
      #would lead to disconnected network
      check <- suppressMessages(suppressWarnings(
        check.indirect.drops(data, comp=c(as.numeric(comparisons[i,1]),
                                       as.numeric(comparisons[i,2])))
      ))

      if (!is.null(check)) {
        # Identify the path made by the indirect evidence
        path <- as.numeric(igraph::shortest_paths(igraph::as.undirected(g),
                                                  comparisons[i,1], comparisons[i,2],
                                                  weights=NA
        )[["vpath"]][[1]])

        loop <- sort(path)

        splits1 <- append(splits1, comparisons[["t1"]][i])
        splits2 <- append(splits2, comparisons[["t2"]][i])
        paths <- append(paths, paste(path, collapse="->"))
        loops <- append(loops, paste(loop, collapse="->"))
      }
    }
  }

  splits <- data.frame("t1"=splits1, "t2"=splits2, "path"=paths, "loops"=loops)

  # Ensures only one comparison given per inconsistent loop
  splits <- splits[seq(dim(splits)[1],1),]
  splits <- splits[duplicated(splits[["loops"]])==FALSE, 1:3]

  if (nrow(splits)==0 | (nrow(splits)==1 & any(is.na(splits$t1)))) {
    stop("No closed loops of treatments arising from independent sources of evidence are present in the data - testing for consistency is not possible in this network")
  }

  return(splits)
}






#' Identify unique contrasts relative to each study reference within a network.
#' Repetitions of the same treatment comparison are grouped together.
#'
#' @param data A data frame containing variables `studyID` and `treatment` (as
#'   numeric codes) that indicate which treatments are used in which studies.
#'
#' @return A data frame of unique comparisons in which each row represents a
#'   different comparison. `t1` and `t2` indicate the treatment codes that make
#'   up the comparison. `nr` indicates the number of times the given comparison
#'   is made within the network.
#'
#'   If there is only a single observation for each study within the dataset
#'   (i.e. as for standard network meta-analysis) `nr` will represent the number
#'   of studies that compare treatments `t1` and `t2`.
#'
#'   If there are multiple observations for each study within the dataset (as in
#'   MBNMAtime) `nr` will represent the number of time points in the
#'   dataset in which treatments `t1` and `t2` are compared.
#'
#' @examples
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' ref.comparisons(data)
ref.comparisons <- function(data)
{
  # Assert checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data, add=argcheck)
  checkmate::assertNames(names(data), must.include=c("studyID", "treatment"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  checkmate::assert(
    checkmate::checkFactor(data[["treatment"]]),
    checkmate::checkNumeric(data[["treatment"]])
  )

  # if (all(names(data) %in% c("studyID", "treatment") != TRUE)) {
  #   stop("data must contain variables 'studyID' and 'treatment'")
  # }
  #
  # if (!(is.factor(data[["treatment"]]) | is.numeric(data[["treatment"]]))) {
  #   stop("`treatment` must be either factor or numeric")
  # }

  data <- dplyr::arrange(data, studyID, treatment)

  t1 <- vector()
  t2 <- vector()
  for (i in seq_along(unique(data[["studyID"]]))) {
    subset <- subset(data, studyID==unique(data[["studyID"]])[i])
    for (k in 2:nrow(subset)) {
      t1 <- append(t1, subset[["treatment"]][1])
      t2 <- append(t2, subset[["treatment"]][k])
      if (is.na(subset[["treatment"]][k])) {
        stop()
      }
    }
  }

  comparisons <- data.frame(t1 = t1, t2 = t2)
  comparisons <- comparisons %>% dplyr::group_by(t1, t2) %>%
    dplyr::mutate(nr = n())
  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)
}







#' Drop treatments from multi-arm (>2) studies for node-splitting
#'
#' Drops arms in a way which preserves connectivity and equally removes
#' data from each treatment in a nodesplit comparison (so as to maximise precision)
drop.comp <- function(ind.df, drops, comp, start=rbinom(1,1,0.5)) {
  index <- start
  #print(index)
  for (i in seq_along(drops)) {
    #index <- rbinom(1,1,0.5)
    #print(index)

    switchloop <- FALSE
    temp.df <- ind.df[!(ind.df$studyID %in% drops[i] &
                          ind.df$treatment==comp[index+1]),]

    temp.net <- suppressMessages(plot(MBNMA.network(temp.df), doseparam = 1000))

    connectcheck <- is.finite(igraph::shortest.paths(igraph::as.undirected(temp.net),
                                                     to=comp[index+1])[
                                                       c(comp[1], comp[2])
                                                       ])

    if (!(comp[index+1] %in% temp.df$treatment[!(temp.df$studyID %in% drops[i])] &
          all(connectcheck==TRUE))) {

      index <- !index
      temp.df <- ind.df[!(ind.df$studyID %in% drops[i] &
                            ind.df$treatment==comp[index+1]),]
    }

    ind.df <- temp.df
    index <- !index
  }
  return(ind.df)
}





#' Ensures indirect evidence can be estimated for comparisons identified
#' within inconsistency.loops
check.indirect.drops <- function(data=data, comp) {

  # Drop studies/comparisons that compare comps
  dropID <- vector()
  dropcomp <- vector()
  studies <- unique(data$studyID)
  for (study in seq_along(studies)) {
    subset <- data[data$studyID==studies[study],]
    if (all(comp %in% subset$treatment)) {
      if (subset$narm[1]<=2) {
        dropID <- append(dropID, subset$studyID[1])
      } else if (subset$narm[1]>2) {
        dropcomp <- append(dropcomp, subset$studyID[1])
      }
    }
  }

  # Drop studies
  data <- data[!(data$studyID %in% dropID),]

  # Drop comparisons from studies
  stoploop <- FALSE
  count <- 1
  while(stoploop==FALSE) {
    temp <- drop.comp(data, drops=dropcomp, comp=comp)
    temp.net <- MBNMA.network(temp)
    nt <- length(temp.net$treatments)
    if (nt==length(unique(data$treatment))) {
      g <- plot(temp.net, doseparam=1000)
      connectcheck <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                                       to=1)[
                                                         c(comp[1], comp[2])
                                                         ])
      if (all(connectcheck==TRUE)) {
        data <- temp
        stoploop <- TRUE
      }
    }
    #print("Restarting drop.comp")
    count <- count+1
    if (count==5) {
      break()
    }
  }
  if (count<5) {
    return(data)
  } else {
    return(NULL)
  }
}
