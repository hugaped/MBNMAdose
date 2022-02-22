# Functions for assessing inconsistency in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-30

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))

#' Node-splitting model for testing consistency at the treatment-level
#'
#' Splits contributions for a given set of treatment comparisons into direct and
#' indirect evidence. A discrepancy between the two suggests that the consistency
#' assumption required for NMA (and subsequently MBNMA) may violated.
#'
#' @param drop.discon A boolean object that indicates whether to drop treatments
#' that are disconnected at the treatment level. Default is `TRUE`. If set to `FALSE` then
#' this could lead to identification of nodesplit comparisons that are not connected
#' to the network reference treatment, or lead to errors in running the nodesplit models, though it
#' can be useful for error checking.
#' @param comparisons A matrix specifying the comparisons to be split (one row per comparison).
#' The matrix must have two columns indicating each treatment for each comparison. Values can
#' either be character (corresponding to the treatment names given in `network`) or
#' numeric (corresponding to treatment codes within the `network` - note that these
#' may change if `drop.discon = TRUE`).
#' @inheritParams mbnma.run
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(triptans)
#'
#' split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
#'   method="common")
#'
#'
#'
#' #### To perform nodesplit on selected comparisons ####
#'
#' # Check for closed loops of treatments with independent evidence sources
#' loops <- inconsistency.loops(network$data.ab)
#'
#' # This...
#' single.split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
#'              method="random", comparisons=rbind(c("sumatriptan_1", "almotriptan_1")))
#'
#' #...is the same as...
#' single.split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
#'              method="random", comparisons=rbind(c(6, 12)))
#'
#'
#' # Plot results
#' plot(split, plot.type="density") # Plot density plots of posterior densities
#' plot(split, plot.type="forest") # Plot forest plots of direct and indirect evidence
#'
#' # Print and summarise results
#' print(split)
#' summary(split) # Generate a data frame of summary results
#' }
#' @export
nma.nodesplit <- function(network, likelihood=NULL, link=NULL, method="common",
                            comparisons=NULL, drop.discon=TRUE,
                            ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::assertLogical(drop.discon, add=argcheck)
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

    agents <- unique(sapply(trt.labs, function(x) strsplit(x, "_")[[1]][1]))
    data.ab$agent <- factor(data.ab$agent, labels = agents)
  } else if (drop.discon==FALSE) {
    data.ab <- network$data.ab
    trt.labs <- network$treatments
    data.ab$agent <- factor(data.ab$agent, labels = network$agents)
  }

  # Identify closed loops of treatments
  if (is.null(comparisons)) {
    comparisons <- inconsistency.loops(data.ab)
  } else {
    comparisons <- check.nodesplit.comparisons(data.ab, network, comparisons, trt.labs)
  }


  ##### Run NMA #####
  nma.net <- suppressMessages(mbnma.network(data.ab))
  nma.jags <- nma.run(nma.net, method=method,
                      likelihood=likelihood, link=link,
                      warn.rhat=FALSE, drop.discon = FALSE, ...)

  nodesplit.result <- list()
  for (split in seq_along(comparisons[,1])) {

    comp <- as.numeric(comparisons[split,1:2])
    print(paste0("Calculating nodesplit for: ",
                 paste0(trt.labs[comp[2]], " vs ", trt.labs[comp[1]])))

    ##### Estimate NMA (overall) relative effect #####
    nma.res <- nma.jags$jagsresult$BUGSoutput$sims.list$d[,comp[2]] -
      nma.jags$jagsresult$BUGSoutput$sims.list$d[,comp[1]]

    ##### Estimate Indirect evidence #####

    ind.df <- data.ab

    # Drop studies/comparisons that compare comps
    dropID <- vector()
    dropcomp <- vector()
    studies <- unique(ind.df$studyID)
    for (study in seq_along(studies)) {
      subset <- ind.df[ind.df$studyID==studies[study],]
      if (all(comp %in% subset$treatment)) {
        if (subset$narm[1]<=2) {
          if (is.factor(studies)) {
            dropID <- append(dropID, as.character(subset$studyID[1]))
          } else {
            dropID <- append(dropID, subset$studyID[1])
          }
        } else if (subset$narm[1]>2) {
          if (is.factor(studies)) {
            dropcomp <- append(dropcomp, as.character(subset$studyID[1]))
          } else {
            dropcomp <- append(dropcomp, subset$studyID[1])
          }
        }
      }
    }

    # Drop studies
    ind.df <- ind.df[!(ind.df$studyID %in% dropID),]

    # Drop comparisons from studies
    ind.df <- suppressWarnings(drop.comp(ind.df, drops=dropcomp, comp=comp))

    # Run NMA of indirect evidence
    ind.net <- suppressMessages(mbnma.network(ind.df))
    ind.jags <- nma.run(ind.net, method=method, likelihood=likelihood, link=link,
                        warn.rhat=FALSE, drop.discon = FALSE, ...)
    ind.res <- ind.jags$jagsresult$BUGSoutput$sims.list$d[,comp[2]] -
      ind.jags$jagsresult$BUGSoutput$sims.list$d[,comp[1]]


    ##### Estimate Direct Evidence #####

    dir.net <- suppressMessages(change.netref(mbnma.network(data.ab), ref=comp[1]))

    # Run NMA of direct evidence (using unrelated mean effects)
    dir.jags <- nma.run(dir.net, method=method, likelihood=likelihood, link=link,
                        warn.rhat=FALSE, drop.discon=FALSE, UME=TRUE, ...)
    dir.res <- dir.jags$jagsresult$BUGSoutput$sims.matrix[
      ,colnames(dir.jags$jagsresult$BUGSoutput$sims.matrix)==paste0("d[", comp[2],",1]")
    ]


    ##### Generate plots/results #####

    # Overlaps
    overlap.mat <- list("direct"=dir.res, "indirect"=ind.res)
    overlap <- overlapping::overlap(overlap.mat, plot=FALSE)
    diff <- sum((dir.res-ind.res)>0) / length(dir.res)
    p.values <- overlap$OV

    # Quantiles
    quantile_dif <- stats::quantile(ind.res - dir.res, c(0.025, 0.5, 0.975))
    quantile_dir <- stats::quantile(dir.res, c(0.025, 0.5, 0.975))
    quantile_ind <- stats::quantile(ind.res, c(0.025, 0.5, 0.975))
    quantile_nma <- stats::quantile(nma.res, c(0.025, 0.5, 0.975))
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
      ggplot2::ggplot(data=plotdata, ggplot2::aes(x=plotdata$source, y=plotdata$med, ymin=plotdata$l95, ymax=plotdata$u95)) +
      ggplot2::geom_pointrange() +
      ggplot2::coord_flip() +
      ggplot2::xlab("") + ggplot2::ylab("Treatment effect (95% CrI)") + ggplot2::ggtitle(title) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=15),
                     axis.title = ggplot2::element_text(size=18),
                     title=ggplot2::element_text(size=18)) +
      ggplot2::theme(plot.margin=ggplot2::unit(c(1,1,1,1),"cm")) +
      ggplot2::theme_bw()

    # Density plots (with shaded area of overlap)
    molten <- data.frame(ind.res, dir.res)
    molten <- suppressMessages(reshape2::melt(molten))
    names(molten) <- c("Estimate", "value")
    linetypes <- c("solid", "dash")
    levels(molten$Estimate) <- c("Indirect", "Direct")

    dens <- ggplot2::ggplot(molten, ggplot2::aes(x=molten$value, linetype=molten$Estimate, fill=molten$Estimate)) +
      ggplot2::geom_density(alpha=0.2) +
      ggplot2::xlab(title) +
      ggplot2::ylab("Posterior density") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size=12)) +
      ggplot2::theme(axis.text = ggplot2::element_text(size=12),
                     axis.title = ggplot2::element_text(size=14)) +
      ggplot2::theme_bw()

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

  class(nodesplit.result) <- "nodesplit"

  return(nodesplit.result)
}






#' Identify comparisons in loops that fulfill criteria for node-splitting
#'
#' Identify comparisons informed by both direct and indirect evidence from
#' independent sources, which therefore fulfill the criteria for testing for
#' inconsistency via node-splitting.
#'
#' @param df A data frame containing variables `studyID` and `treatment` (as
#'   numeric codes) that indicate which treatments are used in which studies. If `checkindirect = TRUE`
#'   then variables `agent` and `dose` are also required.
#' @param checkindirect A boolean object to indicate whether or not to perform an additional
#'   check to ensure network remains connected even after dropping direct evidence on a comparison.
#'   Default is `TRUE` and should be kept as `TRUE` if working with dose-response data, though this requires
#'   further computational iterations to confirm. If set to `FALSE`, additional comparisons may be identified, though computation will be much more
#'   rapid.
#' @param incldr A boolean object indicating whether or not to allow for indirect evidence contributions via
#' the dose-response relationship. This can be used when node-splitting in dose-response MBNMA to allow
#' for a greater number of potential loops in which to check for consistency.
#'
#' @details Similar to `gemtc::mtc.nodesplit.comparisons()` but uses a fixed
#'   reference treatment and therefore identifies fewer loops in which to test for
#'   inconsistency. Heterogeneity can also be parameterised as inconsistency and
#'   so testing for inconsistency in additional loops whilst changing the
#'   reference treatment would also be identifying heterogeneity. Depends on
#'   \code{\link[igraph]{igraph}}.
#'
#' @return A data frame of comparisons that are informed by direct and indirect
#'   evidence from independent sources. Each row of the data frame is a
#'   different treatment comparison. Numerical codes in `t1` and `t2` correspond
#'   to treatment codes. `path` indicates the treatment codes that connect the
#'   shortest path of indirect evidence.
#'
#'   If `incldr=TRUE` then `path` may indicate `doseresp` for some comparisons.
#'   These are comparisons for which indirect evidence is only available via the
#'   dose-response relationship. The two numbers given after (e.g. `3 2`) indicate the
#'   number of doses available in the indirect evidence with which to estimate the
#'   dose-response function for the treatments in `t1` and `t2` respectively/
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' \donttest{
#' # Identify comparisons informed by direct and indirect evidence
#' #in triptans dataset
#' network <- mbnma.network(triptans)
#' inconsistency.loops(network$data.ab)
#'
#' # Include indirect evidence via dose-response relationship
#' inconsistency.loops(network$data.ab, incldr=TRUE)
#' }
#'
#'
#' # Do not perform additional connectivity check on data
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'             treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'             )
#' inconsistency.loops(data, checkindirect=FALSE)
#' @export
inconsistency.loops <- function(df, checkindirect=TRUE, incldr=FALSE)
{
  # Assert checks
  checkmate::assertDataFrame(df)
  checkmate::assertLogical(checkindirect)
  checkmate::assertLogical(incldr)

  treatments <- factor(unique(df$treatment))

  df <- df %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(design=list(as.numeric(df$treatment)))

  df <- df %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(narm=dplyr::n())


  comparisons <- ref.comparisons(df)

  if (incldr==TRUE) {
    # Add columns for agent
    lookup <- unique(df[,c("treatment", "agent")])
    comparisons$a1 <-lookup$agent[match(comparisons$t1, lookup$treatment)]
    comparisons$a2 <-lookup$agent[match(comparisons$t2, lookup$treatment)]
  }

  splits1 <- vector()
  splits2 <- vector()
  paths <- vector()
  loops <- vector()

  # Initiate progress bar
  pb <- utils::txtProgressBar(0, nrow(comparisons), style = 3)
  for (i in 1:nrow(comparisons)) {
    utils::setTxtProgressBar(pb, i)

    added <- FALSE
    drops <- comparisons[-i,]

    # Alternative graph create (non-gemtc)
    g <- igraph::graph.empty()
    g <- g + igraph::vertex(levels(treatments))

    ed <- t(matrix(c(drops[["t1"]], drops[["t2"]]), ncol = 2))
    edges <- igraph::edges(as.vector(ed), arrow.mode=0)
    g <- g + edges

    # Check whether there is still an indirect connection once direct evidence studies are removed
    if (as.logical(is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                                    comparisons[i,1], comparisons[i,2]))) == TRUE) {

      # Check if dropping 2-arm studies with both treatments and then either arm from multi-arm
      #would lead to disconnected network
      check <- 1
      if (checkindirect==TRUE) {
        check <- suppressMessages(suppressWarnings(
          check.indirect.drops(df, comp=c(as.numeric(comparisons[i,1]),
                                            as.numeric(comparisons[i,2])))
        ))
      }

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

        added <- TRUE
      }
    }

    # Check for additional indirect loops via the dose-response relationship
    if (added==FALSE & incldr==TRUE) {
      if (all(comparisons[i,c("a1", "a2")] %in% c(drops$a1, drops$a2))) {
        splits1 <- append(splits1, comparisons[["t1"]][i])
        splits2 <- append(splits2, comparisons[["t2"]][i])

        # Count number of doses available to estimate indirect
        lookup <- data.frame("a"=c(drops$a1, drops$a2), "t"=c(drops$t1, drops$t2))
        lookup <- lookup %>% dplyr::group_by(a) %>% dplyr::mutate(count= dplyr::n_distinct(t))

        param1 <- lookup$count[lookup$a %in% comparisons[i,"a1"]][1]
        param2 <- lookup$count[lookup$a %in% comparisons[i,"a2"]][1]

        # Add param if placebo is in dataset
        if (all(df$dose[df$treatment==1] == 0)) { # if placebo is in dataset
          if (1 %in% lookup$t) { # if placebo is in drops
            param1 <- param1 + 1
            param2 <- param2 + 1
          }

          # Assign placebo the max dose-response params
          maxparam <- max(lookup$count)
          if (comparisons[i,"a1"] == 1) {
            param1 <- maxparam
          }
          if (comparisons[i,"a2"] == 1) {
            param2 <- maxparam
          }
        }

        paths <- append(paths, paste("drparams", param1, param2, sep=" "))
        loops <- append(loops, paste(c(drops$t1, "dresp", drops$t2), collapse="->"))
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
#' @param df A data frame containing variables `studyID` and `treatment` (as
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
#'   `MBNMAtime`) `nr` will represent the number of time points in the
#'   dataset in which treatments `t1` and `t2` are compared.
#'
#' @noRd
#' @examples
#' df <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'             treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'             )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' ref.comparisons(df)
ref.comparisons <- function(df)
{
  # Assert checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(df, add=argcheck)
  checkmate::assertNames(names(df), must.include=c("studyID", "treatment"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  checkmate::assert(
    checkmate::checkFactor(df[["treatment"]]),
    checkmate::checkNumeric(df[["treatment"]])
  )

  df <- dplyr::arrange(df, df$studyID, df$treatment)

  t1 <- vector()
  t2 <- vector()
  for (i in seq_along(unique(df[["studyID"]]))) {
    subset <- subset(df, studyID==unique(df[["studyID"]])[i])
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
    dplyr::mutate(nr = dplyr::n())
  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)
}







#' Drop treatments from multi-arm (>2) studies for node-splitting
#'
#' Drops arms in a way which preserves connectivity and equally removes
#' data from each treatment in a nodesplit comparison (so as to maximise precision)
#'
#' @param ind.df A data frame in long format (one arm per row) from which to drop treatments
#' @param drops A vector of study identifiers from which to drop treatments
#' @param comp A numeric vector of length 2 that contains treatment codes corresponding to the comparison
#' for node-splitting
#' @param start Can take either `0` or `1` to indicate whether to drop the treatment
#' in `comp[1]` (`0`) or `comp[2]` (`1`)
#'
drop.comp <- function(ind.df, drops, comp, start=1) {
  index <- start

  for (i in seq_along(drops)) {

    switchloop <- FALSE
    temp.df <- ind.df[!(ind.df$studyID %in% drops[i] &
                          ind.df$treatment==comp[index+1]),]

    if (all(comp %in% temp.df$treatment)) {
      temp.net <- suppressMessages(plot.invisible(mbnma.network(temp.df), doseparam = 1000))

      connectcheck <- is.finite(igraph::shortest.paths(igraph::as.undirected(temp.net),
                                                       to=comp[index+1])[
                                                         c(comp[1], comp[2])
                                                         ])
    } else {
      connectcheck <- FALSE
    }



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
#'
#' Requires repeatedly compiling `mbnma.network` objects to identify whether
#' nodes remain connected by indirect evidence.
#'
#' @param df A data frame containing arm data (one arm per row)
#' @param comp The comparison in which the function will check that both treatments are connected by studies
#' in `df`
#' @noRd
check.indirect.drops <- function(df, comp) {

  # Drop studies/comparisons that compare comps
  dropID <- vector()
  dropcomp <- vector()
  studies <- unique(df$studyID)
  for (study in seq_along(studies)) {
    subset <- df[df$studyID==studies[study],]
    if (all(comp %in% subset$treatment)) {
      if (subset$narm[1]<=2) {
        dropID <- append(dropID, subset$studyID[1])
      } else if (subset$narm[1]>2) {
        dropcomp <- append(dropcomp, subset$studyID[1])
      }
    }
  }

  # Drop studies
  df <- df[!(df$studyID %in% dropID),]

  # Drop comparisons from studies
  stoploop <- FALSE
  count <- 1
  while(stoploop==FALSE) {
    temp <- drop.comp(df, drops=dropcomp, comp=comp)
    temp.net <- mbnma.network(temp)
    nt <- length(temp.net$treatments)
    if (nt==length(unique(df$treatment))) {
      g <- plot.invisible(temp.net, doseparam=1000)
      connectcheck <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                                       to=1)[
                                                         c(comp[1], comp[2])
                                                         ])
      if (all(connectcheck==TRUE)) {
        df <- temp
        stoploop <- TRUE
      }
    }

    count <- count+1
    if (count==5) {
      break()
    }
  }
  if (count<5) {
    return(df)
  } else {
    return(NULL)
  }
}







#' Node-splitting model for testing consistency at the treatment level using MBNMA
#'
#' Splits contributions for a given set of treatment comparisons into direct and
#' indirect evidence. A discrepancy between the two suggests that the consistency
#' assumption required for NMA and MBNMA may violated.
#'
#' @param comparisons A matrix specifying the comparisons to be split (one row per comparison).
#' The matrix must have two columns indicating each treatment for each comparison. Values can
#' either be character (corresponding to the treatment names given in `network`) or
#' numeric (corresponding to treatment codes within the `network` - note that these
#' may change if `drop.discon = TRUE`).
#' @param ... Arguments to be sent to `mbnma.run()`
#' @inheritParams mbnma.run
#' @inheritParams inconsistency.loops
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(triptans)
#'
#' split <- mbnma.nodesplit(network, fun=demax(), likelihood = "binomial", link="logit",
#'   method="common")
#'
#'
#' #### To perform nodesplit on selected comparisons ####
#'
#' # Check for closed loops of treatments with independent evidence sources
#' # Including indirect evidence via the dose-response relationship
#' loops <- inconsistency.loops(network$data.ab, incldr=TRUE)
#'
#' # This...
#' single.split <- mbnma.nodesplit(network, fun=dexp(), likelihood = "binomial", link="logit",
#'              method="random", comparisons=rbind(c("sumatriptan_1", "almotriptan_1")))
#'
#' #...is the same as...
#' single.split <- mbnma.nodesplit(network, fun=dexp(), likelihood = "binomial", link="logit",
#'              method="random", comparisons=rbind(c(6, 12)))
#'
#'
#' # Plot results
#' plot(split, plot.type="density") # Plot density plots of posterior densities
#' plot(split, txt_gp=forestplot::fpTxtGp(cex=0.5)) # Plot forest plots (with smaller label size)
#'
#' # Print and summarise results
#' print(split)
#' summary(split) # Generate a data frame of summary results
#' }
#' @export
mbnma.nodesplit <- function(network, fun=dloglin(),
                            method="common",
                            comparisons=NULL,
                            incldr=TRUE,
                            beta.1="rel", beta.2="rel", beta.3="rel", beta.4="rel",
                            user.fun=NULL,
                            ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check fun
  fun <- check.fun(fun=fun, network=network, beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                   user.fun=user.fun)
  if ("nonparam" %in% fun$name) {
    stop("Node-splitting cannot currently be performed for non-parametric models")
    # But this could be added
  }

  # Load data
  data.ab <- network$data.ab
  trt.labs <- network$treatments


  # Identify closed loops of treatments
  if (is.null(comparisons)) {
    comparisons <- inconsistency.loops(data.ab, incldr = incldr)
  } else {
    comparisons <- check.nodesplit.comparisons(data.ab, network, comparisons, trt.labs)
  }

  # Set default for pd="pv" unless otherwise specified for faster running
  args <- list(...)
  if (!"pd" %in% names(args)) {
    args[["pd"]] <- "pv"
  }


  ##### Run MBNMA #####
  mbnma.jags <- do.call(mbnma.run, c(args, list(network=network,
                                                method = method,
                                                fun = fun,
                                                warn.rhat=FALSE))
  )
  # mbnma.jags <- mbnma.run(network, method=method, fun=fun, warn.rhat=FALSE, pd=pd, ...)


  # Get comparison treatment names (in format for get.relative)
  compnames <- mbnma.jags$network$treatments[unique(c(comparisons[,1], comparisons[,2]))]

  temp <- sapply(compnames,
                 FUN=function(x) {
                   y <- strsplit(x, split="_")[[1]]
                   list(y[1], as.numeric(y[2]))
                 })

  comp.list <- list()
  for (i in 1:ncol(temp)) {
    if (!temp[,i][[1]] %in% names(comp.list)) {
      comp.list[[temp[,i][[1]]]] <- temp[,i][[2]]
    } else {
      comp.list[[temp[,i][[1]]]] <- append(comp.list[[temp[,i][[1]]]], temp[,i][[2]])
      comp.list[[temp[,i][[1]]]] <- sort(comp.list[[temp[,i][[1]]]])
    }
  }


  # Calculate relative effects for MBNMA
  mbnma.rel <- get.relative(mbnma.jags, treatments = comp.list)$relarray

  nodesplit.result <- list()
  for (split in seq_along(comparisons[,1])) {
    print(paste0("Comparison ", split,"/",nrow(comparisons)))

    comp <- as.numeric(comparisons[split,1:2])
    print(paste0("Calculating nodesplit for: ",
                 paste0(trt.labs[comp[2]], " vs ", trt.labs[comp[1]])))

    ##### Store MBNMA relative effects #####
    compnames <- mbnma.jags$network$treatments[c(comp[1], comp[2])]
    mbnma.res <- mbnma.rel[compnames[2], compnames[1], ]


    ####### Estimate Indirect and Direct in same model #########

    ind.net <- suppressMessages(change.netref(mbnma.jags$network, ref=comp[1]))

    ind.jags <- do.call(mbnma.run, c(args, list(network=ind.net,
                                                method = method,
                                                fun = fun,
                                                warn.rhat=FALSE,
                                                nodesplit=c(1, comp[2])))
    )

    # ind.jags <- mbnma.run(ind.net, method=method, fun=fun,
    #                       warn.rhat=FALSE, nodesplit=c(1, comp[2]), ...)

    # Get indirect
    ind.res <- get.relative(ind.jags, treatments = comp.list)$relarray[compnames[2],
                                                                       compnames[1],
                                                                       ]

    # Get direct
    dir.res <- ind.jags$BUGSoutput$sims.matrix[
      ,colnames(ind.jags$BUGSoutput$sims.matrix)=="direct"
      ]

    ##### Generate plots/results #####

    # Overlaps
    overlap.mat <- list("direct"=dir.res, "indirect"=ind.res)
    overlap <- overlapping::overlap(overlap.mat, plot=FALSE)
    p.values <- overlap$OV

    # Quantiles
    quantile_dif <- stats::quantile(ind.res - dir.res, c(0.025, 0.5, 0.975))
    quantile_dir <- stats::quantile(dir.res, c(0.025, 0.5, 0.975))
    quantile_ind <- stats::quantile(ind.res, c(0.025, 0.5, 0.975))
    quantile_mbnma <- stats::quantile(mbnma.res, c(0.025, 0.5, 0.975))
    quantiles <- list("difference" = quantile_dif, "mbnma"=quantile_mbnma,
                      "direct"=quantile_dir, "indirect"=quantile_ind)


    nodesplit <- list("comparison"= c(trt.labs[comp[2]], trt.labs[comp[1]]),
                      "direct"=dir.res, "indirect"=ind.res, "mbnma"=mbnma.res,
                      "overlap matrix"=overlap.mat,
                      "p.values"=p.values, "quantiles"=quantiles,
                      #"forest.plot"=gg, "density.plot"=dens,
                      "split.model"=ind.jags,
                      "mbnma.model"=mbnma.jags)

    nodesplit.result[[paste0(trt.labs[comp[2]], "v", trt.labs[comp[1]])]] <-
      nodesplit
  }

  class(nodesplit.result) <- "nodesplit"

  return(nodesplit.result)
}






#' Calculates relative effects between treatments in an MBNMA model
#'
#' @param mbnma An object of `class("mbnma")`
#' @param treatments A list whose elements each represent different treatments.
#' Treatment is defined as a combination of agent and dose. Only agents specified in
#' `mbnma` can be included. Each element in `treatments` is named corresponding to the
#' agent and contains a numeric vector of doses. Relative effects will be calculated between
#' all treatments specified in `treatments`. If `treatments` is left empty then the maximum
#' dose for all agents in `mbnma` will be used as the default.
#' @param eform Whether outputted results should be presented in their exponential form (e.g. for
#' models with log or logit link functions)
#' @inheritParams predict.mbnma
#'
#'
#' @return An array of `length(treatments) x length(treatments) x nsims`, where `nsims`
#' is the number of iterations monitored in `mbnma`. The array contains the individual
#' MCMC values for each relative effect calculated between all `treatments` on the link scale
#' specified in the `mbnma` model. The direction of effect is for the row-defined treatment
#' versus the column-defined
#' treatment.
#'
#'
#' @examples
#' # Using the osteoarthritis data
#' network <- mbnma.network(osteopain)
#'
#' expon <- mbnma.run(network, fun=dexp(), method="random")
#'
#' # Calculate relative effects between:
#' # Celebrex 100mg/d, Celebrex 200mg/d, Tramadol 100mg/d
#' rel.eff <- get.relative(expon, treatments=list("Celebrex"=c(100,200), "Tramadol"=100))
#'
#' @export
get.relative <- function(mbnma, treatments=list(), eform=FALSE, lim="cred") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertList(treatments, null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(lim, choices=c("cred", "pred"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Ensure prediction intervals are used where appropriate
  if (lim=="pred" & mbnma$model.arg$method=="random") {
    addsd <- TRUE

    if (!"sd" %in% mbnma$parameters.to.save) {
      stop(crayon::red("'sd' not included in parameters.to.save - cannot calculate prediction intervals"))
    }

  } else {
    addsd <- FALSE
  }

  # If treatments is not specified use the max dose of each agent in the dataset
  if (length(treatments)==0) {
    jagsdat <- mbnma$model.arg$jagsdata
    for (i in seq_along(mbnma$network$agents)) {
      treatments[[mbnma$network$agents[i]]] <- max(jagsdat$dose[jagsdat$agent==i], na.rm=TRUE)
    }
  }

  if (length(treatments)<2) {
    if (length(treatments[[1]])<2) {
      stop("`treatments` must have at least two elements to estimate relative\neffects between them")
    }

    if (!all(names(treatments) %in% mbnma$network$agents)) {
      stop("names(treatments) are not all in mbnma$network$agents")
    }
  }

  # Change `placebo` to dose=0 of an agent
  if ("Placebo" %in% names(treatments)) {
    treatments[[mbnma$network$agents[2]]] <- c(0, treatments[[mbnma$network$agents[2]]])
    treatments$Placebo <- NULL
  }

  # Identify dose-response function used in mbnma
  DR <- suppressMessages(
    write.dose.fun(fun=mbnma$model.arg$fun, effect="abs")[[1]])
  DR <- gsub("(^.+<-)(.+)", "\\2", DR)
  DR <- gsub("s\\.", "", DR) # Could remove if needed
  if (length(mbnma$model.arg$fun$name)>1) {
    DR <- DR[-1]
  }


  # Get dose-response parameter values
  betaparams <- get.model.vals(mbnma)


  trtnew <- treatments

  # Generate spline basis matrix if required
  splineopt <- c("rcs", "bs", "ns", "ls")
  fun <- mbnma$model.arg$fun

  # Get indices of non-placebo agents
  index <- match(names(trtnew), mbnma$network$agents[!mbnma$network$agents %in% "Placebo"])

  # If there are multiple DR functions
  if ("posvec" %in% names(fun)) {
    posvec <- fun$posvec

    # Remove 1st element if Placebo in network
    if ("Placebo" %in% mbnma$network$agents) {
      posvec <- posvec[-1]
    }
  } else {
    posvec <- rep(1, length(index))
  }
  for (i in seq_along(index)) {

    if (fun$name[posvec[index[i]]] %in% splineopt) {
      trtnew[[i]] <- genspline(trtnew[[i]],
                               spline = fun$name[posvec[index[i]]],
                               knots=fun$knots[[posvec[index[i]]]],
                               degree = fun$degree[posvec[index[i]]],
                               max.dose=max(mbnma$network$data.ab$dose[mbnma$network$data.ab$agent==which(names(trtnew)[i] == mbnma$network$agents)]))
    }
  }

  # Create list of treatments with doses
  trtlist <- list()
  for (i in seq_along(trtnew)) {
    if (is.matrix(trtnew[[i]])) { # Allows for splines
      for (k in 1:nrow(trtnew[[i]])) {
        trtlist[[length(trtlist)+1]] <- list(names(trtnew)[i], as.vector(trtnew[[i]][k,]))
      }
    } else {
      for (k in seq_along(trtnew[[i]])) {
        trtlist[[length(trtlist)+1]] <- list(names(trtnew)[i], trtnew[[i]][k])
      }
    }
  }


  # Generate array of relative effects between all treatments in network
  rel <- array(dim=c(length(trtlist), length(trtlist), mbnma$BUGSoutput$n.sims))
  for (i in 1:length(trtlist)) {
    for (k in 1:length(trtlist)) {

      agnum.i <- which(mbnma$network$agents %in% trtlist[[i]][[1]])
      agnum.k <- which(mbnma$network$agents %in% trtlist[[k]][[1]])

      # Account for lack of placebo
      if (!"Placebo_0" %in% mbnma$network$treatments) {
        agnum.i <- agnum.i+1
        agnum.k <- agnum.k+1
      }
      agnum <- c(agnum.i, agnum.k)

      # For agent-specific dose-response functions
      if (length(mbnma$model.arg$fun$name)>1) {
        posi <- mbnma$model.arg$fun$posvec[agnum.i]
        tempDR1 <- DR[posi]

        posk <- mbnma$model.arg$fun$posvec[agnum.k]
        tempDR2 <- DR[posk]

        pos <- c(posi, posk)
      } else {
        tempDR1 <- DR
        tempDR2 <- DR
      }

      DR1 <- gsub("agent\\[i,k\\]", paste0(",", agnum.i-1), tempDR1)
      DR1 <- gsub("dose\\[i,k\\]", trtlist[[i]][[2]][1], DR1)

      DR2 <- gsub("agent\\[i,k\\]", paste0(",", agnum.k-1), tempDR2)
      DR2 <- gsub("dose\\[i,k\\]", trtlist[[k]][[2]][1], DR2)

      # Replace spline
      for (m in 1:length(betaparams)) {
        DR1 <- gsub(paste0("spline\\[i,k,", m, "\\]"), trtlist[[i]][[2]][m], DR1)
        DR2 <- gsub(paste0("spline\\[i,k,", m, "\\]"), trtlist[[k]][[2]][m], DR2)
      }

      for (beta in seq_along(betaparams)) {
        assign(names(betaparams)[beta], betaparams[[beta]])

        if (length(mbnma$model.arg$fun$name)==1) {
          if (!is.matrix(betaparams[[beta]])) {
            DR1 <- gsub(paste0("(",names(betaparams)[beta], ")(\\[,[0-9]+\\])"), "\\1", DR1)
            DR2 <- gsub(paste0("(",names(betaparams)[beta], ")(\\[,[0-9]+\\])"), "\\1", DR2)
          }
        } else {

          if (is.matrix(betaparams[[beta]])) {

            DRcomb <- c(DR1, DR2)
            for (m in 1:2) {
              # Look for correct column index for each beta param
              veci <- mbnma$model.arg$fun$posvec[1:agnum[m]]
              veci <- table(veci)[names(table(veci))==pos[m]]

              if (names(veci)=="1") {
                # Check if placebo in dataset
                if (mbnma$network$agents[1]=="Placebo") {
                  veci <- veci - 1
                }
              }

              # Swap index in DR1 for veci
              DRcomb[m] <- gsub(paste0("(", names(betaparams)[beta], "\\[,)([0-9]+\\])"),
                          paste0("\\1", veci, "]"), DRcomb[m])
            }
            DR1 <- DRcomb[1]
            DR2 <- DRcomb[2]

          }
        }
      }

      chunk <- eval(parse(text=paste0(DR1, " - ", DR2)))

      if (length(rel)<=1) {stop("length(rel)<=1")}

      # Incorporate between-study SD
      if (addsd==TRUE) {
        mat <- matrix(nrow=length(chunk), ncol=2)
        mat[,1] <- chunk
        mat[,2] <- mbnma$BUGSoutput$sims.list[["sd"]]
        chunk <- apply(mat, MARGIN=1, FUN=function(x) stats::rnorm(1, x[1], x[2]))
      }

      rel[i,k,] <- chunk
      rel[k,i,] <- -chunk
      rel[i,i,] <- 0
      rel[k,k,] <- 0

    }
  }

  trtnames <- vector()
  for (i in seq_along(treatments)) {
    for (k in seq_along(treatments[[i]])) {
      # Change dose=0 back to `placebo` if needed
      if (treatments[[i]][k]==0) {
        trtnames <- append(trtnames, paste("Placebo", 0, sep="_"))
      } else {
        trtnames <- append(trtnames, paste(names(treatments)[i], treatments[[i]][k], sep="_"))
      }
    }
  }

  rownames(rel) <- trtnames
  colnames(rel) <- trtnames

  if (eform==FALSE) {
    outmat <- rel
  } else {
    outmat <- exp(rel)
  }


  ######### Summary matrixes ######

  xmat <- outmat

  meanmat <- matrix(nrow=nrow(xmat), ncol=ncol(xmat))
  semat <- meanmat
  medmat <- meanmat
  l95mat <- medmat
  u95mat <- medmat

  for (i in 1:nrow(xmat)) {
    for (k in 1:ncol(xmat)) {
      if (!is.na(xmat[i,k,1])) {
        meanmat[i,k] <- mean(xmat[i,k,])
        semat[i,k] <- stats::sd(xmat[i,k,])
        medmat[i,k] <- stats::median(xmat[i,k,])
        l95mat[i,k] <- stats::quantile(xmat[i,k,], probs = 0.025)
        u95mat[i,k] <- stats::quantile(xmat[i,k,], probs = 0.975)
      }
    }
  }

  sumlist <- list("mean"=meanmat, "se"=semat, "median"=medmat, "lower95"=l95mat, "upper95"=u95mat)

  for (i in seq_along(sumlist)) {
    dimnames(sumlist[[i]])[[1]] <- dimnames(xmat)[[1]]
    dimnames(sumlist[[i]])[[2]] <- dimnames(xmat)[[2]]
  }

  out <- list("relarray"=outmat)
  out <- c(out, sumlist)
  out$lim <- lim

  class(out) <- "relative.array"

  return(out)
}







#' Check validity of object supplied to `comparisons` in nodesplit
#'
#' @param data.ab A data frame with data for which to nodesplit
#' @param trt.labs A character vector of treatment labels corresponding to treatment codes in `network`
#' @inheritParams nma.nodesplit
#' @noRd
check.nodesplit.comparisons <- function(data.ab, network, comparisons, trt.labs) {
  if (!is.data.frame(comparisons) & !is.matrix(comparisons)) {
    stop("`comparisons` must be either a matrix or a data frame of comparisons on which to nodesplit")
  }
  if (is.data.frame(comparisons)) {
    if (all(c("t1", "t2") %in% names(comparisons))) {
      comparisons <- data.frame(comparisons$t1, comparisons$t2)
    }
    comparisons <- as.matrix(comparisons)
  }
  if (ncol(comparisons)!=2) {
    stop("`comparisons` must be a matrix of comparisons on which to split containing exactly two columns")
  }
  if (is.numeric(comparisons)) {
    # To ensure numbers correspond to original treatment numbers even if treatments are dropped from the network
    comparisons <- apply(comparisons, MARGIN=2, FUN=function(x) (network$treatments[x]))
  }
  if (is.character(comparisons)) {
    if (!all(comparisons %in% trt.labs)) {
      stop("Treatment names given in `comparisons` do not match those within `network` or they match treatments that have been dropped from the network due to being disconnected")
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
  fullcomp <- inconsistency.loops(data.ab, incldr=TRUE)
  match <- match(check, paste(fullcomp[,1], fullcomp[,2], sep="_"))
  if (any(is.na(match))) {
    out <- comparisons[is.na(match),]
    out <- matrix(unlist(lapply(out, FUN=function(x) {trt.labs[x]})), ncol=2)
    printout <- c()
    for (i in 1:nrow(out)) {
      printout <- paste(printout, paste(out[i,], collapse=" "), sep="\n")
    }
    stop(cat(paste0("\nThe following `comparisons` are not part of closed loops of treatments informed by direct and indirect evidence from independent sources:\n",
                    printout, "\n\n")))
  }

  return(comparisons)
}
