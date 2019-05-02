# Functions for assessing inconsistency in MBNMAdose
# Author: Hugo Pedder
# Date created: 2019-04-30


#' Node-splitting model for testing consistency
#'
#' Uses GeMTC
#'
#' @param ... Arguments to be sent to `gemtc::mtc.nodesplit`
MBNMA.nodesplit <- function(network, likelihood="binomial", link="logit", method="common",
                            drop.discon=TRUE, comparisons=NULL,
                            level="treatment", ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "MBNMA.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::assertLogical(drop.discon, add=argcheck)
  checkmate::assertDataFrame(comparisons, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  if (link=="probit") {
    stop("Node-splitting does not work with probit link function...  :-(")
  }

  if (method=="common") {method <- "fixed"}

  if (likelihood=="binomial") {likelihood <- "binom"}

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

  # Change names to allow for GeMTC
  varlist <- list(
    c("studyID", "study"),
    c("y", "mean"),
    c("se", "std.err"),
    c("r", "responders"),
    c("N", "sampleSize"),
    c("E", "exposure")
    )

  for (i in seq_along(varlist)) {
    if (varlist[[i]][1] %in% names(data.ab)) {
      names(data.ab)[names(data.ab)==varlist[[i]][1]] <- varlist[[i]][2]
    }
  }

  # Convert to GeMTC
  mtc <- gemtc::mtc.network(data.ab)



  # Identify closed loops of treatments
  if (is.null(comparisons)) {
    comparisons <- inconsistency.loops(network$data.ab)
  }

  # Nodesplit
  nodesplit <- gemtc::mtc.nodesplit(mtc, likelihood=likelihood, link=link,
                          comparisons=comparisons[1:2], linearModel=method, ...
                          )

  return(nodesplit)
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

  for (i in 1:nrow(comparisons)) {

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
