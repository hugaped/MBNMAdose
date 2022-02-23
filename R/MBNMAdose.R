#' MBNMAdose for dose-response Model-Based Network Meta-Analysis
#'
#' @description
#' `MBNMAdose` provides a collection of useful commands that allow users to run dose-response
#' Model-Based Network Meta-Analyses (MBNMA).
#'
#' @section Introduction:
#' `MBNMAdose` allows meta-analysis of studies that compare multiple doses of different agents in a way that can
#' account for the dose-response relationship.
#'
#' Whilst making use of all the available evidence in a statistically robust and biologically plausible framework,
#' this also can help connect networks at the agent level that may otherwise be disconnected at the dose/treatment
#' level, and help improve precision of estimates \insertCite{pedder2021}{MBNMAdose}. The modelling framework is based on synthesising relative effects
#' which avoids the necessity to adjust for baseline predictors, thereby making fewer assumptions than in typical
#' Model-Based Meta-Analysis.
#'
#' By modelling the dose-response, MBNMA avoids heterogeneity and inconsistency that can arise from "lumping" different
#' doses together (a technique sometimes done in Network Meta-Analysis). All models and analyses are implemented
#' in a Bayesian framework, following an extension of the standard NMA methodology presented by
#' \insertCite{lu2004;textual}{MBNMAdose} and are run in \insertCite{jags;textual}{MBNMAdose}. For full details of
#' dose-response MBNMA methodology see \insertCite{mawdsley2016;textual}{MBNMAdose}. Within this package we
#' refer to a **treatment** as a specific **dose** or a specific **agent**.
#'
#' @section Workflow:
#' Functions within `MBNMAdose` follow a clear pattern of use:
#'
#' 1. Load your data into the correct format using `mbnma.network()`
#' 2. Analyse your data using `mbnma.run()` with a wide range of dose-response functions
#' 3. Examine model results using forest plots and treatment rankings
#' 4. Check model fit and test for consistency using functions like `mbnma.nodesplit()`
#' 5. Use your model to predict responses using `predict()`
#'
#' At each of these stages there are a number of informative plots that can be generated to help understand
#' the data and to make decisions regarding model fitting.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # Generate an "mbnma.network" object that stores data in the correct format
#' network <- mbnma.network(triptans)
#'
#' # Generate a network plot at the dose/treatment level
#' plot(network, level="treatment")
#'
#' # Generate a network plot at the agent level
#' plot(network, level="agent", remove.loops=TRUE)
#'
#' \donttest{
#' # Perform "split" NMA to examine dose-response relationship
#' nma <- nma.run(network)
#' plot(nma)
#'
#' # Analyse data using mbnma.run() with an Emax dose-response function
#' # and common treatment effects
#' result <- mbnma.run(network, fun=demax(),
#'   method="common")
#'
#' # Generate forest plots for model results
#' plot(result)
#'
#' # Rank results and plot rankograms
#' ranks <- rank(result)
#' plot(ranks, params="emax")
#'
#' # Predict responses
#' pred <- predict(result, E0=0.2)
#'
#' # Plot predicted response with "split" NMA results displayed
#' plot(pred, overlay.split=TRUE)
#' }
#'
#' @keywords internal
"_PACKAGE"

