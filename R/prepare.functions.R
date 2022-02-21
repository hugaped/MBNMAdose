# Functions for manipulating/preparing MBNMA datasets
# Author: Hugo Pedder
# Date created: 2019-04-07

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))

#' Create an mbnma.network object
#'
#' Creates an object of `class("mbnma.network")`. Various MBNMA functions can subsequently be applied
#' to this object.
#'
#' @param data.ab A data frame of arm-level data in "long" format containing the columns:
#' * `studyID` Study identifiers
#' * `dose` Numeric data indicating the dose (must take positive values)
#' * `agent` Agent identifiers (can be numeric, factor or character)
#' * `y` Numeric data indicating the aggregate response for a continuous outcome. Required for
#' continuous data.
#' * `se` Numeric data indicating the standard error for a given observation. Required for
#' continuous data.
#' * `r` Numeric data indicating the number of responders within a study arm. Required for
#' binomial or poisson data.
#' * `n` Numeric data indicating the total number of participants within a study arm. Required for
#' binomial data or when modelling Standardised Mean Differences
#' * `E` Numeric data indicating the total exposure time for participants within a study arm. Required
#' for poisson data.
#' * `class` An optional column indicating a particular class code. Agents with the same identifier
#' must also have the same class code.
#' @param description Optional. Short description of the network.
#'
#' @details Agents/classes for arms that have dose = 0 will be relabelled as `"Placebo"`.
#' Missing values (`NA`) cannot be included in the dataset. Single arm studies cannot
#' be included.
#'
#' @return `mbnma.network()`: An object of `class("mbnma.network")` which is a list containing:
#' * `description` A short description of the network
#' * `data.ab` A data frame containing the arm-level network data (treatment identifiers will have
#' been recoded to a sequential numeric code)
#' * `studyID` A character vector with the IDs of included studies
#' * `agents` A character vector indicating the agent identifiers that correspond to the
#' new agent codes.
#' * `treatments` A character vector indicating the treatment identifiers that correspond
#' to the new treatment codes.
#' * `classes` A character vector indicating the class identifiers (if included in the original data)
#' that correspond to the new class codes.
#'
#' @examples
#' # Using the triptans headache dataset
#' print(triptans)
#'
#' # Define network
#' network <- mbnma.network(triptans, description="Example network")
#'
#' summary(network)
#' plot(network)
#'
#' @export
mbnma.network <- function(data.ab, description="Network") {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(description, len=1, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  mbnma.validate.data(data.ab)

  # Add indices for arms and narms and change agent/class codes
  index.data <- add_index(data.ab)

  network <- index.data
  network <- c(list("description" = description), network)

  class(network) <- "mbnma.network"
  return(network)
}



#' Validates that a dataset fulfills requirements for MBNMA
#'
#' @inheritParams mbnma.network
#' @param single.arm A boolean object to indicate whether to allow single arm studies in the dataset (`TRUE`)
#' or not (`FALSE`)
#'
#' @details Checks done within the validation:
#' * Checks data.ab has required column names
#' * Checks there are no NAs
#' * Checks that all SEs are >0 (if variables are included in dataset)
#' * Checks that all doses are >=0
#' * Checks that all r and n are positive (if variables are included in dataset)
#' * Checks that all y, se, r, n and E are numeric
#' * Checks that class codes are consistent within each agent
#' * Checks that agent/class names do not contain restricted characters
#' * Checks that studies have at least two arms (if `single.arm = FALSE`)
#' * Checks that each study includes at least two treatments
#' * Checks that agent names do not inlcude underscores
#'
#' @return An error if checks are not passed. Runs silently if checks are passed
mbnma.validate.data <- function(data.ab, single.arm=FALSE) {

  varnames <- c("studyID", "agent", "dose")
  var_norm <- c("y", "se")
  var_bin <- c("r", "n")
  var_pois <- c("r", "E")
  data.ab <- dplyr::arrange(data.ab, data.ab$studyID, data.ab$agent, data.ab$dose)

  # Check data.ab has required column names
  msg <- "Required variable names are: 'studyID', 'agent', 'dose' and either `y` and `se` for data with a normal likelihood, `r` and `n` for data with a binomial likelihood, or `r` and `E` for data with a poisson likelihood"
  if (all(varnames %in% names(data.ab))==FALSE) {
    if ("time" %in% names(data.ab)) {
      message(paste(
        "`time` is included as a variable in the dataset but required variables for dose-response",
        "MBNMA are not. Are you trying to run time-course MBNMA?",
        "If so use the MBNMAtime package rather than MBNMAdose.",
        sep="\n"))
    }
    stop(msg)
  }
  if (all(var_norm %in% names(data.ab))==FALSE &
      all(var_bin %in% names(data.ab))==FALSE &
      all(var_pois %in% names(data.ab))==FALSE) {
    stop(msg)
  }

  # Check there are no NAs
  na.vars <- vector()
  for (i in seq_along(varnames)) {
    if (anyNA(data.ab[[varnames[i]]])) {
      na.vars <- append(na.vars, varnames[i])
    }
  }

  numeric.error <- vector()
  if (!is.numeric(data.ab$dose)) {
    numeric.error <- append(numeric.error, "dose")
  }

  # Check for NA values in data variables
  for (i in 1:2) {
    if (all(var_norm %in% names(data.ab))) {
      if (anyNA(data.ab[[var_norm[i]]])) {
        na.vars <- append(na.vars, var_norm[i])
      }
      if (!is.numeric(data.ab[[var_norm[i]]])) {
        numeric.error <- append(numeric.error, var_norm[i])
      }
    }
    if (all(var_bin %in% names(data.ab))) {
      if (anyNA(data.ab[[var_bin[i]]])) {
        na.vars <- append(na.vars, var_bin[i])
      }
      if (!is.numeric(data.ab[[var_bin[i]]])) {
        numeric.error <- append(numeric.error, var_bin[i])
      }
    }
    if (all(var_pois %in% names(data.ab))) {
      if (anyNA(data.ab[[var_pois[i]]])) {
        na.vars <- append(na.vars, var_pois[i])
      }
      if (!is.numeric(data.ab[[var_pois[i]]])) {
        numeric.error <- append(numeric.error, var_pois[i])
      }
    }
  }

  if (length(na.vars)>0) {
    stop(paste0("NA values in:\n", paste(na.vars, collapse="\n")))
  }

  if (length(numeric.error)>0) {
    stop(paste0("Data must be numeric for:\n", paste(numeric.error, collapse="\n")))
  }


  # Check that all SEs are >0
  var_tmp <- c("se", var_bin, "E")
  for (i in seq_along(var_tmp)) {
    if (var_tmp[i] %in% names(data.ab)) {
      if (!all(data.ab[[var_tmp[i]]]>=0)) {
        stop(paste("All values for", var_tmp[i], "must be >=0", sep=" "))
      }
    }
  }

  # Check that all doses are positive
  if (!all(data.ab$dose>=0)) {
    stop("All values for `dose` must be >=0")
  }

  # Check that if numeric, agent codes are >0
  if (is.numeric(data.ab$agent)) {
    if (!all(data.ab$agent>0)) {
      stop("Agent codes in dataset must be numbered sequentially from 1")
    }
  }

  # Check that there are no underscores in agent names
  if (any(grepl("_", data.ab$agent))) {
    stop("Underscores (_) cannot be used in agent names")
  }


  # Generate narms index for checking if studies are only single-arm
  if (single.arm==FALSE) {
    data.ab <- data.ab %>%
      dplyr::group_by(studyID) %>%
      dplyr::mutate(narms = dplyr::n())
  }

  singlearm.studyID <- vector()
  singletreat.studyID <- vector()
  for (i in seq_along(unique(data.ab$studyID))) {
    subset <- data.ab[data.ab$studyID==unique(data.ab$studyID)[i],]

    # Check that no studies are single arm
    if (all(subset$narms<2)) {
      singlearm.studyID <- append(singlearm.studyID, as.character(subset$studyID[1]))
    }

    # Check that no studies include only the same treatment
    checkdoses <- paste(subset$agent, subset$dose, sep="_")
    if (!(length(unique(checkdoses))>1)) {
      singletreat.studyID <- append(singletreat.studyID, as.character(subset$studyID[1]))
    }
  }

  if (length(singlearm.studyID) >0) {
    stop(paste0("The following studies do not contain more than a single study arm:\n",
                paste(unique(singlearm.studyID), collapse="\n")))
  }
  if (length(singletreat.studyID) >0) {
    stop(paste0("The following studies only include comparison(s) of the same agent at the same dose:\n",
                paste(unique(singletreat.studyID), collapse="\n")))
  }


  # Check that agent/class names do not contain forbidden characters
  agents <- unique(data.ab$agent)
  if (any(grepl("_", agents))) {
    stop("Agent names cannot contain `_` (underscore) character")
  }
  if ("class" %in% names(data.ab)) {
    classes <- unique(data.ab$class)
    if (any(grepl("_", classes))) {
      stop("Class names cannot contain `_` (underscore) character")
    }
  }


  # Check that class codes are consistent within each agent
  if ("class" %in% names(data.ab)) {
    class.mismatch <- vector()
    for (i in seq_along(unique(data.ab$agent))) {
      match <- data.ab$class[data.ab$agent==data.ab$agent[i]]
      if (length(unique(match)) > 1) {
        class.mismatch <- append(class.mismatch, data.ab$agent[i])
      }
    }
    if (length(unique(class.mismatch))>0) {
      stop(paste0("Class codes are different within the same agent for the following treatments:\n",
                  paste(unique(class.mismatch), collapse="\n")))
    }

    # Check that if numeric, class codes are >0
    if (is.numeric(data.ab$class)) {
      if (!all(data.ab$class>0)) {
        stop("Class codes in dataset must be numbered sequentially from 1")
      }
    }
  }

}




#' Add arm indices and agent identifiers to a dataset
#'
#' Adds arm indices (`arms`, `narms`) to a dataset and adds numeric identifiers for
#' agent and class (if included in the data).
#'
#' @inheritParams mbnma.network
#' @param agents A character string of agent names used to force a particular agent ordering.
#' Default is `NULL`, which automatically orders `Placebo` (dose=0) as agent `1` and then
#' subsequent agents by the order given in `data.ab`
#' @param treatments A character string of treatment names used to force a particular treatment ordering.
#' Default is `NULL`, which automatically orders `Placebo` (dose=0) as treatment `1` and then
#' subsequent treatments by the order of agents and doses (smallest to highest) given in `data.ab`
#'
#' @return A data frame similar to `data.ab` but with additional columns:
#' * `arm` Arm identifiers coded for each study
#' * `narm` The total number of arms in each study
#'
#' If `agent` or `class` are non-numeric or non-sequential (i.e. with missing numeric codes),
#' agents/classes in the returned data frame will be numbered and recoded to enforce sequential
#' numbering (a warning will be shown stating this).
#'
add_index <- function(data.ab, agents=NULL, treatments=NULL) {

  # Run Checks
  checkmate::assertDataFrame(data.ab)

  if ("agent" %in% names(data.ab)) {
    if (is.null(agents)) {
      recoded <- recode.agent(data.ab, level = "agent")
      agents <- recoded[["lvlnames"]]
      data.ab <- recoded[["data.ab"]]

      # Generate treatment labels
      treatments.df <- recoded[["data.ab"]]
      treatments.df$agent.fac <- factor(treatments.df$agent, labels=agents)
      treatments.df <- dplyr::arrange(treatments.df, treatments.df$agent.fac, treatments.df$dose)
      treatments <- unique(paste(treatments.df$agent.fac, treatments.df$dose, sep="_"))

      # Generate treatment variable
      data.ab$treatment <- as.numeric(factor(paste(data.ab$agent,
                                                   data.ab$dose,
                                                   sep="_"),
                                             labels=treatments,
                                             levels=unique(paste(treatments.df$agent,
                                                                 treatments.df$dose,
                                                                 sep="_"))
      ))

      data.ab <- dplyr::arrange(data.ab, studyID, agent, dose)

    }
  }


  #### Add indices ####

  # Do not run this function with pylr loaded!!
  data.ab <- data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(arm = sequence(dplyr::n()))

  data.ab <- data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(narm=dplyr::n())


  # Reorder columns in data.ab
  ord <- c("agent", "dose", "treatment", "class", "narm", "arm", "y", "se", "r", "E", "n")
  newdat <- data.frame("studyID"=data.ab$studyID)
  for (i in seq_along(ord)) {
    if (ord[i] %in% names(data.ab)) {
      newdat <- cbind(newdat, data.ab[,which(names(data.ab)==ord[i])])
    }
  }
  olddat <- data.ab[,!(names(data.ab) %in% c("studyID", ord))]
  newdat <- cbind(newdat, olddat)

  newdat <- dplyr::arrange(newdat, dplyr::desc(newdat$narm), newdat$studyID, newdat$arm)

  output <- list("data.ab"=newdat,
                 "studyID"=as.character(unique(newdat$studyID)))

  if ("agent" %in% names(data.ab)) {
    output[["agents"]] <- agents
    output[["treatments"]] <- treatments
  }


  # Store class labels and recode (if they exist in data.ab)
  if ("class" %in% names(newdat)) {

    recoded <- recode.agent(newdat, level = "class")
    classes <- recoded[["lvlnames"]]

    # Generate class key
    classdata <- recoded$data.ab[, names(recoded$data.ab) %in% c("agent", "class")]
    classkey <- unique(classdata)
    classkey$agent <- factor(classkey$agent, labels=agents)
    classkey$class <- factor(classkey$class, labels=classes)

    output[["data.ab"]] <- recoded[["data.ab"]]
    output[["classes"]] <- classes
    output[["classkey"]] <- classkey
  }

  return(output)
}



#' Assigns agent or class variables numeric identifiers
#'
#' @param level Can take either `"agent"` or `"class"`
#' @inheritParams add_index
#'
#' @details Also relabels the agent for any arms in which dose = 0 to "Placebo_0"
#'
#' @return A list containing a data frame with recoded agent/class identifiers and
#'   a character vector of original agent/class names
#'
recode.agent <- function(data.ab, level="agent") {
  # Run Checks
  checkmate::assertDataFrame(data.ab)
  checkmate::assertChoice(level, choices = c("agent", "class"))

  # Check that there are no NA values
  if (any(is.na(data.ab[[level]]))) {
    stop(paste0("NA values not allowed for ", level))
  }

  print.msg <- FALSE

  lvls <- as.character(sort(unique(data.ab[[level]])))

  # Convert agent/class to numeric codes
  if (is.factor(data.ab[[level]])) {
    data.ab[[level]] <- as.numeric(data.ab[[level]])
  } else if (is.character(data.ab[[level]])) {
    data.ab[[level]] <- as.numeric(factor(data.ab[[level]], labels=lvls))
  }

  agent.seq <- sort(unique(data.ab[[level]]))

  for (i in seq_along(agent.seq)) {

    # If all doses of a particular agent/class = 0 then recode to "Placebo"
    allzero <- FALSE
    if (all(is.na(data.ab$dose[data.ab[[level]]==agent.seq[i]]))) {
      allzero <- TRUE
    } else if (all(data.ab$dose[data.ab[[level]]==agent.seq[i]]==0)) {
      allzero <- TRUE
    }
    if (allzero==TRUE) {
      if (lvls[1]!="Placebo") {
        # Swap current lvls for "Placebo"
        lvls <- c("Placebo", lvls[-agent.seq[i]])
      }
      data.ab[[level]][data.ab[[level]]==agent.seq[i]] <- 0
      print.msg <- TRUE
    }

    # Else if agent/class contains any dose=0, convert those doses to "Placebo"
    anyzero <- FALSE
    if (any(is.na(data.ab$dose[data.ab[[level]]==agent.seq[i]]))) {
      anyzero <- TRUE
    } else if (any(data.ab$dose[data.ab[[level]]==agent.seq[i]]==0)) {
      anyzero <- TRUE
    }
    if (anyzero==TRUE) {
      if (lvls[1]!="Placebo") {
        # Add "Placebo"
        lvls <- c("Placebo", lvls)
      }
      data.ab[[level]][data.ab[[level]]==agent.seq[i] & data.ab$dose==0] <- 0
      print.msg <- TRUE
    }
  }

  # Messages
  if (print.msg==TRUE) {
    message(paste0("Values for `", level, "` with dose = 0 have been recoded to `Placebo`"))
  }
  if (!identical(1:max(data.ab[[level]]), sort(unique(data.ab[[level]]))[-1])) {
    message(paste0(level, " is being recoded to enforce sequential numbering and allow inclusion of `Placebo`"))
  }

  # Reorder by number sequentially (meaning that "Placebo" now is 1)
  data.ab[[level]] <- as.numeric(factor(data.ab[[level]], labels=lvls))

  return(list("data.ab"=data.ab, "lvlnames"=lvls))
}








#' Prepares data for JAGS
#'
#' Converts MBNMA data frame to a list for use in JAGS model
#'
#' @inheritParams mbnma.run
#' @inheritParams mbnma.network
#' @param class A boolean object indicating whether or not `data.ab` contains
#'   information on different classes of treatments
#' @param level Can take either `"agent"` to indicate that data should be at the agent-
#'   level (for MBNMA) or `"treatment"` to indicate that data should be at the treatment-
#'   level (for NMA)
#' @param nodesplit A numeric vector of length 2 containing treatment codes on which to perform
#'   an MBNMA nodesplit (see \code{\link{mbnma.nodesplit}}).
#'
#' @return A named list of numbers, vector, matrices and arrays to be sent to
#'   JAGS. List elements are:
#'   * If `likelihood="normal"`:
#'     - `y` An array of mean responses for each arm within each study
#'     - `se` An array of standard errors for each arm within each study
#'   * If `likelihood="binomial"`:
#'     - `r` An array of the number of responses/count for each each arm within each study
#'     - `n` An array of the number of participants for each arm within each study
#'   * If `likelihood="poisson"`:
#'     - `r` An array of the number of responses/count for each each arm within each study
#'     - `E` An array of the total exposure time for each arm within each study
#'   * `dose` A matrix of doses for each arm within each study (if `level="agent"`)
#'   * `narm` A numeric vector with the number of arms per study
#'   * `NS` The total number of studies in the dataset
#'   * `Nagent` The total number of agents in the dataset (if `level="agent"`)
#'   * `agent` A matrix of agent codes within each study (if `level="agent"`)
#'   * `NT` The total number of treatment in the dataset (if `level="treatment"`)
#'   * `treatment` A matrix of treatment codes within each study (if `level="treatment"`)
#'   * `Nclass` Optional. The total number of classes in the dataset
#'   * `class` Optional. A matrix of class codes within each study
#'   * `classkey` Optional. A vector of class codes that correspond to agent codes.
#'   Same length as the number of agent codes.
#'   * `split.ind` Optional. A matrix indicating whether a specific arm contributes evidence
#'    to a nodesplit comparison.
#'
#' @examples
#' # Using the triptans headache dataset
#' network <- mbnma.network(triptans)
#' jagsdat <- getjagsdata(network$data.ab, likelihood="binomial", link="logit")
#'
#'
#' # Get JAGS data with class
#' netclass <- mbnma.network(osteopain)
#' jagsdat <- getjagsdata(netclass$data.ab, class=TRUE)
#'
#'
#' # Get JAGS data at the treatment level for split Network Meta-Analysis
#' network <- mbnma.network(triptans)
#' jagsdat <- getjagsdata(network$data.ab, level="treatment")
#'
#' @export
getjagsdata <- function(data.ab, class=FALSE,
                        likelihood=check.likelink(data.ab)$likelihood,
                        link=check.likelink(data.ab)$link,
                        level="agent", fun=NULL, nodesplit=NULL) {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(level, choices=c("agent", "treatment"), null.ok=FALSE, add=argcheck)
  checkmate::assertClass(fun, "dosefun", null.ok=TRUE, add=argcheck)
  checkmate::assertNumeric(nodesplit, len=2, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  df <- data.ab

  # Create vector of variable names for which data will be required
  varnames <- c("studyID", "arm", "narm")

  if (level=="agent") {
    varnames <- append(varnames, c("dose", "agent"))
  } else if (level=="treatment") {
    varnames <- append(varnames, c("treatment"))
  }

  if (class==TRUE) {
    varnames <- append(varnames, "class")
  }

  # Add variables depending on likelihood
  if (likelihood == "binomial") {
    datavars <- c("r", "n")
  } else if (likelihood == "poisson") {
    datavars <- c("r", "E")
  } else if (likelihood=="normal") {
    datavars <- c("y", "se")
  } else {
    stop("`likelihood` can be either `binomial`, `poisson`, or `normal`")
  }
  if (link=="smd") {
    datavars <- append(datavars, "n")
  }
  varnames <- append(varnames, datavars)

  # Check required variables are in df
  if (!all(varnames %in% names(df))) {
    msg <- paste0("Variables are missing from dataset:\n",
                  paste(varnames[!(varnames %in% names(df))], collapse="\n"))
    stop(msg)
  }

  df <- dplyr::arrange(df, dplyr::desc(df$narm), df$studyID, df$arm)

  df$studynam <- df$studyID
  df <- transform(df, studyID=as.numeric(factor(studyID, levels=as.character(unique(df$studyID)))))

  # Create a separate object for each datavars
  for (i in seq_along(datavars)) {
    assign(datavars[i], array(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                              dim=c(max(as.numeric(df$studyID)),
                                    max(df$arm)
                              ))
           )
  }

  narm <- vector()
  NS <- max(as.numeric(df$studyID))

  # Generate list in which to store individual data variables
  datalist <- list(narm=narm, NS=NS, studyID=vector())
  for (i in seq_along(datavars)) {
    datalist[[datavars[i]]] <- get(datavars[i])
  }

  # datalist <- list(get(datavars[1]), get(datavars[2]),
  #                  narm=narm, NS=NS, studyID=vector())
  # names(datalist)[1:2] <- datavars

  if (level=="agent") {
    datalist[["Nagent"]] <- max(df$agent)
    datalist[["agent"]] <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                                  nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
                                  )
    datalist[["dose"]] <- datalist[["agent"]]


    # Generate empty dmult matrix
    if (length(fun$name)>1) {
      datalist[["f"]] <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                                nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
      )
    }

    # Generate empty spline matrix
    splineopt <- c("rcs", "ns", "bs", "ls")
    if (any(splineopt %in% fun$name)) {

      doses <- df[, colnames(df) %in% c("agent", "dose")]
      doses <- unique(dplyr::arrange(doses, agent, dose))

      # If there are multiple spline functions
      if (dplyr::n_distinct(fun$name[fun$name %in% splineopt])>1) {
        doses$splinefun <- fun$name[fun$posvec][doses$agent]
      } else {
        doses$splinefun <- unique(fun$name[fun$name %in% splineopt])
      }

      # Remove non-spline functions
      doses <- subset(doses, splinefun %in% splineopt)

      # If there are multiple spline degrees
      degcheck <- fun$degree[fun$posvec][fun$name[fun$posvec] %in% splineopt]
      if (dplyr::n_distinct(degcheck)>1) {
        doses$degree <- fun$degree[fun$posvec][doses$agent]
      } else {
        doses$degree <- unique(fun$degree[fun$name %in% splineopt])
      }

      # If there are multiple spline knots
      #if (dplyr::n_distinct(fun$knots[[fun$name %in% splineopt]])>1) {
      if (length(unique(fun$knots))>1) {
        nknot <- max(unlist(lapply(fun$knots, length)), na.rm = TRUE)

        # If this dose not work unhash section below and hash this instead
        knotposvec <- fun$knots[fun$posvec]

        temp <- matrix(nrow=nrow(doses), ncol=nknot)
        for (r in 1:nrow(temp)) {
          temp[r,1:length(knotposvec[[doses$agent[r]]])] <- knotposvec[[doses$agent[r]]]
        }
        doses$knots <- temp

        # if (nknot>1) {
        #   temp <- fun$knots[fun$posvec]
        #
        #   doses$knots <- matrix(nrow=nrow(doses), ncol=nknot)
        #   print(temp)
        #   for (r in seq_along(doses$agent)) {
        #     doses$knots[r] <- temp[[doses$agent[r]]]
        #   }
        #
        # } else {
        #   doses$knots <- fun$knots[[fun$posvec]][doses$agent]
        # }

      } else {
        nknot <- length(unique(fun$knots[[fun$name %in% splineopt]]))

        # If this dose not work unhash section below and hash this instead
        temp <- matrix(rep(unique(fun$knots[[fun$name %in% splineopt]]), nrow(doses)),
                       byrow = TRUE, nrow=nrow(doses))
        doses$knots <- temp
        # if (nknot>1) {
        #   temp <- matrix(rep(unique(fun$knots[[fun$name %in% splineopt]]), nrow(doses)),
        #                  byrow = TRUE, nrow=nrow(doses))
        #   doses$knots <- temp
        # } else {
        #   doses$knots <- unique(fun$knots[[fun$name %in% splineopt]])
        # }
      }

      # Fill non-spline splinefun with any spline function
      doses$splinefun[!doses$splinefun %in% splineopt] <- doses$splinefun[doses$splinefun %in% splineopt][1]


      # Fill NAs with 1 for non-spline functions
      # doses$degree[is.na(doses$degree)] <- 1
      # doses$knots[is.na(doses$knots)] <- 1


      # Generate spline basis matrix
      # Run for placebo separately to avoid matrix dimension problems
      if (0 %in% doses$dose[doses$agent==1] & length(doses$dose[doses$agent==1])==1) {

        dosespline <- doses[doses$agent!=1,]

        uniag <- unique(dosespline$agent)
        splinemat <- matrix(nrow=nrow(dosespline), ncol=100)
        for (i in seq_along(uniag)) {
          sub <- subset(dosespline, agent==uniag[i])

          subspline <- genspline(sub$dose,
                                 spline=unique(sub$splinefun),
                                 #knots=as.vector(unique(sub$knots)),
                                 knots=unique(sub$knots),
                                 degree=unique(sub$degree))

          splinemat[which(dosespline$agent==uniag[i]),1:ncol(subspline)] <- subspline
        }
        # Remove excess columns
        splinemat <- splinemat[,1:(min(which(apply(splinemat, MARGIN=2, FUN=function(x) all(is.na(x)))))-1)]

        dosespline$spline <- splinemat

        # dosespline <- doses[doses$agent!=1,] %>% dplyr::group_by(agent) %>%
        #   dplyr::mutate(spline=genspline(dose, spline=splinefun, knots=knots, degree=degree))

        matsize <- ncol(dosespline$spline)

        pldose <- doses[doses$agent==1,] %>% dplyr::mutate(spline=matrix(rep(0,matsize), ncol=matsize))
        dosespline <- rbind(dosespline, pldose)

      } else {
        dosespline <- doses

        uniag <- unique(dosespline$agent)
        splinemat <- matrix(nrow=nrow(dosespline), ncol=100)
        for (i in seq_along(uniag)) {
          sub <- subset(dosespline, agent==uniag[i])

          subspline <- genspline(sub$dose,
                                  spline=unique(sub$splinefun),
                                  knots=unique(sub$knots),
                                  degree=unique(sub$degree))

          splinemat[which(dosespline$agent==uniag[i]),1:ncol(subspline)] <- subspline
        }
        # Remove excess columns
        splinemat <- splinemat[,1:(min(which(apply(splinemat, MARGIN=2, FUN=function(x) all(is.na(x)))))-1)]

        dosespline$spline <- splinemat

        # dosespline <- doses %>% dplyr::group_by(agent) %>%
        #   dplyr::mutate(spline=genspline(dose,
        #                                              spline=unique(splinefun),
        #                                              knots=unique(knots),
        #                                              degree=unique(degree)))

        matsize <- ncol(dosespline$spline)
      }

      df <- suppressMessages(dplyr::left_join(df, dosespline))

      matsize <- ncol(dosespline$spline)

      datalist[["spline"]] <- array(dim=c(nrow(datalist[["dose"]]),
                                          ncol(datalist[["dose"]]),
                                          matsize))

    }

  } else if (level=="treatment") {
    datalist[["NT"]] <- max(df$treatment)

    datalist[["treatment"]] <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                                      nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
                                      )
  }


  # Add empty list element for class data if classes are to be modelled
  if (class==TRUE) {
    Nclass <- max(df$class)

    codes <- data.frame(df$agent, df$class)
    codes <- dplyr::arrange(codes, codes$df.agent)
    classcode <- unique(codes)$df.class

    datalist[["Nclass"]] <- Nclass
    datalist[["class"]] <- classcode
  }

  # Add empty matrix indicating which data points contribute to direct/indirect in node-split model
  if (!is.null(nodesplit)) {
    datalist[["split.ind"]] <- datalist[["agent"]]
  }


  # Add data to datalist elements
  for (i in 1:max(as.numeric(df$studyID))) {
    for (k in 1:max(df$arm[df$studyID==i])) {
      datalist[["studyID"]] <- append(datalist[["studyID"]], df$studynam[as.numeric(df$studyID)==i &
                                                       df$arm==k])
      for (m in seq_along(datavars)) {
        datalist[[datavars[m]]][i,k] <- df[[datavars[m]]][as.numeric(df$studyID)==i &
                              df$arm==k]
      }

      if (level=="agent") {
        datalist[["agent"]][i,k] <- max(df$agent[as.numeric(df$studyID)==i &
                                                   df$arm==k])
        datalist[["dose"]][i,k] <- max(df$dose[as.numeric(df$studyID)==i &
                                                 df$arm==k])

        # Add spline matrix
        if (any(c("rcs", "ns", "bs", "ls") %in% fun$name)) {
          datalist[["spline"]][i,k,] <- df[as.numeric(df$studyID)==i &
                                             df$arm==k,
                                           grepl("spline$", colnames(df))]
        }
      } else if (level=="treatment") {
        datalist[["treatment"]][i,k] <- max(df$treatment[as.numeric(df$studyID)==i &
                                                   df$arm==k])
      }

      if (!is.null(nodesplit)) {
        if (df$treatment[as.numeric(df$studyID)==i & df$arm==1] == nodesplit[1] &
            df$treatment[as.numeric(df$studyID)==i & df$arm==k] == nodesplit[2]
        ) {
          datalist[["split.ind"]][i,k] <- 1
        } else {
          datalist[["split.ind"]][i,k] <- 0
        }
      }

    }

    datalist[["narm"]] <- append(datalist[["narm"]], max(df$arm[as.numeric(df$studyID)==i]))
  }

  # Add maxdose for nonparametric dose-response functions
  if ("nonparam" %in% fun$name) {
    data.ab <- data.ab %>% dplyr::group_by(agent) %>%
      dplyr::mutate(maxdose=max(dose, na.rm=TRUE)) %>%
      dplyr::slice_head(temp)

    datalist[["maxdose"]] <- data.ab$maxdose
  }

  if ("spline" %in% names(datalist)) {
    datalist[["spline"]][is.na(datalist[["spline"]])] <- 0
  }

  # Add agent-specific index matrix f
  if (length(fun$name)>1) {
    datalist[["f"]] <- matrix(fun$posvec[datalist$agent],
                              nrow=nrow(datalist$agent),
                              ncol=ncol(datalist$agent))
  }

  return(datalist)

}







#' Identify unique comparisons within a network
#'
#' Identify unique contrasts within a network that make up all the head-to-head comparisons. Repetitions
#' of the same treatment comparison are grouped together.
#'
#' @param df A data frame containing variables `studyID` and `treatment` (as numeric codes) that
#' indicate which treatments are used in which studies.
#'
#' @return A data frame of unique comparisons in which each row represents a different comparison.
#' `t1` and `t2` indicate the treatment codes that make up the comparison. `nr` indicates the number
#' of times the given comparison is made within the network.
#'
#' If there is only a single follow-up observation for each study within the dataset (i.e. as for standard
#' network meta-analysis) `nr` will represent the number of studies that compare treatments `t1` and
#' `t2`.
#'
#' If there are multiple observations for each study within the dataset (as in time-course MBNMA)
#' `nr` will represent the number of time points in the dataset in which treatments `t1` and `t2` are
#' compared.
#'
#' @examples
#' df <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify unique comparisons within the data
#' mbnma.comparisons(df)
#'
#'
#' # Using the triptans headache dataset
#' network <- mbnma.network(triptans) # Adds treatment identifiers
#' mbnma.comparisons(network$data.ab)
#'
#' @export
mbnma.comparisons <- function(df)
{
  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(df, add=argcheck)
  checkmate::assertNames(names(df), must.include = c("studyID", "treatment"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  df <- dplyr::arrange(df, df$studyID, df$treatment)

  t1 <- vector()
  t2 <- vector()

  # Generate vectors of treatment comparisons
  for (i in seq_along(df[["studyID"]])) {

    k <- i+1

    while (k<=nrow(df) &
           df[["studyID"]][k] == df[["studyID"]][i] &
           !is.null(df[["studyID"]][k])) {

      # Ensures ordering of t1 to t2 is lowest to highest
      t <- sort(c(df[["treatment"]][i], df[["treatment"]][k]))

      t1 <- append(t1, t[1])
      t2 <- append(t2, t[2])

      k <- k+1
    }

  }

  comparisons <- data.frame("t1"=t1, "t2"=t2)

  # Count number of repeats for each comparison
  comparisons <- comparisons %>%
    dplyr::group_by(t1, t2) %>%
    dplyr::mutate(nr=dplyr::n())

  # Make unique set of comparisons
  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)

}







#' Drop studies that are not connected to the network reference treatment
#'
#' @param connect.dose A boolean object to indicate whether treatments should be
#' kept in the network if they connect via the simplest possible dose-response
#' relationship (`TRUE`) or not (`FALSE`). Simplest possible dose-response relationship
#' is any function with a single dose-response parameter (e.g. linear, exponential)
#' @inheritParams mbnma.run
#'
#' @return A list containing a single row per arm data frame containing only studies that are
#' connected to the network reference treatment, and a character vector of treatment labels
#'
#' @examples
#' # Using the triptans headache dataset
#' network <- mbnma.network(triptans)
#' drops <- drop.disconnected(network)
#'
#' # No studies have been dropped since network is fully connected
#' length(unique(network$data.ab$studyID))==length(unique(drops$data.ab$studyID))
#'
#'
#' # Make data with no placebo
#' noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
#' net.noplac <- mbnma.network(noplac.df)
#'
#' # Studies are dropped as some only connect via the dose-response function
#' drops <- drop.disconnected(net.noplac, connect.dose=FALSE)
#' length(unique(net.noplac$data.ab$studyID))==length(unique(drops$data.ab$studyID))
#'
#' # Studies are not dropped if they connect via the dose-response function
#' drops <- drop.disconnected(net.noplac, connect.dose=TRUE)
#' length(unique(net.noplac$data.ab$studyID))==length(unique(drops$data.ab$studyID))
#'
#' @export
drop.disconnected <- function(network, connect.dose=FALSE) {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertLogical(connect.dose, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Set number of parameters for minimum connection
  if (connect.dose==FALSE) {
    doselink <- 10000 # Impossibly large number of doses (i.e. forces disconnectedness)
  } else {doselink <- 1}

  trt.labs <- network$treatments

  # Check connectivity
  discon <- suppressMessages(suppressWarnings(check.network(plot.invisible(network, level="treatment", v.color = "connect", doselink=doselink))))

  data.ab <- network$data.ab

  data.ab$treatment <- as.character(factor(data.ab$treatment, labels=network$treatments))

  # Identify disconnected studies via disconnected treatments dropped from discon
  drops <- vector()
  studies <- unique(data.ab$studyID)
  for (i in seq_along(studies)) {
    if (any(discon %in% data.ab$treatment[data.ab$studyID==studies[i]])) {
      drops <- append(drops, studies[i])
    }
  }

  # Drop drops studies from data
  data.ab <- data.ab[!(data.ab$studyID %in% drops),]
  trt.labs <- network$treatments[!(network$treatments %in% discon)]

  data.ab$treatment <- as.numeric(factor(data.ab$treatment, levels = trt.labs))

  return(list("data.ab"=data.ab, "trt.labs"=trt.labs))
}







#' Replace doses with indices of doses in order
#'
#' @inheritParams add_index
#' @noRd
index.dose <- function(data.ab) {
  agents <- sort(unique(data.ab$agent))
  maxdose <- vector()
  for (i in seq_along(agents)) {
    df <- data.ab[data.ab$agent==agents[i],]
    doses <- sort(unique(df$dose))
    maxdose <- append(maxdose, max(doses))

    for (k in seq_along(doses)) {

      if (doses[k]==0) {
        data.ab$dose[data.ab$agent==agents[i] &
                       data.ab$dose==doses[k]] <- -1
      }

      data.ab$dose[data.ab$agent==agents[i] &
                     data.ab$dose==doses[k]] <- -k -1
    }
  }
  data.ab$dose <- -data.ab$dose
  return(list("data.ab"=data.ab, "maxdose"=maxdose))
}





#' Adds placebo comparisons for dose-response relationship
#'
#' Function adds additional rows to a data.frame of comparisons in a network that account
#' for the relationship between placebo and other agents via the dose-response
#' relationship.
#'
#' @param data.ab A data frame stored in an `mbnma.network` object (`mbnma.network$data.ab`)
#' @param level A character that can take either `"treatment"` or `"agent"` to indicate the level of the
#' network for which to identify dose-response
#' @inheritParams mbnma.network
#'
DR.comparisons <- function(data.ab, level="treatment", doselink=NULL) {
  t1 <- vector()
  t2 <- vector()

  studies <- unique(data.ab$studyID)
  for (i in seq_along(studies)) {
    subset <- data.ab[data.ab$studyID==studies[i],]
    subset <- subset %>%
      dplyr::group_by(agent) %>%
      dplyr::mutate(nagent=dplyr::n())

    if (any(subset$nagent>=doselink)) {
      for (k in 1:nrow(subset)) {
        t1 <- append(t1, 0)
        t2 <- append(t2, subset[[level]][k])
      }
    }
  }

  comparisons <- data.frame("t1"=t1, "t2"=t2)

  comparisons <- comparisons %>%
    dplyr::group_by(t1, t2) %>%
    dplyr::mutate(nr=dplyr::n())

  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)
}



#' Change the network reference treatment
#'
#' @param ref A positive integer indicating the *treatment* code of the new reference
#' treatment to use
#'
#' @return An object of `class("mbnma.network")` that has a new reference treatment.
#' The new object is only really used as an intermediate in other package functions
#' and it should not be used separately, as some characteristics of the dataset may
#' not be properly encoded.
#'
#' @inheritParams mbnma.network
change.netref <- function(network, ref=1) {
  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertIntegerish(ref, len=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  data.ab <- network$data.ab

  if (!(ref %in% data.ab$treatment)) {
    stop("`ref` does not match any of the treatment codes in `network`")
  }

  trtcodes <- data.ab$treatment
  trtcodes[trtcodes==ref] <- 0
  trtcodes <- as.numeric(factor(trtcodes))

  trtnames <- network$treatments
  trtnames <- c(trtnames[ref], trtnames[-ref])

  if ("agents" %in% names(network)) {

    # Recode by agent
    data.ab$treatment <- trtcodes
    data.ab <- dplyr::arrange(data.ab, data.ab$studyID, data.ab$treatment)
    data.ab <- add_index(data.ab, agents = network$agents, treatments=trtnames)

    network$data.ab <- data.ab$data.ab
    network$treatments <- trtnames

  } else {

    # Recode by treatment
    data.ab$treatment <- trtcodes
    data.ab$agent <- NULL
    data.ab <- dplyr::arrange(data.ab, data.ab$studyID, data.ab$treatment)
    data.ab <- add_index(data.ab, agents = NULL, treatments=NULL)

    network$data.ab <- data.ab$data.ab
    network$treatments <- trtnames
    network$agents <- NULL
  }

  return(network)
}






#' Used to remove superfluous outputs from rjags object if agent-specific dose-response functions are fitted
#'
#' Function not currently in use
#'
#' @noRd
cutjags <- function(jagsresult) {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(jagsresult, "rjags", add=argcheck)
  checkmate::reportAssertions(argcheck)

  fun <- jagsresult$model.arg$fun

  if (length(fun$name)>1) {

    # Drop first element in posvec if placebo is included
    posvec <- fun$posvec
    if (jagsresult$network$agents[1]=="Placebo") {
      posvec <- posvec[-1]
    }

    apool <- fun$paramlist

    univec <- unique(posvec)
    for (i in seq_along(univec)) {
      index <- which(posvec==univec[i])

      for (k in seq_along(apool[[i]])) {
        if ("rel" %in% apool[[i]][k]) {

          tag <- gsub("\\.", "\\\\.", names(apool[[i]])[k])

          # Cut sims.array
          dropi <- grep(paste0(tag, "\\["), dimnames(jagsresult$BUGSoutput$sims.array)[[3]])
          dropi <- dropi[-index]
          jagsresult$BUGSoutput$sims.array <- jagsresult$BUGSoutput$sims.array[,,-dropi]

          # Cut sims.list
          jagsresult$BUGSoutput$sims.list[[names(apool[[i]])[k]]] <-
            jagsresult$BUGSoutput$sims.list[[names(apool[[i]])[k]]][,index]

          # Cut sims.matrix
          dropi <- grep(paste0(tag, "\\["), colnames(jagsresult$BUGSoutput$sims.matrix))
          dropi <- dropi[-index]
          jagsresult$BUGSoutput$sims.matrix <- jagsresult$BUGSoutput$sims.matrix[,-dropi]

          # Cut summary
          dropi <- grep(paste0(tag, "\\["), rownames(jagsresult$BUGSoutput$summary))
          dropi <- dropi[-index]
          jagsresult$BUGSoutput$summary <- jagsresult$BUGSoutput$summary[-dropi,]

          # Cut mean, SD and median
          jagsresult$BUGSoutput$mean[[names(apool[[i]])[k]]] <- jagsresult$BUGSoutput$mean[[names(apool[[i]])[k]]][index]
          jagsresult$BUGSoutput$sd[[names(apool[[i]])[k]]] <- jagsresult$BUGSoutput$sd[[names(apool[[i]])[k]]][index]
          jagsresult$BUGSoutput$median[[names(apool[[i]])[k]]] <- jagsresult$BUGSoutput$median[[names(apool[[i]])[k]]][index]

        }
      }
    }

    # Yet to be changed
    #jagsresult$BUGSoutput$long.short
    #jagsresult$BUGSoutput$indexes.short

  }

  return(jagsresult)
}




#' Assigns different parameters for agent-specific (multiple) dose-response functions in a model
#'
#' FUNCTION IS NOW DEPRACATED
#'
#' For a given  function (or set of functions), it indicates which parameter corresponds
#' to which function (including a parameter name) and which agents are modelled by that parameter
#' (and therefore that function)
#'
#' @noRd
assignfuns <- function(fun, agents, user.fun, wrapper=FALSE, knots=3) {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(fun, add=argcheck)
  checkmate::assertCharacter(agents, add=argcheck)
  checkmate::assertFormula(user.fun, null.ok = TRUE, add=argcheck)
  checkmate::assertNumeric(knots, null.ok = FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Convert user.fun to string
  user.str <- as.character(user.fun[2])

  # Ensure placebo isn't contained within the function
  if ("Placebo" %in% agents) {
    agents <- agents[agents!="Placebo"]
    if (length(fun)>1) {
      fun <- fun[-1]
    }
  }


  if (length(fun)==1) {
    fun <- rep(fun, length(agents))
  }

  funlist <- list("linear"="slope", "exponential"="lambda",
                  "emax"=c("emax", "ed50"), "emax.hill"=c("emax", "ed50", "hill"),
                  "rcs"=paste0("beta.", 1:knots))

  betas <- list()
  count <- 0
  if ("user" %in% fun) {
    for (i in 1:4) {
      if (grepl(paste0("beta.", i), user.str)) {
        betas[[paste0("beta.", i)]]$fun <- "user"
        betas[[paste0("beta.", i)]]$betaname <- paste0("beta.", i)
        betas[[paste0("beta.", i)]]$param <- paste0("beta.", i)
        betas[[paste0("beta.", i)]]$agents <- which(fun %in% "user")
        count <- i
      }
    }
  }
  for (i in seq_along(funlist)) {
    if (names(funlist)[i] %in% fun) {
      for (k in seq_along(funlist[[i]])) {
        count <- count+1
        betas[[paste0("beta.", count)]]$fun <- names(funlist)[i]
        betas[[paste0("beta.", count)]]$param <- funlist[[i]][k]
        betas[[paste0("beta.", count)]]$agents <- which(fun %in% names(funlist)[i])

        if (wrapper==FALSE) {
          betas[[paste0("beta.", count)]]$betaname <- paste0("beta.", count)
        } else if (wrapper==TRUE) {
          betas[[paste0("beta.", count)]]$betaname <- funlist[[i]][k]
        }
      }
    }
  }

  return(betas)
}






#' Check if all nodes in the network are connected (identical to function in `MBNMAtime`)
#'
#' @param g An network plot of `class("igraph")`
#' @param reference A numeric value indicating which treatment code to use as the reference treatment for
#' testing that all other treatments connect to it
#'
check.network <- function(g, reference=1) {

  connects <- is.finite(igraph::shortest.paths(igraph::as.undirected(g),
                                               to=reference))
  treats <- rownames(connects)[connects==FALSE]

  if (length(treats>0)) {
    warning(paste0("The following treatments/agents are not connected\nto the network reference:\n",
                   paste(treats, collapse = "\n")))
  }
  return(treats)
}





#' Generates spline basis matrices for fitting to dose-response function
#'
#' @param x A numeric vector indicating all time points available in the dataset
#' @param spline Indicates the type of spline function. Can be either a piecewise linear spline (`"ls"`),
#' natural cubic spline (`"ns"`), or B-spline (`"bs"`).
#' @param degree a positive integer giving the degree of the polynomial from which the spline function is composed
#'  (e.g. `degree=3` represents a cubic spline).
#' @param max.dose A number indicating the maximum dose between which to calculate the spline function.
#' @param knots The number/location of knots. If a single integer is given it indicates the number of knots (they will
#'   be equally spaced across the range of doses *for each agent*). If a numeric vector is given it indicates the quantiles of the knots as
#'   a proportion of the maximum dose in the dataset. For example, if the maximum dose in the dataset
#'   is 100mg/d, `knots=c(0.1,0.5)` would indicate knots should be fitted at 10mg/d and 50mg/d.
#'
#' @return A spline basis matrix with number of rows equal to `length(x)` and the number of columns equal to the number
#' of coefficients in the spline.
#'
#' @examples
#' x <- 0:100
#'
#' genspline(x)
#'
#' # Generate a quadratic B-spline with 1 equally spaced internal knot
#' genspline(x, spline="bs", knots=2, degree=2)
#'
#' # Generate a natural cubic spline with 3 knots at selected quantiles
#' genspline(x, spline="ns", knots=c(0.1, 0.5, 0.7))
#'
#' # Generate a piecewise linear spline with 3 equally spaced knots
#' genspline(x, spline="ls", knots=3)
#'
#' @export
genspline <- function(x, spline="bs", knots=1, degree=1, max.dose=max(x)){

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertNumeric(knots, add=argcheck, lower=0)
  checkmate::assertIntegerish(degree, add=argcheck)
  checkmate::assertNumeric(max.dose, null.ok = FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Remove NA values
  knots <- knots[!is.na(knots)]

  # Check knot specification
  # if (spline=="rcs") {
  #   err <- "Minimum number of knots for 'rcs' is 3"
  #   if (length(knots)==1) {
  #     if (knots<3) {
  #       stop(err)
  #     }
  #   } else if (length(knots)>1) {
  #     if (length(knots)<3){
  #       stop(err)
  #     }
  #   }
  # }

  # Add 0 (for placebo) if not in original data to ensure spline incorporates x=0
  if (x[1]==0 & length(unique(x))==1) {

    if (length(knots)==1) {
      if (knots[1]>1) {
        return(matrix(rep(0,knots-1), nrow=1))
      } else {
        return(matrix(rep(0,1), nrow=1))
      }
    } else if (length(knots)>1 | knots[1]<1) {
      return(matrix(rep(0,length(knots)-1), nrow=1))
    }

  } else {

    x0 <- x
    if (!0 %in% x0) {
      x0 <- c(0,x0)
    }

    # Calculate quantiles for knots
    if (length(knots)==1 & knots[1]>=1) {
      p <- seq(0,1,1/(knots+1))
      p <- p[-c(1,length(p))]
      knots <- stats::quantile(0:max.dose, probs = p)
      names(knots) <- NULL
    } else if (length(knots)>1 | knots[1]<1) {
      knots <- stats::quantile(0:max.dose, probs = knots)
    }

    # Generate spline basis matrix
    if (spline=="bs") {
      splinedesign <- splines::bs(x=x0, knots=knots, degree=degree)
    # } else if (spline=="rcs") {
    #   splinedesign <- Hmisc::rcspline.eval(x0, knots = knots, inclx = TRUE)
    } else if (spline=="ns") {
      splinedesign <- splines::ns(x=x0, knots=knots)

    } else if (spline=="ls") {
      splinedesign <- lspline::lspline(x=x0, knots=knots, marginal = FALSE)
    }
    rownames(splinedesign) <- x0

    # Drop 0 if it was originally added to vector to ensure returned matrix has same size as x
    splinedesign <- splinedesign[rownames(splinedesign) %in% x,]

    if (!is.matrix(splinedesign)) {
      splinedesign <- matrix(splinedesign, nrow=1)
    }

    if (ncol(splinedesign)>4) {
      stop("splines of this complexity cannot currently be modelled using 'dspline()'...\nand your data is unlikely to be able to support it!")
    }


    return(splinedesign)

  }
}






mult2agent <- function(fun, param) {
  ind <- which(fun$params %in% param)
  listlen <- lengths(fun$paramlist)
  count <- 0
  k <- 0
  while (count<ind) {
    count <- count + listlen[k+1]
    k <- k+1
  }
  subcodes <- which(fun$posvec==k)
  return(subcodes)
}





#' Guesses the upper limit for absolute and relative effects on the link scale
#'
#' Used to identify an upper bound for SD priors
#'
#' @inheritParams getjagsdata
#' @param buffer a scaling factor by which the the limits should be multiplied by - e.g. 1.2 (the default)
#' indicates that 20% should be added to the limit values to ensure that the limits do not constrain the posterior
#' distribution
#'
#' @return a list containing two numeric elements. `"rel"` contains the maximum limit on the relative scale, `"abs"` contains the
#' maximum limit on the absolute scale.
#'
#' @noRd
calcom <- function(data.ab, link, likelihood, buffer=1.2) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit", "log", "smd"), null.ok=FALSE, add=argcheck)
  checkmate::assertNumeric(buffer, lower=0, null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  df <- data.ab

  if (likelihood=="binomial") {
    df$x <- df$r / df$n
    df$x[df$x==1 | df$x==0] <- (df$r[df$x==1 | df$x==0]+0.5) / (df$n[df$x==1 | df$x==0]+1)
  } else if (likelihood=="normal") {
    df$x <- df$y
  } else if (likelihood=="poisson") {
    df$x <- df$r / df$E
    df$x[df$x==0] <- (df$r[df$x==0]+0.5) / df$E[df$x==0]
  }

  df$xlink <- rescale.link(df$x, direction="link", link=link)

  if (link=="smd") {
    df$xlink <- df$y / (df$se * (df$n^0.5))
  }

  rel <- max(df$xlink) - min(df$xlink)
  abs <- max(df$xlink)

  # Add 20% of value to ensure prior is not restrictive
  rel <- rel*buffer
  abs <- abs*buffer

  return(list(rel=round(rel, 3), abs=round(abs,3)))
}
