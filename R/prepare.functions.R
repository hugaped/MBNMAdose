# Functions for manipulating/preparing MBNMA datasets
# Author: Hugo Pedder
# Date created: 2019-04-07

#' Create an MBNMA.network object
#'
#' Creates an object of class `MBNMA.network`. Various MBNMA functions can subsequently be applied
#' to this object.
#'
#' @param data.ab A data frame of arm-level data in "long" format containing the columns:
#' * `studyID` Study identifiers
#' * `dose` Numeric data indicating the dose
#' * `agent` Agent identifiers (can be numeric, factor or character)
#' * `y` Numeric data indicating the aggregate response for a continuous outcome. Required for
#' continuous data.
#' * `se` Numeric data indicating the standard error for a given observation. Required for
#' continuous data.
#' * `r` Numeric data indicating the number of responders within a study arm. Required for
#' binomial or poisson data.
#' * `N` Numeric data indicating the total number of participants within a study arm. Required for
#' binomial data
#' * `E` Numeric data indicating the total exposure time for participants within a study arm. Required
#' for poisson data.
#' * `class` An optional column indicating a particular class code. Agents with the same identifier
#' must also have the same class code.
#' @param reference A number or character (depending on the format of `treatment` within `data.ab`)
#' indicating the reference treatment in the network (i.e. those for which estimated relative treatment
#' effects estimated by the model will be compared to).
#' @param description Optional. Short description of the network.
#'
#' @details Missing values (`NA`) cannot be included in the dataset. Single arm studies cannot
#' be included.
#'
#' @return An object of class `MBNMA.network` which is a list containing:
#' * `description` A short description of the network
#' * `data.ab` A data frame containing the arm-level network data (treatment identifiers will have
#' been recoded to a sequential numeric code)
#' * `agents` A character vector indicating the agent identifiers that correspond to the
#' new agent codes.
#' * `classes` A character vector indicating the class identifiers (if included in the original data)
#' that correspond to the new class codes.
#'
#' @examples
#' # Using the triptans headache dataset
#' print(HF2PPITT)
#'
#' # Define network
#' network <- MBNMA.network(HF2PPITT, description="Example")
#'
#' # Define network with different network reference agent
#' network <- MBNMA.network(HF2PPITT, reference="Ce_200", description="Example")
#' @export
MBNMA.network <- function(data.ab, description="Network") {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(description, len=1, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  MBNMA.validate.data(data.ab)

  # Add indices for arms and narms and change agent/class codes
  index.data <- add_index(data.ab)

  network <- index.data
  network <- c(list("description" = description), network)

  class(network) <- "MBNMA.network"
  return(network)
}



#' Validates that a dataset fulfills requirements for MBNMA
#'
#' @inheritParams MBNMA.network
#'
#' @details Checks done within the validation:
#' * Checks data.ab has required column names
#' * Checks there are no NAs
#' * Checks that all SEs are positive (if variables are included in dataset)
#' * Checks that all r and N are positive (if variables are included in dataset)
#' * Checks that class codes are consistent within each agent
#' * Checks that studies have at least two arms (if `single.arm = FALSE`)
#' * Checks that each study includes at least two treatments
#'
#' @return An error if checks are not passed. Runs silently if checks are passed
#'
MBNMA.validate.data <- function(data.ab, single.arm=FALSE) {
  # data.ab must have columns c("studyID", "agent", "dose", and either "y" and "se" or "r" and "N")
  # optional column of class

  varnames <- c("studyID", "agent", "dose")
  var_norm <- c("y", "se")
  var_bin <- c("r", "N")
  var_pois <- c("r", "E")
  data.ab <- dplyr::arrange(data.ab, studyID, agent, dose)

  # Check data.ab has required column names
  msg <- "Required variable names are: 'studyID', 'agent', 'dose' and either `y` and `se` for data with a normal likelihood, `r` and `N` for data with a binomial likelihood, or `r` and `E` for data with a poisson likelihood"
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
  if (all(var_norm %in% names(data.ab))) {
    if (anyNA(data.ab[[var_norm[i]]])) {
      na.vars <- append(na.vars, var_norm[i])
    }
  }
  if (all(var_bin %in% names(data.ab))) {
    if (anyNA(data.ab[[var_bin[i]]])) {
      na.vars <- append(na.vars, var_bin[i])
    }
  }
  if (all(var_pois %in% names(data.ab))) {
    if (anyNA(data.ab[[var_pois[i]]])) {
      na.vars <- append(na.vars, var_pois[i])
    }
  }
  if (length(na.vars)>0) {
    stop(paste0("NA values in:\n", paste(na.vars, collapse="\n")))
  }


  # Check that all SEs are positive
  var_tmp <- c("se", var_bin, "E")
  for (i in seq_along(var_tmp)) {
    if (var_tmp[i] %in% names(data.ab)) {
      if (!all(data.ab[[var_tmp[i]]]>0)) {
        stop(paste("All values for", var_tmp[i], "must be >0", sep=" "))
      }
    }
  }


  # Generate narms index for checking if studies are only single-arm
  if (single.arm==FALSE) {
    data.ab <- data.ab %>%
      dplyr::group_by(studyID) %>%
      dplyr::mutate(narms = n())
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
  }

}




#' Add arm indices and agent identifiers to a dataset
#'
#' Adds arm (`arms`, `narms`) indices to a dataset and adds numeric identifiers for
#' agent and class (if included in the data).
#'
#' @inheritParams MBNMA.network
#'
#' @return A data frame similar to `data.ab` but with additional columns:
#' * `arm` Arm identifiers coded for each study
#' * `narm` The total number of arms in each study
#'
#' If `agent` or `class` are non-numeric or non-sequential (i.e. with missing numeric codes),
#' agents/classes in the returned data frame will be numbered and recoded to enforce sequential
#' numbering (a warning will be shown stating this).
#'
#' @examples
#' # Add indices to triptans headache dataset
#' data.ab <- add_index(HF2PPITT)
#'
#' @export
add_index <- function(data.ab) {

  # Run Checks
  checkmate::assertDataFrame(data.ab)

  if ("agent" %in% names(data.ab)) {
    recoded <- recode.agent(data.ab, level = "agent")
    treatments.df <- dplyr::arrange(data.ab, agent, dose)
    treatments <- unique(paste(treatments.df$agent, treatments.df$dose, sep="_"))
    agents <- recoded[["lvlnames"]]
    data.ab <- recoded[["data.ab"]]

    data.ab$treatment <- as.numeric(factor(paste(data.ab$agent,
                                                 data.ab$dose,
                                                 sep="_"),
                                           labels=treatments
    ))

    data.ab <- dplyr::arrange(data.ab, studyID, agent, dose)
  }


  #### Add indices

  # Do not run this function with pylr loaded!!
  data.ab <- data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(arm = sequence(n()))

  data.ab <- data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(narm=n())


  # Reorder columns in data.ab
  ord <- c("agent", "dose", "treatment", "class", "narm", "arm", "y", "se", "r", "E", "N")
  newdat <- data.frame("studyID"=data.ab$studyID)
  for (i in seq_along(ord)) {
    if (ord[i] %in% names(data.ab)) {
      newdat <- cbind(newdat, data.ab[,which(names(data.ab)==ord[i])])
    }
  }
  olddat <- data.ab[,!(names(data.ab) %in% c("studyID", ord))]
  newdat <- cbind(newdat, olddat)

  output <- list("data.ab"=newdat)

  if ("agent" %in% names(data.ab)) {
    output[["agents"]] <- agents
    output[["treatments"]] <- treatments
  }


  # Store class labels and recode (if they exist in data.ab)
  if ("class" %in% names(data.ab)) {

    recoded <- recode.agent(data.ab, level = "class")
    classes <- recoded[["lvlnames"]]

    # Generate class key
    classdata <- data.ab[, names(data.ab) %in% c("agent", "class")]
    classkey <- unique(classdata)
    classkey$agent <- factor(classkey$agent, labels=agents)
    classkey$class <- factor(classkey$class, labels=classes)

    output[["classes"]] <- classes
    output[["classkey"]] <- classkey
  }

  return(output)
}



#' Assigns agent and class variables numeric identifiers
#'
#' @param level Can take either `"agent"` or `"class"`
recode.agent <- function(data.ab, level="agent") {
  # Run Checks
  checkmate::assertDataFrame(data.ab)
  checkmate::assertChoice(level, choices = c("agent", "class"))

  if (is.factor(data.ab[[level]])) {
    agents <- levels(data.ab[[level]])
    match <- match(agents, as.character(data.ab[[level]]))
    if (any(is.na(match))) {
      agents <- agents[!is.na(match)]
    }
  }

  if (is.numeric(data.ab[[level]])) {
    if (max(data.ab[[level]]) != length(unique(data.ab[[level]])) |
        !all.equal(data.ab[[level]], as.integer(data.ab[[level]]))
    ) {
      print(paste0(level, " is being recoded to enforce sequential numbering"))
    }
    agents <- sort(unique(data.ab[[level]]))
  }

  if (is.character(data.ab[[level]])) {
    agents <- sort(unique(data.ab[[level]]))
  }

  # Numeric data must be checked that sequence is consistent for sequential numbering
  # Factor data must be allocated codes based on factor levels
  # Character data must be allocated codes automatically (alphabetically)

  # Must be numeric for MBNMA.run
  data.ab[[level]] <- as.numeric(factor(data.ab[[level]],
                                     levels=agents)) # provide factor for sorting so that reference is as given by user


  return(list("data.ab"=data.ab, "lvlnames"=agents))
}








#' Prepares data for JAGS
#'
#' Converts MBNMA data frame to a list for use in JAGS model
#'
#' @inheritParams MBNMA.run
#' @inheritParams MBNMA.network
#' @param class A boolean object indicating whether or not `data.ab` contains
#'   information on different classes of treatments
#' @param level Can take either `"agent"` to indicate that data should be at the agent-
#'   level (for MBNMA) or `"treatment"` to indicate that data should be at the treatent-
#'   level (for NMA)
#'
#' @return A named list of numbers, vector, matrices and arrays to be sent to
#'   JAGS. List elements are:
#'   * If `likelihood="normal"`:
#'     - `y` An array of mean responses for each arm within each study
#'     - `se` An array of standard errors for each arm within each study
#'   * If `likelihood="binomial"`:
#'     - `r` An array of the number of responses/count for each each arm within each study
#'     - `N` An array of the number of participants for each arm within each study
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
#' @export
getjagsdata <- function(data.ab, class=FALSE, likelihood="binomial", link="logit",
                        level="agent") {

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(level, choices=c("agent", "treatment"), null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  df <- data.ab

  varnames <- c("studyID", "arm", "narm")

  if (level=="agent") {
    varnames <- append(varnames, c("dose", "agent"))
  } else if (level=="treatment") {
    varnames <- append(varnames, c("treatment"))
  }

  if (class==TRUE) {
    varnames <- append(varnames, "class")
  }

  if (likelihood == "binomial") {
    datavars <- c("r", "N")
  } else if (likelihood == "poisson") {
    datavars <- c("r", "E")
  } else if (likelihood=="normal") {
    datavars <- c("y", "se")
  } else {
    stop("`likelihood` can be either `binomial`, `poisson`, or `normal`")
  }
  varnames <- append(varnames, datavars)

  if (!all(varnames %in% names(df))) {
    msg <- paste0("Variables are missing from dataset:\n",
                  paste(varnames[!(varnames %in% names(df))], collapse="\n"))
    stop(msg)
  }

  df <- dplyr::arrange(df, dplyr::desc(narm), studyID, arm)

  df <- transform(df,studyID=as.numeric(factor(studyID, levels=as.character(unique(df$studyID)))))

  for (i in seq_along(datavars)) {
    assign(datavars[i], array(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                              dim=c(max(as.numeric(df$studyID)),
                                    max(df$arm)
                              ))
           )
  }

  narm <- vector()
  NS <- max(as.numeric(df$studyID))

  datalist <- list(get(datavars[1]), get(datavars[2]),
                   narm=narm, NS=NS)
  names(datalist)[1:2] <- datavars

  if (level=="agent") {
    datalist[["Nagent"]] <- max(df$agent)
    datalist[["agent"]] <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                                  nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
                                  )
    datalist[["dose"]] <- datalist[["agent"]]

  } else if (level=="treatment") {
    datalist[["NT"]] <- max(df$treatment)

    datalist[["treatment"]] <- matrix(rep(NA, max(as.numeric(df$studyID))*max(df$arm)),
                                      nrow = max(as.numeric(df$studyID)), ncol = max(df$arm)
                                      )
  }


  if (class==TRUE) {
    Nclass <- max(df$class)

    codes <- data.frame(df$agent, df$class)
    codes <- dplyr::arrange(codes, df.agent)
    classcode <- unique(codes)$df.class

    datalist[["Nclass"]] <- Nclass
    datalist[["classcode"]] <- classcode
  }


  for (i in 1:max(as.numeric(df$studyID))) {
    for (k in 1:max(df$arm[df$studyID==i])) {
      for (m in seq_along(datavars)) {
        datalist[[m]][i,k] <- df[[datavars[m]]][as.numeric(df$studyID)==i &
                              df$arm==k]
      }

      if (level=="agent") {
        datalist[["agent"]][i,k] <- max(df$agent[as.numeric(df$studyID)==i &
                                                   df$arm==k])
        datalist[["dose"]][i,k] <- max(df$dose[as.numeric(df$studyID)==i &
                                                 df$arm==k])
      } else if (level=="treatment") {
        datalist[["treatment"]][i,k] <- max(df$treatment[as.numeric(df$studyID)==i &
                                                   df$arm==k])
      }

    }

    datalist[["narm"]] <- append(datalist[["narm"]], max(df$arm[as.numeric(df$studyID)==i]))
  }

  return(datalist)

}







#' Identify unique comparisons within a network (identical to MBNMAtime)
#'
#' Identify unique contrasts within a network that make up all the head-to-head comparisons. Repetitions
#' of the same treatment comparison are grouped together.
#'
#' @param data A data frame containing variables `studyID` and `treatment` (as numeric codes) that
#' indicate which treatments are used in which studies.
#'
#' @return A data frame of unique comparisons in which each row represents a different comparison.
#' `t1` and `t2` indicate the treatment codes that make up the comparison. `nr` indicates the number
#' of times the given comparison is made within the network.
#'
#' If there is only a single observation for each study within the dataset (i.e. as for standard
#' network meta-analysis) `nr` will represent the number of studies that compare treatments `t1` and
#' `t2`.
#'
#' If there are multiple observations for each study within the dataset (as in time-course MBNMA)
#' `nr` will represent the number of time points in the dataset in which treatments `t1` and `t2` are
#' compared.
#'
#' @examples
#' data <- data.frame(studyID=c(1,1,2,2,3,3,4,4,5,5,5),
#'   treatment=c(1,2,1,3,2,3,3,4,1,2,4)
#'   )
#'
#' # Identify comparisons informed by direct and indirect evidence
#' MBNMA.comparisons(data)
#' @export
MBNMA.comparisons <- function(data)
{
  # Assert checks
  checkmate::assertDataFrame(data)

  if (all(names(data) %in% c("studyID", "treatment") != TRUE)) {
    stop("data must contain variables 'studyID' and 'treatment'")
  }

  t1 <- vector()
  t2 <- vector()

  for (i in seq_along(data[["studyID"]])) {

    k <- i+1

    while (k<=nrow(data) &
           data[["studyID"]][k] == data[["studyID"]][i] &
           !is.null(data[["studyID"]][k])) {

      # Ensures ordering of t1 to t2 is lowest to highest
      t <- sort(c(data[["treatment"]][i], data[["treatment"]][k]))

      t1 <- append(t1, t[1])
      t2 <- append(t2, t[2])

      k <- k+1
    }

  }

  comparisons <- data.frame("t1"=t1, "t2"=t2)

  comparisons <- comparisons %>%
    dplyr::group_by(t1, t2) %>%
    dplyr::mutate(nr=n())

  comparisons <- unique(comparisons)
  comparisons <- dplyr::arrange(comparisons, t1, t2)

  return(comparisons)

}







#' Drop studies that are not connected to the network reference treatment
#'
#' @return A single row per arm data frame containing only studies that are
#' connected to the network reference treatment
#'
#' @export
drop.disconnected <- function(network) {

  trt.labs <- network$treatments
  #discon <- suppressWarnings(check.network(plot(network, level="treatment", v.color = "connect")))
  png("NUL")
  discon <- suppressWarnings(check.network(plot(network, level="treatment", v.color = "connect")))
  dev.off()

  data.ab <- network$data.ab

  data.ab$treatment <- as.character(factor(data.ab$treatment, labels=network$treatments))

  drops <- vector()
  studies <- unique(data.ab$studyID)
  for (i in seq_along(studies)) {
    if (any(discon %in% data.ab$treatment[data.ab$studyID==studies[i]])) {
      drops <- append(drops, studies[i])
    }
  }

  #data.ab <- data.ab[!(data.ab$treatment %in% discon),]
  data.ab <- data.ab[!(data.ab$studyID %in% drops),]
  trt.labs <- network$treatments[!(network$treatments %in% discon)]

  data.ab$treatment <- as.numeric(factor(data.ab$treatment, levels = trt.labs))

  return(list("data.ab"=data.ab, "trt.labs"=trt.labs))
}
