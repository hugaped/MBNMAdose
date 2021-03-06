# Functions for manipulating/preparing MBNMA datasets
# Author: Hugo Pedder
# Date created: 2019-04-07

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose"))

#' Create an mbnma.network object
#'
#' Creates an object of class `mbnma.network`. Various MBNMA functions can subsequently be applied
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
#' * `N` Numeric data indicating the total number of participants within a study arm. Required for
#' binomial data
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
#' * `agents` A character vector indicating the agent identifiers that correspond to the
#' new agent codes.
#' * `treatments` A character vector indicating the treatment identifiers that correspond
#' to the new treatment codes.
#' * `classes` A character vector indicating the class identifiers (if included in the original data)
#' that correspond to the new class codes.
#'
#' @examples
#' # Using the triptans headache dataset
#' print(HF2PPITT)
#'
#' # Define network
#' network <- mbnma.network(HF2PPITT, description="Example")
#'
#' # Plot network
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



#' Validates that a dataset fulfils requirements for MBNMA
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
#' * Checks that all r and N are positive (if variables are included in dataset)
#' * Checks that all y, se, r, N and E are numeric
#' * Checks that class codes are consistent within each agent
#' * Checks that agent/class names do not contain restricted characters
#' * Checks that studies have at least two arms (if `single.arm = FALSE`)
#' * Checks that each study includes at least two treatments
#'
#' @return An error if checks are not passed. Runs silently if checks are passed
mbnma.validate.data <- function(data.ab, single.arm=FALSE) {
  # data.ab must have columns c("studyID", "agent", "dose", and either "y" and "se" or "r" and "N")
  # optional column of class

  varnames <- c("studyID", "agent", "dose")
  var_norm <- c("y", "se")
  var_bin <- c("r", "N")
  var_pois <- c("r", "E")
  data.ab <- dplyr::arrange(data.ab, data.ab$studyID, data.ab$agent, data.ab$dose)

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

  numeric.error <- vector()
  if (!is.numeric(data.ab$dose)) {
    numeric.error <- append(numeric.error, "dose")
  }

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
      if (!all(data.ab[[var_tmp[i]]]>0)) {
        stop(paste("All values for", var_tmp[i], "must be >0", sep=" "))
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


  #### Add indices

  # Do not run this function with pylr loaded!!
  data.ab <- data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(arm = sequence(dplyr::n()))

  data.ab <- data.ab %>%
    dplyr::group_by(studyID) %>%
    dplyr::mutate(narm=dplyr::n())


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
  if ("class" %in% names(newdat)) {

    recoded <- recode.agent(newdat, level = "class")
    classes <- recoded[["lvlnames"]]

    # Generate class key
    #classdata <- data.ab[, names(data.ab) %in% c("agent", "class")]
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
#' @noRd
recode.agent <- function(data.ab, level="agent") {
  # Run Checks
  checkmate::assertDataFrame(data.ab)
  checkmate::assertChoice(level, choices = c("agent", "class"))

  # Check for consistency across all dose=0
  #dose.df <- data.ab[, names(data.ab) %in% c(level, "dose")]

  # Check that there are no NA values
  if (any(is.na(data.ab[[level]]))) {
    stop(paste0("NA values not allowed for ", level))
  }

  print.msg <- FALSE

  lvls <- as.character(sort(unique(data.ab[[level]])))

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
#'   level (for MBNMA) or `"treatment"` to indicate that data should be at the treatemnt-
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
#'
#' @examples
#' # Using the triptans headache dataset
#' network <- mbnma.network(HF2PPITT)
#'
#' jagsdat <- getjagsdata(network$data.ab, likelihood="binomial", link="logit")
#'
#'
#' # Get JAGS data with class
#' df <- HF2PPITT
#' df$class <- ifelse(df$agent=="placebo", "placebo", "active")
#' netclass <- mbnma.network(df)
#'
#' jagsdat <- getjagsdata(netclass$data.ab, class=TRUE)
#'
#'
#' # Get JAGS data at the treatment level for Network Meta-Analysis
#' network <- mbnma.network(HF2PPITT)
#'
#' jagsdat <- getjagsdata(network$data.ab, level="treatment")
#'
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

  df <- dplyr::arrange(df, dplyr::desc(df$narm), df$studyID, df$arm)

  df <- transform(df, studyID=as.numeric(factor(studyID, levels=as.character(unique(df$studyID)))))

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
    codes <- dplyr::arrange(codes, codes$df.agent)
    classcode <- unique(codes)$df.class

    datalist[["Nclass"]] <- Nclass
    #datalist[["classcode"]] <- classcode
    datalist[["class"]] <- classcode
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







#' Identify unique comparisons within a network
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
#' If there is only a single follow-up observation for each study within the dataset (i.e. as for standard
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
#' # Identify unique comparisons within the data
#' mbnma.comparisons(data)
#'
#'
#' # Using the triptans headache dataset
#' network <- mbnma.network(HF2PPITT) # Adds treatment identifiers
#' mbnma.comparisons(network$data.ab)
#'
#' @export
mbnma.comparisons <- function(data)
{
  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data, add=argcheck)
  checkmate::assertNames(names(data), must.include = c("studyID", "treatment"), add=argcheck)
  #checkmate::assertInt(doselink, null.ok = TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  data <- dplyr::arrange(data, data$studyID, data$treatment)

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
    dplyr::mutate(nr=dplyr::n())

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
#' network <- mbnma.network(HF2PPITT)
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

  if (connect.dose==FALSE) {
    doselink <- 10000
  } else {doselink <- 1}

  trt.labs <- network$treatments
  discon <- suppressMessages(suppressWarnings(check.network(plot.invisible(network, level="treatment", v.color = "connect", doselink=doselink))))

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
    maxdose <- append(maxdose, length(doses))
    for (k in seq_along(doses)) {
      data.ab$dose[data.ab$agent==agents[i] &
                     data.ab$dose==doses[k]] <- -k
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
#' @param data.ab A data frame stored in `mbnma.network` object (`mbnma.network$data.ab`)
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
      # temp <- subset[subset$nagent>=doselink,]
      # for (k in 1:nrow(temp)) {
      #   t1 <- append(t1, 0)
      #   t2 <- append(t2, temp[[level]][k])
      # }
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
#' @noRd
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

  data.ab$treatment <- trtcodes
  data.ab$agent <- NULL
  data.ab <- dplyr::arrange(data.ab, data.ab$studyID, data.ab$treatment)
  data.ab <- add_index(data.ab)

  network$data.ab <- data.ab$data.ab
  network$treatments <- trtnames
  network$agents <- NULL

  return(network)
}
