# Functions for running MBNMA models
# Author: Hugo Pedder
# Date created: 2019-04-18


#' Run MBNMA dose-response models
#'
#' Fits a Bayesian dose-response for model-based network meta-analysis
#' (MBNMA) that can account for multiple doses of different agents by
#' applying a desired dose-response function. Follows the methods
#' of Mawdsley (REF).
#'
#' @param network An object of class `MBNMA.network`.
#' @param parameters.to.save A character vector containing names of parameters
#'   to monitor in JAGS
#' @param fun is a character specifying a functional form to be assigned to the
#'   dose-response. Options are given in `details`.
#' @param user.fun A character specifying any relationship including `dose` and
#'   one/several of: `beta.1`, `beta.2`, `beta.3`.
#' @param model.file A JAGS model written as a character object that can be used
#'   to overwrite the JAGS model that is automatically written based on the
#'   specified options. Useful when ammending priors using replace.prior()
#'
#' @param beta.1 Refers to dose-parameter(s) specified within the dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.2 Refers to dose-parameter(s) specified within the dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.3 Refers to dose-parameter(s) specified within the dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @param method Can take either `"common"` or `"random"` to indicate whether relative effects
#'   should be modelled with between-study heterogeneity or not (see LINK `Details`).
#' @param class.effect A list of named strings that determines which dose-response
#'   parameters to model with a class effect and what that effect should be
#'   (`"common"` or `"random"`). For example: `list("beta.2"="fixed", "beta.3"="random")`.
#' @param likelihood A string indicating the likelihood to use in the model. Can take either `"binomial"`,
#'   `"normal"` or `"poisson"`. If left as `NULL` the likelihood will be inferred from the data.
#' @param link A string indicating the link function to use in the model. Can take any link function
#'   defined within JAGS (`"logit"`, `"log"`, `"probit"`, `"cloglog"`) or be assigned the value `"identity"` for
#'   and identity link function. If left as `NULL` the link function will be automatically assigned based
#'   on the likelihood.
#' @param var.scale A numeric vector indicating the relative scale of variances between
#' correlated dose-response parameters when relative effects are modelled on more than
#' one dose-response parameter and `method="random"` (see Details LINK). Each element of
#' the vector refers to the relative scale of each of the dose-response parameters that is
#' modelled using relative effects.
#'
#' @param pd Can take either:
#'   * `pv` only pV will be reported (as automatically outputted by R2jags).
#'   * `plugin` calculates pD by the plug-in
#'   method REF. It is faster, but may output negative
#'   non-sensical values, due to skewed deviances that can arise with non-linear models.
#'   * `pd.kl` calculates pD by the Kullback–Leibler divergence REF. This
#'   will require running the model for additional iterations but
#'   will always produce a positive result.
#'   * `popt` calculates pD using an optimism adjustment which allows for calculation
#'   of the penalized expected deviance REF.
#' @param parallel A boolean value that indicates whether JAGS should be run in
#'   parallel (`TRUE`) or not (`FALSE`). If `TRUE` then the number of cores to
#'   use is automatically calculated.
#' @param arg.params Contains a list of arguments sent to `MBNMA.run()` by time-course
#' specific wrapper functions
#' @param ... Arguments to be sent to R2jags.
#'
#' @inheritParams replace.prior
#'
#' @details When relative effects are modelled on more than one dose-response parameter and
#' `method="random"`, correlation between the dose-response parameters is automatically
#' estimated using a vague Wishart prior. This prior can be made slightly more informative
#' by specifying the relative scale of variances between the dose-response parameters using
#' `var.scale`.
#'
#' @return An object of S3 class `c("MBNMA", "rjags")`` containing parameter
#'   results from the model. Can be summarized by `print()` and can check
#'   traceplots using `traceplot()` (from the `R2jags` package).
#'
#'   Nodes that are automatically monitored (if present in the model) have the
#'   following interpretation. They will have an additional suffix that relates
#'   to the name/number of the dose-response paramter to which they correspond
#'   (e.g. `d.et50` or `d.1`):
#'   TO ADD FROM MBNMAtime
#'
#'
#'   If there are errors in the JAGS model code then the object will be a list
#'   consisting of two elements - an error message from JAGS that can help with
#'   debugging and `model.arg`, a list of arguments provided to `MBNMA.run()`
#'   which includes `jagscode`, the JAGS code for the model that can help
#'   users identify the source of the error.
#'
#' @section Dose-response parameters:
#' EXPLAIN DIFFERENCE BETWEEN "rel", "common", "random", and numeric value
#' * `"common"` implies that all studies estimate the same true effect
#' (akin to a "fixed effect" meta-analysis)
#' * `"random"` implies that all studies estimate a separate true effect, but that each
#' of these true effects vary randomly around a true mean effect. This approach allows
#' for modelling of between-study heterogeneity.
#' * `numeric()` Assigned a numeric value - this can only be used if `pool="const"`. It indicates that
#' this dose-response parameter should not be estimated from the data but should be assigned
#' the numeric value determined by the user. This can be useful for fixing specific dose-response
#' parameters (e.g. Hill parameters in Emax functions).
#'
#'
#' @section Dose-response function:
#'   Several general dose-response functions are provided, but a
#'   user-defined dose-response relationship can instead be used.
#'
#'   Built-in time-course functions are:
#'   * `"linear"`: `beta.1` refers to the gradient
#'   * `"exponential"`: `beta.1` refers to the rate of gain/decay
#'   * `"emax"` (emax without a Hill parameter): `beta.1` refers to
#'   Emax parameter, `beta.2` refers to ET50 parameter
#'   * `"emax.hill"` (emax with a Hill parameter): `beta.1` refers to Emax parameter, `beta.2` refers
#'   to ET50 parameter, `beta.3` refers to Hill parameter
#'   * `"nonparam.up"` (monotonically increasing non-parametric dose-response relationship following
#'   the method of OWENS REF)
#'   * `"nonparam.down"` (monotonically decreasing non-parametric dose-response relationship following
#'   the method of OWENS REF)
#'   * `"user"` (user-defined function: `user.fun` must be specified in arguments)
#'   * `"none"` (no dose-response relationship).
#'
#' @importFrom Rdpack reprompt
#' @importFrom magrittr "%>%"
#'
#' @references
#'   \insertAllCited
#'
#'
#'
#' @examples
#' @export
MBNMA.run <- function(network, parameters.to.save=NULL,
                      fun="linear", user.fun=NULL,
                      model.file=NULL,
                      beta.1="rel",
                      beta.2=NULL, beta.3=NULL,
                      method="common",
                      class.effect=list(),
                      var.scale=NULL,
                      pd="pv", parallel=TRUE,
                      likelihood=NULL, link=NULL,
                      priors=NULL,
                      arg.params=NULL, ...
) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "MBNMA.network", add=argcheck)
  checkmate::assertCharacter(model.file, len=1, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(pd, choices=c("pv", "pd.kl", "plugin", "popt"), null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertList(arg.params, unique=TRUE, null.ok=TRUE, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]


  if (!is.null(arg.params)) {
    if (!all((names(arg.params)) %in% c("wrap.params", "run.params"))) {
      stop("arg.params has been incorrectly specified")
    }
    wrap.params <- arg.params$wrap.params
    run.params <- arg.params$run.params

    fun.params <- names(class.effect)
    for (k in seq_along(fun.params)) {
      fun.params[k] <- run.params[which(fun.params[k]==wrap.params)]
    }
    names(class.effect) <- fun.params

  } else if (is.null(arg.params)) {
    wrap.params <- list(beta.1, beta.2, beta.3)
    wrap.params <- which(sapply(wrap.params,
                                is.character))
  }

  if (is.null(model.file)) {
    model <- MBNMA.write(fun=fun, user.fun=user.fun,
                         beta.1=beta.1, beta.2=beta.2, beta.3=beta.3,
                         method=method,
                         class.effect=class.effect,
                         var.scale=var.scale,
                         likelihood=likelihood, link=link
    )

    # Change beta.1 and beta.2 to emax and et50, etc. if necessary
    # Change beta.1 and beta.2 to emax and et50, etc. if necessary
    if (!is.null(arg.params)) {
      code.params <- c("d", "beta", "sd", "tau", "D", "sd.D")
      for (i in seq_along(wrap.params)) {
        for (k in seq_along(code.params)) {
          model <- gsub(paste(code.params[k], strsplit(run.params[i], split="[.]")[[1]][2], sep="."),
                        paste(code.params[k], wrap.params[i], sep="."), model)
        }
      }

      wrap.params <- wrap.params[which(sapply(list(beta.1, beta.2, beta.3),
                                              is.character))]



      # This bit is messy...remove if problematic though it will cause issues with parameters.to.save
      # param.list <- list("beta.1"=beta.1.str, "beta.2"=beta.2.str, "beta.3"=beta.3.str, "beta.4"=beta.4.str)
      # temp <- names(unlist(sapply(param.list, function(x) {
      #   if(!is.null(x)) {
      #     if (x %in% c("arm.common", "arm.random", "const.common", "const.random")) {
      #       x
      #     }
      #   }
      # })))
      # temp <- unlist(sapply(temp, FUN=function(x) strsplit(x, split="[.]")[[1]][2]))
      # wrap.params <- append(wrap.params, temp)

    } else {
      wrap.params <- which(sapply(list(beta.1, beta.2, beta.3),
                                  is.character))
    }

    if (!is.null(priors)) {
      model <- replace.prior(priors=priors, model=model)
    }

  } else {
    warning("All parameter specifications (dose-response parameters, class effects, priors, etc.) are being overwritten by `model.file`")
    model <- model.file
  }

  assigned.parameters.to.save <- parameters.to.save
  if (is.null(parameters.to.save)) {
    parameters.to.save <-
      gen.parameters.to.save(model.params=wrap.params, model=model)
  }

  # Add nodes to monitor to calculate plugin pd
  if (pd=="plugin") {
    pluginvars <- c("theta", "resdev")
    for (param in seq_along(pluginvars)) {
      if (!(pluginvars[param] %in% parameters.to.save)) {
        parameters.to.save <- append(parameters.to.save, pluginvars[param])
      }
    }
    message("The following parameters have been monitored to allow pD plugin calculation: ",
            paste(pluginvars, collapse=", "))
  }

  if (length(class.effect)>0) {
    class <- TRUE
  } else {class <- FALSE}


  #### Run jags model ####

  if (fun %in% c("nonparam.up", "nonparam.down")) {
    # Change doses to dose indices
    data.ab <- index.dose(network[["data.ab"]])[["data.ab"]]
  } else {
    data.ab <- network[["data.ab"]]
  }

  result.jags <- MBNMA.jags(data.ab, model,
                            class=class,
                            parameters.to.save=parameters.to.save,
                            likelihood=likelihood, link=link,
                            ...)
  result <- result.jags[["jagsoutput"]]
  jagsdata <- result.jags[["jagsdata"]]

  if (pd == "pd.kl" | pd == "popt") {
    if (pd=="pd.kl") {
      temp <- rjags::dic.samples(result$model, n.iter=1000, type="pD")
    } else if (pd=="popt") {
      temp <- rjags::dic.samples(result$model, n.iter=1000, type="popt")
    }
    result$BUGSoutput$pD <- sum(temp$penalty)

    # temp <-
    #   runjags::run.jags(model.file=model,
    #            monitor=c("deviance", "dic", "full.pd"),
    #            data=jagsdata, n.chains=3,
    #            burnin=10000, adapt=5000, sample=2000,
    #            jags.refresh=0.1, method="rjags")
    # result$BUGSoutput$pD <- temp$deviance.sum[2]
  } else if (pd == "plugin") {
    # plugin method
    warning("Plugin method only works for normal likelihood")
    result$BUGSoutput$pD <- pDcalc(y=jagsdata[["y"]], se=jagsdata[["se"]], fups=jagsdata[["fups"]], narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                                   theta.result=result$BUGSoutput$mean$theta, resdev.result=result$BUGSoutput$mean$resdev)
  }

  # Recalculate DIC so it is adjusted for choice of pD
  result$DIC <- result$BUGSoutput$pD + result$BUGSoutput$median$deviance

  # Add variables for other key model characteristics (for predict and plot functions)
  model.arg <- list("parameters.to.save"=assigned.parameters.to.save,
                    "fun"=fun, "user.fun"=user.fun,
                    "jagscode"=model,
                    "beta.1"=beta.1, "beta.2"=beta.2,
                    "beta.3"=beta.3,
                    "method"=method,
                    "likelihood"=likelihood, "link"=link,
                    "class.effect"=class.effect,
                    "parallel"=parallel, "pd"=pd,
                    "priors"=get.prior(model), "arg.params"=arg.params)
  result[["model.arg"]] <- model.arg
  result[["type"]] <- "dose"

  if (!("error" %in% names(result))) {
    class(result) <- c("MBNMA", class(result))
  }

  return(result)

}





MBNMA.jags <- function(data.ab, model,
                       class=FALSE, rho=NULL, covar=NULL,
                       parameters.to.save=parameters.to.save,
                       likelihood=NULL, link=NULL,
                       warn.rhat=FALSE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(model, any.missing=FALSE, len=1, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertCharacter(parameters.to.save, any.missing=FALSE, unique=TRUE,
                             null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  if (is.null(likelihood) & is.null(link)) {
    # For MBNMAtime
    jagsdata <- getjagsdata(data.ab, class=class, rho=rho, covstruct=covar) # get data into jags correct format (list("fups", "NT", "NS", "narm", "y", "se", "treat", "time"))
  } else if (is.null(rho) & is.null(covar)) {
    # For MBNMAdose
    jagsdata <- getjagsdata(data.ab, class=class,
                            likelihood=likelihood, link=link) # get data into jags correct format
  }


  # Add variable for maxtime to jagsdata if required
  if (grepl("maxtime", model)) {
    jagsdata[["maxtime"]] <- max(data.ab$time)
  } else if (grepl("maxdose", model)) {
    jagsdata[["maxdose"]] <- index.dose(network[["data.ab"]])[["maxdose"]]
  }

  # Put data from jagsdata into separate R objects
  for (i in seq_along(jagsdata)) {
    ##first extract the object value
    temp <- jagsdata[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(jagsdata)[[i]],"<- temp")))
  }

  # Take names of variables in jagsdata for use in rjags
  jagsvars <- list()
  for (i in seq_along(names(jagsdata))) {
    jagsvars[[i]] <- names(jagsdata)[i]
  }

  # Create a temporary model file
  tmpf=tempfile()
  tmps=file(tmpf,"w")
  cat(model,file=tmps)
  close(tmps)

  out <- tryCatch({
    result <- R2jags::jags(data=jagsvars, model.file=tmpf,
                           parameters.to.save=parameters.to.save,
                           ...
    )
  },
  error=function(cond) {
    message(cond)
    return(list("error"=cond))
  }
  )

  # Gives warning if any rhat values > 1.02
  if (warn.rhat==TRUE) {
    if (!("error" %in% names(out))) {
      rhat.warning(out)
    }
  }

  return(list("jagsoutput"=out, "jagsdata"=jagsdata))
}





#' Automatically generate parameters to save for a dose-response MBNMA model
#'
#' Identical to `gen.parameters.to.save()` in MBNMAtime
gen.parameters.to.save <- function(model.params, model) {
  # model.params is a vector (numeric/character) of the names of the dose-response parameters in the model
  #e.g. c(1, 2, 3) or c("emax", "et50")
  # model is a JAGS model written as a character object

  checkmate::assertCharacter(model, len=1)

  model.params <- as.character(model.params)

  # Set some automatic parameters based on the model code
  parameters.to.save <- vector()
  for (i in seq_along(model.params)) {
    if (grepl(paste0("\\\nd\\.", model.params[i], "\\[(c,)?k\\] ~"), model)==TRUE |
        grepl(paste0("\\\nd\\.", model.params[i], "\\[k\\] <- mult\\["), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("d.", model.params[i]))
    } else if (grepl(paste0("\\\nd\\.", model.params[i], "\\[k\\] ~"), model)==FALSE) {
      if (grepl(paste0("\\\nbeta\\.", model.params[i], "(\\[k\\])? ~"), model)==TRUE) {
        parameters.to.save <- append(parameters.to.save, paste0("beta.", model.params[i]))
      }
    }
    if (grepl(paste0("\\\nsd\\.", model.params[i], " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", model.params[i]))
    }
    if (grepl(paste0("\\\nsd\\.beta.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd\\.beta\\.", model.params[i]))
    }
    if (grepl(paste0("\\\nD\\.", model.params[i], " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("D.", model.params[i]))
    }
    if (grepl(paste0("\\\nsd\\.D\\.", model.params[i], " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.D.", model.params[i]))
    }
    if (grepl(paste0("\\\nBETA\\.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("BETA.", model.params[i]))
    }
    if (grepl(paste0("\\\nsd\\.BETA\\.", model.params[i]), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.BETA.", model.params[i]))
    }
  }

  for (i in 1:4) {
    if (grepl(paste0("\\\nd\\.", i, " ~"), model)==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("d.", i))
    }
  }

  # For MBNMAtime
  if (grepl("rho", model)==TRUE) {
    parameters.to.save <- append(parameters.to.save, "rho")
  } else {
    parameters.to.save <- append(parameters.to.save, c("totresdev"))
  }

  # For MBNMAdose
  if (grepl("\\\nsd ~", model)==TRUE) {
    parameters.to.save <- append(parameters.to.save, "sd")
  }

  return(unique(parameters.to.save))

}






#' Run an NMA model
#'
#' Used for calculating split NMA results, either when comparing models that do not
#' account for dose-response relationship, or to estimate split results for `overlay.split`.
#' Results can also be compared between consistency (`UME=FALSE`) and inconsistency
#' (`UME=TRUE`) models to test the validity of the consistency assumption.
#'
#' @inheritParams MBNMA.run
#' @param drop.discon A boolean object that indicates whether or not to drop disconnected
#'   studies from the network.
#' @param UME A boolean object to indicate whether to fit an Unrelated Mean Effects model
#'   that does not assume consistency and so can be used to test if the consistency
#'   assumption is valid.
#'
#' @examples
#' @export
NMA.run <- function(network, method="common", likelihood="binomial", link="logit",
                    warn.rhat=TRUE, n.iter=5000, drop.discon=TRUE, UME=FALSE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "MBNMA.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::assertLogical(warn.rhat, add=argcheck)
  checkmate::assertIntegerish(n.iter, null.ok = TRUE, add=argcheck)
  checkmate::assertLogical(drop.discon, add=argcheck)
  checkmate::assertLogical(UME, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  #### Write model for NMA ####
  model <- write.NMA(method=method, likelihood=likelihood, link=link, UME=UME)


  #### Parameters ####
  parameters.to.save <- c("d", "totresdev")
  if (method=="random") {
    parameters.to.save <- append(parameters.to.save, "sd")
  }

  #### Prepare data ####
  data.ab <- network$data.ab
  trt.labs <- network$treatments

  # data.ab$treatment <- paste(as.character(factor(data.ab$agent, labels=network$agents)),
  #                            data.ab$dose,
  #                            sep="_"
  # )

  # Check treatments that are not connected and remove if not
  if (drop.discon==TRUE) {
    connect <- drop.disconnected(network)
    data.ab <- connect[["data.ab"]]
    trt.labs <- connect[["trt.labs"]]
  }

  jagsdata <- getjagsdata(data.ab, likelihood = likelihood, link=link, level="treatment")


  #### Run JAGS model ####
  # Put data from jagsdata into separate R objects
  for (i in seq_along(jagsdata)) {
    ##first extract the object value
    temp <- jagsdata[[i]]
    ##now create a new variable with the original name of the list item
    eval(parse(text=paste(names(jagsdata)[[i]],"<- temp")))
  }

  # Take names of variables in jagsdata for use in rjags
  jagsvars <- list()
  for (i in seq_along(names(jagsdata))) {
    jagsvars[[i]] <- names(jagsdata)[i]
  }

  # Create a temporary model file
  tmpf=tempfile()
  tmps=file(tmpf,"w")
  cat(model,file=tmps)
  close(tmps)

  out <- tryCatch({
    result <- R2jags::jags(data=jagsvars, model.file=tmpf,
                           parameters.to.save=parameters.to.save,
                           n.iter=n.iter,
                           ...
    )
  },
  error=function(cond) {
    message(cond)
    return(list("error"=cond))
  }
  )

  # Gives warning if any rhat values > 1.02
  if (warn.rhat==TRUE) {
    if (!("error" %in% names(out))) {
      rhat.warning(out)
    }
  }

  output <- list("jagsresult"=out, "trt.labs"=trt.labs)
  class(output) <- "NMA"
  return(output)

}







#' Check likelihood and link function
#'
#' Checks that likelihood and link function is provided and confirm that the correct
#' form of data is provided.
#'
#' @export
check.likelink <- function(data.ab, likelihood=NULL, link=NULL) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit"), null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  if (is.null(likelihood)) {
    if (all(c("r", "N") %in% names(data.ab))) {
      likelihood <- "binomial"
      message("`likelihood` not given by user - set to `binomial` based on data provided")
    } else if (all(c("y", "se") %in% names(data.ab))) {
      likelihood <- "normal"
      message("`likelihood` not given by user - set to `normal` based on data provided")
    } else if (all(c("r", "E") %in% names(data.ab))) {
      likelihood <- "poisson"
      message("`likelihood` not given by user - set to `poisson` based on data provided")
    }
  }
  if (is.null(link)) {
    if (likelihood=="binomial") {
      link <- "logit"
      message("`link` not given by user - set to `logit` based on assigned value for `likelihood`")
    } else if (likelihood=="normal") {
      link <- "identity"
      message("`link` not given by user - set to `identity` based on assigned value for `likelihood`")
    } else if (likelihood=="poisson") {
      link <- "log"
      message("`link` not given by user - set to `log` based on assigned value for `likelihood`")
    }
  }

  # Check valid likelihood is used
  if (likelihood=="binomial" & !all(c("r", "N") %in% names(data.ab))) {
    stop("Binomial likelihood - columns `r` and `N` must be included in `data.ab`")
  } else if (likelihood=="poisson" & !all(c("E", "N") %in% names(data.ab))) {
    stop("Poisson likelihood - columns `E` and `N` must be included in `data.ab`")
  } else if (likelihood=="normal" & !all(c("y", "se") %in% names(data.ab))) {
    stop("Normal likelihood - columns `y` and `se` must be included in `data.ab`")
  }

  return(list("likelihood"=likelihood, "link"=link))
}










#######################################################
#########   MBNMA.run Wrapper Functions   #############
#######################################################


#' Run MBNMA model with a linear dose-response function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' dose-response function. Follows the methods
#' of MAWDSLEY REF. This function acts as a wrapper for `MBNMA.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams MBNMA.run
#' @inherit MBNMA.run return references
#' @param slope Refers to the slope parameter of the linear dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection MBNMA.run Dose-response parameters
#'
#' @examples
#'
#' @export
MBNMA.linear <- function(network, parameters.to.save=NULL,
                         slope="rel",
                         method="common",
                         class.effect=list(),
                         var.scale=NULL,
                         pd="pv", parallel=TRUE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("slope"),
    run.params=c("beta.1")
  )

  result <- MBNMA.run(network=network, parameters.to.save=parameters.to.save,
                      fun="linear", user.fun=NULL,
                      model.file=NULL,
                      beta.1=slope,
                      method=method,
                      class.effect=class.effect,
                      var.scale=var.scale,
                      pd=pd, parallel=parallel,
                      likelihood=likelihood, link=link,
                      priors=priors,
                      arg.params=arg.params, ...)

  return(result)
}





#' Run MBNMA model with a exponential dose-response function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' dose-response function. Follows the methods
#' of MAWDSLEY REF. This function acts as a wrapper for `MBNMA.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams MBNMA.run
#' @inherit MBNMA.run return references
#' @param lambda Refers to the rate of growth/decay of the exponential dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection MBNMA.run Dose-response parameters
#'
#' @examples
#'
#' @export
MBNMA.exponential <- function(network, parameters.to.save=NULL,
                         lambda="rel",
                         method="common",
                         class.effect=list(),
                         var.scale=NULL,
                         pd="pv", parallel=TRUE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("lambda"),
    run.params=c("beta.1")
  )

  result <- MBNMA.run(network=network, parameters.to.save=parameters.to.save,
                      fun="linear", user.fun=NULL,
                      model.file=NULL,
                      beta.1=lambda,
                      method=method,
                      class.effect=class.effect,
                      var.scale=var.scale,
                      pd=pd, parallel=parallel,
                      likelihood=likelihood, link=link,
                      priors=priors,
                      arg.params=arg.params, ...)

  return(result)
}






#' Run MBNMA model with an Emax dose-response function (without Hill parameter)
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' dose-response function. Follows the methods
#' of MAWDSLEY REF. This function acts as a wrapper for `MBNMA.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams MBNMA.run
#' @inherit MBNMA.run return references
#' @param emax Refers to the Emax parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param ed50 Refers to the ED50 parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection MBNMA.run Dose-response parameters
#'
#' @examples
#'
#' @export
MBNMA.emax <- function(network, parameters.to.save=NULL,
                         emax="rel",
                         ed50="rel",
                         method="common",
                         class.effect=list(),
                         var.scale=NULL,
                         pd="pv", parallel=TRUE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("emax", "ed50"),
    run.params=c("beta.1", "beta.2")
  )

  result <- MBNMA.run(network=network, parameters.to.save=parameters.to.save,
                      fun="emax", user.fun=NULL,
                      model.file=NULL,
                      beta.1=emax,
                      beta.2=ed50,
                      method=method,
                      class.effect=class.effect,
                      var.scale=var.scale,
                      pd=pd, parallel=parallel,
                      likelihood=likelihood, link=link,
                      priors=priors,
                      arg.params=arg.params, ...)

  return(result)
}






#' Run MBNMA model with an Emax dose-response function (with a Hill parameter)
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' dose-response function. Follows the methods
#' of MAWDSLEY REF. This function acts as a wrapper for `MBNMA.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams MBNMA.run
#' @inherit MBNMA.run return references
#' @param emax Refers to the Emax parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param ed50 Refers to the ED50 parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param hill Refers to the Hill parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection MBNMA.run Dose-response parameters
#'
#' @examples
#'
#' @export
MBNMA.emax.hill <- function(network, parameters.to.save=NULL,
                       emax="rel",
                       ed50="rel",
                       hill="common",
                       method="common",
                       class.effect=list(),
                       var.scale=NULL,
                       pd="pv", parallel=TRUE,
                       likelihood=NULL, link=NULL,
                       priors=NULL,
                       arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("emax", "ed50", "hill"),
    run.params=c("beta.1", "beta.2", "beta.3")
  )

  result <- MBNMA.run(network=network, parameters.to.save=parameters.to.save,
                      fun="emax.hill", user.fun=NULL,
                      model.file=NULL,
                      beta.1=emax,
                      beta.2=ed50,
                      beta.3=hill,
                      method=method,
                      class.effect=class.effect,
                      var.scale=var.scale,
                      pd=pd, parallel=parallel,
                      likelihood=likelihood, link=link,
                      priors=priors,
                      arg.params=arg.params, ...)

  return(result)
}
