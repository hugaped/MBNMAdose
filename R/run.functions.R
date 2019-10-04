# Functions for running MBNMA models
# Author: Hugo Pedder
# Date created: 2019-04-18

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Run MBNMA dose-response models
#'
#' Fits a Bayesian dose-response for model-based network meta-analysis
#' (MBNMA) that can account for multiple doses of different agents by
#' applying a desired dose-response function. Follows the methods
#' of \insertCite{mawdsley2016;textual}{MBNMAdose}.
#'
#' @param network An object of class `mbnma.network`.
#' @param parameters.to.save A character vector containing names of parameters
#'   to monitor in JAGS
#' @param fun A string specifying a functional form to be assigned to the
#'   dose-response. Options are given in `details`.
#' @param user.fun A string specifying any relationship including `dose` and
#'   one/several of: `beta.1`, `beta.2`, `beta.3`.
#' @param model.file A JAGS model written as a character object that can be used
#'   to overwrite the JAGS model that is automatically written based on the
#'   specified options.
#'
#' @param beta.1 Refers to dose-parameter(s) specified within the dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.2 Refers to dose-parameter(s) specified within the dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.3 Refers to dose-parameter(s) specified within the dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @param method Can take either `"common"` or `"random"` to indicate whether relative effects
#'   should be modelled with between-study heterogeneity or not (see details).
#' @param class.effect A list of named strings that determines which dose-response
#'   parameters to model with a class effect and what that effect should be
#'   (`"common"` or `"random"`). Element names should match dose-response parameter names (which will therefore
#'   depend on whether or not a wrapper function has been used for `mbnma.run()`).
#'   For example: `list("beta.2"="fixed", "beta.3"="random")` when using
#'   `mbnma.run()` or `list("ed50"="fixed", "hill"="random")` when using `mbnma.emax.hill()`.
#' @param likelihood A string indicating the likelihood to use in the model. Can take either `"binomial"`,
#'   `"normal"` or `"poisson"`. If left as `NULL` the likelihood will be inferred from the data.
#' @param link A string indicating the link function to use in the model. Can take any link function
#'   defined within JAGS (e.g. `"logit"`, `"log"`, `"probit"`, `"cloglog"`) or be assigned the value `"identity"` for
#'   and identity link function. If left as `NULL` the link function will be automatically assigned based
#'   on the likelihood.
#' @param cor A boolean object that indicates whether correlation should be modelled
#' between relative effect dose-response parameters (`TRUE`) or not (`FALSE`). This is
#' automatically set to `FALSE` if class effects are modelled.
#' @param var.scale A numeric vector indicating the relative scale of variances between
#' correlated dose-response parameters when relative effects are modelled on more than
#' one dose-response parameter and `cor=TRUE` (see details). Each element of
#' the vector refers to the relative scale of each of the dose-response parameters that is
#' modelled using relative effects.
#' @param priors A named list of parameter values (without indices) and
#'   replacement prior distribution values given as strings
#'   **using distributions as specified in JAGS syntax** (see examples).
#'
#' @param pd Can take either:
#'   * `pv` only pV will be reported (as automatically outputted by `R2jags`).
#'   * `plugin` calculates pD by the plug-in
#'   method \insertCite{spiegelhalter2002}{MBNMAdose}. It is faster, but may output negative
#'   non-sensical values, due to skewed deviances that can arise with non-linear models.
#'   * `pd.kl` calculates pD by the Kullback-Leibler divergence \insertCite{plummer2008}{MBNMAdose}. This
#'   will require running the model for additional iterations but
#'   will always produce a positive result.
#'   * `popt` calculates pD using an optimism adjustment which allows for calculation
#'   of the penalized expected deviance \insertCite{plummer2008}{MBNMAdose}
#' @param parallel A boolean value that indicates whether JAGS should be run in
#'   parallel (`TRUE`) or not (`FALSE`). If `TRUE` then the number of cores to
#'   use is automatically calculated.
#' @param arg.params Contains a list of arguments sent to `mbnma.run()` by dose-response
#' specific wrapper functions
#' @param n.iter number of total iterations per chain (including burn in; default: 15000)
#' @param n.thin thinning rate. Must be a positive integer. Set `n.thin > 1`` to save memory
#' and computation time if n.iter is large. Default is
#' `max(1, floor(n.chains * (n.iter-n.burnin) / 1000))`` which will only thin if there are at least 2000
#' simulations.
#' @param n.chains number of Markov chains (default: 3)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the
#' beginning. Default is `n.iter/2``, that is, discarding the first half of the
#' simulations. If n.burnin is 0, jags() will run 100 iterations for adaption.
#' @param ... Arguments to be sent to R2jags.
#'
#'
#' @details When relative effects are modelled on more than one dose-response parameter and
#' `cor = TRUE`, correlation between the dose-response parameters is automatically
#' estimated using a vague Wishart prior. This prior can be made slightly more informative
#' by specifying the relative scale of variances between the dose-response parameters using
#' `var.scale`.
#'
#' @return An object of S3 `class(c("mbnma", "rjags"))` containing parameter
#'   results from the model. Can be summarized by `print()` and can check
#'   traceplots using `R2jags::traceplot()` or various functions from the package `mcmcplots`.
#'
#'   Nodes that are automatically monitored (if present in the model) have the
#'   following interpretation. These will have an additional suffix that relates
#'   to the name/number of the dose-response parameter to which they correspond
#'   (e.g. `d.ed50` or `d.1`):
#'   * `d` The pooled effect for each agent for a given dose-response
#'   parameter. Will be estimated by the model if dose-response parameters (`beta.1`,
#'   `beta.2`, `beta.3`) are set to `"rel"`.
#'   * `sd` (without a suffix) - the between-study SD (heterogeneity) for relative effects, reported if
#'   `method="random"`.
#'   * `D` The class effect for each class for a given dose-response
#'   parameter. Will be estimated by the model if specified in `class.effect`.
#'   * `sd.D` The within-class SD for different agents within the same class. Will
#'   be estimated by the model if any dose-response parameter in `class.effect` is
#'   set to `"random"`.
#'   * `beta` The absolute value of a given dose-response parameter across the whole
#'   network (does not vary by agent/class). Will be estimated by the model if
#'   dose-response parameters (`beta.1`, `beta.2`, `beta.3`) are set to `"common"`
#'   or `"random"`.
#'   * `sd` (with a suffix) - the between-study SD (heterogeneity) for absolute dose-response
#'   parameters, reported if `beta.1`, `beta.2` or `beta.3` are set to `"random"`
#'   * `totresdev` The residual deviance of the model
#'   * `deviance` The deviance of the model
#'
#'
#'   If there are errors in the JAGS model code then the object will be a list
#'   consisting of two elements - an error message from JAGS that can help with
#'   debugging and `model.arg`, a list of arguments provided to `mbnma.run()`
#'   which includes `jagscode`, the JAGS code for the model that can help
#'   users identify the source of the error.
#'
#' @section Dose-response parameters:
#' * `"rel"` implies that relative effects should be pooled for this dose-response
#' parameter, that vary by agent.
#' * `"common"` implies that all studies estimate the same true absolute effect
#' (akin to a "fixed effects" meta-analysis) across the whole network
#' * `"random"` implies that all studies estimate a separate true absolute effect, but
#' that each of these true effects vary randomly around a true mean effect. This
#' approach allows for modelling of between-study heterogeneity.
#' * `numeric()` Assigned a numeric value. It indicates that
#' this dose-response parameter should not be estimated from the data but should be
#' assigned the numeric value determined by the user. This can be useful for fixing
#' specific dose-response parameters (e.g. Hill parameters in Emax functions) to a value.
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
#'   the method of \insertCite{owen2015;textual}{MBNMAdose})
#'   * `"nonparam.down"` (monotonically decreasing non-parametric dose-response relationship following
#'   the method of \insertCite{owen2015;textual}{MBNMAdose})
#'   * `"user"` (user-defined function: `user.fun` must be specified in arguments)
#'
#' @importFrom Rdpack reprompt
#' @importFrom magrittr "%>%"
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#'
#' ######## Dose-response functions ########
#'
#' # Fit a dose-response MBNMA with a linear function and common treatment effects
#' result <- mbnma.run(network, fun="linear", beta.1="rel", method="common")
#'
#' # Fit a dose-response MBNMA with an exponential function and random treatment effects
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random")
#'
#' # Fit a user-defined function (quadratic)
#' fun.def <- "(beta.1 * dose) + (beta.2 * (dose^2))"
#' result <- mbnma.run(network, fun="user", user.fun=fun.def,
#'               beta.1="rel", beta.2="rel", method="common")
#'
#' # Fit an Emax function with a single random (exchangeable) parameter estimated
#' #for ED50 and common treatment effects on relative Emax effects
#' result <- mbnma.run(network, fun="emax",
#'               beta.1="rel", beta.2="random", method="common")
#'
#' # Fit an Emax function with a Hill parameter, with a fixed value for the Hill parameter
#' #provided to the model and random relative effects on Emax and ED50 (which will
#' #therefore be modelled with a correlation between them).
#' result <- mbnma.run(network, fun="emax.hill",
#'               beta.1="rel", beta.2="rel", beta.3=5, method="random")
#'
#'
#' ########## Class effects ##########
#'
#' # Generate a dataset with one class for active treatments and one for placebo
#' class.df <- HF2PPITT
#' class.df$class <- ifelse(class.df$agent=="placebo", "placebo", "active")
#' netclass <- mbnma.network(class.df)
#'
#' # Fit an Emax function with common relative effects on Emax and ED50 and
#' #a random class effect on ED50.
#' result <- mbnma.run(netclass, fun="emax",
#'               beta.1="rel", beta.2="rel", method="common",
#'               class.effect=list(beta.2="random"))
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from an Emax function with random relative effects on Emax and ED50
#' result <- mbnma.run(network, fun="emax",
#'               beta.1="rel", beta.2="rel", method="random")
#' print(result$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)",
#'                  inv.R = "dwish(Omega[,],100)")
#'
#' result <- mbnma.run(network, fun="emax",
#'               beta.1="rel", beta.2="rel", method="random",
#'               priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
#'               n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
#'               pd="plugin")
#'
#' # Calculate effective number of parameters via Kullback-Leibler method
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
#'               pd="pd.kl")
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(result)
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "d.1")
#'
#'
#' ########## Output ###########
#'
#' # Print R2jags output and summary
#' print(result)
#' summary(result)
#'
#' # Plot forest plot of results
#' plot(result)
#' }
#'
#' @export
mbnma.run <- function(network,
                      fun="linear",
                      beta.1="rel",
                      beta.2=NULL, beta.3=NULL,
                      method="common",
                      class.effect=list(),
                      cor=TRUE,
                      var.scale=NULL,
                      user.fun=NULL,
                      parameters.to.save=NULL,
                      pd="pv", parallel=TRUE,
                      likelihood=NULL, link=NULL,
                      priors=NULL,
                      model.file=NULL,
                      n.iter=10000, n.chains=3,
                      n.burnin=floor(n.iter/2), n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                      arg.params=NULL, ...
) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertCharacter(model.file, len=1, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(pd, choices=c("pv", "pd.kl", "plugin", "popt"), null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(cor, len=1, add=argcheck)
  checkmate::assertList(arg.params, unique=TRUE, null.ok=TRUE, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  # Ensure rjags parameters make sense
  if (n.iter<=n.burnin) {
    stop(paste0("`n.iter` must be greater than `n.burnin`. `n.burnin` = ", n.burnin))
  }

  # Check if placebo has been included
  if (network$agents[1]=="Placebo" & network$treatments[1]=="Placebo_0") {
    plac.incl <- TRUE
  } else {
    plac.incl <- FALSE
    if (fun %in% c("nonparam.up", "nonparam.down")) {
      stop("Placebo (or an agent with dose=0) must be included in the network to model a nonparametric dose-response relationship")
    }
  }


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
    model <- mbnma.write(fun=fun, user.fun=user.fun,
                         beta.1=beta.1, beta.2=beta.2, beta.3=beta.3,
                         method=method,
                         class.effect=class.effect,
                         cor=cor, var.scale=var.scale,
                         likelihood=likelihood, link=link
    )

    # Change code for if plac not included in network
    if (plac.incl==FALSE) {
      model <- gsub("\\\nfor \\(k in 2:Nagent\\)\\{ # Priors on relative treatment effects\\\n",
                    "for (k in 1:Nagent){ # Priors on relative treatment effects\n",
                    model)
      model <- gsub("\\\nfor \\(k in 2:Nclass\\)\\{ # Priors on relative class effects\\\n",
                    "for (k in 1:Nclass){ # Priors on relative class effects\n",
                    model)

      model <- gsub("s\\.beta\\.[1-3]\\[1\\] <- 0", "", model)
    }

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
    pluginvars <- c("psi", "resdev")
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

  result.jags <- mbnma.jags(data.ab, model,
                            class=class,
                            parameters.to.save=parameters.to.save,
                            likelihood=likelihood, link=link,
                            n.iter=n.iter,
                            n.thin=n.thin,
                            n.chains=n.chains,
                            n.burnin=n.burnin,
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

  } else if (pd == "plugin") {
    # plugin method
    if (likelihood=="normal") {
      obs1 <- jagsdata[["y"]]
      obs2 <- jagsdata[["se"]]
    } else if (likelihood=="binomial") {
      obs1 <- jagsdata[["r"]]
      obs2 <- jagsdata[["N"]]
    } else if (likelihood=="poisson") {
      obs1 <- jagsdata[["r"]]
      obs2 <- jagsdata[["E"]]
    }
    result$BUGSoutput$pD <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                                   theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                                   likelihood=likelihood, type="dose")
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
                    "cor"=cor,
                    "var.scale"=var.scale,
                    "parallel"=parallel, "pd"=pd,
                    "priors"=get.prior(model), "arg.params"=arg.params)
  result[["model.arg"]] <- model.arg
  result[["type"]] <- "dose"
  result[["agents"]] <- network[["agents"]]
  result[["treatments"]] <- network[["treatments"]]
  if (length(class.effect)>0) {
    result[["classes"]] <- network[["classes"]]
  }

  if (!("error" %in% names(result))) {
    class(result) <- c("mbnma", class(result))
  }

  return(result)

}





mbnma.jags <- function(data.ab, model,
                       class=FALSE,
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
    #jagsdata <- getjagsdata(data.ab, class=class, rho=rho, covstruct=covar) # get data into jags correct format (list("fups", "NT", "NS", "narm", "y", "se", "treat", "time"))
    jagsdata <- getjagsdata(data.ab, class=class) # get data into jags correct format (list("fups", "NT", "NS", "narm", "y", "se", "treat", "time"))
  } else {
    # For MBNMAdose
    jagsdata <- getjagsdata(data.ab, class=class,
                            likelihood=likelihood, link=link) # get data into jags correct format
  }


  # Add variable for maxtime to jagsdata if required
  if (grepl("maxtime", model)) {
    jagsdata[["maxtime"]] <- max(data.ab$time)
  } else if (grepl("maxdose", model)) {
    #jagsdata[["maxdose"]] <- index.dose(network[["data.ab"]])[["maxdose"]]
    jagsdata[["maxdose"]] <- index.dose(data.ab)[["maxdose"]]
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
#' Identical to `gen.parameters.to.save()` in `MBNMAtime`
#'
#' @param model.params A character or numeric vector containing the names of the
#' dose-response parameters in the model
#' @param model A JAGS model written as a character object
#' @noRd
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
      parameters.to.save <- append(parameters.to.save, paste0("sd.beta.", model.params[i]))
    }
    if (grepl(paste0("\\\nD\\.", model.params[i], "(\\[k\\])? ~"), model)==TRUE) {
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
#' @inheritParams mbnma.run
#' @param warn.rhat A boolean object to indicate whether to return a warning if Rhat values
#' for any monitored parameter are >1.02 (suggestive of non-convergence).
#' @param drop.discon A boolean object that indicates whether or not to drop disconnected
#'   studies from the network.
#' @param UME A boolean object to indicate whether to fit an Unrelated Mean Effects model
#'   that does not assume consistency and so can be used to test if the consistency
#'   assumption is valid.
#' @param n.iter number of total iterations per chain (including burn in; default: 10000)
#'
#' @examples
#' \donttest{
#' # Run random effects NMA on the alogliptin dataset
#' network <- mbnma.network(alog_pcfb)
#' nma <- nma.run(network, method="random")
#' print(nma)
#' plot(nma)
#'
#' # Run common effects NMA keeping treatments that are disconnected in the NMA
#' network <- mbnma.network(GoutSUA_2wkCFB)
#' nma <- nma.run(network, method="common", drop.discon=FALSE)
#'
#' # Run an Unrelated Mean Effects (UME) inconsistency model on triptans dataset
#' network <- mbnma.network(HF2PPITT)
#' ume <- nma.run(network, method="random", UME=TRUE)
#' }
#'
#' @export
nma.run <- function(network, method="common", likelihood=NULL, link=NULL,
                    warn.rhat=TRUE, n.iter=10000, drop.discon=TRUE, UME=FALSE, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
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
  model <- write.nma(method=method, likelihood=likelihood, link=link, UME=UME)


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
  class(output) <- "nma"
  return(output)

}







#' Check likelihood and link function
#'
#' Checks that likelihood and link function is provided and confirm that the correct
#' form of data is provided.
#'
#' @inheritParams mbnma.run
#' @inheritParams mbnma.network
#'
#' @export
check.likelink <- function(data.ab, likelihood=NULL, link=NULL) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit", "log"), null.ok=TRUE, add=argcheck)
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
#########   mbnma.run Wrapper Functions   #############
#######################################################


#' Run MBNMA model with a linear dose-response function
#'
#' Fits a Bayesian model-based network meta-analysis (MBNMA) with a defined
#' dose-response function. Follows the methods
#' of \insertCite{mawdsley2016;textual}{MBNMAdose}. This function acts as a wrapper for `mbnma.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams mbnma.run
#' @inherit mbnma.run return references
#' @param slope Refers to the slope parameter of the linear dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection mbnma.run Dose-response parameters
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Fit a linear dose-response MBNMA with random treatment effects
#' linear <- mbnma.linear(network, slope="rel", method="random")
#'
#' # Fit a linear dose-response MBNMA using a cloglog link function
#' linear <- mbnma.linear(network, slope="rel", link="cloglog")
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from linear dose-response MBNMA
#' linear <- mbnma.linear(network, slope="rel", method="random")
#' print(linear$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)")
#'
#' linear <- mbnma.linear(network, slope="rel", method="random",
#'               priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' linear <- mbnma.linear(network, slope="rel", method="random",
#'               n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' linear <- mbnma.linear(network, slope="rel", method="random",
#'               pd="plugin")
#'
#' # Calculate effective number of parameters via Kullback-Leibler method
#' linear <- mbnma.linear(network, slope="rel", method="random",
#'               pd="pd.kl")
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(linear)
#'
#' # Traceplots
#' mcmcplots::traplot(linear)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(linear, "d.slope")
#'
#'
#' ########## Output ###########
#'
#' # Print R2jags output and summary
#' print(linear)
#' summary(linear)
#'
#' # Plot forest plot of results
#' plot(linear)
#' }
#'
#' @export
mbnma.linear <- function(network,
                         slope="rel",
                         method="common",
                         class.effect=list(),
                         cor=TRUE,
                         var.scale=NULL,
                         parameters.to.save=NULL,
                         pd="pv", parallel=TRUE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("slope"),
    run.params=c("beta.1")
  )

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="linear", user.fun=NULL,
                      model.file=NULL,
                      beta.1=slope,
                      method=method,
                      class.effect=class.effect,
                      cor=cor, var.scale=var.scale,
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
#' of \insertCite{mawdsley2016;textual}{MBNMAdose}. This function acts as a wrapper for `mbnma.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams mbnma.run
#' @inherit mbnma.run return references
#' @param lambda Refers to the rate of growth/decay of the exponential dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection mbnma.run Dose-response parameters
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Fit a exponential dose-response MBNMA with random treatment effects
#' exponential <- mbnma.exponential(network, lambda="rel", method="random")
#'
#' # Fit a exponential dose-response MBNMA using a cloglog link function
#' exponential <- mbnma.exponential(network, lambda="rel", link="cloglog")
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from exponential dose-response MBNMA
#' exponential <- mbnma.exponential(network, lambda="rel", method="random")
#' print(exponential$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)")
#'
#' exponential <- mbnma.exponential(network, lambda="rel", method="random",
#'                    priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' exponential <- mbnma.exponential(network, lambda="rel", method="random",
#'                    n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' exponential <- mbnma.exponential(network, lambda="rel", method="random",
#'                    pd="plugin")
#'
#' # Calculate effective number of parameters via Kullback-Leibler method
#' exponential <- mbnma.exponential(network, lambda="rel", method="random",
#'                    pd="pd.kl")
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(exponential)
#'
#' # Traceplots
#' mcmcplots::traplot(exponential)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(exponential, "d.lambda")
#'
#'
#' ########## Output ###########
#'
#' # Print R2jags output and summary
#' print(exponential)
#' summary(exponential)
#'
#' # Plot forest plot of results
#' plot(exponential)
#' }
#'
#' @export
mbnma.exponential <- function(network,
                         lambda="rel",
                         method="common",
                         class.effect=list(),
                         cor=TRUE,
                         var.scale=NULL,
                         parameters.to.save=NULL,
                         pd="pv", parallel=TRUE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("lambda"),
    run.params=c("beta.1")
  )

  # Increase precision of prior for beta.1 to prevent errors
  # if (is.null(priors)) {
  #   maxdose <- max(network$data.ab$dose)
  #   lim <- 700/maxdose
  #   prec <- 1/((lim/2)^2)
  #   dist <- gsub("prec", signif(prec,5), "dnorm(0,prec)")
  #   priors <- list("d.lambda"=dist)
  # }

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="exponential", user.fun=NULL,
                      model.file=NULL,
                      beta.1=lambda,
                      method=method,
                      class.effect=class.effect,
                      cor=cor, var.scale=var.scale,
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
#' of \insertCite{mawdsley2016;textual}{MBNMAdose}. This function acts as a wrapper for `mbnma.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams mbnma.run
#' @inherit mbnma.run return references
#' @param emax Refers to the Emax parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param ed50 Refers to the ED50 parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection mbnma.run Dose-response parameters
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Fit an Emax dose-response MBNMA with random treatment effects on Emax and ED50
#' emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#'
#' # Fit an Emax dose-response MBNMA with common treatment effects on Emax and
#' #a single common parameter estimated for ED50
#' emax <- mbnma.emax(network, emax="rel", ed50="common", method="common")
#'
#'
#' ########## Class effects ##########
#'
#' # Generate a dataset with one class for active treatments and one for placebo
#' class.df <- HF2PPITT
#' class.df$class <- ifelse(class.df$agent=="placebo", "placebo", "active")
#' netclass <- mbnma.network(class.df)
#'
#' # Fit an Emax function with common relative effects on Emax and ED50 and
#' #a random class effect on ED50.
#' emax <- mbnma.emax(netclass,
#'             emax="rel", ed50="rel", method="common",
#'             class.effect=list(ed50="random"))
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from an Emax function with random relative effects on Emax and ED50
#' emax <- mbnma.emax(network,
#'             emax="rel", ed50="rel", method="random")
#' print(emax$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)",
#'                  inv.R = "dwish(Omega[,],100)")
#'
#' emax <- mbnma.emax(network,
#'             emax="rel", ed50="rel", method="random",
#'             priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' emax <- mbnma.emax(network, emax="rel", ed50="rel",
#'             n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' emax <- mbnma.emax(network, emax="rel", ed50="rel",
#'             pd="plugin")
#'
#' # Calculate effective number of parameters via Kullback-Leibler method
#' emax <- mbnma.emax(network, emax="rel", ed50="rel",
#'             pd="pd.kl")
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(emax)
#'
#' # Traceplots
#' mcmcplots::traplot(emax)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(emax, "d.emax")
#'
#'
#' ########## Output ###########
#'
#' # Print R2jags output and summary
#' print(emax)
#' summary(emax)
#'
#' # Plot forest plot of results
#' plot(emax)
#' }
#'
#' @export
mbnma.emax <- function(network,
                         emax="rel",
                         ed50="rel",
                         method="common",
                         class.effect=list(),
                         cor=TRUE,
                         var.scale=NULL,
                         parameters.to.save=NULL,
                         pd="pv", parallel=TRUE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("emax", "ed50"),
    run.params=c("beta.1", "beta.2")
  )

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="emax", user.fun=NULL,
                      model.file=NULL,
                      beta.1=emax,
                      beta.2=ed50,
                      method=method,
                      class.effect=class.effect,
                      cor=cor, var.scale=var.scale,
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
#' of \insertCite{mawdsley2016;textual}{MBNMAdose}. This function acts as a wrapper for `mbnma.run()` that
#' uses more clearly defined parameter names.
#'
#' @inheritParams mbnma.run
#' @inherit mbnma.run return references
#' @param emax Refers to the Emax parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param ed50 Refers to the ED50 parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param hill Refers to the Hill parameter of the Emax dose-response function.
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#'
#' @inheritSection mbnma.run Dose-response parameters
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Fit an Emax (with Hill parameter) dose-response MBNMA with random treatment
#' #effects on Emax, ED50 and Hill
#' emax.hill <- mbnma.emax.hill(network, emax="rel", ed50="rel", hill="rel",
#'                  method="random")
#'
#' # Fit an Emax (with Hill parameter) dose-response MBNMA with common treatment
#' #effects on Emax, a single random parameter estimated for ED50
#' #and a single common parameter estimated for Hill
#' emax.hill <- mbnma.emax.hill(network, emax="rel", ed50="random", hill="common",
#'                  method="common")
#'
#' # Assign a specific numerical value for Hill parameter
#' emax.hill <- mbnma.emax.hill(network, emax="rel", ed50="rel", hill=5)
#'
#'
#' ########## Class effects ##########
#'
#' # Generate a dataset with one class for active treatments and one for placebo
#' class.df <- HF2PPITT
#' class.df$class <- ifelse(class.df$agent=="placebo", "placebo", "active")
#' netclass <- mbnma.network(class.df)
#'
#' # Fit an Emax (with Hill parameter) function with common relative effects on
#' #all parameters and common class effects on ED50 and Hill.
#' emax.hill <- mbnma.emax.hill(netclass,
#'                  emax="rel", ed50="rel", hill="rel", method="common",
#'                  class.effect=list(ed50="common", hill="common"))
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from an Emax (with Hill parameter) function with
#' #relative effects on Emax and ED50 and a single common parameter for Hill
#' emax.hill <- mbnma.emax.hill(network,
#'                  emax="rel", ed50="rel", hill="common", method="common")
#' print(emax.hill$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(beta.hill = "dnorm(0,0.5) T(,0)")
#'
#' emax.hill <- mbnma.emax.hill(network,
#'                  emax="rel", ed50="rel", hill="common", method="common",
#'                  priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' emax.hill <- mbnma.emax.hill(network, emax="rel", ed50="rel", hill=2,
#'                  n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' emax.hill <- mbnma.emax.hill(network, emax="rel", ed50="rel", hill=2,
#'                  pd="plugin")
#'
#' # Calculate effective number of parameters via Kullback-Leibler method
#' emax.hill <- mbnma.emax.hill(network, emax="rel", ed50="rel", hill=2,
#'                  pd="pd.kl")
#'
#'
#' ####### Examine MCMC diagnostics (using mcmcplots package) #######
#'
#' # Density plots
#' mcmcplots::denplot(emax.hill)
#'
#' # Traceplots
#' mcmcplots::traplot(emax.hill)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(emax.hill, "d.emax")
#'
#'
#' ########## Output ###########
#'
#' # Print R2jags output and summary
#' print(emax.hill)
#' summary(emax.hill)
#'
#' # Plot forest plot of results
#' plot(emax.hill)
#' }
#'
#' @export
mbnma.emax.hill <- function(network,
                       emax="rel",
                       ed50="rel",
                       hill="common",
                       method="common",
                       class.effect=list(),
                       cor=TRUE,
                       var.scale=NULL,
                       parameters.to.save=NULL,
                       pd="pv", parallel=TRUE,
                       likelihood=NULL, link=NULL,
                       priors=NULL,
                       arg.params=NULL, ...)
{

  arg.params <- list(
    wrap.params=c("emax", "ed50", "hill"),
    run.params=c("beta.1", "beta.2", "beta.3")
  )

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="emax.hill", user.fun=NULL,
                      model.file=NULL,
                      beta.1=emax,
                      beta.2=ed50,
                      beta.3=hill,
                      method=method,
                      class.effect=class.effect,
                      cor=cor, var.scale=var.scale,
                      pd=pd, parallel=parallel,
                      likelihood=likelihood, link=link,
                      priors=priors,
                      arg.params=arg.params, ...)

  return(result)
}





#' Calculate plugin pD from a JAGS model with univariate likelihood for studies
#' with repeated measurements
#'
#' Uses results from MBNMA JAGS models to calculate pD via the
#' plugin method \insertCite{spiegelhalter2002}{MBNMAdose}. Can only be used for models with known
#' standard errors or covariance matrices (typically univariate likelihoods).
#'
#' @param obs1 A matrix (study x arm) or array (study x arm x time point) containing
#'   observed data for `y` (normal likelihood) or `r` (binomial or poisson likelihood)
#'   in each arm of each study. This will be the same array
#'   used as data for the JAGS model.
#' @param obs2 A matrix (study x arm) or array (study x arm x time point) containing
#'   observed data for `se` (normal likelihood), `N` (binomial likelihood) or `E` (poisson likelihood)
#'   in each arm of each study. This will be the same array
#'   used as data for the JAGS model.
#' @param fups A numeric vector of length equal to the number of studies,
#'   containing the number of follow-up mean responses reported in each study. Required for
#'   time-course MBNMA models (if `type="time"`)
#' @param narm A numeric vector of length equal to the number of studies,
#'   containing the number of arms in each study.
#' @param NS A single number equal to the number of studies in the dataset.
#' @param theta.result A matrix (study x arm) or array (study x arm x time point)
#'   containing the posterior mean predicted means/probabilities/rate in each arm of each
#'   study. This will be estimated by the JAGS model.
#' @param resdev.result A matrix (study x arm) or array (study x arm x time point)
#'   containing the posterior mean residual deviance contributions in each arm of each
#'   study. This will be estimated by the JAGS model.
#'
#' @param likelihood A character object of any of the following likelihoods:
#' * `normal`
#' * `binomial` (does not work with time-course MBNMA models)
#' * `poisson` (does not work with time-course MBNMA models)
#' @param type The type of MBNMA model fitted. Can be either `"time"` or `"dose"`
#'
#' @return A single numeric value for pD calculated via the plugin method.
#'
#' @details Method for calculating pD via the plugin method proposed by
#'   Spiegelhalter \insertCite{spiegelhalter2002}{MBNMAdose}. Standard errors / covariance matrices must be assumed
#'   to be known. To obtain values for `theta.result` and `resdev.result` these
#'   parameters must be monitored when running the MBNMA model (using `parameters.to.save`).
#'
#'   For non-linear time-course MBNMA models residual deviance contributions may be skewed, which
#'   can lead to non-sensical results when calculating pD via the plugin method.
#'   Alternative approaches are to use pV as an approximation or
#'   pD calculated by Kullback-Leibler divergence \insertCite{plummer2008}{MBNMAdose}.
#'
#' @references
#'   \insertAllCited
#'
#' @inherit mbnma.run references
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Fit a dose-response MBNMA, monitoring "psi" and "resdev"
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
#'               parameters.to.save=c("psi", "resdev"))
#'
#'
#' #### Calculate pD for binomial data ####
#'
#' # Prepare data for pD calculation
#' r <- result$model$data()$r
#' N <- result$model$data()$N
#' narm <- result$model$data()$narm
#' NS <- result$model$data()$NS
#'
#' psi <- result$BUGSoutput$median$psi
#' resdevs <- result$BUGSoutput$median$resdev
#'
#' # Calculate pD via plugin method
#' pD <- pDcalc(obs1=r, obs2=N, narm=narm, NS=NS,
#'           theta.result=psi, resdev.result=resdevs,
#'           likelihood="binomial", type="dose")
#' }
#'
#' @export
pDcalc <- function(obs1, obs2, fups=NULL, narm, NS, theta.result, resdev.result,
                   likelihood="normal", type="time") {
  # For univariate models only!!

  # likelihood could in theory be c("normal", "multivar.normal", "binomial")
  # theta.result = model$BUGSoutput$mean$theta
  # resdev.result = model$BUGSoutput$mean$totresdev

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertArray(obs1, add=argcheck)
  checkmate::assertArray(obs2, add=argcheck)
  checkmate::assertArray(theta.result, add=argcheck)
  checkmate::assertArray(resdev.result, add=argcheck)
  checkmate::assertNumeric(fups, null.ok=TRUE, add=argcheck)
  checkmate::assertNumeric(narm, add=argcheck)
  checkmate::assertNumeric(NS, add=argcheck)
  checkmate::assertChoice(likelihood, choices=c("normal", "binomial", "poisson"), null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(type, choices=c("dose", "time"), null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (type=="time") {
    if (is.null(fups)) {
      stop("`fups` cannot be NA in pDcalc for time-course MBNMA")
    }
    dev.post <- array(dim=c(NS,max(narm),max(fups)))
    pD <- array(dim=c(NS,max(narm),max(fups)))
  } else if (type=="dose") {
    dev.post <- matrix(nrow=NS, ncol=max(narm))
    pD <- matrix(nrow=NS, ncol=max(narm))
    rhat <- matrix(nrow=NS, ncol=max(narm))
  }

  for (i in 1:NS) {
    for (k in 1:narm[i]) {

      if (type=="time") {
        for (m in 1:fups[i]) {
          # Need to use formula for residual deviance as plugin
          if (likelihood=="normal") {
            dev.post[i,k,m] <- ((obs1[i,k,m] - theta.result[i,k,m])/obs2[i,k,m])^2
            pD[i,k,m] <- resdev.result[i,k,m] - dev.post[i,k,m]
          } else {
            stop("pD cannot be calculated via `plugin` method for time-course MBNMA models without data following a normal likelihood")
          }
        }
      } else if (type=="dose") {
        if (likelihood=="normal") {
          dev.post[i,k] <- ((obs1[i,k] - theta.result[i,k])/obs2[i,k])^2

        } else if (likelihood=="binomial") {
          rhat[i,k] <- theta.result[i,k] * obs2[i,k]
          dev.post[i,k] <- 2*(obs1[i,k] * (log(obs1[i,k])-log(rhat[i,k]))  +
                                (obs2[i,k]-obs1[i,k]) * (log(obs2[i,k]-obs1[i,k]) -
                                                           log(obs2[i,k]-rhat[i,k])))
        } else if (likelihood=="poisson") {
          rhat[i,k] <- theta.result[i,k] * obs2[i,k]
          dev.post[i,k] <- 2*((rhat[i,k]-obs1[i,k]) + (obs1[i,k] * (log(obs1[i,k]/rhat[i,k]))))
        }

        pD[i,k] <- resdev.result[i,k] - dev.post[i,k]

      }

    }
  }

  pD <- sum(pD, na.rm=TRUE)

  return(pD)
}






#' Update MBNMA to monitor deviance nodes in the model
#'
#' Useful for obtaining deviance contributions or fitted values
#'
#' @param mbnma An S3 object of class `"mbnma"` generated by running
#'   a dose-response MBNMA model
#' @param param Used to indicate which node to monitor in the model. Can be any parameter
#' in the model code that varies by all arms within all studies. These are some typical
#' parameters that it might be of interest to monitor, provided they are in the original
#' model code:
#'  * `"theta"` for fitted values
#'  * `"psi"` for fitted values on natural scale (e.g. probabilities)
#'  * `"dev"` for deviance contributions
#'  * `"resdev"` for residual deviance contributions
#'  * `"delta"` for within-study relative effects versus the study reference treatment
#' @inheritParams R2jags::jags
#'
#' @return A data frame containing the posterior mean of the updates by arm and study,
#' with arm and study identifiers.
#'
#' For MBNMAdose:
#' * `facet` indicates the agent identifier in the given arm of a study
#' * `fupdose` indicates the dose in the given arm of a study
#'
#' For MBNMAtime:
#' * `facet` indicates the treatment identifier in the given arm of the study
#' * `fupdose` indicates the follow-up time at the given observation in the given
#' arm of the study
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Fit a dose-response MBNMA, monitoring "psi" and "resdev"
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
#'   parameters.to.save=c("psi", "resdev"))
#'
#' mbnma.update(result, param="theta") # monitor theta
#'
#' mbnma.update(result, param="rhat") # monitor rhat
#'
#' mbnma.update(result, param="delta") # monitor delta
#' }
#'
#' @export
mbnma.update <- function(mbnma, param="theta",
                         n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", add=argcheck)
  checkmate::assertCharacter(param, len = 1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  modelcode <- mbnma$model.arg$jagscode
  # Ensure param is in model code
  if (grepl(paste0("\\\n", param), mbnma$model.arg$jagscode)==FALSE &
      grepl(paste0("\\(", param), mbnma$model.arg$jagscode)==FALSE) {
    stop(paste0(param, " not in model code"))
  }

  # Ensure param varies by study and arm (could add by observation for MBNMAtime)
  str <- sub("(.+)(Run through all arms within a study.+?[}])(.+)", "\\2", modelcode)
  paramstr <- paste0(param, "\\[i,k\\]")
  if (grepl(paramstr, str)==FALSE) {
    stop(paste0(param, " does not vary by study and arm and so array will not be the correct size"))
  }


  result <- rjags::jags.samples(mbnma$model, variable.names = param,
                                n.iter=n.iter, n.thin=n.thin)

  # Take means of posteriors and convert to data.frame with indices
  if (mbnma$type=="time") {
    update.mat <- apply(result[[param]], c(1,2,3), function(x) mean(x, na.rm=TRUE))
    update.df <- reshape2::melt(update.mat)
    names(update.df) <- c("study", "arm", "fupdose", "mean")

    # Remove missing values
    update.df <- update.df[stats::complete.cases(update.df),]

    # Treatment as facet
    temp <- replicate(max(update.df$fupdose), mbnma$model$data()$treatment)
    update.df$facet <- as.vector(temp)[
      stats::complete.cases(as.vector(temp))
      ]

    # Studyarm as group
    update.df$groupvar <- paste(as.numeric(update.df$study), as.numeric(update.df$arm), sep="_")

  } else if (mbnma$type=="dose") {
    update.mat <- apply(result[[param]], c(1,2), function(x) mean(x, na.rm=TRUE))
    update.df <- reshape2::melt(update.mat)
    names(update.df) <- c("study", "arm", "mean")

    # Remove missing values
    update.df <- update.df[stats::complete.cases(update.df),]

    # Agent as facet
    update.df$facet <- as.vector(mbnma$model$data()$agent)[
      stats::complete.cases(as.vector(mbnma$model$data()$agent))
      ]

    update.df$fupdose <- as.vector(mbnma$model$data()$dose)[
      stats::complete.cases(as.vector(mbnma$model$data()$dose))
      ]

    # Study as group
    update.df$groupvar <- as.numeric(update.df$study)
  }

  return(update.df)
}
