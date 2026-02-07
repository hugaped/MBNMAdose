# Functions for running MBNMA models
# Author: Hugo Pedder
# Date created: 2019-04-18

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
#' @param fun An object of `class("dosefun")` that specifies a functional form to be assigned to the
#'   dose-response. See Details.
#' @param method Can take either `"common"` or `"random"` to indicate whether relative effects
#'   should be modelled with between-study heterogeneity or not (see details).
#' @param regress A formula of effect modifiers (variables that
#'  interact with the treatment effect) to incorporate using Network Meta-Regression
#'  (E.g. `~ Population + Age`). All variables in the formula are modelled as interacting
#'  with the treatment effect (i.e. prognostic variables cannot be included in this way).
#'  Effects modifiers must be named variables in `network$data.ab` and must be identical
#'  within a study. Factor and character effect modifiers will be converted to a series of
#'  named dummy variables.
#' @param regress.effect Indicates whether effect modification should be assumed to be
#'  `"common"` (assumed to be equal versus Placebo throughout the network),
#'  `"random"` (assumed to be exchangeable versus Placebo throughout the network),
#'  `"agent"` (assumed to be equal versus Placebo within each agent), or
#'  `"class"` (assumed to be equal versus Placebo within each class).
#' @param class.effect A list of named strings that determines which dose-response
#'   parameters to model with a class effect and what that effect should be
#'   (`"common"` or `"random"`). Element names should match dose-response parameter names.
#'   Note that assuming class effects on some dose-response parameters may be unreasonable if
#'   the range of doses differ substantially across agents within a class.
#' @param UME A boolean object to indicate whether to fit an Unrelated Mean Effects model
#'   that does not assume consistency and so can be used to test if the consistency
#'   assumption is valid.
#' @param likelihood A string indicating the likelihood to use in the model. Can take either `"binomial"`,
#'   `"normal"` or `"poisson"`. If left as `NULL` the likelihood will be inferred from the data.
#' @param link A string indicating the link function to use in the model. Can take any link function
#'   defined within JAGS (e.g. `"logit"`, `"log"`, `"probit"`, `"cloglog"`), be assigned the value `"identity"` for
#'   an identity link function, or be assigned the value `"smd"` for modelling Standardised Mean Differences using an
#'   identity link function. If left as `NULL` the link function will be automatically assigned based
#'   on the likelihood.
#' @param sdscale Logical object to indicate whether to write a model that specifies a reference SD
#'  for standardising when modelling using Standardised Mean Differences. Specifying `sdscale=TRUE`
#'  will therefore only modify the model if link function is set to SMD (`link="smd"`).
#' @param cor A boolean object that indicates whether correlation should be modelled
#' between relative effect dose-response parameters. This is
#' automatically set to `FALSE` if class effects are modelled or if multiple dose-response
#' functions are fitted.
#' @param omega A scale matrix for the inverse-Wishart prior for the covariance matrix used
#' to model the correlation between dose-response parameters (see Details for dose-response functions). `omega` must
#' be a symmetric positive definite matrix with dimensions equal to the number of dose-response parameters modelled using
#' relative effects (`"rel"`). If left as `NULL` (the default) a diagonal matrix with elements equal to 100
#' is used.
#' @param priors A named list of parameter values (without indices) and
#'   replacement prior distribution values given as strings
#'   **using distributions as specified in JAGS syntax** (see \insertCite{jagsmanual;textual}{MBNMAdose}). Note
#'   that normal distributions in JAGS are specified as \deqn{N(\mu, prec)}, where \deqn{prec = 1 / {\sigma^2}}.
#'
#' @param pD logical; if `TRUE` (the default) then adds the computation of pD, using the method
#' of \insertCite{plummer2008}{MBNMAdose}. If `FALSE` then uses the
#' approximation of `pD=var(deviance) / 2` (often referred to as pV).
#' @param n.iter number of total iterations per chain (including burn in; default: 20000)
#' @param n.thin thinning rate. Must be a positive integer. Set `n.thin > 1`` to save memory
#' and computation time if n.iter is large. Default is
#' `max(1, floor(n.chains * (n.iter-n.burnin) / 1000))`` which will only thin if there are at least 2000
#' simulations.
#' @param n.chains number of Markov chains (default: 3)
#' @param n.burnin length of burn in, i.e. number of iterations to discard at the
#' beginning. Default is `n.iter/2``, that is, discarding the first half of the
#' simulations. If n.burnin is 0, jags() will run 100 iterations for adaption.
#' @param autojags A boolean value that indicates whether the model should be continually updated until
#' it has converged below a specific cutoff of `Rhat`
#' @param Rhat A cutoff value for the Gelman-Rubin convergence diagnostic\insertCite{gelmanrubin}{MBNMAdose}.
#' Unless all parameters have Rhat values lower than this the model will continue to sequentially update up
#' to a maximum of `n.update`. Default is `1.05`.
#' @param n.update The maximum number of updates. Each update is run for 1000 iterations, after which the
#' Rhat values of all parameters are checked against `Rhat`. Default maximum updates
#' is `10` (i.e. 10,000 additional iterations in total).
#'
#' @param model.file The file path to a JAGS model (.jags file) that can be used
#'   to overwrite the JAGS model that is automatically written based on the
#'   specified options in `MBNMAdose`. Useful for adding further model flexibility.
#' @param jagsdata A named list of the data objects to be used in the JAGS model. Only
#'   required if users are defining their own JAGS model using `model.file`. Format
#'   should match that of standard models fitted in `MBNMAdose`
#'   (see `mbnma$model.arg$jagsdata`)
#'
#' @param ... Arguments to be sent to R2jags.
#'
#' @details When relative effects are modelled on more than one dose-response parameter and
#' `cor = TRUE`, correlation between the dose-response parameters is automatically
#' estimated using a vague Wishart prior. This prior can be made slightly more informative
#' by specifying the relative scale of variances between the dose-response parameters using
#' `omega`. `cor` will automatically be set to `FALSE` if class effects are modelled.
#'
#' @return An object of S3 `class(c("mbnma", "rjags"))` containing parameter
#'   results from the model. Can be summarized by `print()` and can check
#'   traceplots using `R2jags::traceplot()` or various functions from the package `mcmcplots`.
#'
#'   Nodes that are automatically monitored (if present in the model) have the
#'   following interpretation:
#'
#'   | \strong{Parameters(s)/Parameter Prefix}              | \strong{Interpretation} |
#'   | -------------------------- | -------------- |
#'   | `<named dose-response parameter>` (e.g. `emax`) | The pooled effect for each dose-response parameter, as defined in dose-response functions. Will vary by agent if pooling is specified as `"rel"` in the dose-response function. |
#'   | `sd` | The between-study SD (heterogeneity) for relative effects, reported if `method="random"` |
#'   | `sd.<named dose-response parameter>` (e.g. `sd.emax`) |  Between-study SD (heterogeneity) for absolute dose-response parameters specified as `"random"`. |
#'   | `<named capitalized dose-response parameter>` (e.g. `EMAX`) | The class effect within each class for a given dose-response parameter. These will be estimated by the model if specified in `class.effects` for a given dose-response parameter. |
#'   | `sd.<named capitalized dose-response parameter>` (e.g. `sd.EMAX`) | The within-class SD for different agents within the same class. Will be estimated by the model if any dose-response parameter in `class.effect` is set to `"random"`. |
#'   | `totresdev` | The residual deviance of the model |
#'   | `deviance` | The deviance of the model |
#'
#'
#'   If there are errors in the JAGS model code then the object will be a list
#'   consisting of two elements - an error message from JAGS that can help with
#'   debugging and `model.arg`, a list of arguments provided to `mbnma.run()`
#'   which includes `jagscode`, the JAGS code for the model that can help
#'   users identify the source of the error.
#'
#' @section Dose-response parameter arguments:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#'
#'
#' @section Dose-response function:
#'   Several general dose-response functions are provided, but a
#'   user-defined dose-response relationship can instead be used.
#'
#'   As of version 0.4.0 dose-response functions are specified as an object of `class("dosefun")`. See
#'   help details for each of the functions below for the interpretation of specific dose-response parameters.
#'
#'   Built-in dose-response functions are:
#'   * `dpoly()`: polynomial (e.g. for a linear model - `dpoly(degree=1)`)
#'   * `dloglin()`: log-linear
#'   * `dexp()`: exponential
#'   * `demax()`: (emax with/without a Hill parameter)
#'   * `dspline()`: splines (can fit B-splines (`type="bs"`), restricted cubic splines (`type="rcs"`), natural splines (`type="ns"`), or
#'   piecewise linear splines (`type="ls"`))
#'   * `dfpoly()`: fractional polynomials
#'   * `dnonparam()`: Non-parametric monotonic function (`direction` can be either `"increasing"` or `"decreasing"`) following the method
#'   of \insertCite{owen2015;textual}{MBNMAdose}
#'   * `duser()`: user-defined function
#'   * `dmulti()`: allows agent-specific dose-response functions to be fitted. A separate function must be provided for each agent
#'   in the network.
#'
#'
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
#' network <- mbnma.network(triptans)
#'
#'
#' ######## Dose-response functions ########
#'
#' # Fit a dose-response MBNMA with a linear function
#' # with common treatment effects
#' result <- mbnma.run(network, fun=dpoly(degree=1), method="common")
#'
#' # Fit a dose-response MBNMA with a log-linear function
#' # with random treatment effects
#' result <- mbnma.run(network, fun=dloglin(), method="random")
#'
#' # Fit a dose-response MBNMA with a fractional polynomial function
#' # with random treatment effects
#' # with a probit link function
#' result <- mbnma.run(network, fun=dfpoly(), method="random", link="probit")
#'
#' # Fit a user-defined function (quadratic)
#' fun.def <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
#' result <- mbnma.run(network, fun=duser(fun=fun.def), method="common")
#'
#' # Fit an Emax function
#' # with a single random (exchangeable) parameter for ED50
#' # with common treatment effects
#' result <- mbnma.run(network, fun=demax(emax="rel", ed50="random"),
#'               method="common")
#'
#' # Fit an Emax function with a Hill parameter
#' # with a fixed value of 5 for the Hill parameter
#' # with random relative effects
#' result <- mbnma.run(network, fun=demax(hill=5), method="random")
#'
#' # Fit a model with natural cubic splines
#' # with 3 knots at 10% 30% and 60% quartiles of dose ranges
#' depnet <- mbnma.network(ssri) # Using the sSRI depression dataset
#' result <- mbnma.run(depnet, fun=dspline(type="ns", knots=c(0.1,0.3,0.6)))
#'
#' # Fit a model with different dose-response functions for each agent
#' multifun <- dmulti(list(dloglin(), # for placebo (can be any function)
#'                        demax(), # for eletriptan
#'                        demax(), # for sumatriptan
#'                        dloglin(), # for frovatriptan
#'                        demax(), # for almotriptan
#'                        demax(), # for zolmitriptan
#'                        dloglin(), # for naratriptan
#'                        demax())) # for rizatriptan
#'
#' multidose <- mbnma.run(network, fun=multifun)
#'
#'
#' ########## Class effects ##########
#'
#'  # Using the osteoarthritis dataset
#'  pain.df <- osteopain
#'
#'  # Set a shared class (NSAID) only for Naproxcinod and Naproxen
#'  pain.df <- pain.df %>% dplyr::mutate(
#'                 class = dplyr::case_when(agent %in% c("Naproxcinod", "Naproxen") ~ "NSAID",
#'                         !agent %in% c("Naproxcinod", "Naproxen") ~ agent
#'                         )
#'                 )
#'
#'  # Run an Emax MBNMA with a common class effect on emax
#'  painnet <- mbnma.network(pain.df)
#'  result <- mbnma.run(painnet, fun = demax(),
#'                 class.effect = list(emax = "common"))
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from a fractional polynomial function
#' result <- mbnma.run(network, fun=dfpoly(degree=1), method="random")
#' print(result$model.arg$priors)
#'
#' # Change the prior distribution for the power
#' newpriors <- list(power.1 = "dnorm(0,0.001) T(0,)")
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)")
#'
#' result <- mbnma.run(network, fun=dfpoly(degree=1), method="random",
#'               priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' result <- mbnma.run(network, fun=dloglin(), method="random",
#'               n.iter=5000, n.thin=5, n.chains=4)
#'
#' ####### Examine MCMC diagnostics (using mcmcplots or coda packages) #######
#'
#' # Density plots
#' mcmcplots::denplot(result)
#'
#' # Traceplots
#' mcmcplots::traplot(result)
#'
#' # Caterpillar plots
#' mcmcplots::caterplot(result, "rate")
#'
#' # Autocorrelation plots (using the coda package)
#' coda::autocorr.plot(coda::as.mcmc(result))
#'
#'####### Automatically run jags until convergence is reached #########
#'
#' # Rhat of 1.08 is set as the criteria for convergence
#' #on all monitored parameters
#' conv.res <- mbnma.run(network, fun=demax(),
#'               method="random",
#'               n.iter=10000, n.burnin=9000,
#'               autojags=TRUE, Rhat=1.08, n.update=8)
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
                      fun=dpoly(degree=1),
                      method="common",
                      regress=NULL,
                      regress.effect="common",
                      class.effect=list(), UME=FALSE,
                      sdscale=FALSE,
                      cor=FALSE,
                      omega=NULL,
                      parameters.to.save=NULL,
                      pD=TRUE,
                      likelihood=NULL, link=NULL,
                      priors=NULL,
                      n.iter=20000, n.chains=3,
                      n.burnin=floor(n.iter/2), n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                      autojags=FALSE, Rhat=1.05, n.update=10,
                      model.file=NULL, jagsdata=NULL, ...
) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertCharacter(model.file, len=1, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(pD, null.ok=FALSE, add=argcheck)
  #checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(cor, len=1, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  args <- as.list(environment())

  # Check fun
  fun <- check.fun(fun=fun, network=network)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link, warnings=TRUE)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  # Check regression
  if (!is.null(regress)) {
    regress.mat <- check.regress(network=network, regress=regress)
  } else {
    regress.mat <- NULL
  }

  # Check sdscale
  if (sdscale==TRUE) {
    if (link!="smd") {
      sdscale <- FALSE
    }
  }

  # Reduce n.burnin by 1 to avoid JAGS error if n.burnin=n.iter
  if (n.iter==n.burnin) {
    n.burnin <- n.burnin - 1
  }

  # Ensure pd.kl not run with parallel
  parallel <- FALSE
  if (parallel==TRUE & pD==TRUE) {
    warning("pD cannot be calculated using Kullback-Leibler divergence for\nmodels run in parallel. Defaulting to pD=FALSE")
    pD <- FALSE
  }

  # Ensure cor set to FALSE if multiple dose-response functions are modelled
  # if (cor==TRUE & length(fun)>1) {cor <- FALSE}
  if (cor==TRUE & length(fun$name)>1) {cor <- FALSE}

  # If multiple dose-response parameters are relative effects then add omega default
  if (cor==TRUE & is.null(omega)) {
    relparam <- fun$apool %in% "rel" & !names(fun$apool) %in% names(class.effect)
    if (sum(relparam)>1) {
      omega <- diag(rep(100,sum(relparam)))
    }
  }

  # Ensure rjags parameters make sense
  if (n.iter<=n.burnin) {
    stop(paste0("`n.iter` must be greater than `n.burnin`. `n.burnin` = ", n.burnin))
  }

  # Check if placebo has been included
  if (network$agents[1]=="Placebo" & ("Placebo_0" %in% network$treatments)) {
    plac.incl <- TRUE
  } else {
    plac.incl <- FALSE

    if ("nonparam" %in% fun$name) {
      stop("Non-parametric dose-response models must currently include placebo - due to be updated in subsequent versions")
    }
  }

  if (is.null(model.file)) {

    # Write JAGS model code
    model <- mbnma.write(fun=fun,
                         method=method,
                         regress.mat=regress.mat, regress.effect=regress.effect,
                         class.effect=class.effect, UME=UME,
                         sdscale=sdscale,
                         cor=cor, omega=omega,
                         om=calcom(data.ab=network$data.ab, link=link, likelihood=likelihood),
                         likelihood=likelihood, link=link
    )

    # Change code for if plac not included in network
    if (plac.incl==FALSE) {
      model <- gsub("for \\(k in 2:Nagent\\)\\{ # Priors on relative agent effects",
                    "for (k in 1:Nagent){ # Priors on relative agent effects",
                    model)
      model <- gsub("for \\(k in 2:Nclass\\)\\{ # Priors on relative class effects",
                    "for (k in 1:Nclass){ # Priors on relative class effects",
                    model)

      drop <- grep("s\\.beta\\.[(0-9)+]\\[1\\] <- \\.?[0-9]+", model)
      if (length(drop)>0) {
        model <- model[-drop]
      }
    }

    # Add user-defined priors to the model
    if (!is.null(priors)) {
      model <- replace.prior(priors=priors, model=model)
    }

  } else {
    warning("All parameter specifications (dose-response parameters, class effects, priors, etc.) are being overwritten by `model.file`")
    model <- readLines(model.file)
  }

  # Generate default parameters to monitor
  assigned.parameters.to.save <- parameters.to.save
  if (is.null(parameters.to.save)) {
    parameters.to.save <-
      gen.parameters.to.save(fun=fun, model=model, regress.mat = regress.mat)
  }

  # Set boolean for presence of class effects in model
  class <- ifelse(length(class.effect)>0 | "class" %in% regress.effect, TRUE, FALSE)

  # Change doses to dose indices for non-parametric models
  if ("nonparam" %in% fun$name) {
    data.ab <- index.dose(network[["data.ab"]])[["data.ab"]]
  } else {
    data.ab <- network[["data.ab"]]
  }


  #### Run jags model ####

  result.jags <- mbnma.jags(data.ab, model,
                            regress=regress, regress.effect=regress.effect,
                            class=class, omega=omega, sdscale=sdscale,
                            parameters.to.save=parameters.to.save,
                            likelihood=likelihood, link=link, fun=fun,
                            jagsdata=jagsdata,
                            n.iter=n.iter,
                            n.thin=n.thin,
                            n.chains=n.chains,
                            n.burnin=n.burnin,
                            parallel=parallel, pD=pD,
                            autojags=autojags, Rhat=Rhat, n.update=n.update,
                            UME=UME, class.effect=class.effect,
                            ...)
  result <- result.jags[["jagsoutput"]]
  jagsdata <- result.jags[["jagsdata"]]


  # Calculate model fit statistics (using differnt pD as specified)
  # if (!"error" %in% names(result)) {
  #   fitstats <- changepd(model=result, jagsdata=jagsdata, pd=pd, likelihood=likelihood, type="dose")
  #   result$BUGSoutput$pD <- fitstats$pd
  #   result$BUGSoutput$DIC <- fitstats$dic
  # }

  # Define model arguments
  model.arg <- list("parameters.to.save"=assigned.parameters.to.save,
                    "fun"=fun,
                    "jagscode"=result.jags$model,
                    "jagsdata"=jagsdata,
                    "method"=method,
                    "likelihood"=likelihood, "link"=link,
                    "regress"=regress, "regress.mat"=regress.mat, "regress.effect"=regress.effect,
                    "class.effect"=class.effect,
                    "cor"=cor,
                    "omega"=omega,
                    "UME"=UME,
                    "sdscale"=sdscale,
                    #"parallel"=parallel,
                    "pD"=pD,
                    "priors"=get.prior(model))

  result[["model.arg"]] <- model.arg
  result[["type"]] <- "dose"
  result[["network"]] <- network

  if (length(class.effect)>0) {
    result[["classes"]] <- network[["classes"]]
  }

  if (!("error" %in% names(result))) {

    result <- cutjags(result)

    class(result) <- c("mbnma", class(result))
  }

  return(result)

}





mbnma.jags <- function(data.ab, model,
                       class=FALSE, sdscale=FALSE,
                       regress=NULL, regress.effect="common",
                       omega=NULL,
                       likelihood=NULL, link=NULL, fun=NULL,
                       nodesplit=NULL, jagsdata=NULL,
                       warn.rhat=FALSE, parallel=FALSE,
                       autojags=FALSE, Rhat=1.1, n.update=10,
                       class.effect=NULL, UME=NULL,
                       ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(model, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertClass(fun, "dosefun", null.ok=TRUE, add=argcheck)
  checkmate::assertNumeric(nodesplit, len=2, null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(autojags, null.ok=FALSE, add=argcheck)
  checkmate::assertNumeric(Rhat, lower=1, add=argcheck)
  checkmate::assertNumeric(n.update, lower=1, add=argcheck)
  checkmate::assertList(jagsdata, null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(sdscale, len = 1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  args <- list(...)

  # For MBNMAdose
  if (is.null(jagsdata)) {
    jagsdata <- getjagsdata(data.ab, class=class, sdscale=sdscale,
                            regress=regress, regress.effect=regress.effect,
                            likelihood=likelihood, link=link, fun=fun,
                            nodesplit=nodesplit)

    if (!is.null(omega) & class==FALSE) {
      jagsdata[["omega"]] <- omega
    }

    # Add parameter to monitor for direct evidence in node-splitting
    if (!is.null(nodesplit)) {
      model <- add.nodesplit(model=model)

      if ("parameters.to.save" %in% names(args)) {
        args[["parameters.to.save"]] <- append(args[["parameters.to.save"]], "direct")
      }
    }
  }

  # Drop dose from jagsdata in spline models
  dosedat <- jagsdata[["dose"]]
  if (all(fun$name %in% c("rcs", "ns", "bs", "ls"))) {
    jagsdata[["dose"]] <- NULL
  }


  # Remove studyID from jagsdata (not used in model)
  tempjags <- jagsdata
  tempjags[["studyID"]] <- NULL

  # Put data from jagsdata into separate R objects
  for (i in seq_along(tempjags)) {
    # extract the object value
    temp <- tempjags[[i]]
    # create a new variable with the original name of the list item
    eval(parse(text=paste(names(tempjags)[[i]],"<- temp")))
  }

  # Take names of variables in jagsdata for use in rjags
  jagsvars <- list()
  for (i in seq_along(names(tempjags))) {
    jagsvars[[i]] <- names(tempjags)[i]
  }

  # Create a temporary model file
  tmpf=tempfile()
  tmps=file(tmpf,"w")
  if (length(model)==1) {
    cat(model,file=tmps)
  } else if (length(model)>1) {
    cat(paste(model, collapse="\n"),file=tmps)
  }
  close(tmps)

  # Generate default initial values if not supplied
  if (!("inits" %in% names(args)) && !is.null(class.effect) && !is.null(UME)) {
    args$inits <- gen.inits(jagsdata, fun=fun, class.effect=class.effect, UME=UME, n.chains=args$n.chains)
  }

  out <- tryCatch({
    withCallingHandlers({
      if (parallel == FALSE) {
        result <- do.call(R2jags::jags, c(args, list(data = jagsvars, model.file = tmpf)))

        # Automatically update
        if (autojags == TRUE) {
          result <- R2jags::autojags(result, Rhat = Rhat, n.update = n.update, n.iter = 1000, refresh = 100)
        } else if (autojags == FALSE) {
          result <- result
        }
      } else if (parallel == TRUE) {
        if (autojags == TRUE) {
          stop("autojags=TRUE cannot be used with parallel=TRUE")
        }

        # Run jags in parallel
        result <- do.call(R2jags::jags.parallel, c(args, list(data = jagsvars, model.file = tmpf)))
      }
    }, warning = function(w) {
      if (grepl("missing in parameter", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
  },
  error = function(cond) {
    message(cond)
    return(list(error = cond))
  })

  # Gives warning if any rhat values > 1.02
  if (warn.rhat==TRUE) {
    if (!("error" %in% names(out))) {
      rhat.warning(out)
    }
  }

  jagsdata[["dose"]] <- dosedat # add dose back to jagsdata if spline models were used
  names(model) <- NULL

  return(list("jagsoutput"=out, "jagsdata"=jagsdata, "model"=model))
}




#' Generate initial values for non-parametric dose-response functions
#'
#' Ensures model runs properly
#'
#' @noRd
gen.inits <- function(jagsdata, fun, class.effect, UME, n.chains)
{
    Nclass <- jagsdata$Nclass
    Nagent <- jagsdata$Nagent
    NS <- jagsdata$NS

    inits <- vector("list", n.chains)
    for (i in 1:n.chains) {

        inits.i <- list(
            mu = rep(0, jagsdata$NS)
        )

        if (fun$name == "nonparam") {
            inits.i[["d1"]] <- gen.init(jagsdata, fun)
        }
        else {
            for (j in seq_along(fun$apool)) {

                pname <- names(fun$apool[j])

                truncated <- (fun$p.expon == FALSE) && ((pname %in% c("ed50", "onset")) || (pname=="rate" && "itp"==fun$name))

                if (fun$apool[j] == "rel") {

                    if (pname %in% names(class.effect)) {
                        ## If class effects are modelled for this parameter

                        ## Initialize truncated parameters to a small positive value away from the boundary
                        mu0 <- if (truncated) 0.001 else 0

                        inits.i[[toupper(pname)]] <-  c(NA, rep(mu0, Nclass-1))

                    } else {

                        if (fun$apool[j] %in% c("common", "random")) {
                            ## Scalar prior
                            init.pname <- if (truncated) 0 else 0.001
                        } else if (UME==FALSE) {
                            ## Vector prior
                            init.pname <- c(NA, rep(0, Nagent - 1))
                        }
                        else {
                            ## Matrix prior
                            init.pname <- matrix(NA, Nagent, Nagent-1)
                            for (c in 1:(Nagent-1)) {
                                for (k in (c+1):Nagent) { # UME priors
                                    init.pname[k,c] <- 0
                                }
                            }
                        }
                        inits.i[[pname]] <- init.pname
                    }
                }
            }
        }
        inits[[i]] <- inits.i
    }
    return (inits)
}

gen.init <- function(jagsdata, fun) {
  checkmate::assert_set_equal(fun$name, "nonparam")

  sapply(jagsdata$maxdose, FUN=function(x, direction=fun$direction) {
   val <- x-1

   if (direction=="increasing") {
     start <- stats::runif(1,0,3)
     c(NA, seq(start, start+stats::runif(1,0,3), length.out=val), rep(NA, max(jagsdata$maxdose)-(val+1)))
   } else if (direction=="decreasing") {
     start <- stats::runif(1,-3,0)
     c(NA, seq(start, start-stats::runif(1,0,3), length.out=val), rep(NA, max(jagsdata$maxdose)-(val+1)))
   }
  })
}




#' Automatically generate parameters to save for a dose-response MBNMA model
#'
#' @inheritParams mbnma.run
#' @inheritParams mbnma.write
#' @param model A JAGS model written as a character object
gen.parameters.to.save <- function(fun, model, regress.mat=NULL) {
  # model.params is a vector (numeric/character) of the names of the dose-response parameters in the model
  #e.g. c(1, 2, 3) or c("emax", "et50")
  # model is a JAGS model written as a character object

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(model, min.len = 10, add=argcheck)
  checkmate::assertClass(fun, "dosefun", add=argcheck)
  checkmate::assertMatrix(regress.mat, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  # Set some automatic parameters based on the model code
  parameters.to.save <- vector()

  for (i in seq_along(fun$params)) {
    # For unnamed parameters
    if (any(grepl(paste0("^d\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("d.", i))
    }
    if (any(grepl(paste0("^D\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("D.", i))
    }
    if (any(grepl(paste0("^sd\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", i))
    }
    if (any(grepl(paste0("^sd\\.D\\.", i), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.D.", i))
    }

    # For named parameters
    if (any(grepl(fun$params[i], model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, fun$params[i])
    }
    if (any(grepl(toupper(fun$params[i]), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, toupper(fun$params[i]))
    }
    if (any(grepl(paste0("^sd\\.", fun$params[i]), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", fun$params[i]))
    }
    if (any(grepl(paste0("^sd\\.", toupper(fun$params[i])), model))==TRUE) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.", toupper(fun$params[i])))
    }

    # Remove beta if both d and beta are in for any parameter
    if (paste0("d.",i) %in% parameters.to.save & paste0("beta.",i) %in% parameters.to.save) {
      parameters.to.save <- parameters.to.save[!parameters.to.save %in% paste0("beta.",i)]
    }
  }

  for (i in seq_along(colnames(regress.mat))) {
    # Save B
    parameters.to.save <- append(parameters.to.save, paste0("B.", colnames(regress.mat)[i]))

    # save sd.B
    if (any(grepl("sd\\.B\\.", model))) {
      parameters.to.save <- append(parameters.to.save, paste0("sd.B.", colnames(regress.mat)[i]))
    }
  }


  # Include nonparametric
  if ("nonparam" %in% fun$name) {
    parameters.to.save <- append(parameters.to.save, "d.1")
  }
  if (any(grepl("totresdev", model))) {
    parameters.to.save <- append(parameters.to.save, c("totresdev"))
  }
  if (any(grepl("^sd ", model))) {
    parameters.to.save <- append(parameters.to.save, c("sd"))
  }

  return(unique(parameters.to.save))
}







#' Run an NMA model
#'
#' Used for calculating treatment-level NMA results, either when comparing MBNMA models to models that
#' make no assumptions regarding dose-response , or to estimate split results for `overlay.split`.
#' Results can also be compared between consistency (`UME=FALSE`) and inconsistency
#' (`UME=TRUE`) models to test the validity of the consistency assumption at the treatment-level.
#'
#' @inheritParams mbnma.run
#' @param warn.rhat A boolean object to indicate whether to return a warning if Rhat values
#' for any monitored parameter are >1.02 (suggestive of non-convergence).
#' @param drop.discon A boolean object that indicates whether or not to drop disconnected
#'   studies from the network.
#' @param n.iter number of total iterations per chain (including burn in; default: 20000)
#'
#' @examples
#' \donttest{
#' # Run random effects NMA on the alogliptin dataset
#' alognet <- mbnma.network(alog_pcfb)
#' nma <- nma.run(alognet, method="random")
#' print(nma)
#' plot(nma)
#'
#' # Run common effects NMA keeping treatments that are disconnected in the NMA
#' goutnet <- mbnma.network(gout)
#' nma <- nma.run(goutnet, method="common", drop.discon=FALSE)
#'
#' # Run an Unrelated Mean Effects (UME) inconsistency model on triptans dataset
#' tripnet <- mbnma.network(triptans)
#' ume <- nma.run(tripnet, method="random", UME=TRUE)
#' }
#'
#' @export
nma.run <- function(network, method="common", likelihood=NULL, link=NULL, priors=NULL,
                    sdscale=FALSE,
                    warn.rhat=TRUE, n.iter=20000, drop.discon=TRUE, UME=FALSE, pD=TRUE,
                    parameters.to.save=NULL, ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::assertLogical(warn.rhat, add=argcheck)
  checkmate::assertIntegerish(n.iter, null.ok = TRUE, add=argcheck)
  checkmate::assertLogical(drop.discon, add=argcheck)
  checkmate::assertLogical(UME, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(pD, null.ok=FALSE, add=argcheck)
  checkmate::assertCharacter(parameters.to.save, null.ok=TRUE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  args <- list(...)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  # Check sdscale
  if (sdscale==TRUE) {
    if (link!="smd") {
      sdscale <- FALSE
    }
  }

  #### Write model for NMA ####
  model <- write.nma(method=method, likelihood=likelihood, link=link, UME=UME,
                     sdscale=sdscale,
                     om=calcom(data.ab=network$data.ab, link=link, likelihood=likelihood))

  #### Add priors ####
  if (!is.null(priors)) {
    # Check that order is correct if names in priors match network$agents
    for (i in seq_along(priors)) {
      if (length(priors[[i]])>1) {
        if (!is.null(names(priors[[i]]))) { # If treatment-specific priors are named
          if (length(class.effect)>0) {
            warning("MBNMAdose defaults to agent-specific priors if multiple priors are specified for a parameter")
          }

          priornam <- names(priors[[i]])

          # If at least 2 names match those in network$agents then sort priors to match network order
          if (sum(priornam %in% network$agents)>=2) {
            priors[[i]] <- priors[[i]][match(network$agents, priornam, nomatch=0)]
          }
        }
      }
      prior <- priors[[i]]
    }

    model <- replace.prior(priors=priors, model=model)
  }

  #### Parameters ####
  if (is.null(parameters.to.save)){
    parameters.to.save <- c("d", "totresdev")
    if (method=="random") {
      parameters.to.save <- append(parameters.to.save, "sd")
    }
  }


  # # Add nodes to monitor to calculate plugin pd
  # if (pd=="plugin") {
  #   pluginvars <- c("psi", "resdev")
  #   for (param in seq_along(pluginvars)) {
  #     if (!(pluginvars[param] %in% parameters.to.save)) {
  #       parameters.to.save <- append(parameters.to.save, pluginvars[param])
  #     }
  #   }
  #   message("The following parameters have been monitored to allow pD plugin calculation: ",
  #           paste(pluginvars, collapse=", "))
  # }

  #### Prepare data ####
  data.ab <- network$data.ab
  trt.labs <- network$treatments

  # Check treatments that are not connected and remove if not
  if (drop.discon==TRUE) {
    connect <- drop.disconnected(network)
    data.ab <- connect[["data.ab"]]
    trt.labs <- connect[["trt.labs"]]
  }

  jagsdata <- getjagsdata(data.ab, likelihood = likelihood, link=link, sdscale=sdscale,
                          level="treatment")


  #### Run JAGS model ####
  # Put data from jagsdata into separate R objects
  for (i in seq_along(jagsdata)) {
    # extract the object value
    temp <- jagsdata[[i]]
    # create a new variable with the original name of the list item
    eval(parse(text=paste(names(jagsdata)[[i]],"<- temp")))
  }

  # Take names of variables in jagsdata for use in rjags
  jagsvars <- list()
  tempjags <- jagsdata
  tempjags[["studyID"]] <- NULL # Remove studyID from jagsdata (not used in model)
  for (i in seq_along(names(tempjags))) {
    jagsvars[[i]] <- names(tempjags)[i]
  }

  # Create a temporary model file
  tmpf=tempfile()
  tmps=file(tmpf,"w")
  if (length(model)==1) {
    cat(model,file=tmps)
  } else if (length(model)>1) {
    cat(paste(model, collapse="\n"),file=tmps)
  }
  close(tmps)

  out <- tryCatch({
    withCallingHandlers({
      result <- do.call(R2jags::jags, c(args, list(data=jagsvars, model.file=tmpf, n.iter=n.iter,
                                                   parameters.to.save=parameters.to.save)))
    }, warning = function(w) {
      if (grepl("missing in parameter", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
  },
  error = function(cond) {
    message(cond)
    return(list(error = cond))
  })

  # Gives warning if any rhat values > 1.02
  if (warn.rhat==TRUE) {
    if (!("error" %in% names(out))) {
      rhat.warning(out)
    }
  }

  # Calculate model fit statistics (using differnt pD as specified)
  # fitstats <- changepd(model=result, jagsdata=jagsdata, pd=pd, likelihood=likelihood, type="dose")
  # out$BUGSoutput$pD <- fitstats$pd
  # out$BUGSoutput$DIC <- fitstats$dic

  output <- list("jagsresult"=out, "trt.labs"=trt.labs, "UME"=UME)
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
#' @noRd
check.likelink <- function(data.ab, likelihood=NULL, link=NULL, warnings=FALSE) {

  # Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit", "log", "smd"), null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(warnings, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (is.null(likelihood)) {
    if (all(c("r", "n") %in% names(data.ab))) {
      likelihood <- "binomial"

      if (warnings==TRUE) {
        message("`likelihood` not given by user - set to `binomial` based on data provided")
      }
    } else if (all(c("y", "se") %in% names(data.ab))) {
      likelihood <- "normal"

      if (warnings==TRUE) {
        message("`likelihood` not given by user - set to `normal` based on data provided")
      }
    } else if (all(c("r", "E") %in% names(data.ab))) {
      likelihood <- "poisson"

      if (warnings==TRUE) {
        message("`likelihood` not given by user - set to `poisson` based on data provided")
      }
    }
  }
  if (is.null(link)) {
    if (likelihood=="binomial") {
      link <- "logit"

      if (warnings==TRUE) {
        message("`link` not given by user - set to `logit` based on assigned value for `likelihood`")
      }
    } else if (likelihood=="normal") {
      link <- "identity"

      if (warnings==TRUE) {
        message("`link` not given by user - set to `identity` based on assigned value for `likelihood`")
      }
    } else if (likelihood=="poisson") {
      link <- "log"

      if (warnings==TRUE) {
        message("`link` not given by user - set to `log` based on assigned value for `likelihood`")
      }
    }
  }

  # Check valid likelihood is used
  if (likelihood=="binomial" & !all(c("r", "n") %in% names(data.ab))) {
    stop("Binomial likelihood - columns `r` and `n` must be included in `data.ab`")
  } else if (likelihood=="poisson" & !all(c("E", "n") %in% names(data.ab))) {
    stop("Poisson likelihood - columns `E` and `n` must be included in `data.ab`")
  } else if (likelihood=="normal" & !all(c("y", "se") %in% names(data.ab))) {
    stop("Normal likelihood - columns `y` and `se` must be included in `data.ab`")
  }

  # Check valid data available for SMD
  if (link=="smd" & !all("n" %in% names(data.ab))) {
    stop("Treatment effects modelled as Standardised Mean Differences (SMD)\n columns 'n' must be included in data.ab")
  }

  return(list("likelihood"=likelihood, "link"=link))
}





#' Check dose-response functions
#'
#' Checks correct specification of dose-response function and and converts from previous version's syntax.
#' Can be deprecated eventually
#'
#' @inheritParams mbnma.run
#' @inheritParams mbnma.network
#'
#' @noRd
check.fun <- function(fun, network) {

  checkmate::assertClass(network, "mbnma.network")
  checkmate::assertClass(fun, "dosefun")

  if (length(fun[["name"]])>1) {
    if (length(fun[["posvec"]])!=length(network$agents)) {
      stop("Number of agent-specific dose-response functions in dmulti() must be equal to the number of agents in network$agents")
    }
    if (!is.null(fun[["agents"]])) {
      err <- which(is.na(match(fun$agents, network$agents)))
      if (length(err)>1) {
        stop(paste0("Agent names specified in dmulti() are not in network object:\n",
                    paste(fun$agents[err], collapse="\t")))
      }

      # Order of agents must be same as in network
      err <- match(fun$agents, network$agents)
      if (!all(err==1:length(fun$agents))) {
        stop("Agent names specified in dmulti() must be ordered the same as agents in network object")
      }
    }
  }

  return(fun)
}






#' Calculate plugin pD from a JAGS model with univariate likelihood for studies
#' with repeated measurements
#'
#' Uses results from MBNMA JAGS models to calculate pD via the
#' plugin method \insertCite{spiegelhalter2002}{MBNMAdose}. Can only be used for models with known
#' standard errors or covariance matrices. Currently only functions with univariate likelihoods. Function
#' is identical in MBNMAdose and MBNMAtime packages.
#'
#' @param obs1 A matrix (study x arm) or array (study x arm x time point) containing
#'   observed data for `y` (normal likelihood) or `r` (binomial or poisson likelihood)
#'   in each arm of each study. This will be the same array
#'   used as data for the JAGS model.
#' @param obs2 A matrix (study x arm) or array (study x arm x time point) containing
#'   observed data for `se` (normal likelihood), `n` (binomial likelihood) or `E` (poisson likelihood)
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
#' network <- mbnma.network(triptans)
#'
#' # Fit a dose-response MBNMA, monitoring "psi" and "resdev"
#' result <- mbnma.run(network, fun=dloglin(), method="random",
#'               parameters.to.save=c("psi", "resdev"))
#'
#'
#' #### Calculate pD for binomial data ####
#'
#' # Prepare data for pD calculation
#' r <- result$model$data()$r
#' n <- result$model$data()$n
#' narm <- result$model$data()$narm
#' NS <- result$model$data()$NS
#'
#' psi <- result$BUGSoutput$median$psi
#' resdevs <- result$BUGSoutput$median$resdev
#'
#' # Calculate pD via plugin method
#' pD <- pDcalc(obs1=r, obs2=n, narm=narm, NS=NS,
#'           theta.result=psi, resdev.result=resdevs,
#'           likelihood="binomial", type="dose")
#' }
#'
#' @export
pDcalc <- function(obs1, obs2, fups=NULL, narm, NS, theta.result, resdev.result,
                   likelihood="normal", type="time") {

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

  # If function is run in MBNMAtime package
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

      # For MBNMAtime
      if (type=="time") {
        for (m in 1:fups[i]) {
          # Use formula for residual deviance as plugin
          if (likelihood=="normal") {
            dev.post[i,k,m] <- ((obs1[i,k,m] - theta.result[i,k,m])/obs2[i,k,m])^2
            pD[i,k,m] <- resdev.result[i,k,m] - dev.post[i,k,m]
          } else {
            stop("pD cannot be calculated via `plugin` method for time-course MBNMA models without data following a normal likelihood")
          }
        }

        # For MBNMAdose
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

        # Calculate pd contribution for each data point
        pD[i,k] <- resdev.result[i,k] - dev.post[i,k]

      }

    }
  }

  # Sum individual pd contributions
  pD <- sum(pD, na.rm=TRUE)

  return(pD)
}






#' Update MBNMA to monitor deviance nodes in the model
#'
#' Useful for obtaining deviance contributions or fitted values. Same function used in MBNMAdose and MBNMAtime packages.
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
#' @param armdat Include raw arm-level data for each data point (agent, dose, study grouping)
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
#' network <- mbnma.network(triptans)
#'
#' # Fit a dose-response MBNMA, monitoring "psi" and "resdev"
#' result <- mbnma.run(network, fun=dloglin(), method="random",
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
mbnma.update <- function(mbnma, param="theta", armdat=TRUE,
                         n.iter=mbnma$BUGSoutput$n.iter, n.thin=mbnma$BUGSoutput$n.thin) {
  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "rjags", add=argcheck)
  checkmate::assertCharacter(param, len = 1, add=argcheck)
  checkmate::assertLogical(armdat, add=argcheck)
  checkmate::reportAssertions(argcheck)

  modelcode <- mbnma$model$model()

  # Ensure param is in model code
  if (all(grepl(paste0("^", param, "\\[i\\,k\\]"), modelcode)==FALSE)) {
    stop(paste0(param, " not in model code or does not vary by arm"))
  }


  # Update JAGS model for additional iterations
  result <- rjags::jags.samples(mbnma$model, variable.names = param,
                                n.iter=n.iter, n.thin=n.thin)

  # Take means of posteriors and convert to data.frame with indices
  update.mat <- apply(result[[param]], c(1,2), function(x) mean(x, na.rm=TRUE))
  update.df <- reshape2::melt(update.mat)
  names(update.df) <- c("study", "arm", "mean")

  # Remove missing values
  update.df <- update.df[stats::complete.cases(update.df),]

  if (armdat==TRUE) {
    # Agent as facet
    update.df$facet <- as.vector(mbnma$model.arg$jagsdata$agent)[
      stats::complete.cases(as.vector(mbnma$model.arg$jagsdata$agent))
    ]

    update.df$fupdose <- as.vector(mbnma$model.arg$jagsdata$dose)[
      stats::complete.cases(as.vector(mbnma$model.arg$jagsdata$dose))
    ]

    # Study as group
    update.df$groupvar <- as.numeric(update.df$study)
  }

  return(update.df)
}



