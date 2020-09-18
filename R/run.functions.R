# Functions for running MBNMA models
# Author: Hugo Pedder
# Date created: 2019-04-18

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))

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
#' @param fun A character vector specifying a functional form to be assigned to the
#'   dose-response. Options are given in `details`.
#' @param user.fun A formula specifying any relationship including `dose` and
#'   one/several of: `beta.1`, `beta.2`, `beta.3`, `beta.4`.
#' @param model.file A JAGS model written as a character object that can be used
#'   to overwrite the JAGS model that is automatically written based on the
#'   specified options.
#'
#' @param beta.1 Refers to dose-parameter(s) specified within the dose-response function(s).
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.2 Refers to dose-parameter(s) specified within the dose-response function(s).
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.3 Refers to dose-parameter(s) specified within the dose-response function(s).
#' Can take either `"rel"`, `"common"`, `"random"`, or be assigned a numeric value (see details).
#' @param beta.4 Refers to dose-parameter(s) specified within the dose-response function(s).
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
#' @param UME A boolean object to indicate whether to fit an Unrelated Mean Effects model
#'   that does not assume consistency and so can be used to test if the consistency
#'   assumption is valid.
#' @param likelihood A string indicating the likelihood to use in the model. Can take either `"binomial"`,
#'   `"normal"` or `"poisson"`. If left as `NULL` the likelihood will be inferred from the data.
#' @param link A string indicating the link function to use in the model. Can take any link function
#'   defined within JAGS (e.g. `"logit"`, `"log"`, `"probit"`, `"cloglog"`) or be assigned the value `"identity"` for
#'   and identity link function. If left as `NULL` the link function will be automatically assigned based
#'   on the likelihood.
#' @param knots The number/location of knots if a restricted cubic spline dose-response function is fitted (`fun="rcs"`).
#' If a single number is given it indicates the number of knots (they will
#'   be equally spaced across the range of doses). If a numeric vector is given it indicates the location of the knots.
#'   Minimum number of knots is 3.
#' @param cor A boolean object that indicates whether correlation should be modelled
#' between relative effect dose-response parameters (`TRUE`) or not (`FALSE`). This is
#' automatically set to `FALSE` if class effects are modelled or if multiple dose-response
#' functions are fitted.
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
#'   will require running the model for additional iterations but is a more robust calculation for the effective
#'   number of parameters in non-linear models.
#'   * `popt` calculates pD using an optimism adjustment which allows for calculation
#'   of the penalized expected deviance \insertCite{plummer2008}{MBNMAdose}
#' @param parallel A boolean value that indicates whether JAGS should be run in
#'   parallel (`TRUE`) or not (`FALSE`). If `TRUE` then the number of cores to
#'   use is automatically calculated.
#' @param arg.params Contains a list of arguments sent to `mbnma.run()` by dose-response
#' specific wrapper functions
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
#' to a maximum of `n.update`. Default is `1.1`
#' @param n.update The maximum number of updates. Each update is run for 1000 iterations, after which the
#' Rhat values of all parameters are checked against `Rhat`. Default maximum updates
#' is `10` (i.e. 10,000 additional iterations in total).
#' @param ... Arguments to be sent to R2jags.
#'
#'
#' @details When relative effects are modelled on more than one dose-response parameter and
#' `cor = TRUE`, correlation between the dose-response parameters is automatically
#' estimated using a vague Wishart prior. This prior can be made slightly more informative
#' by specifying the relative scale of variances between the dose-response parameters using
#' `var.scale`. `cor` will automatically be set to `FALSE` if class effects are modelled.
#'
#' @return An object of S3 `class(c("mbnma", "rjags"))` containing parameter
#'   results from the model. Can be summarized by `print()` and can check
#'   traceplots using `R2jags::traceplot()` or various functions from the package `mcmcplots`.
#'
#'   Nodes that are automatically monitored (if present in the model) have the
#'   following interpretation:
#'
#'   | \strong{Parameters(s)}              | \strong{Interpretation} |
#'   | -------------------------- | -------------- |
#'   | `d.1`, `d.2`, `d.3`, `d.4` | The pooled effect for each agent for a given dose-response parameter. These will be estimated by the model for dose-response parameters (`beta.1`, `beta.2`, `beta.3`, `beta.4`) specified as `"rel"` (e.g.`mbnma.run(beta.1="rel")`) |
#'   | `sd` | The between-study SD (heterogeneity) for relative effects, reported if `method="random"` |
#'   | `D.1`, `D.2`, `D.3`, `D.4` | The class effect for each class for a given dose-response parameter. These will be estimated by the model if specified in `class.effects` for a given dose-response parameter. |
#'   | `sd.D.1`, `sd.D.2`, `sd.D.3`, `sd.D.4` | The within-class SD for different agents within the same class. Will be estimated by the model if any dose-response parameter in `class.effect` is set to `"random"`. |
#'   | `beta.1`, `beta.2`, `beta.3`, `beta.4` | The absolute value of a given dose-response parameter across the whole network (i.e. does not vary by agent/class). Will be estimated by the model for dose-response parameters (`beta.1`, `beta.2`, `beta.3`, `beta.4`) specified as `"common"` or `"random"`. |
#'   | `sd.1`, `sd.2`, `sd.3`, `sd.4` |  Between-study SD (heterogeneity) for absolute dose-response parameters (`beta.1`, `beta.2`, `beta.3`, `beta.4`) specified as `"random"`. |
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
#' | `"common"` | Implies that all studies estimate the same true \emph{absolute} effect for this dose-response parameter across the whole network. |
#' | `"random"` | Implies that all studies estimate separate true \emph{absolute} effects for this dose-response parameter, but that each of these true effects are exchangeable around a true mean effect. This approach allows for modelling of between-study heterogeneity of absolute effects. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#'
#'
#' @section Dose-response function:
#'   Several general dose-response functions are provided, but a
#'   user-defined dose-response relationship can instead be used.
#'
#'   Built-in dose-response functions are:
#'   * `"linear"`: `beta.1` refers to the gradient
#'   * `"exponential"`: `beta.1` refers to the rate of gain/decay
#'   * `"emax"` (emax without a Hill parameter): `beta.1` refers to
#'   Emax parameter, `beta.2` refers to ET50 parameter
#'   * `"emax.hill"` (emax with a Hill parameter): `beta.1` refers to Emax parameter, `beta.2` refers
#'   to ET50 parameter, `beta.3` refers to Hill parameter
#'   * `"rcs"` restricted cubic splines with knot number/location defined by `knot`.`beta.1` refers to the
#'   first spline coeffficient, `beta.2` to the second coefficient, etc. Follows the method of \insertCite{hamza2020;textual}{MBNMAdose}
#'   * `"nonparam.up"` (monotonically increasing non-parametric dose-response relationship following
#'   the method of \insertCite{owen2015;textual}{MBNMAdose})
#'   * `"nonparam.down"` (monotonically decreasing non-parametric dose-response relationship following
#'   the method of \insertCite{owen2015;textual}{MBNMAdose})
#'   * `"user"` (user-defined function: `user.fun` must be specified in arguments)
#'
#'
#'   As of version 0.2.5, separate dose-response functions can be specified for
#'   different agents in the network by passing a character vector with multiple elements to `fun`.
#'   Each agent in `network` is assigned the dose-response function in the corresponding element in `fun`.
#'   `fun` must therefore be the same length as the number of agents in `network`. Dose-response parameters
#'   `beta.1`, `beta.2`, `beta.3` and `beta.4` refer to the corresponding dose-response parameters across
#'   the multiple functions in the following order: `"user"`, `"linear"`, `"exponential"`, `"emax"`, `"emax.hill"`, `"rcs"`.
#'
#'   This would mean that if `fun` included `"linear"`, `"exponential"` and `"emax"` within it then for the
#'   corresponding agents:
#'
#'    * `beta.1` would refer to linear slope parameters
#'    * `beta.2` to exponential rate of growth/decay parameters
#'    * `beta.3` to Emax parameters
#'    * `beta.4` to ED50 parameters
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
#' fun.def <- ~ (beta.1 * dose) + (beta.2 * (dose^2))
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
#'# Fit a model with restricted cubic splines and 3 knots
#' #at 10% 30% and 60% quartiles of dose ranges
#'depnet <- mbnma.network(ssri) # Using the sSRI depression dataset
#'result <- mbnma.run(depnet, fun="rcs", knots=c(0.1,0.3,0.6))
#'
#' # Fit a model with different dose-response functions for each agent
#' multidose <- mbnma.run(network, fun=c("emax", "emax", "emax", "exponential",
#'                  "emax", "emax", "exponential", "emax"))
#'
#'
#' ########## Class effects ##########
#'
#'  # Using the osteoarthritis dataset
#'  pain.df <- osteopain_2wkabs
#'
#'  Set a shared class (NSAID) only for Naproxcinod and Naproxen
#'  pain.df <- pain.df %>% mutate(
#'                 class = case_when(agent %in% c("Naproxcinod", "Naproxen") ~ "NSAID",
#'                         !agent %in% c("Naproxcinod", "Naproxen") ~ agent
#'                         )
#'                 )
#'
#'  # Run a restricted cubic spline MBNMA with a common class effect on beta.1
#'  painnet <- mbnma.network(pain.df)
#'  splines <- mbnma.run(painnet, fun = "rcs",
#'                 class.effect = list(beta.1 = "common"))
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
#' # Calculate effective number of parameters using penalized expected deviance
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
#'               pd="popt")
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
#'####### Automatically run jags until convergence is reached #########
#'
#' # Rhat of 1.08 is set as the criteria for convergence
#' #on all monitored parameters
#' conv.res <- mbnma.run(network, fun="emax",
#'               beta.1="rel", beta.2="rel", method="random",
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
                      fun="linear",
                      beta.1="rel",
                      beta.2="rel", beta.3="rel", beta.4="rel",
                      method="common",
                      class.effect=list(), UME=FALSE,
                      knots=3,
                      cor=TRUE,
                      var.scale=NULL,
                      user.fun=NULL,
                      parameters.to.save=NULL,
                      pd="pd.kl", parallel=FALSE,
                      likelihood=NULL, link=NULL,
                      priors=NULL,
                      model.file=NULL,
                      n.iter=20000, n.chains=3,
                      n.burnin=floor(n.iter/2), n.thin=max(1, floor((n.iter - n.burnin) / 1000)),
                      autojags=FALSE, Rhat=1.1, n.update=10,
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

  # Ensure pd.kl or popt not run with parallel
  if (parallel==TRUE & pd %in% c("pd.kl", "popt")) {
    warning("pd cannot be calculated using Kullback-Leibler divergence (pd=`pk.kl` or pd=`popt`) for\nmodels run in parallel. Defaulting to pd=`pv`")
    pd <- "pv"
  }

  # Ensure cor set to FALSE if multiple dose-response functions are modelled
  if (cor==TRUE & length(fun)>1) {cor <- FALSE}

  # Ensure rjags parameters make sense
  if (n.iter<=n.burnin) {
    stop(paste0("`n.iter` must be greater than `n.burnin`. `n.burnin` = ", n.burnin))
  }

  # Check if placebo has been included
  if (network$agents[1]=="Placebo" & ("Placebo_0" %in% network$treatments)) {
    plac.incl <- TRUE
  } else {
    plac.incl <- FALSE
  }


  # Switch beta parameters for wrapper parameters
  assigned.class <- class.effect
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
    wrap.params <- list(beta.1, beta.2, beta.3, beta.4)
    wrap.params <- which(sapply(wrap.params,
                                is.character))
  }

  if (is.null(model.file)) {

    # Write JAGS model code
    model <- mbnma.write(fun=fun, user.fun=user.fun,
                         beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4,
                         method=method, knots=knots,
                         class.effect=class.effect, UME=UME,
                         cor=cor, var.scale=var.scale,
                         likelihood=likelihood, link=link
    )

    # Edit beta parameters if they aren't in dose-response function
    for (i in 1:4) {
      if (!(grepl(paste0("beta\\.", i), model) | grepl(paste0("d\\.", i), model))) {
        assign(paste0("beta.", i), NULL)
      }
    }

    # Change code for if plac not included in network
    if (plac.incl==FALSE) {
      model <- gsub("\\\nfor \\(k in 2:Nagent\\)\\{ # Priors on relative treatment effects\\\n",
                    "for (k in 1:Nagent){ # Priors on relative treatment effects\n",
                    model)
      model <- gsub("\\\nfor \\(k in 2:Nclass\\)\\{ # Priors on relative class effects\\\n",
                    "for (k in 1:Nclass){ # Priors on relative class effects\n",
                    model)

      model <- gsub("s\\.beta\\.[1-4]\\[1\\] <- 0", "", model)
    }

    # Change beta.1 and beta.2 to wrapper parameters (e.g. emax, et50) if necessary
    if (!is.null(arg.params)) {
      code.params <- c("d", "beta", "sd", "tau", "D", "sd.D")
      for (i in seq_along(wrap.params)) {
        for (k in seq_along(code.params)) {
          model <- gsub(paste(code.params[k], strsplit(run.params[i], split="[.]")[[1]][2], sep="."),
                        paste(code.params[k], wrap.params[i], sep="."), model)
        }
      }

      wrap.params <- wrap.params[which(sapply(list(beta.1, beta.2, beta.3, beta.4),
                                              is.character))]

    } else {
      wrap.params <- which(sapply(list(beta.1, beta.2, beta.3, beta.4),
                                  is.character))
    }

    # Add user-defined priors to the model
    if (!is.null(priors)) {
      model <- replace.prior(priors=priors, model=model)
    }

  } else {
    warning("All parameter specifications (dose-response parameters, class effects, priors, etc.) are being overwritten by `model.file`")
    model <- model.file
  }

  # Generate default parameters to monitor
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

  # Set boolean for presence of class effects in model
  class <- ifelse(length(class.effect)>0, TRUE, FALSE)

  # Change doses to dose indices for non-parametric models
  if (any(c("nonparam.up", "nonparam.down") %in% fun)) {
    data.ab <- index.dose(network[["data.ab"]])[["data.ab"]]
  } else {
    data.ab <- network[["data.ab"]]
  }

  #### Run jags model ####

  result.jags <- mbnma.jags(data.ab, model,
                            class=class,
                            parameters.to.save=parameters.to.save,
                            likelihood=likelihood, link=link, fun=fun, knots=knots,
                            n.iter=n.iter,
                            n.thin=n.thin,
                            n.chains=n.chains,
                            n.burnin=n.burnin,
                            parallel=parallel,
                            autojags=autojags, Rhat=Rhat, n.update=n.update,
                            ...)
  result <- result.jags[["jagsoutput"]]
  jagsdata <- result.jags[["jagsdata"]]


  # Calculate model fit statistics (using differnt pD as specified)
  fitstats <- changepd(model=result, jagsdata=jagsdata, pd=pd, likelihood=likelihood, type="dose")
  result$BUGSoutput$pD <- fitstats$pd
  result$BUGSoutput$DIC <- fitstats$dic

  # Add variables for other key model characteristics (for predict and plot functions)
  model.arg <- list("parameters.to.save"=assigned.parameters.to.save,
                    "fun"=fun, "user.fun"=user.fun,
                    "jagscode"=result.jags$model,
                    "jagsdata"=jagsdata,
                    "beta.1"=beta.1, "beta.2"=beta.2,
                    "beta.3"=beta.3,
                    "beta.4"=beta.4,
                    "method"=method,
                    "likelihood"=likelihood, "link"=link,
                    "class.effect"=assigned.class,
                    "knots"=knots,
                    "cor"=cor,
                    "var.scale"=var.scale,
                    "parallel"=parallel, "pd"=pd,
                    "priors"=get.prior(model), "arg.params"=arg.params)
  result[["model.arg"]] <- model.arg
  result[["type"]] <- "dose"
  result[["network"]] <- network

  if (length(class.effect)>0) {
    result[["classes"]] <- network[["classes"]]
  }

  # Remove dose-response parameters for agents if multi dr functions are used
  if (!is.list(result)) {
    result <- cutjags(result)
  }

  if (!("error" %in% names(result))) {
    class(result) <- c("mbnma", class(result))
  }

  return(result)

}





mbnma.jags <- function(data.ab, model,
                       class=FALSE,
                       likelihood=NULL, link=NULL, fun=NULL, knots=3,
                       nodesplit=NULL,
                       warn.rhat=FALSE, parallel=FALSE,
                       autojags=FALSE, Rhat=1.1, n.update=10,
                       ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertDataFrame(data.ab, add=argcheck)
  checkmate::assertCharacter(model, any.missing=FALSE, len=1, add=argcheck)
  checkmate::assertLogical(parallel, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertLogical(class, len=1, null.ok=FALSE, any.missing=FALSE, add=argcheck)
  checkmate::assertCharacter(fun, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertNumeric(nodesplit, len=2, null.ok=TRUE, add=argcheck)
  checkmate::assertLogical(autojags, null.ok=FALSE, add=argcheck)
  checkmate::assertNumeric(Rhat, lower=1, add=argcheck)
  checkmate::assertNumeric(n.update, lower=1, add=argcheck)
  checkmate::reportAssertions(argcheck)

  args <- list(...)

  # Get data into JAGS list format
  if (is.null(likelihood) & is.null(link)) {
    # For MBNMAtime (function is equivalent in both packages)
    jagsdata <- getjagsdata(data.ab, class=class)
  } else {
    # For MBNMAdose
    jagsdata <- getjagsdata(data.ab, class=class,
                            likelihood=likelihood, link=link, fun=fun, knots=knots,
                            nodesplit=nodesplit)
  }


  # Add variable for maxtime\maxdose to jagsdata if required for non-parametric functions
  if (grepl("maxtime", model)) {
    jagsdata[["maxtime"]] <- max(data.ab$time)
  } else if (grepl("maxdose", model)) {
    jagsdata[["maxdose"]] <- index.dose(data.ab)[["maxdose"]]

    # Generate monotonically increasing/decreasing initial values
    # Check for user-defined initial values
    if (!("inits" %in% names(args))) {
      if (grepl("T\\(d\\.1\\[c-1,k\\],\\)", model)) {
        fun="nonparam.up"
      } else if (grepl("T\\(,d\\.1\\[c-1,k\\]\\)", model)) {
        fun=="nonparam.down"
      }

      args$inits <- gen.inits(jagsdata, fun=fun, n.chains=args$n.chains)
    }
  }

  # Add parameter to monitor for direct evidence in node-splitting
  if (!is.null(nodesplit)) {
    model <- add.nodesplit(model=model)

    if ("parameters.to.save" %in% names(args)) {
      args[["parameters.to.save"]] <- append(args[["parameters.to.save"]], "direct")
    }
  }

  # Drop dose from jagsdata in spline models
  dosedat <- jagsdata[["dose"]]
  if (length(fun)==1) {
    if (fun %in% c("rcs", "ns", "bs")) {
      jagsdata[["dose"]] <- NULL
    }
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
  cat(model,file=tmps)
  close(tmps)

  out <- tryCatch({
    if (parallel==FALSE) {
      result <- do.call(R2jags::jags, c(args, list(data = jagsvars,
                                                   model.file = tmpf)))

      # AUtomatically update
      if (autojags==TRUE) {
        result <- R2jags::autojags(result, Rhat=Rhat, n.update=n.update, n.iter=1000, refresh=100)
      } else if (autojags==FALSE) {
          result <- result
        }
    } else if (parallel==TRUE) {
      if (autojags==TRUE) {
        stop("autojags=TRUE cannot be used with parallel=TRUE")
      }

      # Run jags in parallel
      result <- do.call(R2jags::jags.parallel, c(args, list(data=jagsvars, model.file=tmpf)))
    }
  },
  error=function(cond) {
    message(cond)
    return(list(error=cond))
  }
  )

  # Gives warning if any rhat values > 1.02
  if (warn.rhat==TRUE) {
    if (!("error" %in% names(out))) {
      rhat.warning(out)
    }
  }

  jagsdata[["dose"]] <- dosedat # add dose back to jagsdata if spline models were used

  return(list("jagsoutput"=out, "jagsdata"=jagsdata, "model"=model))
}




#' Generate initial values for non-parametric dose-response functions
#'
#' Ensures model runs properly
#'
gen.inits <- function(jagsdata, fun, n.chains) {
  if ("maxdose" %in% names(jagsdata)) {
    inits <- list()
    for (i in 1:n.chains) {
      inits[[length(inits)+1]] <- list("d.1"=gen.init(jagsdata, fun))
    }
  }
  return(inits)
}





gen.init <- function(jagsdata, fun) {
  sapply(jagsdata$maxdose, FUN=function(x, direction=fun) {
   val <- x-1

   if (direction=="nonparam.up") {
     start <- stats::runif(1,0,3)
     c(NA, seq(start, start+stats::runif(1,0,3), length.out=val), rep(NA, max(jagsdata$maxdose)-(val+1)))
   } else if (direction=="nonparam.down") {
     start <- stats::runif(1,-3,0)
     c(NA, seq(start, start-stats::runif(1,0,3), length.out=val), rep(NA, max(jagsdata$maxdose)-(val+1)))
   }
  })
}




#' Automatically generate parameters to save for a dose-response MBNMA model
#'
#' Identical to `gen.parameters.to.save()` in `MBNMAtime`
#'
#' @param model.params A character or numeric vector containing the names of the
#' dose-response parameters in the model
#' @param model A JAGS model written as a character object
gen.parameters.to.save <- function(model.params, model) {
  # model.params is a vector (numeric/character) of the names of the dose-response parameters in the model
  #e.g. c(1, 2, 3) or c("emax", "et50")
  # model is a JAGS model written as a character object

  checkmate::assertCharacter(model, len=1)

  model.params <- as.character(model.params)

  # Set some automatic parameters based on the model code
  parameters.to.save <- vector()
  for (i in seq_along(model.params)) {
    if (grepl(paste0("\\\nd\\.", model.params[i], "\\[k(,c)?\\] ~"), model)==TRUE |
        grepl(paste0("\\\nd\\.", model.params[i], "\\[k(,c)?\\] <- mult\\["), model)==TRUE) {
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

  # For non-monotonic DR functions
  if (grepl(paste0("\\\nd\\.1\\[c,k\\] ~"), model)==TRUE) {
    parameters.to.save <- append(parameters.to.save, "d.1")
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
#' goutnet <- mbnma.network(GoutSUA_2wkCFB)
#' nma <- nma.run(goutnet, method="common", drop.discon=FALSE)
#'
#' # Run an Unrelated Mean Effects (UME) inconsistency model on triptans dataset
#' tripnet <- mbnma.network(HF2PPITT)
#' ume <- nma.run(tripnet, method="random", UME=TRUE)
#' }
#'
#' @export
nma.run <- function(network, method="common", likelihood=NULL, link=NULL, priors=NULL,
                    warn.rhat=TRUE, n.iter=20000, drop.discon=TRUE, UME=FALSE, pd="pd.kl", ...) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(network, "mbnma.network", add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), add=argcheck)
  checkmate::assertLogical(warn.rhat, add=argcheck)
  checkmate::assertIntegerish(n.iter, null.ok = TRUE, add=argcheck)
  checkmate::assertLogical(drop.discon, add=argcheck)
  checkmate::assertLogical(UME, add=argcheck)
  checkmate::assertList(priors, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(pd, choices=c("pv", "pd.kl", "plugin", "popt"), null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  args <- list(...)

  # Check/assign link and likelihood
  likelink <- check.likelink(network$data.ab, likelihood=likelihood, link=link)
  likelihood <- likelink[["likelihood"]]
  link <- likelink[["link"]]

  #### Write model for NMA ####
  model <- write.nma(method=method, likelihood=likelihood, link=link, UME=UME)

  #### Add priors ####
  if (!is.null(priors)) {
    model <- replace.prior(priors=priors, model=model)
  }

  #### Parameters ####
  parameters.to.save <- c("d", "totresdev")
  if (method=="random") {
    parameters.to.save <- append(parameters.to.save, "sd")
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

  #### Prepare data ####
  data.ab <- network$data.ab
  trt.labs <- network$treatments

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
  cat(model,file=tmps)
  close(tmps)

  out <- tryCatch({
    result <- do.call(R2jags::jags, c(args, list(data=jagsvars, model.file=tmpf, n.iter=n.iter,
                                                 parameters.to.save=parameters.to.save)))
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

  # Calculate model fit statistics (using differnt pD as specified)
  fitstats <- changepd(model=result, jagsdata=jagsdata, pd=pd, likelihood=likelihood, type="dose")
  out$BUGSoutput$pD <- fitstats$pd
  out$BUGSoutput$DIC <- fitstats$dic

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
#' @inheritSection mbnma.run Dose-response parameter arguments
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' tripnet <- mbnma.network(HF2PPITT)
#'
#' # Fit a linear dose-response MBNMA with random treatment effects
#' linear <- mbnma.linear(tripnet, slope="rel", method="random")
#'
#' # Fit a linear dose-response MBNMA using a cloglog link function
#' linear <- mbnma.linear(tripnet, slope="rel", link="cloglog")
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from linear dose-response MBNMA
#' linear <- mbnma.linear(tripnet, slope="rel", method="random")
#' print(linear$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)")
#'
#' linear <- mbnma.linear(tripnet, slope="rel", method="random",
#'               priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' linear <- mbnma.linear(tripnet, slope="rel", method="random",
#'               n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' linear <- mbnma.linear(tripnet, slope="rel", method="random",
#'               pd="plugin")
#'
#' # Calculate effective number of parameters via penalized expected deviance
#' linear <- mbnma.linear(tripnet, slope="rel", method="random",
#'               pd="popt")
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
                         class.effect=list(), UME=FALSE,
                         cor=TRUE,
                         var.scale=NULL,
                         parameters.to.save=NULL,
                         pd="pd.kl", parallel=FALSE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  # Assign corresponding run and wrapper parameters
  arg.params <- list(
    wrap.params=c("slope"),
    run.params=c("beta.1")
  )

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="linear", user.fun=NULL,
                      model.file=NULL,
                      beta.1=slope,
                      method=method,
                      class.effect=class.effect, UME=UME,
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
#' @inheritSection mbnma.run Dose-response parameter arguments
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' tripnet <- mbnma.network(HF2PPITT)
#'
#' # Fit a exponential dose-response MBNMA with random treatment effects
#' exponential <- mbnma.exponential(tripnet, lambda="rel", method="random")
#'
#' # Fit a exponential dose-response MBNMA using a cloglog link function
#' exponential <- mbnma.exponential(tripnet, lambda="rel", link="cloglog")
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from exponential dose-response MBNMA
#' exponential <- mbnma.exponential(tripnet, lambda="rel", method="random")
#' print(exponential$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)")
#'
#' exponential <- mbnma.exponential(tripnet, lambda="rel", method="random",
#'                    priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' exponential <- mbnma.exponential(tripnet, lambda="rel", method="random",
#'                    n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' exponential <- mbnma.exponential(tripnet, lambda="rel", method="random",
#'                    pd="plugin")
#'
#' # Calculate effective number of parameters via penalized expected deviance
#' exponential <- mbnma.exponential(tripnet, lambda="rel", method="random",
#'                    pd="popt")
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
                         class.effect=list(), UME=FALSE,
                         cor=TRUE,
                         var.scale=NULL,
                         parameters.to.save=NULL,
                         pd="pd.kl", parallel=FALSE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  # Assign corresponding run and wrapper parameters
  arg.params <- list(
    wrap.params=c("lambda"),
    run.params=c("beta.1")
  )

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="exponential", user.fun=NULL,
                      model.file=NULL,
                      beta.1=lambda,
                      method=method,
                      class.effect=class.effect, UME=UME,
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
#' @inheritSection mbnma.run Dose-response parameter arguments
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' tripnet <- mbnma.network(HF2PPITT)
#'
#' # Fit an Emax dose-response MBNMA with random treatment effects on Emax and ED50
#' emax <- mbnma.emax(tripnet, emax="rel", ed50="rel", method="random")
#'
#' # Fit an Emax dose-response MBNMA with common treatment effects on Emax and
#' #a single common parameter estimated for ED50
#' emax <- mbnma.emax(tripnet, emax="rel", ed50="common", method="common")
#'
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from an Emax function with random relative effects on Emax and ED50
#' emax <- mbnma.emax(tripnet,
#'             emax="rel", ed50="rel", method="random")
#' print(emax$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(sd = "dnorm(0,0.5) T(0,)",
#'                  inv.R = "dwish(Omega[,],100)")
#'
#' emax <- mbnma.emax(tripnet,
#'             emax="rel", ed50="rel", method="random",
#'             priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' emax <- mbnma.emax(tripnet, emax="rel", ed50="rel",
#'             n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' emax <- mbnma.emax(tripnet, emax="rel", ed50="rel",
#'             pd="plugin")
#'
#' # Calculate effective number of parameters via penalized expected deviance
#' emax <- mbnma.emax(tripnet, emax="rel", ed50="rel",
#'             pd="popt")
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
                         class.effect=list(), UME=FALSE,
                         cor=TRUE,
                         var.scale=NULL,
                         parameters.to.save=NULL,
                         pd="pd.kl", parallel=FALSE,
                         likelihood=NULL, link=NULL,
                         priors=NULL,
                         arg.params=NULL, ...)
{

  # Assign corresponding run and wrapper parameters
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
                      class.effect=class.effect, UME=UME,
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
#' @inheritSection mbnma.run Dose-response parameter arguments
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' \donttest{
#' # Using the triptans data
#' tripnet <- mbnma.network(HF2PPITT)
#'
#' # Fit an Emax (with Hill parameter) dose-response MBNMA with random treatment
#' #effects on Emax, ED50 and Hill
#' emax.hill <- mbnma.emax.hill(tripnet, emax="rel", ed50="rel", hill="rel",
#'                  method="random")
#'
#' # Fit an Emax (with Hill parameter) dose-response MBNMA with common treatment
#' #effects on Emax, a single random parameter estimated for ED50
#' #and a single common parameter estimated for Hill
#' emax.hill <- mbnma.emax.hill(tripnet, emax="rel", ed50="random", hill="common",
#'                  method="common")
#'
#' # Assign a specific numerical value for Hill parameter
#' emax.hill <- mbnma.emax.hill(tripnet, emax="rel", ed50="rel", hill=5)
#'
#'
#'
#' ####### Priors #######
#'
#' # Obtain priors from an Emax (with Hill parameter) function with
#' #relative effects on Emax and ED50 and a single common parameter for Hill
#' emax.hill <- mbnma.emax.hill(tripnet,
#'                  emax="rel", ed50="rel", hill="common", method="common")
#' print(emax.hill$model.arg$priors)
#'
#' # Set new more informative prior distributions
#' newpriors <- list(beta.hill = "dnorm(0,0.5) T(,0)")
#'
#' emax.hill <- mbnma.emax.hill(tripnet,
#'                  emax="rel", ed50="rel", hill="common", method="common",
#'                  priors=newpriors)
#'
#'
#' ########## Sampler options ##########
#'
#' # Change the number of MCMC iterations, the number of chains, and the thin
#' emax.hill <- mbnma.emax.hill(tripnet, emax="rel", ed50="rel", hill=2,
#'                  n.iter=5000, n.thin=5, n.chains=4)
#'
#' # Calculate effective number of parameters via plugin method
#' emax.hill <- mbnma.emax.hill(tripnet, emax="rel", ed50="rel", hill=2,
#'                  pd="plugin")
#'
#' # Calculate effective number of parameters via penalized expected deviance
#' emax.hill <- mbnma.emax.hill(tripnet, emax="rel", ed50="rel", hill=2,
#'                  pd="popt")
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
                       class.effect=list(), UME=FALSE,
                       cor=TRUE,
                       var.scale=NULL,
                       parameters.to.save=NULL,
                       pd="pd.kl", parallel=FALSE,
                       likelihood=NULL, link=NULL,
                       priors=NULL,
                       arg.params=NULL, ...)
{

  # Assign corresponding run and wrapper parameters
  arg.params <- list(
    wrap.params=c("emax", "ed50", "hill"),
    run.params=c("beta.1", "beta.2", "beta.3", "beta.4")
  )

  result <- mbnma.run(network=network, parameters.to.save=parameters.to.save,
                      fun="emax.hill", user.fun=NULL,
                      model.file=NULL,
                      beta.1=emax,
                      beta.2=ed50,
                      beta.3=hill,
                      method=method,
                      class.effect=class.effect, UME=UME,
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
#' standard errors or covariance matrices. Currently only functions with univariate likelihoods. Function
#' is identical in MBNMAdose and MBNMAtime packages.
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

  # Update JAGS model for additional iterations
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








#' Update model fit statistics depending on calculation for pD
#'
#' @param model A model object of class `"rjags"`
#' @param jagsdata A list object containing data used to estimate `model`
#' @param type Can take either `"dose"` for a dose-response MBNMA or `"time"` for a
#' time-course MBNMA (this accounts for multiple observations within an arm)
#'
#' @return A list containing `pd` (effective number of parameters calculated using the method
#' specified in arguments), `deviance` (the posterior mdedian of the total residual deviance)
#' and `dic` (the model DIC)
#'
#' @inheritParams mbnma.run
changepd <- function(model, jagsdata=NULL, pd="pv", likelihood=NULL, type="dose") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertList(jagsdata, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(pd, choices=c("pv", "pd.kl", "plugin", "popt"), null.ok=FALSE, add=argcheck)
  checkmate::assertCharacter(likelihood, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(type, choices=c("dose", "time"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  pdout <- model$BUGSoutput$pD

  # pd by MCMC sampling methods (pd.kl or popt)
  if (pd == "pd.kl" | pd == "popt") {
    if (pd=="pd.kl") {
      temp <- rjags::dic.samples(model$model, n.iter=1000, type="pD")
    } else if (pd=="popt") {
      temp <- rjags::dic.samples(model$model, n.iter=1000, type="popt")
    }
    pdout <- sum(temp$penalty)

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
    pdout <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                    theta.result=model$BUGSoutput$mean$psi, resdev.result=model$BUGSoutput$mean$resdev,
                    likelihood=likelihood, type=type)
  }

  # Recalculate DIC so it is adjusted for choice of pD
  dicout <- pdout + model$BUGSoutput$median$deviance
  devout <- model$BUGSoutput$median$deviance

  return(list("pd"=pdout, "deviance"=devout, "dic"=dicout))
}
