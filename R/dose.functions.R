######################################
###### Dose-response functions #######
######################################


#' Exponential dose-response function
#'
#' Similar parameterisation to the Emax model but with non-asymptotic maximal effect (Emax). Can fit
#' a 1-parameter (Emax only) or 2-parameter model (includes onset parameter that controls the curvature of
#' the dose-response relationship)
#'
#' 1-parameter model:
#' \eqn{emax\times{(1-exp(-x))}}
#'
#' 2-parameter model:
#' \eqn{emax\times{(1-exp(onset*-x))}}
#'
#' where emax is the maximum efficacy of an agent and rate is the speed
#'
#' @param emax Pooling for Emax  parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param onset Pooling for onset parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @inheritParams demax
#'
#' @return An object of `class("dosefun")`
#'
#' @details
#'
#' Dose-response parameter arguments:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Single parameter exponential function is default
#' dexp()
#'
#' # Two parameter exponential function
#' dexp(onset="rel")
#'
#' @export
dexp <- function(emax="rel", onset=NULL, p.expon=FALSE) {

  # Run checks
  params <- list(emax=emax, onset=onset)
  for (i in seq_along(params)) {
    err <- TRUE
    if (!is.null(params[[i]])) {
      if (length(params[[i]])==1) {
        if (any(c("rel", "common", "random") %in% params[[i]])) {
          err <- FALSE
        } else if (is.numeric(params[[i]])) {
          err <- FALSE
        }
      }
    } else if (names(params)[i]=="onset") {
      err <- FALSE
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'rel', 'common', 'random' or be assigned a numeric value"))
    }
  }


  # Define function
  if (is.null(onset)) {
    fun <- ~ emax * (1 - exp(-dose))
    jags <- "s.beta.1 * (1 - exp(- dose[i,k]))"

  } else {
    if (p.expon==TRUE) {
      fun <- ~ emax * (1 - exp(exp(onset) * -dose))
      jags <- "s.beta.1 * (1 - exp(exp(s.beta.2) * - dose[i,k]))"
    } else if (p.expon==FALSE) {
      fun <- ~ emax * (1 - exp(onset * -dose))
      jags <- "s.beta.1 * (1 - exp(s.beta.2 * - dose[i,k]))"
    }
  }

  for (i in seq_along(params)) {
    jags <- gsub(paste0("s\\.beta\\.", i), paste0("s.beta.",i,"[agent[i,k]]"), jags)
  }

  if (p.expon==TRUE) {
    f <- function(dose, beta.1, beta.2=0) {
      y <- beta.1 * (1-exp(exp(onset)*-dose))
      return(y)
    }
  } else if (p.expon==FALSE) {
    f <- function(dose, beta.1, beta.2=0) {
      y <- beta.1 * (1-exp(onset*-dose))
      return(y)
    }
  }



  # Generate output values
  paramnames <- c("emax")
  apool <- c(emax)

  if (!is.null(onset)) {
    paramnames <- append(paramnames, "onset")
    apool <- append(apool, onset)
  }
  nparam <- length(paramnames)
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="exp", fun=fun, p.expon=p.expon,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname)
  class(out) <- "dosefun"

  if (p.expon==TRUE) {
    message("'onset' parameters are on exponential scale to ensure they take positive values on the natural scale")
  }

  return(out)
}





#' Integrated Two-Component Prediction (ITP) function
#'
#' Similar parameterisation to the Emax model but with non-asymptotic maximal effect (Emax). Proposed
#' by proposed by \insertCite{fumanner;textual}{MBNMAdose}
#'
#' @param emax Pooling for Emax  parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param rate Pooling for Rate parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @inheritParams demax
#'
#' @return An object of `class("dosefun")`
#'
#' @details
#'
#' Emax represents the maximum response.
#' Rate represents the rate at which a change in the dose of the drug leads to
#' a change in the effect
#'
#' \deqn{{E_{max}}\times\frac{(1-exp(-{rate}\times{x}))}{(1-exp(-{rate}\times{max(x)}))}}
#'
#'
#' @section Dose-response parameters:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#' When relative effects are modelled on more than one dose-response parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#' This prior can be made slightly more informative by specifying the scale matrix `omega`
#' and by changing the degrees of freedom of the inverse-Wishart prior
#' using the `priors` argument in `mbnma.run()`.
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Model a common effect on rate
#' ditp(emax="rel", rate="common")
#'
#' @export
ditp <- function(emax="rel", rate="rel", p.expon=FALSE) {

  params <- list(emax=emax, rate=rate)
  for (i in seq_along(params)) {
    if (!is.null(params[[i]])) {
      err <- TRUE
      if (length(params[[i]])==1) {
        if (any(c("rel", "common", "random") %in% params[[i]])) {
          err <- FALSE
        } else if (is.numeric(params[[i]])) {
          err <- FALSE
        }
      }
      if (err) {
        stop(paste0(names(params)[i], " must take either 'rel', 'common', 'random' or be assigned a numeric value"))
      }
    }
  }

  # Define dose-response function
  if (p.expon==TRUE) {
    fun <- ~ emax * (1 - exp(-exp(rate)*dose)) / (1 - exp(-exp(rate)*max(dose)))
    jags <- "s.beta.1 * ((1-exp(-exp(s.beta.2)*dose[i,k])) / (1-exp(-exp(s.beta.2)*maxdose)))"

  } else if (p.expon==FALSE) {
    fun <- ~ emax * (1 - exp(-rate*dose)) / (1 - exp(-rate*max(dose)))
    jags <- "s.beta.1 * ((1-exp(-s.beta.2*dose[i,k])) / (1-exp(-s.beta.2*maxdose)))"
  }

  for (i in seq_along(params)) {
    jags <- gsub(paste0("s\\.beta\\.", i), paste0("s.beta.",i,"[agent[i,k]]"), jags)
  }


  f <- function(dose, beta.1, beta.2) {
    y <- beta.1 * (1-exp(-beta.2*dose)) / (1-exp(-beta.2*max(dose)))
    return(y)
  }

  # if (emax=="rel") {
  #   #jags <- gsub("beta\\.1", "beta.1[i,k]", jags)
  # } else if (any(c("common", "random") %in% emax))) {
  #   jags <- gsub("beta\\.1", "i.beta.1[i,k]", jags)
  # }
  # if (rate=="rel") {
  #   #jags <- gsub("beta\\.2", "beta.2[i,k]", jags)
  # } else if (any(c("common", "random") %in% rate)) {
  #   jags <- gsub("beta\\.2", "i.beta.2[i,k]", jags)
  # }


  # Generate output values
  paramnames <- c("emax", "rate")
  nparam <- 2

  apool <- c(emax, rate)
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="itp", fun=fun,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname, p.expon=p.expon)
  class(out) <- "dosefun"

  if (p.expon==TRUE) {
    message("'ed50' parameters are on exponential scale to ensure they take positive values on the natural scale")
  }

  return(out)
}






#' Log-linear (exponential) dose-response function
#'
#' Modelled assuming relative effects (`"rel"`)
#'
#' \eqn{rate\times{log(x + 1)}}
#'
#' @return An object of `class("dosefun")`
#'
#' @details
#'
#' Dose-response parameter arguments:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' dloglin()
#'
#' @export
dloglin <- function() {
  rate <- "rel"

  # Define function
  fun <- ~ rate * log(dose + 1)
  jags <- paste0("s.beta.1[agent[i,k]] * log(dose[i,k] + 1)")


  f <- function(dose, beta.1) {
    y <- beta.1 * log(dose + 1)
    return(y)
  }

  # Generate output values
  paramnames <- "rate"
  nparam <- 1

  apool <- rate
  names(apool) <- paramnames
  bname <- paste0("beta.", 1:nparam)
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="loglin", fun=fun,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname)

  class(out) <- "dosefun"

  return(out)

}





#' Emax dose-response function
#'
#' @param emax Pooling for Emax  parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param ed50 Pooling for ED50 parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param hill Pooling for Hill parameter. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param p.expon A logical object to indicate whether `ed50` and `hill` parameters should be
#'  expressed within the dose-response function on an exponential scale
#'
#' @return An object of `class("dosefun")`
#'
#' @details
#'
#' Emax represents the maximum response.
#' exp(ED50) represents the dose at which 50% of the maximum response is achieved.
#' exp(Hill) is the Hill parameter, which allows for a sigmoidal function.
#'
#' Without Hill parameter:
#' \deqn{\frac{E_{max}\times{x}}{ET_{50}+x}}
#'
#' With Hill parameter:
#' \deqn{\frac{E_{max}\times{x^{hill}}}{ET_{50}\times{hill}}+x^{hill}}
#'
#'
#' @section Dose-response parameters:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#' When relative effects are modelled on more than one dose-response parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#' This prior can be made slightly more informative by specifying the scale matrix `omega`
#' and by changing the degrees of freedom of the inverse-Wishart prior
#' using the `priors` argument in `mbnma.run()`.
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Model without a Hill parameter
#' demax(emax="rel", ed50="common")
#'
#' # Model including a Hill parameter and defaults for Emax and ED50 parameters
#' demax(hill="common")
#'
#' @export
demax <- function(emax="rel", ed50="rel", hill=NULL, p.expon=FALSE) {

  # Run checks
  params <- list(emax=emax, ed50=ed50, hill=hill)
  for (i in seq_along(params)) {
    err <- TRUE
    if (!is.null(params[[i]])) {
      if (length(params[[i]])==1) {
        if (any(c("rel", "common", "random") %in% params[[i]])) {
          err <- FALSE
        } else if (is.numeric(params[[i]])) {
          err <- FALSE
        }
      }
    } else if (names(params)[i]=="hill") {
      err <- FALSE
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'rel', 'common', 'random' or be assigned a numeric value"))
    }
  }

  if (p.expon==TRUE) {
    if (!is.null(hill)) {
      fun <- ~ (emax * (dose ^ hill)) / ((exp(ed50) ^ hill) + (dose ^ hill))
      jags <- "(s.beta.1 * (dose[i,k] ^ exp(s.beta.3))) / ((exp(s.beta.2) ^ exp(s.beta.3)) + (dose[i,k] ^ exp(s.beta.3)))"
    } else if (is.null(hill)) {
      fun <- ~ (emax * dose) / (exp(ed50) + dose)
      jags <- "(s.beta.1 * dose[i,k]) / (exp(s.beta.2) + dose[i,k])"
    }
  } else {
    if (!is.null(hill)) {
      fun <- ~ (emax * (dose ^ hill)) / ((ed50 ^ hill) + (dose ^ hill))
      jags <- "(s.beta.1 * (dose[i,k] ^ s.beta.3)) / ((s.beta.2 ^ s.beta.3) + (dose[i,k] ^ s.beta.3))"
    } else if (is.null(hill)) {
      fun <- ~ (emax * dose) / (ed50 + dose)
      jags <- "(s.beta.1 * dose[i,k]) / (s.beta.2 + dose[i,k])"
    }
  }

  for (i in seq_along(params)) {
    jags <- gsub(paste0("s\\.beta\\.", i), paste0("s.beta.",i,"[agent[i,k]]"), jags)
  }


  f <- function(dose, beta.1, beta.2, beta.3) {
    y <- (beta.1 * (dose ^ beta.3) ) / ((beta.2 ^ beta.3) + (dose ^ beta.3))
    return(y)
  }


  # Generate output values
  paramnames <- c("emax", "ed50")
  nparam <- 2

  apool <- c(emax, ed50)

  if (!is.null(hill)) {
    paramnames <- append(paramnames, "hill")
    nparam <- 3
    apool <- append(apool, hill)
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="emax", fun=fun, p.expon=p.expon,
              params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname)
  class(out) <- "dosefun"

  if (p.expon==TRUE) {
    message("'ed50' parameters are on exponential scale to ensure they take positive values on the natural scale")

    if (!is.null(hill)) {
      message("'hill' parameters are on exponential scale to ensure they take positive values on the natural scale")
    }
  }

  return(out)
}






#' Polynomial dose-response function
#'
#' @param degree The degree of the polynomial - e.g. `degree=1` for linear, `degree=2` for quadratic, `degree=3` for cubic.
#' @param beta.1 Pooling for the 1st polynomial coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.2 Pooling for the 2nd polynomial coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.3 Pooling for the 3rd polynomial coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.4 Pooling for the 4th polynomial coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#'
#' @return An object of `class("dosefun")`
#'
#' @details
#' * \eqn{\beta_1} represents the 1st coefficient.
#' * \eqn{\beta_2} represents the 2nd coefficient.
#' * \eqn{\beta_3} represents the 3rd coefficient.
#' * \eqn{\beta_4} represents the 4th coefficient.
#'
#' Linear model:
#' \deqn{\beta_1{x}}
#'
#' Quadratic model:
#' \deqn{\beta_1{x} + \beta_2{x^2}}
#'
#' Cubic model:
#' \deqn{\beta_1{x} + \beta_2{x^2} + \beta_3{x^3}}
#'
#' Quartic model:
#' \deqn{\beta_1{x} + \beta_2{x^2} + \beta_3{x^3} + \beta_4{x^4}}
#'
#'
#' @section Dose-response parameters:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#' When relative effects are modelled on more than one dose-response parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#' This prior can be made slightly more informative by specifying the scale matrix `omega`
#' and by changing the degrees of freedom of the inverse-Wishart prior
#' using the `priors` argument in `mbnma.run()`.
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Linear model with random effects
#' dpoly(beta.1="rel")
#'
#' # Quadratic model dose-response function
#' # with an exchangeable (random) absolute parameter estimated for the 2nd coefficient
#' dpoly(beta.1="rel", beta.2="random")
#'
#' @export
dpoly <- function(degree=1, beta.1="rel", beta.2="rel",
                  beta.3="rel", beta.4="rel") {

  # Run checks
  checkmate::assertIntegerish(degree, lower=1, upper = 4, add=argcheck)

  params <- list(beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4)
  for (i in 1:degree) {
    err <- TRUE
    if (length(params[[i]])==1) {
      if (any(c("rel", "common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'rel', 'common', 'random' or be assigned a numeric value"))
    }
  }


  # Define dose-response function
  fun <- "beta.1 * dose"
  for (i in 1:3) {
    if (degree>i) {
      fun <- paste0(fun, " + beta.", i+1, " * (dose^", i+1, ")")
    }
  }
  jags <- gsub("dose", "dose[i,k]", fun)
  jags <- gsub("beta", "s.beta", jags)
  fun <- stats::as.formula(paste0("~", fun))

  for (i in seq_along(params)) {
    jags <- gsub(paste0("s\\.beta\\.", i), paste0("s.beta.",i,"[agent[i,k]]"), jags)
  }


  # Generate output values
  paramnames <- paste0("beta.", 1:degree)
  nparam <- degree

  apool <- vector()
  for (i in 1:nparam) {
    apool <- append(apool, get(paste0("beta.",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="poly", fun=fun, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname)
  class(out) <- "dosefun"

  return(out)
}






#' Fractional polynomial dose-response function
#'
#' @param degree The degree of the fractional polynomial as defined in  \insertCite{royston1994;textual}{MBNMAdose}
#' @param beta.1 Pooling for the 1st fractional polynomial coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.2 Pooling for the 2nd fractional polynomial coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param power.1 Value for the 1st fractional polynomial power (\eqn{\gamma_1}). Must take any numeric value in the set `-2, -1, -0.5, 0, 0.5, 1, 2, 3`.
#' @param power.2 Value for the 2nd fractional polynomial power (\eqn{\gamma_2}). Must take any numeric value in the set `-2, -1, -0.5, 0, 0.5, 1, 2, 3`.
#'
#' @return An object of `class("dosefun")`
#'
#' @details
#' * \eqn{\beta_1} represents the 1st coefficient.
#' * \eqn{\beta_2} represents the 2nd coefficient.
#' * \eqn{\gamma_1} represents the 1st fractional polynomial power
#' * \eqn{\gamma_2} represents the 2nd fractional polynomial power
#'
#' For a polynomial of `degree=1`:
#' \deqn{{\beta_1}x^{\gamma_1}}
#'
#' For a polynomial of `degree=2`:
#' \deqn{{\beta_1}x^{\gamma_1}+{\beta_2}x^{\gamma_2}}
#'
#' \eqn{x^{\gamma}} is a regular power except where \eqn{\gamma=0}, where \eqn{x^{(0)}=ln(x)}.
#' If a fractional polynomial power \eqn{\gamma} repeats within the function it is multiplied by another \eqn{ln(x)}.
#'
#'
#' @section Dose-response parameters:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#' When relative effects are modelled on more than one dose-response parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#' This prior can be made slightly more informative by specifying the scale matrix `omega`
#' and by changing the degrees of freedom of the inverse-Wishart prior
#' using the `priors` argument in `mbnma.run()`.
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # 1st order fractional polynomial a value of 0.5 for the power
#' dfpoly(beta.1="rel", power.1=0.5)
#'
#' # 2nd order fractional polynomial with relative effects for coefficients
#' # and a value of -0.5 and 2 for the 1st and 2nd powers respectively
#' dfpoly(degree=2, beta.1="rel", beta.2="rel",
#'   power.1=-0.5, power.2=2)
#'
#' @export
dfpoly <- function(degree=1, beta.1="rel", beta.2="rel",
                   power.1=0, power.2=0) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertIntegerish(degree, lower=1, upper = 2, add=argcheck)
  for (i in 1:2) {
    checkmate::assertChoice(get(paste0("power.", i)), choices=c(-2,-1,-0.5,0,0.5,1,2,3), add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  paramscoef <- list(beta.1=beta.1, beta.2=beta.2)
  paramspower <- list(power.1=power.1, beta.4=power.2)
  for (i in 1:degree) {
    err <- TRUE
    if (length(paramscoef[[i]])==1) {
      if (any(c("rel", "common", "random") %in% paramscoef[[i]])) {
        err <- FALSE
      } else if (is.numeric(paramscoef[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(paramscoef)[i], " must take either 'rel', 'common', 'random' or be assigned a numeric value"))
    }

    err <- TRUE
    if (length(paramspower[[i]])==1) {
      if (any(c("common", "random") %in% paramspower[[i]])) {
        err <- FALSE
      } else if (is.numeric(paramspower[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(paramspower)[i], " must take either 'common', 'random' or be assigned a numeric value"))
    }
  }


  # Define dose-response function
  if (degree==1) {
    fun <- ~ beta.1 * ifelse(dose>0, ifelse(power.1==0, log(dose), dose^power.1), 0)
    jags <- "s.beta.1 * ifelse(dose[i,k]>0, ifelse(s.beta.2==0, log(dose[i,k]), dose[i,k]^s.beta.2), 0)"
  } else if (degree==2) {
    fun <- ~ beta.1 * ifelse(dose>0, ifelse(beta.3==0, log(dose), dose^beta.3), 0) + (beta.2 * ifelse(beta.4==beta.3, ifelse(dose>0, ifelse(beta.4==0, log(dose)^2, (dose^beta.4) * log(dose)), 0), ifelse(dose>0, ifelse(beta.4==0, log(dose), dose^beta.4), 0)))
    jags <- "s.beta.1 * ifelse(dose[i,k]>0, ifelse(s.beta.3==0, log(dose[i,k]), dose[i,k]^s.beta.3), 0) + (s.beta.2 * ifelse(s.beta.4==s.beta.3, ifelse(dose[i,k]>0, ifelse(s.beta.4==0, log(dose[i,k])^2, (dose[i,k]^s.beta.4) * log(dose[i,k])), 0), ifelse(dose[i,k]>0, ifelse(s.beta.4==0, log(dose[i,k]), dose[i,k]^s.beta.4), 0)))"
  }

  # Set parameters
  jags <- gsub("(s\\.beta\\.[0-9]+)", "\\1[agent[i,k]]", jags)

  # Write function
  f1 <- function(dose, beta.1, beta.2) {
    if (dose>0) {
      if (beta.2==0) {
        y <- log(dose)
      } else {
        y <- dose^beta.2
      }
    } else {
      y <- 0
    }
    return(beta.1 * y)
  }

  f2 <- function(dose, beta.1, beta.2, beta.3, beta.4) {
    if (dose>0) {
      if (beta.3==0) {
        y1 <- log(dose)
      } else {
        y1 <- dose^beta.3
      }
    } else {
      y1 <- 0
    }
    y1 <- beta.1 * y1

    if (beta.4==beta.3) {
      if (dose>0) {
        if (beta.4==0) {
          y2 <- log(dose)^2
        } else {
          y2 <- dose^beta.4 * log(dose)
        }
      } else {
        if (dose >0) {
          if (beta.4==0) {
            y2 <- log(dose)
          } else {
            y2 <- dose^beta.4
          }
        }
      }
    } else {
      y2 <- 0
    }
    return(y1 + y2)
  }

  if (degree==1) {
    f <- f1
  } else if (degree==2) {
    f <- f2
  }


  # Generate output values
  paramnames <- c(paste0("beta.", 1:degree), paste0("power.", 1:degree))
  nparam <- degree*2

  apool <- vector()
  for (i in 1:degree) {
    apool <- append(apool, get(paste0("beta.",i)))
  }
  for (i in 1:degree) {
    apool <- append(apool, get(paste0("power.",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="fpoly", fun=fun,
               params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname)
  class(out) <- "dosefun"

  return(out)

}













#' Spline dose-response functions
#'
#' Used to fit B-splines, natural cubic splines, and
#' piecewise linear splines\insertCite{perperoglu2019}{MBNMAdose}.
#'
#' @param type The type of spline. Can take `"bs"` (\href{https://mathworld.wolfram.com/B-Spline.html}{B-spline}),
#'   `"ns"` (\href{https://mathworld.wolfram.com/CubicSpline.html}{natural cubic spline}),
#'   or `"ls"` (piecewise linear spline)
#' @param knots The number/location of spline internal knots. If a single number is given it indicates the number of knots (they will
#'   be equally spaced across the range of doses *for each agent*). If a numeric vector is given it indicates the location of the knots.
#' @param degree The degree of the piecewise B-spline polynomial - e.g. `degree=1` for linear, `degree=2` for quadratic, `degree=3` for cubic.
#' @param betas A vector of beta parameters corresponding to each spline coefficient. The length must be equal
#' to the number of dose-response parameters generated by the spline function. If a single value is given then this
#' specification will be applied to all beta paramters in the model. Can take `"rel"`, `"common"`, `"random"` or be
#' assigned a numeric value (see details).
#'
#' @return An object of `class("dosefun")`
#'
#'
#' @section Dose-response parameters:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#' When relative effects are modelled on more than one dose-response parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#' This prior can be made slightly more informative by specifying the scale matrix `omega`
#' and by changing the degrees of freedom of the inverse-Wishart prior
#' using the `priors` argument in `mbnma.run()`.
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Second order B spline with 2 knots and random effects on the 2nd coefficient
#' dspline(type="bs", knots=2, degree=2,
#'   betas="rel")
#'
#' # Piecewise linear spline with knots at 0.1 and 0.5 quantiles
#' # Single parameter independent of treatment estimated for 1st coefficient
#' #with random effects
#' dspline(type="ls", knots=c(0.1,0.5),
#'   betas=c("random", "rel", "rel"))
#'
#' @export
dspline <- function(type="bs", knots=1, degree=1,
                    betas="rel") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertIntegerish(degree, lower=1, upper = 4, add=argcheck)
  checkmate::assertChoice(type, choices=c("bs", "ns", "ls", "is"), add=argcheck)
  checkmate::assertNumeric(knots, null.ok=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Check knots and degrees
  x <- c(0:100)
  x <- genspline(x, spline=type, knots = knots, degree=degree)

  nparam <- ncol(x)

  # Assign values to all beta parameters if length(betas)==1
  if (nparam!=length(betas)) {
    if (length(betas)!=1) {
      stop("Different betas have been specified in betas, but length(betas)!=nparam")
    } else if (length(betas)==1) {
      betas <- rep(betas, nparam)
    }
  }

  params <- list()
  for (i in seq_along(betas)) {
    params[[paste0("beta.", i)]] <- betas[i]
    }

  err <- TRUE
  for (i in seq_along(params)) {
    if (length(params[[i]])==1) {
      if (any(c("rel", "common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
  }
  if (err) {
    stop("params must take values that are either 'rel', 'common', 'random' or a numeric value")
  }


  # Define function
  base <- "beta.1 * spline.1"
  jags <- base
  if (nparam>1) {
    for (i in 2:(nparam)) {
      temp <- gsub("1", i, base)
      jags <- paste(jags, "+", temp)
    }
  }
  fun <- stats::as.formula(paste("~", jags))
  jags <- gsub("(spline)\\.([0-9]+)", "\\1[i,k,\\2]", jags)
  jags <- gsub("beta", "s.beta", jags)


  # Define parameters
  for (i in seq_along(params)) {
    jags <- gsub(paste0("s\\.beta\\.", i, " "), paste0("s.beta.",i,"[agent[i,k]] "), jags)
  }


  # Generate output values
  paramnames <- paste0("beta.", 1:nparam)

  apool <- vector()
  for (i in 1:nparam) {
    apool <- append(apool, params[[paste0("beta.",i)]])
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- paramnames

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name=type, fun=fun, params=paramnames,
              nparam=nparam, knots=list(knots), degree=degree, jags=jags,
              apool=apool, bname=bname)
  class(out) <- "dosefun"

  return(out)

}





#' Non-parameteric dose-response functions
#'
#' Used to fit monotonically increasing non-parametric dose-response relationship following
#'   the method of \insertCite{owen2015;textual}{MBNMAdose})
#'
#' @param direction Can take either `"increasing"` or `"decreasing"` to indicate the monotonic direction
#'   of the dose-response relationship
#'
#' @return An object of `class("dosefun")`
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#' # Monotonically increasing dose-response
#' dnonparam(direction="increasing")
#'
#' # Monotonically decreasing dose-response
#' dnonparam(direction="decreasing")
#'
#' @export
dnonparam <- function(direction="increasing") {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(direction, choices=c("increasing", "decreasing"), add=argcheck)
  checkmate::reportAssertions(argcheck)


  # Define function
  fun <- stats::as.formula("~d.1")
  jags <- "d.1[dose[i,k], agent[i,k]]"


  # Generate output values
  out <- list(name="nonparam", direction=direction, fun=fun, jags=jags, apool=NA, bname=NA)
  class(out) <- "dosefun"

  return(out)

}





#' User-defined dose-response function
#'
#' @param fun A formula specifying any relationship including `dose` and
#'   one/several of: `beta.1`, `beta.2`, `beta.3`, `beta.4`.
#' @param beta.1 Pooling for the 1st coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.2 Pooling for the 2nd coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.3 Pooling for the 3rd coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#' @param beta.4 Pooling for the 4th coefficient. Can take `"rel"`, `"common"`, `"random"` or be
#'   assigned a numeric value (see details).
#'
#' @return An object of `class("dosefun")`
#'
#' @section Dose-response parameters:
#'
#' | \strong{Argument} | \strong{Model specification} |
#' | ----------------- | ---------------------------- |
#' | `"rel"` | Implies that \emph{relative} effects should be pooled for this dose-response parameter separately for each agent in the network. |
#' | `"common"` | Implies that all agents share the same common effect for this dose-response parameter. |
#' | `"random"` | Implies that all agents share a similar (exchangeable) effect for this dose-response parameter. This approach allows for modelling of variability between agents. |
#' | `numeric()` | Assigned a numeric value, indicating that this dose-response parameter should not be estimated from the data but should be assigned the numeric value determined by the user. This can be useful for fixing specific dose-response parameters (e.g. Hill parameters in Emax functions) to a single value. |
#'
#' When relative effects are modelled on more than one dose-response parameter,
#' correlation between them is automatically estimated using a vague inverse-Wishart prior.
#' This prior can be made slightly more informative by specifying the scale matrix `omega`
#' and by changing the degrees of freedom of the inverse-Wishart prior
#' using the `priors` argument in `mbnma.run()`.
#'
#'
#' @references
#'   \insertAllCited
#'
#' @examples
#'
#' dr <- ~ beta.1 * (1/(dose+1)) + beta.2 * dose^2
#' duser(fun=dr,
#'   beta.1="common", beta.2="rel")
#'
#' @export
duser <- function(fun, beta.1="rel", beta.2="rel", beta.3="rel", beta.4="rel") {


  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertFormula(fun, add=argcheck)
  checkmate::reportAssertions(argcheck)

  params <- list(beta.1=beta.1, beta.2=beta.2, beta.3=beta.3, beta.4=beta.4)
  for (i in 1:4) {
    err <- TRUE
    if (length(params[[i]])==1) {
      if (any(c("rel", "common", "random") %in% params[[i]])) {
        err <- FALSE
      } else if (is.numeric(params[[i]])) {
        err <- FALSE
      }
    }
    if (err) {
      stop(paste0(names(params)[i], " must take either 'rel', 'common', 'random' or be assigned a numeric value"))
    }
  }

  # Check user function
  user.str <- as.character(fun[2])
  if (grepl("beta\\.2", user.str)==TRUE & grepl("beta\\.1", user.str)==FALSE) {
    stop("'fun' cannot contain beta.2 if beta.1 is not present")
  } else if (grepl("beta\\.3", user.str)==TRUE & grepl("beta\\.2", user.str)==FALSE | grepl("beta\\.1", user.str)==FALSE) {
    stop("'fun' cannot contain beta.3 if beta.2 and beta.1 are not present")
  }
  if (!(grepl("dose", user.str))) {
    stop("'fun' must be a function of beta parameters and dose")
  }
  jags <- gsub("dose", "dose[i,k]", user.str)
  jags <- gsub("beta", "s.beta", jags)


  # Get number of parameters
  nparam <- 1
  if (grepl("beta\\.4", user.str)) {
    nparam <- 4
  } else if (grepl("beta\\.3", user.str)) {
    nparam <- 3
  } else if (grepl("beta\\.2", user.str)) {
    nparam <- 2
  }

  # Define parameters
  for (i in seq_along(params)) {
    jags <- gsub(paste0("s\\.beta\\.", i), paste0("s.beta.",i,"[agent[i,k]]"), jags)
  }

  # Generate output values
  paramnames <- paste0("beta.", 1:nparam)

  apool <- vector()
  for (i in 1:nparam) {
    apool <- append(apool, get(paste0("beta.",i)))
  }
  bname <- paste0("beta.", 1:nparam)

  names(apool) <- paramnames
  names(bname) <- names(bname)

  if (!any("rel" %in% apool)) {
    stop("Dose-response functions must include at least one parameter modelled using relative effects ('rel')")
  }

  out <- list(name="user", fun=fun, params=paramnames, nparam=nparam, jags=jags,
              apool=apool, bname=bname)
  class(out) <- "dosefun"

  return(out)
}




#' Agent-specific dose-response function
#'
#' Function combines different dose-response functions together to create an object containing
#' parameters for multiple dose-response functions.
#'
#' @param funs A list of objects of `class("dosefun")`, each element of which corresponds to
#'   an agent in the dataset to be modelled. The list length must be equal to the number of
#'   agents in `network$agents` used in `mbnma.run()`, and the order of the dose-response
#'   functions in the list is assumed to correspond to the same order of agents in `network$agents`.
#'
#' @return An object of `class("dosefun")`
#'
#' @examples
#'
#' funs <- c(rep(list(demax()),3),
#'           rep(list(dloglin()),2),
#'           rep(list(demax(ed50="common")),3),
#'           rep(list(dexp()),2))
#'
#' dmulti(funs)
#' @export
dmulti <- function(funs=list()) {

  # Run checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertList(funs, add=argcheck)
  for (i in seq_along(funs)) {
    checkmate::assertClass(funs[[i]], "dosefun", add=argcheck)
  }
  checkmate::reportAssertions(argcheck)

  if (any(sapply(funs, FUN=function(x) {x[["name"]]})=="nonparam")) {
    stop("dmulti() cannot currently be used with non-parametric dose-response functions")
  }



  # Identify which functions are unique (check name, degree, apool, knots)
  pos <- 0 # a vector of unique codes
  posvec <- rep(NA, length(funs))
  for (i in seq_along(funs)) {
    if (is.na(posvec[i])) {
      if (pos==0) {
        pos <- 1
      } else {
        pos <- max(posvec, na.rm = TRUE) +1
      }
      posvec[i] <- pos

      # Search for matches and make them pos
      if (i!=length(funs)) {
        for (k in (i+1):length(funs)) {
          if (isTRUE(all.equal(funs[[i]], funs[[k]]))) {
            posvec[k] <- pos
          }
        }
      }
    }
  }
  univec <- unique(posvec)

  # Then, starting from lowest...concatenate dose-response parameters...if they are unnamed then add index
  name <- vector()
  params <- vector()
  paramlist <- list()
  nparam <- 0
  jags <- vector()
  apool <- vector()
  apoollist <- list()
  bname <- vector()
  knots <- list()
  degree <- vector()
  for (i in seq_along(univec)) {
    fun <- funs[[which(posvec==univec[i])[1]]]
    name <- append(name, fun[["name"]])
    j <- fun[["jags"]]

    if ("degree" %in% names(fun)) {
      degree <- append(degree, fun[["degree"]])
    } else {
      degree <- append(degree, NA)
    }

    if ("knots" %in% names(fun)) {
      knots[[length(knots)+1]] <- fun[["knots"]][[1]]
    } else {
      knots[[length(knots)+1]] <- NA
      #knots <- append(knots, NA)
    }

    if ("p.expon" %in% names(fun)) {
      p.expon <- fun[["p.expon"]]
    } else {
      p.expon <- FALSE
    }

    for (k in seq_along(fun[["params"]])){
      j <- gsub(paste0(fun[["bname"]][k], " "), paste0("betaswap.",length(bname)+1, " "), j)
      bname <- append(bname, paste0("beta.",length(bname)+1))

      if (fun[["params"]][k] %in% params & grepl("beta", fun[["params"]][k])) {
        p <- paste0("beta.", length(params)+1)
      } else if (fun[["params"]][k] %in% params) {
        p <- paste0(fun[["params"]][k], ".", length(params)+1)
      } else {
        p <- fun[["params"]][k]
      }
      params <- append(params, p)
    }
    j <- gsub("swap", "", j)
    jags <- append(jags, j)
    apool <- append(apool, fun[["apool"]])

    pooltemp <- fun[["apool"]]
    names(pooltemp) <- params[(length(params)-length(fun[["params"]]) +1) : length(params)]
    apoollist <- c(apoollist, list(pooltemp))
  }
  names(apool) <- params
  names(bname) <- params

  out <- list(name=name, params=params, nparam=length(params), jags=jags,
              apool=apool, paramlist=apoollist, bname=bname, posvec=posvec, knots=knots, degree=degree,
              p.expon=p.expon,
              agents=names(funs))
  class(out) <- "dosefun"
  return(out)

}
