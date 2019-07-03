# Functions for writing MBNMA models
# Author: Hugo Pedder
# Date created: 2019-04-16



#' Write MBNMA dose-response model JAGS code
#'
#' Writes JAGS code for a Bayesian time-course model for model-based network
#' meta-analysis (MBNMA).
#'
#' @inheritParams mbnma.run
#'
#' @return A single long character string containing the JAGS model generated
#'   based on the arguments passed to the function.
#'
#' @inherit mbnma.run details
#'
#' @examples
#' # Write model code for a model with an exponential dose-response function,
#' # relative effects modelled on the rate of growth/decay (beta.1) with a random
#' # effects model
#' model <- mbnma.write(fun="exponential",
#'              beta.1="rel",
#'              method="random",
#'              likelihood="binomial",
#'              link="logit"
#'              )
#' cat(model)
#'
#' # Write model code for a model with an Emax dose-response function,
#' # relative effects modelled on Emax (beta.1) with a random effects model,
#' # a single parameter estimated for ED50 (beta.2) with a common effects model
#' model <- mbnma.write(fun="emax",
#'              beta.1="rel",
#'              beta.2="common",
#'              likelihood="normal",
#'              link="identity"
#'              )
#' cat(model)
#'
#' # Write model code for a model with an Emax dose-response function,
#' # relative effects modelled on Emax (beta.1) and ED50 (beta.2).
#' # Class effects modelled on ED50 with common effects
#' model <- mbnma.write(fun="emax",
#'              beta.1="rel",
#'              beta.2="rel",
#'              likelihood="normal",
#'              link="identity",
#'              class.effect=list("beta.2"="common")
#'              )
#' cat(model)
#'
#' # Write model code for a model with an Emax dose-response function,
#' # relative effects modelled on Emax (beta.1) and ED50 (beta.2) with a
#' # random effects model that automatically models a correlation between
#' # both parameters.
#' model <- mbnma.write(fun="emax",
#'              beta.1="rel",
#'              beta.2="rel",
#'              method="random",
#'              likelihood="normal",
#'              link="identity",
#'              )
#' cat(model)
#' @export
mbnma.write <- function(fun="linear",
                        beta.1="rel",
                        beta.2=NULL, beta.3=NULL,
                        method="common",
                        cor=TRUE, cor.prior="wishart",
                        var.scale=NULL,
                        class.effect=list(),
                        user.fun=NULL,
                        likelihood="binomial", link=NULL
                        ) {


  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertCharacter(user.fun, len=1, any.missing=FALSE, null.ok=TRUE, add=argcheck)
  checkmate::assertList(class.effect, unique=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)

  # Store argument values
  #argList <- as.list(match.call(expand.dots = TRUE)[-1])



  ####### VECTORS #######

  parameters <- c("beta.1", "beta.2", "beta.3")

  write.check(fun=fun, user.fun=user.fun,
              beta.1=beta.1, beta.2=beta.2, beta.3=beta.3,
              method=method, cor.prior=cor.prior,
              var.scale=var.scale,
              class.effect=class.effect)

  model <- write.model()

  # Add dose-response function
  inserts <- write.inserts()
  dosefun <- write.dose.fun(fun=fun, user.fun=user.fun)
  model <- gsub(inserts[["insert.te"]], paste0("\\1\n", dosefun, "\n\\2"), model)

  # Add likelihood
  model <- write.likelihood(model, likelihood=likelihood, link=link)

  # Add treatment delta effects
  model <- write.delta(model, beta.1=beta.1, beta.2=beta.2, beta.3=beta.3,
                      method=method)

  # Add treatment beta effects
  if (fun %in% c("nonparam.up", "nonparam.down")) {
    model <- write.beta.nonparam(model, method=method, fun=fun)
  } else {
    model <- write.beta(model, beta.1=beta.1, beta.2=beta.2, beta.3=beta.3,
                        method=method, class.effect=class.effect)
  }


  # Add correlation between dose-response parameters
  model <- write.cor(model, beta.1=beta.1, beta.2=beta.2, beta.3=beta.3,
                     method=method, cor=cor, cor.prior=cor.prior, var.scale=var.scale,
                     class.effect=class.effect)

  # Remove empty loops
  model <- write.remove.loops(model)

  return(model)
}






#' Write the basic JAGS model code for MBNMA to which other lines of model
#' code can be added
#'
#' @return A character object of JAGS model code
#'
#' @examples
#' model <- write.model()
#' cat(model)
write.model <- function() {
  model <-
    "
model{ 			# Begin Model Code

for(i in 1:NS){ # Run through all NS trials

#DR[i,1] <- 0 # Dose-response model is 0 for baseline arms
mu[i] ~ dnorm(0,0.001)

for (k in 1:narm[i]){ # Run through all arms within a study

}

resstudydev[i] <- sum(resdev[i, 1:narm[i]])

for(k in 2:narm[i]){ # Treatment effects
}
}

for (k in 2:Nagent){ # Priors on relative treatment effects
}

for (k in 2:Nclass){ # Priors on relative class effects
}

for (c in 1:(Nagent-1)) {
for (k in (c+1):Nagent) { # UME priors
}
}

totresdev <- sum(resstudydev[])

# Model ends
}
  "

  return(model)
}






#' Writes insert points for RegEx in MBNMA JAGS code
#'
#' @return A list with named elements containing character strings that match
#'   points in MBNMA JAGS code. These points can therefore be used to insert
#'   other lines of JAGS code into the correct section within the overall model
#'   code.
#'
#' @examples
#' inserts <- write.inserts()
#'
write.inserts <- function() {
  insert.start <- "(+.# Begin Model Code\n)(+.)"
  insert.study <- "(+.# Run through all NS trials\n)(+.)"
  insert.arm <- "(+.# Run through all arms within a study\n)(+.)"
  insert.te <- "(+.# Treatment effects\n)(+.)"
  insert.te.priors <- "(+.# Priors on relative treatment effects\n)(+.)"
  insert.end <- "(.+)(\n# Model ends)"
  insert.class.priors <- "(+.# Priors on relative class effects\n)(+.)"
  insert.ume.priors <- "(+.# UME priors\n)(+.)"

  return(inserts <- list("insert.start"=insert.start,
                         "insert.study"=insert.study,
                         "insert.arm"=insert.arm,
                         "insert.te"=insert.te,
                         "insert.te.priors"=insert.te.priors,
                         "insert.end"=insert.end,
                         "insert.class.priors"=insert.class.priors,
                         "insert.ume.priors"=insert.ume.priors
  ))
}




#' Write dose-response function component of JAGS code for MBNMA dose-response
#' models
#'
#' Writes a single line of JAGS code representing the dose-response function
#' component of an MBNMA dose-response model, returned as a single character
#' string.
#'
#' @inheritParams mbnma.run
#'
#' @return A single character string containing JAGS model representing the
#'   dose-response function component of an MBNMA dose-response model, generated
#'   based on the arguments passed to the function.
#'
#'
#' @examples
#' # Write a linear dose-response function
#' write.dose.fun(fun="linear")
#'
#' # Write an Emax dose-response function without a Hill parameter
#' write.dose.fun(fun="emax")
#'
#' # Write an Emax dose-response function with a Hill parameter
#' write.dose.fun(fun="emax.hill")
#'
#' # Write a user-defined dose-response function
#' doseresp <- "beta.1 + (dose ^ beta.2)"
#' write.dose.fun(fun="user", user.fun=doseresp)
write.dose.fun <- function(fun="linear", user.fun=NULL, effect="rel") {

  if (fun=="linear") {
    #DR <- "DR[i,k] <- (beta.1[agent[i,k]] * dose[i,k]) - (beta.1[agent[i,1]] * dose[i,1])"
    DR.1 <- "(beta.1[agent[i,k]] * dose[i,k])"
  } else if (fun=="exponential") {
    #DR <- "DR[i,k] <- exp(beta.1[agent[i,k]] * dose[i,k]) - exp(beta.1[agent[i,1]] * dose[i,1])"
    DR.1 <- "exp(beta.1[agent[i,k]] * dose[i,k])"
    message("Results for the rate of increase/decrease (`beta.1`) modelled on the exponential scale")
  } else if (fun=="emax") {
    #DR <- "DR[i,k] <- (beta.1[agent[i,k]] * dose[i,k] / (dose[i,k] + exp(beta.2[agent[i,k]]))) - (beta.1[agent[i,1]] * dose[i,1] / (dose[i,1] + exp(beta.2[agent[i,1]])))"
    DR.1 <- "(beta.1[agent[i,k]] * dose[i,k] / (dose[i,k] + exp(beta.2[agent[i,k]])))"
    message("Results for ED50 (`beta.2`) modelled on the exponential scale")
  } else if (fun=="emax.hill") {
    #DR <- "DR[i,k] <- (beta.1[agent[i,k]] * (dose[i,k]^exp(beta.3[agent[i,k]]))) / ((dose[i,k]^exp(beta.3[agent[i,k]])) + exp(beta.2[agent[i,k]])^exp(beta.3[agent[i,k]])) - (beta.1[agent[i,1]] * (dose[i,1]^exp(beta.3[agent[i,1]]))) / ((dose[i,1]^exp(beta.3[agent[i,1]])) + exp(beta.2[agent[i,1]])^exp(beta.3[agent[i,1]]))"
    DR.1 <- "(beta.1[agent[i,k]] * (dose[i,k]^exp(beta.3[agent[i,k]]))) / ((dose[i,k]^exp(beta.3[agent[i,k]])) + exp(beta.2[agent[i,k]])^exp(beta.3[agent[i,k]]))"
    message("Results for ED50 (`beta.2`) and Hill (`beta.3`) modelled on the exponential scale")
  } else if (fun=="nonparam.up" | fun=="nonparam.down") {
    DR.1 <- "d.1[dose[i,k], agent[i,k]]"
    message("Modelling non-parametric dose-response - arguments for dose-response parameters `beta.1`, `beta.2`, `beta.3` will be ignored")
  } else if (fun=="user") {
    DR.1 <- user.fun
    DR.1 <- gsub("(beta\\.[1-3])", "\\1[agent[i,k]]", DR.1)
    DR.1 <- gsub("(dose)", "\\1[i,k]", DR.1)

    # DR.2 <- gsub("k", "1", DR.1)
    #
    # DR <- paste0("(", DR.1, ") - (", DR.2, ")")
  }

  # Add "s." to indicate within-study betas
  DR.1 <- gsub("(beta\\.[1-3])", "s.\\1", DR.1)

  DR.2 <- gsub("k", "1", DR.1)
  DR <- paste0("(", DR.1, ") - (", DR.2, ")")

  if (effect=="rel") {
    return(paste0("DR[i,k] <- ", DR))
  } else if (effect=="abs") {
    return(paste0("DR[i,k] <- ", DR.1))
  }
}






#' Checks validity of arguments for mbnma.write
#'
#' @inheritParams mbnma.run
#'
#' @return Returns an error if any conditions are not met. Otherwise returns `NULL`.
#'
#' @details Used to check if the arguments given to mbnma.write are valid. The
#'   function will return informative errors if arguments are misspecified.
#'
write.check <- function(fun="linear",
                        beta.1=list(pool="rel", method="common"),
                        beta.2=NULL,
                        beta.3=NULL,
                        method="common",
                        UME=FALSE,
                        cor.prior="wishart",
                        var.scale=NULL,
                        user.fun=NULL,
                        class.effect=list()) {
  parameters <- c("beta.1", "beta.2", "beta.3")

  # Check fun
  dosefuns <- c("none", "linear", "exponential", "emax", "emax.hill",
                "nonparam.up", "nonparam.down", "user")
  if (is.null(fun)) {
    stop("`fun` must be assigned dose-response function(s)")
  }
  if (!all(unique(fun) %in% dosefuns)) {
    stop(paste0("`fun` must be selected from the following dose-response function(s):\n",
                paste(dosefuns, collapse=", ")))
  }

  # Run argument checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(fun, choices=c("none", "nonparam.up", "nonparam.down", "linear", "exponential", "emax", "emax.hill", "user"), null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), null.ok=FALSE, add=argcheck)
  if (method=="random") {
    checkmate::assertChoice(cor.prior, choices=c("wishart", "rho"))
    checkmate::assertNumeric(var.scale, null.ok = TRUE)
  }
  checkmate::reportAssertions(argcheck)

  # Check betas
  # Checks that beta parameters have correct format
  for (i in 1:3) {
    betaparam <- get(paste0("beta.", i))
    if (!is.null(betaparam)) {
      if (!(betaparam %in% c("rel", "common", "random") | is.numeric(betaparam))) {
        paste0("beta.", i, " must take either `rel`, `common`, `random` or a numeric value")
      }
      # if (!identical(sort(names(betaparam)), sort(c("pool", "method")))) {
      #   stop("`pool` and `method` must both be specified for each dose-response parameter")
      # }
      #
      # if (!(betaparam$pool %in% c("rel", "const"))) {
      #   stop(paste0("beta.", i, ": `pool` can take either `rel` or `const`"))
      # }
      # if (!(betaparam$method %in% c("common", "random") | is.numeric(betaparam$method))) {
      #   stop(paste0("beta.", i, ": `method` can take either `common`, `random` or a numeric value if `pool`=`const`"))
      # }
      # if (is.numeric(betaparam$method) & betaparam$pool!="const") {
      #   stop(paste0("beta.", i, ": `method` can only take a numeric value if `pool`=`const`"))
      # }
    }
  }


  # Checks for parameter classifications
  if (fun=="emax.hill" & is.null(beta.3)) {
    stop("Hill parameter (beta.3) for emax.hill function has not been specified.")
  }
  if (fun=="emax" & is.null(beta.2)) {
    stop("ET50 parameter (beta.2) for emax function has not been specified.")
  }


  #### Check user.fun ####
  if (fun!="user" & !is.null(user.fun)) {
    warning(paste0("user.fun is only applied if fun=`user`. Dose-response function used for this model will be ", fun))
  }

  if (fun=="user") {
    if (is.null(user.fun)) {
      stop("user.fun must contain a string that includes a combination of beta parameters (e.g. `beta.1`) and `dose`")
    }

    if (grepl("beta.2", user.fun)==TRUE & grepl("beta.1", user.fun==FALSE)) {
      stop("user.fun cannot contain beta.2 if beta.1 is not present")
    } else if (grepl("beta.3", user.fun)==TRUE & (grepl("beta.2", user.fun==FALSE) | grepl("beta.1", user.fun==FALSE))) {
      stop("user.fun cannot contain beta.3 if beta.2 and beta.1 are not present")
    }

    for (i in 1:3) {
      if (grepl(paste0("beta.",i), user.fun)==TRUE) {
        if(is.null(get(paste0("beta.",i)))) {
          msg <- paste0("beta.",i, " has been specified in `user.fun` dose-response function yet no arguments have been given for it")
          stop(msg)
        }
      }
    }
  }

  # Checks for class effects
  if (length(class.effect)>0) {
    # Cannot model class effects with nonparam functions
    if (fun %in% c("nonparam.up, nonparam.down")) {
      stop("Class effects cannot be used with non-parametric dose-response functions")
    }

    inclparams <- vector()
    for (i in 1:3) {
      if (!is.null(get(paste0("beta.",i)))) {
        inclparams <- append(inclparams, paste0("beta.", i))
      }
    }

    if (!all(names(class.effect) %in% inclparams)) {
      stop("`class.effect` must be a list with element names corresponding to beta parameters")
    }
    for (i in seq_along(class.effect)) {
      if (get(names(class.effect)[i])!="rel") {
        stop("Class effects can only be assigned to beta parameters modelled using relative effects (`rel`)")
      }
    }

    if (!all(class.effect %in% c("common", "random"))) {
      stop("`class.effect` elements must be either `common` or `random`")
    }
  }


  # MIGHT NEED TO INCLUDE CHECKS FOR UME

}





#' Adds sections of JAGS code for an MBNMA model that correspond to the
#' likelihood
#'
#' @inheritParams mbnma.run
#' @inheritParams write.beta
#'
#' @return A character object of JAGS MBNMA model code that includes likelihood
#'   components of the model
#'
write.likelihood <- function(model, likelihood="binomial", link=NULL) {

  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit", "log"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (likelihood=="binomial") {
    like <- "r[i,k] ~ dbin(psi[i,k], N[i,k])"
    if (is.null(link)) {link <- "logit"}
  } else if (likelihood=="normal") {
    like <- "y[i,k] ~ dnorm(psi[i,k], prec[i,k])
      prec[i,k] <- pow(se[i,k], -2)"
    if (is.null(link)) {link <- "identity"}
  } else if (likelihood=="poisson") {
    like <- "r[i,k] ~ dpois(lambda[i,k])
      lambda[i,k] <- psi[i,k] * E[i,k]"
    if (is.null(link)) {link <- "log"}
  }

  transfer <- "psi[i,k] <- theta[i,k]"

  glm <- "psi[i,k] <- theta[i,k]\ntheta[i,k] <- mu[i] + delta[i,k]"
  if (link!="identity") {
    glm <- gsub("(psi\\[i,k\\])(.+)", paste0(link, "(\\1)\\2"), glm)
  }

  inserts <- write.inserts()

  model <- gsub(inserts[["insert.arm"]], paste0("\\1", paste(like, glm, sep="\n"), "\\2"), model)


  # Add deviance contributions
  if (likelihood=="binomial") {
    resdevs <- "
rhat[i,k] <- psi[i,k] * N[i,k]
resdev[i,k] <- 2 * (r[i,k] * (log(r[i,k]) - log(rhat[i,k])) + (N[i,k] - r[i,k]) * (log(N[i,k] - r[i,k]) - log(N[i,k] - rhat[i,k])))
"
  } else if (likelihood=="normal") {
    resdevs <- "
resdev[i,k] <- pow((y[i,k] - psi[i,k]),2) * prec[i,k] # residual deviance for normal likelihood
"
  } else if (likelihood=="poisson") {
    resdevs <- "
resdev[i,k] <- 2*((lambda[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/lambda[i,k]))
"
  }

  model <- gsub(inserts[["insert.arm"]], paste0("\\1", resdevs, "\\2"), model)

  return(model)
}





#' Generate objects required for write.beta and write.beta.ref
write.beta.vars <- function() {

  delta.ref <- "\ndelta[i,1] <- 0\n"
  delta.fe <- "\ndelta[i,k] <- DR[i,k]\n"

  delta.re <- paste0("\n",
                     "delta[i,k] ~ dnorm(md[i,k], taud[i,k])\n",
                     "md[i,k] <- DR[i,k] + sw[i,k]\n",
                     "taud[i,k] <- tau * 2*(k-1)/k\n",
                     "w[i,k] <- delta[i,k] - DR[i,k]\n",
                     "sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)\n"
                     )
  multiarm <- "\nw[i,1] <- 0\n"

  sd.prior <- paste0("\nsd ~ dnorm(0,0.0025) T(0,)\n",
                "tau <- pow(sd, -2)\n"
                )


  for (i in 1:3) {

    # s.beta for reference agent
    assign(paste("s.beta", i, "ref", sep="."),
           paste0("\n", "s.beta.", i, "[1] <- 0\n")
    )

    ############# CONSTANT EFFECTS #############

    # Constant effects
    assign(paste("beta", i, sep="."),
           paste0("\nbeta.", i, " ~ dnorm(0,0.001)\n")
    )

    # Fixed constant effects
    assign(paste("beta.fe", i, sep="."),
           paste0("\ns.beta.", i, "[k] <- beta.", i, "\n")
    )

    # Random constant effects
    assign(paste("beta.re", i, sep="."),
           paste0("\ns.beta.", i, "[k] ~ dnorm(beta.", i, ", tau.", i, ")\n")
    )

    # Random constant SD prior
    assign(paste("beta.sd.prior", i, sep="."),
           paste0("\n", "sd.", i, " ~ dnorm(0,0.0025) T(0,)", "\n",
                  "tau.", i, " <- pow(sd.", i, ", -2)", "\n"
           )
    )


    ############# RELATIVE EFFECTS #############

    # Change s.beta to d
    assign(paste("btod", i, sep="."),
           paste0("\n", "s.beta.", i, "[k] <- d.", i, "[k]\n")
    )

    # Prior active treatment relative effect
    assign(paste("d.prior", i, sep="."),
           paste0("\n", "d.", i, "[k] ~ dnorm(0,0.001)", "\n")
    )

    # Class effect on active treatment relative effect - fixed
    assign(paste("d.class.fe", i, sep="."),
           paste0("\n", "d.", i, "[k]  <- D.", i, "[class[k]]", "\n")
    )
    assign(paste("d.class.re", i, sep="."),
           paste0("\n", "d.", i, "[k]  ~ dnorm(D.", i, "[class[k]], tau.D.", i, ")\n")
    )

    # Prior on active treatment class effect
    assign(paste("D.prior", i, sep="."),
           paste0("\n", "D.", i, "[k] ~ dnorm(0,0.001)", "\n")
    )

    # SD prior D class effects
    assign(paste("sd.D.prior", i, sep="."),
           paste0("\n", "sd.D.", i, " ~ dnorm(0,0.0025) T(0,)", "\n",
                  "tau.D.", i, " <- pow(sd.D.", i, ", -2)", "\n"
           )
    )


    ############# NON-PARAMETRIC EFFECTS #############

    d.zero <- "d.1[1,1] <- 0\n"
    ind <- "ind <- 1\n"

#     d.const.up <- "
# d.1[k,1] ~ dnorm(d.1[1,1],0.0001) T(d.1[1,1],)
# for (c in 2:maxdose[k]) {
# d.1[k,c] ~ dnorm(d.1[k,c-1],0.0001) T(d.1[k,c-1],)
# }
# "

    d.const.down <- "
d.1[1,k] ~ dnorm(0,0.0001) T(,d.1[1,1])
for (c in 2:maxdose[k]) {
d.1[c,k] ~ dnorm(d.1[c-1,k],0.0001) T(,d.1[c-1,k])
}
"

  d.const.up <- "
d.1[1,k] ~ dnorm(0,0.0001) T(d.1[1,1],)
for (c in 2:maxdose[k]) {
d.1[c,k] ~ dnorm(d.1[c-1,k],0.0001) T(-d.1[c-1,k],)
}
"

  }

  varnames <- ls()

  vars <- list()
  for (i in seq_along(varnames)) {
    vars[[varnames[i]]] <- get(varnames[i])
  }
  return(vars)
}








#' Adds sections of JAGS code for an MBNMA model that correspond to delta
#' parameters
#'
#' @param model A character object of JAGS MBNMA model code
#'
#' @inheritParams mbnma.run
#'
#' @return A character object of JAGS MBNMA model code that includes delta
#'   parameter components of the model
#'
write.delta <- function(model,
                       beta.1, beta.2=NULL, beta.3=NULL,
                       method="common"
) {

  inserts <- write.inserts()

  # Assign treatment effect segments
  vars <- write.beta.vars()

  # Add to model
  # Add deltas and everything
  if (method %in% c("common", "random")) {
    model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[["delta.ref"]], "\\2"), model)
    if (method=="common") {
      model <- gsub(inserts[["insert.te"]], paste0("\\1", vars[["delta.fe"]], "\\2"), model)
    } else if (method=="random") {
      model <- gsub(inserts[["insert.te"]], paste0("\\1", vars[["delta.re"]], "\\2"), model)
      model <- gsub(inserts[["insert.study"]], paste0("\\1", vars[["multiarm"]], "\\2"), model)
      model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[["sd.prior"]], "\\2"), model)
    } else {
      stop("`method` must take either `common` or `random`")
    }
  }

  return(model)
}





#' Adds sections of JAGS code for an MBNMA model that correspond to beta
#' parameters
#'
#' @param model A character object of JAGS MBNMA model code
#'
#' @inheritParams mbnma.run
#'
#' @return A character object of JAGS MBNMA model code that includes beta
#'   parameter components of the model
#'
write.beta <- function(model,
                       beta.1, beta.2=NULL, beta.3=NULL,
                       method="common",
                       class.effect=list()
) {

  inserts <- write.inserts()

  # Assign treatment effect segments
  vars <- write.beta.vars()


  for (i in 1:3) {
    # If a beta parameter has relative effects for indirect evidence calculation
    if (!is.null(get(paste0("beta.", i)))) {

      betaparam <- get(paste0("beta.", i))
      betaname <- paste0("beta.", i)

      # Add zero for s.beta reference
      model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[[paste("s.beta", i, "ref", sep=".")]], "\\2"), model)


      # RELATIVE BETA PARAMETERS
      if (betaparam=="rel") {
        # Convert s.beta to d.
        model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("btod", i, sep=".")]], "\\2"), model)

        # Add priors / class effects for d
        if (betaname %in% names(class.effect)) {
          model <- gsub(inserts[["insert.class.priors"]], paste0("\\1", vars[[paste("D.prior", i, sep=".")]], "\\2"), model)
          if (class.effect[[betaname]]=="common") {
            model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("d.class.fe", i, sep=".")]], "\\2"), model)
          } else if (class.effect[[betaname]]=="random") {
            model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("d.class.re", i, sep=".")]], "\\2"), model)
            model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("sd.D.prior", i, sep=".")]], "\\2"), model)
          }
        } else {
          model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("d.prior", i, sep=".")]], "\\2"), model)
        }


        # CONSTANT BETA PARAMETERS
      } else if (betaparam %in% c("common", "random")) {
        model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("beta", i, sep=".")]], "\\2"), model)
        if (betaparam=="common") {
          model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("beta.fe", i, sep=".")]], "\\2"), model)
        } else if (betaparam=="random") {
          model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("beta.re", i, sep=".")]], "\\2"), model)
          model <- gsub(inserts[["insert.end"]], paste0("\\1", vars[[paste("beta.sd.prior", i, sep=".")]], "\\2"), model)
        }
      } else if (is.numeric(betaparam)) {
        model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[[paste("beta.fe", i, sep=".")]], "\\2"), model)
        model <- gsub(inserts[["insert.end"]],
                      paste0("\\1\nbeta.", i, " <- ", betaparam, "\n\\2"),
                      model)
      } else {
        stop(paste0("beta.", i, " must take either `rel`, `common`, `random` or a numeric value"))
      }
    }
  }

  return(model)
}






#' Adds sections of JAGS code for a nonparametric MBNMA model that correspond to beta
#' parameters
write.beta.nonparam <- function(model, method="common", fun="nonparam.up") {
  inserts <- write.inserts()

  # Assign treatment effect segments
  vars <- write.beta.vars()

  model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[["d.zero"]], "\\2"), model)
  #model <- gsub(inserts[["insert.start"]], paste0("\\1", vars[["ind"]], "\\2"), model)

  if (fun=="nonparam.up") {
    model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[["d.const.up"]], "\\2"), model)
  } else if (fun=="nonparam.down") {
    model <- gsub(inserts[["insert.te.priors"]], paste0("\\1", vars[["d.const.down"]], "\\2"), model)
  } else {
    stop("`fun` can only take either `nonparam.up` or `nonparam.down` if write.beta.nonparam() has been called")
  }

  return(model)
}






#' Adds correlation between dose-response relative effects
#'
#' This uses a Wishart prior as default for modelling the correlation
#'
#' @inheritParams mbnma.run
#' @inheritParams write.beta
write.cor <- function(model, cor=TRUE, cor.prior="wishart", var.scale=NULL,
                      beta.1="rel", beta.2=NULL, beta.3=NULL,
                      method="random", class.effect=list()) {
  #or be a numeric vector of values to assign to rho to fill correlation
  #matrix between random effects dose-response parameters
  # (i.e. rho[2,1], rho[3,1], rho[3,2])

  #if (length(class.effect)>0 & method=="random" & (!is.null(beta.2) | !is.null(beta.3))) {
  if (length(class.effect)>0 & cor==TRUE & (!is.null(beta.2) | !is.null(beta.3))) {
    warning("Class effects cannot be modelled with correlation between time-course relative effects - correlation will be ignored")
  } else {
    # if (is.numeric(cor) & cor.prior=="wishart") {
    #   stop("Fixed (rather than estimated) values for `cor`` can only be given for `rho`, not for a wishart prior")
    # }
    # if (!is.numeric(cor) & cor!="estimate") {
    #   stop("`cor` can only take the value `estimate` or be assigned numerical value(s) corresponding to `rho`")
    # }

    #if (method=="random") {
    if (cor==TRUE) {
      corparams <- vector()
      for (i in 1:3) {
        if (!is.null(get(paste0("beta.", i)))) {
          if (get(paste0("beta.", i))=="rel") {
            corparams <- append(corparams, paste0("beta.",i))
          }
        }
      }
      sufparams <- sapply(corparams, function(x) strsplit(x, ".", fixed=TRUE)[[1]][2])
      mat.size <- length(corparams)

    } else {
      mat.size <- 0
    }

    if (mat.size>=2) {
      model <- write.cov.mat(model, sufparams=sufparams, cor="estimate", cor.prior=cor.prior,
                             var.scale=var.scale)
    }
  }

  return(model)

}





#' Function for adding covariance matrix for correlation between relative effects
#'
#' @param sufparams A numeric vector of dose-response/time-course parameter suffixes. It
#'  should be the same length as the number of relative effects (i.e. the covariance
#'  matrix size).
write.cov.mat <- function(model, sufparams, cor="estimate", cor.prior="wishart",
                          var.scale=NULL) {

  inserts <- write.inserts()

  jagswish <- "
for (r in 1:mat.size) {
d.prior[r] <- 0
}

inv.R ~ dwish(Omega[,], mat.size)

for (r in 1:(mat.size-1)) {  # Covariance matrix upper/lower triangles
for (c in (r+1):mat.size) {
Omega[r,c] <- 0   # Lower triangle
Omega[c,r] <- 0   # Upper triangle
}
}
"

  jagsrho <- "
for (r in 1:mat.size) {
d.prior[r] <- 0
R[r,r] <- 1000    # Covariance matrix diagonals
}

for (r in 1:(mat.size-1)) {  # Covariance matrix upper/lower triangles
for (c in (r+1):mat.size) {
R[r,c] <- 1000*rho[1]   # Lower triangle
R[c,r] <- 1000*rho[1]   # Upper triangle
}
}
"

  mat.size <- length(sufparams)
  for (i in seq_along(sufparams)) {
    # Change d.1[k] ~ dnorm(0,0.001)  to   d.1[k] <- d.mult[1,k]
    model <- gsub(paste0("d\\.", sufparams[i], "\\[k\\] ~ [a-z]+\\([0-9]+(\\.[0-9]+)?,[0-9]+(\\.?[0-9]+)?\\)\\\n"),
                  paste0("d.", sufparams[i], "[k] <- mult[", i, ",k]\n"),
                  model
    )
  }

  if (cor.prior=="wishart") {
    addcode <- jagswish
    model <- gsub(inserts[["insert.te.priors"]],
                  paste0("\\1mult[1:", mat.size, ",k] ~ dmnorm(d.prior[], inv.R[1:", mat.size, ", 1:", mat.size, "])\\2"),
                  model
    )

    # Check that var.scale has correct length and add omega to code
    if (is.null(var.scale)) {
      var.scale <- rep(1,mat.size)
    } else if (length(var.scale)!=mat.size) {
      stop(paste0("`var.scale` must be a numeric vector whose length is the size of the covariance matrix for dose-response parameters.\nCovariance matrix size = ", mat.size))
    }
    for (i in seq_along(var.scale)) {
      model <- gsub(inserts[["insert.end"]],
                    paste0("\\1Omega[", i, ",", i, "] <- ", var.scale[i], "\n\\2"),
                    model)
    }

  } else if (cor.prior=="rho") {
    addcode <- jagsrho
    model <- gsub(inserts[["insert.te.priors"]],
                  paste0("\\1mult[1:", mat.size, ",k] ~ dmnorm.vcov(d.prior[], R[1:", mat.size, ", 1:", mat.size, "])\\2"),
                  model
    )

    if (cor=="estimate") {
      addcode <- paste(addcode, "
                       for (m in 1:(mat.size-1)) {
                       rho[m] ~ dunif(-1,1)
                       }
                       ")
    } else if (is.numeric(cor)) {
      # Add values for rho assigned by user
      if (length(cor)!=(mat.size)-1) {
        stop("Length of numeric vector assigned to `cor` must equal the size of the correlation matrix - 1")
      }
      for (m in seq_along(cor)) {
        model <- gsub(inserts[["insert.end"]],
                      paste0("\\1rho[", m, "] <- ", cor[m]))
      }
  }
}

  addcode <- gsub("mat\\.size", mat.size, addcode)
  model <- gsub(inserts[["insert.end"]], paste0("\\1", addcode, "\\2"), model)
  return(model)
}






#' Removes any loops from MBNMA model JAGS code that do not contain any
#' expressions
#'
#' @inheritParams write.beta
#'
#' @return A character object of JAGS MBNMA model code that has had empty loops
#'   removed from it
#'
write.remove.loops <- function(model) {
  # Remove empty loops
  empty.loops <- list(
    "for \\(k in 2\\:Nagent\\)\\{ \\# Priors on relative treatment effects\\\n}", # ume.loop
    "for \\(k in 1\\:Nagent\\)\\{ \\# Priors on absolute treatment effects\\\n}", # absolute.te.loop
    "for\\(k in 2\\:narm\\[i\\]\\)\\{ \\# Treatment effects\\\n\\}", # treatment effects loop
    "for \\(k in 2\\:Nclass\\)\\{ \\# Priors on relative class effects\\\n\\}", # rel class loop
    "for \\(k in 1\\:Nclass\\)\\{ \\# Priors on absolute class effects\\\n\\}", # arm class loop
    "for \\(c in 1\\:\\(Nagent-1\\)\\) \\{\\\nfor \\(k in \\(c\\+1\\)\\:Nagent\\) \\{ \\# UME priors\\\n\\}\\\n\\}" # UME loop
  )

  for (i in seq_along(empty.loops)) {
    if (grepl(empty.loops[[i]], model)==TRUE) {
      model <- gsub(paste0("(.+)(", empty.loops[[i]], ")(.+)"),
                    paste0("\\1", "\\3"), model)
    }
  }

  return(model)
}





#' Get current priors from JAGS model code
#'
#' Identical to `get.prior()` in `MBNMAtime` package.
#' This function takes JAGS model presented as a string and identifies what
#' prior values have been used for calculation.
#'
#' @inheritParams write.beta
#'
#' @return A character vector, each element of which is a line of JAGS code
#'   corresponding to a prior in the JAGS code.
#'
#' @details Even if an MBNMA model that has not initialised successfully and
#'   results have not been calculated, the JAGS model for it is saved in
#'   `mbnma$model.arg$jagscode` and therefore priors can still be obtained.
#'   This allows for priors to be changed even in failing models, which may help
#'   solve issues with compiling or updating.
#'
#' @examples
#' # Using the triptans data
#' network <- mbnma.network(HF2PPITT)
#'
#' # Run an Emax dose-response MBNMA
#' result <- mbnma.emax(network, emax="rel", ed50="rel", method="random")
#'
#' # Obtain model prior values
#' print(result$model.arg$priors)
#'
#' # Priors when using mbnma.run with an exponential function
#' result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random")
#' print(result$model.arg$priors)
#'
#' @export
get.prior <- function(model) {

  # Run Checks
  checkmate::assertCharacter(model, len=1)

  #model <- strsplit(mbnma$model.arg$jagscode, split="\n")[[1]]
  model <- strsplit(model, split="\n")[[1]]
  #priors <- model[grep(".+~ [A-z]+\\([-?0-9]", model)]

  priorcode <- model[c(grep("^.+~ [A-z]+\\([-?0-9]", model),
                       grep("^.+~ [A-z]+\\(Omega", model))]

  priorlist <- strsplit(priorcode, split=" +?~ +?")
  priors <- list()
  for (i in seq_along(priorlist)) {
    priorname <- unlist(strsplit(priorlist[[i]][1], split="\\["))[1]
    priors[[priorname]] <- priorlist[[i]][2]
  }

  return(priors)
}






#' Replace original priors in an MBNMA model with new priors
#'
#' Identical to `get.prior()` in `MBNMAtime` package.
#'
#' This function takes new priors, as specified by the user, and adds them to
#' the JAGS code from an MBNMA model. New priors replace old priors in the JAGS
#' model.
#'
#' @inheritParams get.prior
#' @param mbnma An S3 object of class `c("mbnma", "rjags")` generated by running a
#'   dose-response MBNMA model.
#' @param priors A named list of parameter values (without indices) and
#'   replacement prior distribution values given as strings
#'   **using distributions as specified in JAGS syntax** (see examples).
#'
#' @details Values in `priors` can include any JAGS functions/distributions
#'   (e.g. censoring/truncation).
#'
#' @return A character object of JAGS MBNMA model code that includes the new
#'   priors in place of original priors
#'
replace.prior <- function(priors, model=NULL, mbnma=NULL) {
  # priors is a named list of parameter values (without indices) and replacement
  #prior values given as strings USING DISTRIBUTIONS AS SPECIFIED IN JAGS SYNTAX (i.e.
  #dnorm() is specified using mean and precision rather than mean and SD.
  #It can include JAGS functions (e.g. censoring/truncation)
  #e.g. for a half-normal SD prior list("sd.et50"="dnorm(0,0.5) T(0,)")

  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(mbnma, "mbnma", null.ok=TRUE, add=argcheck)
  checkmate::assertCharacter(model, len=1, null.ok=TRUE, add=argcheck)
  checkmate::assertList(priors, add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (!is.null(mbnma) & !is.null(model)) {
    stop("Must provide EITHER an existing MBNMA model (using `mbnma`) OR MBNMA JAGS code (using `model`)")
  }

  if (!is.null(mbnma)) {
    model <- strsplit(mbnma$model.arg$jagscode, split="\n")[[1]]
  } else if (!is.null(model)) {
    model <- strsplit(model, split="\n")[[1]]
  } else {
    stop("Must provide EITHER an existing MBNMA model (using `mbnma`) OR MBNMA JAGS code (using `model`)")
  }

  for (i in seq_along(priors)) {
    # Checks
    if (length(grep(paste0("^( +)?", names(priors)[i]), model))==0) {
      stop("Prior named ", names(priors)[i], " not found in the model code. Check priors currently present in model code using get.prior()")
      # } else if (length(grep(paste0("^( +)?", names(priors)[i]), model))>1) {
      #   stop("Prior named ", names(priors)[i], " has matched on multiple instances in the model code. Check priors currently present in model code using get.prior()")
    }

    #line <- grep(paste0("^( +)?", names(priors)[i]), model)
    line <- grep(paste0("^( +)?", names(priors)[i], ".+~"), model)
    state <- model[line]
    model[line] <- gsub("(^.+~ )(.+$)", paste0("\\1", priors[[i]]), state)
  }

  # Cut irrelevant section from JAGS code
  start <- grep("^model\\{", model)
  end <- grep("# Model ends", model) + 1

  model <- paste(model[start:end], collapse="\n")

  return(model)
}








#' Write E0 synthesis JAGS model
write.E0.synth <- function(synth="fixed", likelihood=NULL, link=NULL) {
  model <-
"
model{ 			# Begin Model Code

for(i in 1:NS){ # Run through all NS trials

for (k in 1:narm[i]){ # Run through all arms within a study

}

resstudydev[i] <- sum(resdev[i, 1:narm[i]])

}

totresdev <- sum(resstudydev[])

m.mu ~ dnorm(0,0.0001)

# Model ends
}
"

  # Add likelihood
  model <- write.likelihood(model, likelihood = likelihood, link=link)

  if (synth=="fixed") {
    mucode <- "mu[i] <- m.mu\n"

    model <- gsub("(.+Run through all NS trials\n)(.+)", paste0("\\1", mucode, "\\2"), model)
  } else if (synth=="random") {
    mucode <- "mu[i] ~ dnorm(m.mu, tau.mu)\n"
    musdcode <- "tau.mu <- pow(sd.mu, -2)\nsd.mu ~ dnorm(0,0.0025) T(0,)\n"

    model <- gsub("(.+Run through all NS trials\n)(.+)", paste0("\\1", mucode, "\\2"), model)
    model <- gsub("(.+\n)(# Model ends)", paste0("\\1", musdcode, "\\2"), model)
  }

  # Remove delta
  model <- gsub("(.+<- mu\\[i\\])( \\+ delta\\[i,k\\])", "\\1", model)

  return(model)
}



#' Write JAGS code for split NMA
write.nma <- function(method="common", likelihood="binomial", link="logit",
                      UME=FALSE) {
  model <- "
model{ 			# Begin Model Code


for(i in 1:NS){ # Run through all NS trials

mu[i] ~ dnorm(0,0.001)
delta[i,1] <- 0

for (k in 1:narm[i]){ # Run through all arms within a study

}

resstudydev[i] <- sum(resdev[i, 1:narm[i]])

for(k in 2:narm[i]){ # Treatment effects
}
}


totresdev <- sum(resstudydev[])

# Model ends
}
"

  # Add likelihood
  model <- write.likelihood(model, likelihood = likelihood, link=link)


  # Add d[1] <- 0
  if (UME==FALSE) {
    model <- gsub("(.+# Begin Model Code\n)(.+)",
                  "\\1\nd[1] <- 0\n\\2",
                  model)
  }


  # Add treatment effects
  if (method=="common") {
  te <- "
delta[i,k] <- md[i,k]
md[i,k] <- d[treatment[i,k]] - d[treatment[i,1]]
"

  } else if (method=="random") {
    # Insert w[1]=0
    model <- gsub("(.+# Run through all NS trials\n)(.+)",
                  "\\1\nw[i,1] <- 0\n\\2",
                  model)

    # Insert SD prior
    model <- gsub("(.+)(\n# Model ends.+)",
                  "\\1\ntau <- pow(sd,-2)\nsd ~ dnorm(0,0.0025) T(0,)\n\\2",
                  model)

    te <- "
delta[i,k] ~ dnorm(md[i,k], taud[i,k])
md[i,k] <- d[treatment[i,k]] - d[treatment[i,1]] + sw[i,k]
taud[i,k] <- tau *2*(k-1)/k
w[i,k] <- (delta[i,k] - d[treatment[i,k]] + d[treatment[i,1]])
sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)
"
  }

  model <- gsub("(.+# Treatment effects\n)(.+)",
                paste0("\\1", te, "\\2"),
                model)


  # Add treatment effect priors and make UME changes
  if (UME==FALSE) {
    te.prior <- "
for (k in 2:NT){ # Priors on relative treatment effects
d[k] ~ dnorm(0,0.0001)
}
"
  } else if (UME==TRUE) {

    model <- gsub("d\\[treatment\\[i,k\\]\\] (\\+|-) d\\[treatment\\[i,1\\]\\]",
                  "d[treatment[i,k],treatment[i,1]]",
                  model
                  )

    te.prior <- "
for (k in 1:NT) { d[k,k] <- 0 }

for (c in 1:(NT-1)) {
for (k in (c+1):NT) {
d[k,c] ~ dnorm(0,0.0001)
}
}
"
  }

  model <- gsub("(.+)(\n# Model ends.+)",
                paste0("\\1", te.prior, "\\2"),
                model)

  return(model)
}
