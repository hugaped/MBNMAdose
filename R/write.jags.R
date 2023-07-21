# Functions for writing MBNMA models in JAGS
# Author: Hugo Pedder
# Date created: 2021-04-20

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "studyID", "agent", "dose", "Var1", "value",
                                                        "Parameter", "fupdose", "groupvar", "y",
                                                        "network", "a", "param", "med", "l95", "u95", "value",
                                                        "Estimate", "2.5%", "50%", "97.5%", "treatment"))



#' Write MBNMA dose-response model JAGS code
#'
#' Writes JAGS code for a Bayesian time-course model for model-based network
#' meta-analysis (MBNMA).
#'
#' @inheritParams mbnma.run
#' @param cor.prior NOT CURRENTLY IN USE - indicates the prior distribution to use for the correlation/covariance
#' between relative effects. Must be kept as `"wishart"`
#' @param om a list with two elements that report the maximum relative (`"rel"`) and maximum absolute (`"abs"`) efficacies
#' on the link scale.
#' @param regress.mat A Nstudy x Ncovariate design matrix of meta-regression covariates
#'
#' @return A single long character string containing the JAGS model generated
#'   based on the arguments passed to the function.
#'
#' @inherit mbnma.run details
#'
#' @examples
#' # Write model code for a model with an exponential dose-response function,
#' # with random treatment effects
#' model <- mbnma.write(fun=dexp(),
#'              method="random",
#'              likelihood="binomial",
#'              link="logit"
#'              )
#' cat(model)
#'
#' # Write model code for a model with an Emax dose-response function,
#' # relative effects modelled on Emax with a random effects model,
#' # a single parameter estimated for ED50 with a common effects model
#' model <- mbnma.write(fun=demax(emax="rel", ed50="common"),
#'              likelihood="normal",
#'              link="identity"
#'              )
#' cat(model)
#'
#' # Write model code for a model with an Emax dose-response function,
#' # relative effects modelled on Emax and ED50.
#' # Class effects modelled on ED50 with common effects
#' model <- mbnma.write(fun=demax(),
#'              likelihood="normal",
#'              link="identity",
#'              class.effect=list("ed50"="common")
#'              )
#' cat(model)
#'
#' # Write model code for a model with an Emax dose-response function,
#' # relative effects modelled on Emax and ED50 with a
#' # random effects model that automatically models a correlation between
#' # both parameters.
#' model <- mbnma.write(fun=demax(),
#'              method="random",
#'              likelihood="normal",
#'              link="identity",
#'              )
#' cat(model)
#' @export
mbnma.write <- function(fun=dloglin(),
                        method="common",
                        regress.mat=NULL, regress.effect="common",
                        sdscale=FALSE,
                        cor=TRUE, cor.prior="wishart",
                        omega=NULL, om=list("rel"=5, "abs"=10),
                        class.effect=list(), UME=FALSE,
                        likelihood="binomial", link=NULL
) {


  # Run Checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertList(class.effect, unique=FALSE, add=argcheck)
  checkmate::reportAssertions(argcheck)


  ####### VECTORS #######

  # parameters <- c("beta.1", "beta.2", "beta.3", "beta.4")

  write.check(fun=fun,
              method=method, cor.prior=cor.prior,
              regress.mat=regress.mat, regress.effect=regress.effect,
              omega=omega, om=om,
              sdscale=sdscale,
              class.effect=class.effect, UME=UME)

  model <- write.model(UME=UME)

  # Add dose-response function
  dosefun <- write.dose.fun(fun=fun, UME=UME)
  model <- model.insert(model, pos=which(names(model)=="arm"), x=dosefun[[1]])

  # if (length(dosefun)==2) {  # For models with multiple DR functions
  #   model <- model.insert(model, pos=which(names(model)=="study"), x=dosefun[[2]])
  # }


  # Add likelihood
  model <- write.likelihood(model, likelihood=likelihood, link=link, sdscale=sdscale)

  # Add treatment delta effects
  model <- write.delta(model, method=method, om=om)

  # Add treatment beta effects
  model <- write.beta(model, fun=fun, method=method, class.effect=class.effect, UME=UME, om=om)

  # Add correlation between dose-response parameters
  model <- write.cor(model, fun=fun, method=method, cor=cor, omega=omega,
                     class.effect=class.effect, UME=UME)

  # Add regression components
  model <- write.regress(model,
                         regress.mat=regress.mat, regress.effect=regress.effect,
                         om=om)

  # Drop empty loops
  model <- remove.loops(model)

  return(model)
}







#' Checks validity of arguments for mbnma.write
#'
#' @inheritParams mbnma.run
#' @inheritParams nma.run
#' @inheritParams mbnma.write
#'
#' @return Returns an error if any conditions are not met. Otherwise returns `NULL`.
#' @noRd
#' @details Used to check if the arguments given to mbnma.write are valid. The
#'   function will return informative errors if arguments are misspecified.
#'
write.check <- function(fun=dloglin(),
                        method="common",
                        regress.mat=NULL,
                        regress.effect="common",
                        UME=FALSE, sdscale=FALSE,
                        cor.prior="wishart",
                        omega=NULL, om=list("rel"=5, "abs"=10),
                        user.fun=NULL,
                        class.effect=list()) {


  # Run argument checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(fun, "dosefun", null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), null.ok=FALSE, add=argcheck)
  checkmate::assertMatrix(regress.mat, null.ok=TRUE, add=argcheck)
  checkmate::assertChoice(regress.effect, choices=c("common", "random", "agent", "class"),
                          null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(UME, null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(sdscale, null.ok=FALSE, add=argcheck)
  if (method=="random") {
    checkmate::assertChoice(cor.prior, choices=c("wishart", "rho"))
  }
  checkmate::assertList(om, len=2, null.ok=FALSE, add=argcheck)
  checkmate::assertMatrix(omega, null.ok = TRUE)
  checkmate::reportAssertions(argcheck)

  # Checks for class effects
  if (length(class.effect)>0) {
    # Cannot model class effects with UME
    if (UME==TRUE) {
      stop("Class effects cannot be modelled with UME")
    }

    # Cannot model class effects with nonparam functions
    if (any(fun$name=="nonparam")) {
      stop("Class effects cannot be used with non-parametric dose-response functions")
    }

    # Cannot model class effects with multiple dose-response functions
    if (length(fun$name)>1) {
      stop("Class effects can only be modelled when using a single dose-response function")
    }

    # if (any(c("rcs", "ns", "bs") %in% fun$name)) {
    #   warning("Class effects applied to spline function parameters may produce\nnon-interpretable results since knot locations will differ between agents")
    # }


    if (!all(names(class.effect) %in% fun$params)) {
      stop("`class.effect` must be a list with element names corresponding to dose-response parameters")
    }
    if (!all(fun$apool[which(fun$params %in% names(class.effect))] == "rel")) {
      stop("Class effects can only be specified for dose-response parameters in 'fun' that have been modelled using relative effects ('rel')")
    }

    if (!all(class.effect %in% c("common", "random"))) {
      stop("`class.effect` elements must be either `common` or `random`")
    }
  }

  if (!is.null(omega)) {
    if (length(class.effect)>0) {
      #warning("Class effects cannot be modelled with correlation between time-course relative effects - 'omega' will be ignored")
    } else {
      err <- FALSE
      nrel <- sum(fun$apool %in% "rel")
      if (!all(dim(omega)==nrel)) {
        err <- TRUE
      }
      if (!isSymmetric(omega)) {
        err <- TRUE
      }
      if (any(eigen(omega)$values <= 0)) {
        err <- TRUE
      }
      if (err==TRUE) {
        stop("omega must be a symmetric positive definite matrix with dimensions equal to the number of\ndose-course parameters modelled using relative effects ('rel')")
      }
    }
  }

}






#' Write the basic JAGS model code for MBnMA to which other lines of model
#' code can be added
#'
#' @return A character object of JAGS model code
#' @noRd
#' @examples
#' model <- write.model()
#' cat(model)
write.model <- function(UME=FALSE) {

  model <- c(
    start= "model{ 			# Begin Model Code",
    study="for(i in 1:NS){ # Run through all NS trials",
    arm="for (k in 1:narm[i]){ # Run through all arms within a study",
    "}",
    "resstudydev[i] <- sum(resdev[i, 1:narm[i]])",
    te="for(k in 2:narm[i]){ # Treatment effects",
    "}",
    "}",
    agent.prior="for (k in 2:Nagent){ # Priors on relative agent effects",
    "}",
    class.prior="for (k in 2:Nclass){ # Priors on relative class effects",
    "}",
    trt.prior="for (k in 2:NT){ # Priors on relative treatment effects",
    "}"
    )

  if (UME==TRUE) {
    model <- append(model, c(ume.prior.ref="for (k in 1:Nagent){ # UME prior ref",
                             "}",
                             "for (c in 1:(Nagent-1)) {",
                             ume.prior="for (k in (c+1):Nagent) { # UME priors",
                             "}",
                             "}"))
  }

  model <- append(model, c("totresdev <- sum(resstudydev[])",
                           end="",
                           "# Model ends",
                           "}"
                           )
  )

return(model)
}








#' Write dose-response function component of JAGS code for MBNMA dose-response
#' models
#'
#' Writes a single line of JAGS code representing the dose-response function
#' component of an MBNMA dose-response model, returned as a single character
#' string.
#'
#' @inheritParams mbnma.run
#' @param effect Can take either `"rel"` for relative effects or `"abs"` for absolute (arm-pooled) effects
#'
#' @return A single character string containing JAGS model representing the
#'   dose-response function component of an MBNMA dose-response model, generated
#'   based on the arguments passed to the function.
#' @noRd
#'
#' @examples
#' # Write a quadratic dose-response function
#' write.dose.fun(fun=dpoly(degree=2, beta.1="rel", beta.2="rel"))
#'
#' # Write an Emax dose-response function without a Hill parameter
#' write.dose.fun(fun=demax(hill=NULL))
#'
#' # Write an Emax dose-response function with a common Hill parameter
#' write.dose.fun(fun=demax(hill="common"))
#'
#' # Write a user-defined dose-response function
#' doseresp <- ~ beta.1 + (dose ^ beta.2)
#' write.dose.fun(fun=duser(fun=doseresp, beta.1="rel", beta.2="rel"))
write.dose.fun <- function(fun=dloglin(), effect="rel", UME=FALSE) {

  DR.1 <- fun$jags

  DR.2 <- gsub("k", "1", DR.1)

  if (UME==FALSE) {
    DR <- paste0("(", DR.1, ") - (", DR.2, ")")
  } else if (UME==TRUE) {
    DR <- DR.1
    DR <- gsub("(agent\\[i\\,k\\])", "\\1, agent[i,1]", DR)
  }

  # Return final DR function depending on number of DR functions and UME
  if (length(DR.1)==1 | UME==TRUE) {
    if (effect=="rel") {
      return(list(paste0("DR[i,k] <- ", DR)))
    } else if (effect=="abs") {
      return(list(paste0("DR[i,k] <- ", DR.1)))
    }
  } else if (length(fun$jags)>1) {
    drmult <- "DR[i,k] <- DRmult[i,k,f[i,k]]"
    for (i in seq_along(fun$jags)) {
      if (effect=="rel") {
        drmult <- append(drmult, paste0("DRmult[i,k,", i, "] <- ", DR.1[i], " - ", DR.2[i]))
      } else if (effect=="abs") {
        drmult <- append(drmult, paste0("DRmult[i,k,", i, "] <- ", DR.1[i]))
      }
    }
    drmult <- list(drmult,
                   paste0("DR2[i] <- ", DR.2)
                   )
    return(drmult)
  }
}






#' Insert element into model vector at desired location
#'
#' @noRd
model.insert <- function(a, pos, x){

  if (pos>length(a)) {
    stop("'pos' cannot be greater than length(a)")
  }
  start <- a[1:pos]
  end <- a[(pos+1):length(a)]
  return(c(start, x, end))
}





#' Adds sections of JAGS code for an MBNMA model that correspond to the
#' likelihood
#'
#' @inheritParams mbnma.run
#' @inheritParams write.beta
#' @noRd
#' @return A character object of JAGS MBNMA model code that includes likelihood
#'   components of the model
#'
write.likelihood <- function(model, likelihood="binomial", link=NULL, sdscale=FALSE) {

  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit", "log", "smd"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (likelihood=="binomial") {
    like <- "r[i,k] ~ dbin(psi[i,k], n[i,k])"
    if (is.null(link)) {link <- "logit"}
  } else if (likelihood=="normal") {
    like <- c("y[i,k] ~ dnorm(psi[i,k], prec[i,k])",
              "prec[i,k] <- pow(se[i,k], -2)")
    if (is.null(link)) {link <- "identity"}
  } else if (likelihood=="poisson") {
    like <- c("r[i,k] ~ dpois(lambda[i,k])",
              "lambda[i,k] <- psi[i,k] * E[i,k]")
    if (is.null(link)) {link <- "log"}
  }

  transfer <- "psi[i,k] <- theta[i,k]"

  glm <- c("psi[i,k] <- theta[i,k]",
           glm="theta[i,k] <- mu[i] + delta[i,k]")
  if (!link %in% c("identity", "smd")) {
    glm <- gsub("(psi\\[i,k\\])(.+)", paste0(link, "(\\1)\\2"), glm)
  }

  if (link=="smd") {
    if (likelihood!="normal") {
      stop("link='smd' can only be used with likelihood='normal'")
    }
    glm[1] <- paste0(glm[1], " * pool.sd[i]")

    # Add SMD components
    if (sdscale==FALSE) {
      smd.sub <- c(
        "sd.study[i,k] <- se[i,k] * pow(n[i,k],0.5)",
        "nvar[i,k] <- (n[i,k]-1) * pow(sd.study[i,k],2)"
      )
      model <- model.insert(model, pos=which(names(model)=="arm"), x=smd.sub)

      pool.sd <- c(
        "df[i] <- sum(n[i,1:narm[i]]) - narm[i]",
        "pool.var[i] <- sum(nvar[i,1:narm[i]])/df[i]",
        "pool.sd[i] <- pow(pool.var[i], 0.5)"
      )
      model <- model.insert(model, pos=which(names(model)=="study"), x=pool.sd)
    }
  }

  model <- model.insert(model, pos=which(names(model)=="arm"), x=c(like, glm))


  # Add deviance contributions
  if (likelihood=="binomial") {
    resdevs <- c(
      "rhat[i,k] <- psi[i,k] * n[i,k]",
      "resdev[i,k] <- 2 * (r[i,k] * (log(r[i,k]) - log(rhat[i,k])) + (n[i,k] - r[i,k]) * (log(n[i,k] - r[i,k]) - log(n[i,k] - rhat[i,k])))"
    )
  } else if (likelihood=="normal") {
    resdevs <- c(
      "resdev[i,k] <- pow((y[i,k] - psi[i,k]),2) * prec[i,k] # residual deviance for normal likelihood"
    )
  } else if (likelihood=="poisson") {
    resdevs <- c(
      "resdev[i,k] <- 2*((lambda[i,k]-r[i,k]) + r[i,k]*log(r[i,k]/lambda[i,k]))"
    )
  }

  model <- model.insert(model, pos=which(names(model)=="arm"), x=resdevs)

  return(model)
}






#' Adds sections of JAGS code for an MBNMA model that correspond to delta
#' parameters
#'
#' @param model A character vector of JAGS MBNMA model code
#'
#' @inheritParams mbnma.run
#' @noRd
#'
#' @return A character object of JAGS MBNMA model code that includes delta
#'   parameter components of the model
#'
write.delta <- function(model, method="common", om) {

  priors <- default.priors(om=om)

  # Add mu
  model <- model.insert(model, pos=which(names(model)=="study"), x=priors[["mu"]])

  # Add to model
  # Add deltas and everything
  if (method %in% c("common", "random")) {
    model <- model.insert(model, pos=which(names(model)=="study"), x="delta[i,1] <- 0")
    if (method=="common") {
      model <- model.insert(model, pos=which(names(model)=="te"), x="delta[i,k] <- DR[i,k]")
    } else if (method=="random") {
      model <- model.insert(model, pos=which(names(model)=="te"), x=c(
        "delta[i,k] ~ dnorm(md[i,k], taud[i,k])",
        "md[i,k] <- DR[i,k] + sw[i,k]",
        "taud[i,k] <- tau * 2*(k-1)/k",
        "w[i,k] <- delta[i,k] - DR[i,k]",
        "sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)"
      ))
      model <- model.insert(model, pos=which(names(model)=="study"), x="w[i,1] <- 0")
      model <- model.insert(model, pos=which(names(model)=="end"), x=c(
        #"sd ~ dnorm(0,0.0025) T(0,)",
        priors[["sd"]],
        "tau <- pow(sd, -2)"
      ))
    } else {
      stop("`method` must take either `common` or `random`")
    }
  }

  return(model)
}





#' Adds sections of JAGS code for an MBNMA model that correspond to beta
#' parameters
#'
#' @inheritParams mbnma.run
#' @inheritParams get.prior
#' @noRd
#'
#' @return A character vector of JAGS MBNMA model code that includes beta
#'   parameter components of the model
#'
write.beta <- function(model, fun=dloglin(), method="common", om,
                       class.effect=list(), UME=FALSE
) {

  priors <- default.priors(fun=fun, UME=UME, om=om)

  if ("nonparam" %in% fun$name) {
    model <- model.insert(model, pos=which(names(model)=="start"), x="d.1[1,1] <- 0")

    if ("increasing" %in% fun$direction) {
      insert <- c("d.1[1,k] <- 0",
                  "for (c in 2:maxdose[k]) {",
                  priors[["nonparam"]],
                  "}")
      model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)

    } else if ("decreasing" %in% fun$direction) {
      insert <- c("d.1[1,k] <- 0",
                  "for (c in 2:maxdose[k]) {",
                  priors[["nonparam"]],
                  "}")
      model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)
    } else {
      stop("direction must be specified as either 'increasing' or 'decreasing' in dnonparam()")
    }

  } else {
    for (i in seq_along(fun$apool)) {
      pname <- names(fun$apool[i])
      isnum <- suppressWarnings(!is.na(as.numeric(fun$apool[i])))

      if (fun$apool[i] %in% "rel") {
        # Convert s.beta to d.
        if (UME==FALSE) {
          if (pname %in% c("ed50", "onset") | (pname %in% "rate" & "ditp" %in% fun$name)) {
            insert <- paste0("s.beta.", i, "[1] <- 0.00001") # To avoid numerical error with non-negative params
          } else {
            insert <- paste0("s.beta.", i, "[1] <- 0")
          }
          model <- model.insert(model, pos=which(names(model)=="start"), x=insert)

          insert <- paste0("s.beta.", i, "[k] <- ", pname, "[k]")
          model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)
        }

        # If class effects are modelled for this parameter
        if (pname %in% names(class.effect)) {
          insert <- priors[[toupper(pname)]]
          model <- model.insert(model, pos=which(names(model)=="class.prior"), x=insert)

          if (class.effect[[pname]]=="common") {
            insert <- paste0(pname, "[k]  <- ", toupper(pname), "[class[k]]")
            model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)

          } else if (class.effect[[pname]]=="random") {
            insert <- paste0(pname, "[k]  ~ dnorm(", toupper(pname), "[class[k]], tau.", toupper(pname), ")")
            model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)

            insert <- c(priors[[paste0("sd.", toupper(pname))]],
                        paste0("tau.", toupper(pname), " <- pow(sd.", toupper(pname), ", -2)"))
            model <- model.insert(model, pos=which(names(model)=="end"), x=insert)
          }
        } else {
          if (UME==FALSE) {
            insert <- priors[[pname]]
            model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)

          } else if (UME==TRUE) {
            insert <- c(paste0("s.beta.", i, "[k,c] <- ", pname, "[k,c]"),
                        priors[[pname]])
            model <- model.insert(model, pos=which(names(model)=="ume.prior"), x=insert)

            insert <- c(paste0("s.beta.", i, "[k,k] <- 0"),
                        paste0(pname, "[k,k] <- 0"))
            model <- model.insert(model, pos=which(names(model)=="ume.prior.ref"), x=insert)
          }
        }
      } else if (fun$apool[i] %in% c("common", "random")) {
        model <- model.insert(model, pos=which(names(model)=="end"), x=priors[[pname]])

        if (UME==FALSE) {
          insert <- paste0("s.beta.", i, "[1] <- 0")
          model <- model.insert(model, pos=which(names(model)=="start"), x=insert)

          if (fun$apool[i] %in% "random") {

            insert <- paste0("s.beta.", i, "[k] ~ dnorm(", pname, ", tau.", pname, ")")
            model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)

            insert <- c(priors[[paste0("sd.", pname)]],
                        paste0("tau.", pname, " <- pow(sd.", pname, ", -2)"))
            model <- model.insert(model, pos=which(names(model)=="end"), x=insert)

          } else if (fun$apool[i] %in% "common") {
            insert <- paste0("s.beta.", i, "[k] <- ", pname)
            model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)
          }

        } else if (UME==TRUE) {
          insert <- c(paste0("s.beta.", i, "[k,k] <- 0"))
          model <- model.insert(model, pos=which(names(model)=="ume.prior.ref"), x=insert)

          if (fun$apool[i] %in% "random") {

            insert <- paste0("s.beta.", i, "[k,c] ~ dnorm(", pname, ", tau.", pname, ")")
            model <- model.insert(model, pos=which(names(model)=="ume.prior"), x=insert)

            insert <- c(priors[[paste0("sd.", toupper(pname))]],
                        paste0("tau.", pname, " <- pow(sd.", pname, ", -2)"))
            model <- model.insert(model, pos=which(names(model)=="end"), x=insert)

          } else if (fun$apool[i] %in% "common") {
            insert <- paste0("s.beta.", i, "[k,c] <- ", pname)
            model <- model.insert(model, pos=which(names(model)=="ume.prior"), x=insert)
          }
        }


      } else if (isnum) {
        insert <- paste0(pname, " <- ", fun$apool[i])
        model <- model.insert(model, pos=which(names(model)=="end"), x=insert)

        if (UME==FALSE) {
          insert <- paste0("s.beta.", i, "[1] <- 0")
          model <- model.insert(model, pos=which(names(model)=="start"), x=insert)

          insert <- paste0("s.beta.", i, "[k] <- ", pname)
          model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=insert)

        } else if (UME==TRUE) {
          insert <- paste0("s.beta.", i, "[k,k] <- 0")
          model <- model.insert(model, pos=which(names(model)=="ume.prior.ref"), x=insert)

          insert <- paste0("s.beta.", i, "[k,c] <- ", pname)
          model <- model.insert(model, pos=which(names(model)=="ume.prior"), x=insert)

        }

      } else {
        stop(paste0(pname, " must take either `rel`, `common`, `random` or a numeric value"))
      }
    }
  }

  return(model)
}





#' Adds sections of JAGS code for meta-regression parameters
#'
#' @inheritParams mbnma.run
#' @inheritParams get.prior
#' @noRd
#'
#' @return A character vector of JAGS MBNMA model code that includes regression
#'   parameter components of the model
#'
write.regress <- function(model, regress.mat, regress.effect, om) {

  if (!is.null(regress.mat)) {
    priors <- default.priors(regress.mat=regress.mat, regress.effect=regress.effect, om=om)

    vars <- colnames(regress.mat)
    glm <- model["glm"]

    if (any(c("common", "agent", "class") %in% regress.effect)) {
      level <- "agent"

    } else if (any(c("random", "independent") %in% regress.effect)) {
      level <- "treatment"
    }

    # Add zero value for b[1]
    bzero <- paste0("b[1,1:nreg] <- rep(0,nreg)")
    model <- model.insert(model, pos=which(names(model)=="start"), x=bzero)

    # Add regression component to GLM
    modstr <- paste0("((b[", level, "[i,k],]-b[", level, "[i,1],]) %*% regress.mat[i,])")
    glm <- paste(glm, "+", modstr)

    for (i in seq_along(vars)) {

      # Add zero value for b[1]
      # bzero <- paste0("b.", vars[i], "[1] <- 0")
      # model <- model.insert(model, pos=which(names(model)=="start"), x=bzero)
      #
      # # Add regression component to GLM
      # modstr <- paste0("((b.", vars[i], "[", level, "[i,k]]-b.", vars[i], "[", level, "[i,1]]) * " , vars[i], "[i])")
      # glm <- paste(glm, "+", modstr)

      # Add assumption for regress effect (no need to add for agent since that is default if level==agent)
      if ("common" %in% regress.effect) {
        beffect <- paste0("b[k,", i, "] <- B.", vars[i])
        model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=beffect)

        model <- model.insert(model, pos=which(names(model)=="end"), x=priors[[paste0("B.", vars[i])]])

      } else if ("random" %in% regress.effect) {
        beffect <- paste0("b[k,", i, "] ~ dnorm(B.", vars[i], ", prec.B.", vars[i], ")")
        model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=beffect)

        model <- model.insert(model, pos=which(names(model)=="end"), x=priors[[paste0("B.", vars[i])]])

      } else if ("independent" %in% regress.effect) {
        beffect <- paste0("b[k,", i, "] <- B.", vars[i], "[k]")
        model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=beffect)

        model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=priors[[paste0("B.", vars[i])]])

      } else if ("agent" %in% regress.effect) {
        beffect <- paste0("b[k,", i, "] <- B.", vars[i], "[k]")
        model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=beffect)

        model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=priors[[paste0("B.", vars[i])]])

      } else if ("class" %in% regress.effect) {
        beffect <- paste0("b[k,", i, "] <- B.", vars[i], "[class[k]]")
        model <- model.insert(model, pos=which(names(model)=="agent.prior"), x=beffect)

        model <- model.insert(model, pos=which(names(model)=="class.prior"), x=priors[[paste0("B.", vars[i])]])

      } else {
        stop(paste0("regress.mat is incorrectly specified"))
      }


      # Add prior for sd.B
      if ("random" %in% regress.effect) {
        sdbprior <- c(priors[[paste0("sd.B.", vars[i])]],
                      paste0("prec.B.", vars[i], " <- pow(sd.B.", vars[i], ", -2)")
        )
        model <- model.insert(model, pos=which(names(model)=="end"), x=sdbprior)
      }
    }

    # Replace GLM
    model["glm"] <- glm
  }

  return(model)
}






#' Adds correlation between dose-response relative effects
#'
#' This uses a Wishart prior as default for modelling the correlation
#'
#' @param corprior Can be either `"wishart"` or `"spherical"` (for spherical decomposition)
#'
#' @inheritParams mbnma.run
#' @inheritParams write.beta
#' @inheritParams mbnma.write
#' @noRd
write.cor <- function(model, fun=dloglin(), cor=TRUE, omega=NULL, corprior="wishart",
                      method="random", class.effect=list(), UME=FALSE) {

  if (length(class.effect)>0 & cor==TRUE) {
    warning("Class effects cannot be modelled with correlation between dose-response relative effects - correlation will be ignored")
  } else {

    # Prepare covariance matrix if cor=TRUE
    if (cor==TRUE) {
      corparams <- names(fun$apool)[fun$apool %in% "rel"]
      mat.size <- length(corparams)

    } else {
      mat.size <- 0
    }

    if (mat.size>=2) {
      # Write covariance matrix in JAGS code
      model <- write.cov.mat(model, sufparams=corparams, omega=omega, UME=UME, corprior=corprior)
    }
  }

  return(model)

}





#' Function for adding covariance matrix for correlation between relative effects
#'
#' @param sufparams A numeric vector of dose-response/time-course parameter suffixes. It
#'  should be the same length as the number of relative effects (i.e. the covariance
#'  matrix size).
#' @inheritParams mbnma.write
#' @inheritParams get.prior
#' @inheritParams write.cor
#' @noRd
write.cov.mat <- function(model, sufparams, corprior="wishart",
                          omega=NULL, UME=FALSE) {

  if (UME==FALSE) {
    priorloc <- "agent.prior"
    index <- "k"
  } else if (UME==TRUE) {
    priorloc <- "ume.prior"
    index <- "k\\,c"
  }

  mat.size <- length(sufparams)
  for (i in seq_along(sufparams)) {
    model <- gsub(paste0(sufparams[i], "\\[", index, "\\] ~ [a-z]+\\([0-9]+(\\.[0-9]+)?,[0-9]+(\\.?[0-9]+)?\\)"),
                  paste0(sufparams[i], "[", index, "] <- mult[", i, ",", index, "]"),
                  model
                  )
  }

  insert <- paste0("mult[1:", mat.size, ",k] ~ dmnorm(d.prior[], inv.R[1:", mat.size, ", 1:", mat.size, "])")
  model <- model.insert(model, pos=which(names(model)==priorloc), x=insert)

  if (UME==TRUE) {
    model <- gsub("(mult\\[1\\:[0-9],k)(\\] ~)", "\\1,c\\2", model)
  }

  if (corprior=="wishart") {
    # Wishart prior
    jagswish <- c(paste0("for (p in 1:", mat.size, ") {"),
                  "d.prior[p] <- 0",
                  "}",
                  "",
                  paste0("inv.R ~ dwish(omega[,], ", length(sufparams), ")")
    )
    inv.R.prior <- jagswish

  } else if (corprior=="spherical") {
    # Spherical decomposition prior
    sphere <- c(paste0("for (p in 1:", mat.size, ") {"),
                "d.prior[p] <- 0",
                "inv.R[p,p] <- 1000",
                "}",
                paste0("for (p in 1:", mat.size-1, ") {"),
                paste0("for (q in (p+1):", mat.size, ") {"),
                "inv.R[p,q] <- (rho[p,q] * 1000)",
                "inv.R[q,p] <- inv.R[p,q]",
                "rho[p,q] <- inprod(g[,p], g[,q])",
                "g[q,p] <- 0",
                "a[p,q] ~ dunif(0, 3.1415)",
                "}",
                "}"
    )

    g <- matrix(ncol=mat.size, nrow=mat.size)
    g[1,1] <- "1"
    for (i in 2:mat.size) {
      g[i,i] <- paste(paste0("sin(a[", 1:(i-1), ",", i, "])"), collapse="*")
      for (k in 1:(i-1)) {
        g[k,i] <- paste0("cos(a[", k, ",", i, "])")

        if (k>1) {
          g[k,i] <- paste(g[k,i], paste0("sin(a[", 1:(k-1), ",", i, "])", collapse="*"), sep="*")
        }
      }
    }
    for (i in 1:mat.size) {
      sphere <- append(sphere, paste0("g[", i, ",", i, "] <- ", g[i,i]))
    }
    for (i in 1:(mat.size-1)) {
      for (k in (i+1):mat.size) {
        sphere <- append(sphere, paste0("g[", i, ",", k, "] <- ", g[i,k]))
      }
    }

    inv.R.prior <- sphere
  }

  model <- model.insert(model, pos=which(names(model)=="end"), x=inv.R.prior)

  return(model)
}







#' Write JAGS code for mbnma.nodesplit
#'
#' @inheritParams nma.run
#' @noRd
add.nodesplit <- function(model) {

  priors <- default.priors()

  # If method=="common"
  if (any(grepl("delta\\[i\\,k\\] <-", model))) {

    match <- grep("^delta\\[i\\,k\\] <- DR", model)
    model[match] <- "delta[i,k] <- md[i,k]"
    model <- model.insert(model, pos=match, x="md[i,k] <- ifelse(split.ind[i,k]==1, direct, DR[i,k])")

    # Else if method=="random
  } else if (any(grepl("delta\\[i\\,k\\] ~", model))) {

    match <- grep("^md\\[i\\,k\\] <- DR", model)
    model[match] <- "md[i,k] <- ifelse(split.ind[i,k]==1, direct, DR[i,k] + sw[i,k])"
  }

  # Add prior for direct
  model <- model.insert(model, pos=which(names(model)=="end"), priors[["direct"]])

  return(model)
}




#' Removes empty loops from JAGS code
#'
#' @noRd
remove.loops <- function(model) {

  segs <- c("agent.prior", "trt.prior", "class.prior")

  for (i in seq_along(segs)) {
    pos <- which(names(model)==segs[i])
    if (model[pos+1] == "}") {
      model <- model[c(1:(pos-1), (pos+2):length(model))]
    }
  }
  return(model)
}






#' Write JAGS code for split NMA
#'
#' @inheritParams nma.run
#' @noRd
write.nma <- function(method="common", likelihood="binomial", link="logit", om,
                      UME=FALSE) {

  priors <- default.priors(UME=UME, om=om)

  model <- c(start="model{ 			# Begin Model Code",
             study="for(i in 1:NS){ # Run through all NS trials",
             priors[["mu"]],
             "delta[i,1] <- 0",
             arm="for (k in 1:narm[i]){ # Run through all arms within a study",
             "}",
             "",
             "resstudydev[i] <- sum(resdev[i, 1:narm[i]])",
             te="for(k in 2:narm[i]){ # Treatment effects",
             "}",
             "}",
             "",
             "totresdev <- sum(resstudydev[])",
             end="",
             "# Model ends",
             "}"
  )

  # Add likelihood
  model <- write.likelihood(model, likelihood = likelihood, link=link)


  # Add d[1] <- 0
  if (UME==FALSE) {
    model <- model.insert(model, pos=which(names(model)=="start"), x="d[1] <- 0")
  }


  # Add treatment effects
  if (method=="common") {
    teinsert <- c("delta[i,k] <- md[i,k]",
                  "md[i,k] <- d[treatment[i,k]] - d[treatment[i,1]]"
    )

  } else if (method=="random") {

    model <- model.insert(model, pos=which(names(model)=="study"), x="w[i,1] <- 0")

    insert <- c(priors[["sd"]],
                "tau <- pow(sd,-2)")
    model <- model.insert(model, pos=which(names(model)=="end"), x=insert)


    teinsert <- c("delta[i,k] ~ dnorm(md[i,k], taud[i,k])",
                  "md[i,k] <- d[treatment[i,k]] - d[treatment[i,1]] + sw[i,k]",
                  "taud[i,k] <- tau *2*(k-1)/k",
                  "w[i,k] <- (delta[i,k] - d[treatment[i,k]] + d[treatment[i,1]])",
                  "sw[i,k] <- sum(w[i,1:(k-1)])/(k-1)")
  }

  model <- model.insert(model, pos=which(names(model)=="te"), x=teinsert)


  # Add treatment effect priors and make UME changes
  if (UME==FALSE) {
    te.prior <- c("for (k in 2:NT){ # Priors on relative treatment effects",
                  priors[["d"]],
                  "}")

  } else if (UME==TRUE) {

    model <- gsub("d\\[treatment\\[i,k\\]\\] (\\+|-) d\\[treatment\\[i,1\\]\\]",
                  "d[treatment[i,k],treatment[i,1]]",
                  model
    )

    te.prior <- c("for (k in 1:NT) { d[k,k] <- 0 }",
                  "for (c in 1:(NT-1)) {",
                  "for (k in (c+1):NT) {",
                  priors[["d"]],
                  "}",
                  "}")
  }

  model <- model.insert(model, pos=which(names(model)=="end"), x=te.prior)

return(model)
}







#' Write E0 synthesis JAGS model
#'
#' @inheritParams predict.mbnma
#' @noRd
write.E0.synth <- function(synth="fixed", likelihood=NULL, link=NULL, om) {
  model <- c(
    start="model{ 			# Begin Model Code",
    study="for(i in 1:NS){ # Run through all NS trials",
    arm="for (k in 1:narm[i]){ # Run through all arms within a study",
    "}",
    "",
    "resstudydev[i] <- sum(resdev[i, 1:narm[i]])",
    "}",
    "totresdev <- sum(resstudydev[])",
    end="m.mu ~ dnorm(0,0.0001)",
    "# Model ends",
    "}"
  )


  # Add likelihood
  model <- write.likelihood(model, likelihood = likelihood, link=link)

  if (synth=="fixed") {
    mucode <- "mu[i] <- m.mu"
    model <- model.insert(model, pos=which(names(model)=="study"), x=mucode)

  } else if (synth=="random") {
    mucode <- "mu[i] ~ dnorm(m.mu, tau.mu)"
    musdcode <- c("tau.mu <- pow(sd.mu, -2)",
                  ifelse(om$abs>0,
                         paste0("sd.mu ~ dunif(0,", om$abs, ")"),
                         paste0("sd.mu ~ dunif(", om$abs, ", 0)")))

    model <- model.insert(model, pos=which(names(model)=="study"), x=mucode)
    model <- model.insert(model, pos=which(names(model)=="end"), x=musdcode)
  }

  # Remove delta
  model["glm"] <- gsub("\\+ delta\\[i,k\\]", "", model["glm"])

  return(model)
}






#' Sets default priors for JAGS model code
#'
#' This function creates JAGS code snippets for default MBNMA model priors.
#'
#' @inheritParams mb.run
#'
#' @return A list, each element of which is a named JAGS snippet
#'   corresponding to a prior in the MBNMA JAGS code.
#'
#' @examples
#' \donttest{
#' default.priors(fun=demax())
#' }
#'
#' @export
default.priors <- function(fun=dloglin(), UME=FALSE,
                           regress.mat=NULL, regress.effect="common",
                           om=list("rel"=5, "abs"=10)) {

  sufparams <- which(fun$apool=="rel")

  priors <- list(
    rho = "rho ~ dunif(0,1)",
    mu = "mu[i] ~ dnorm(0,0.0001)",
    direct = "direct ~ dnorm(0,0.0001)" # For nodesplit
  )

  priors[["sd"]] <- ifelse(om$rel>0,
                           paste0("sd ~ dunif(0, ", om$rel, ")"),
                           paste0("sd ~ dunif(", om$rel, ", 0)")
                           )

  # Non-parametric
  if ("nonparam" %in% fun$name) {
    temp <- "d.1[c,k] ~ dnorm(d.1[c-1,k],0.0001)"
    if (fun$direction=="increasing") {
      temp <- paste(temp, "T(d.1[c-1,k],)", sep=" ")
    } else {
      temp <- paste(temp, "T(,d.1[c-1,k])", sep=" ")
    }
    priors[["nonparam"]] <- temp
  }

  for (i in seq_along(fun$params)) {
    pname <- fun$params[i]

    if (any(c("common", "random") %in% fun$apool[i])) {
      priors[[pname]] <- paste0(pname, " ~ dnorm(0,0.0001)")
    } else {
      priors[[pname]] <- paste0(pname, "[k] ~ dnorm(0,0.0001)")
    }

    priors[[toupper(pname)]] <- paste0(toupper(pname), "[k] ~ dnorm(0,0.0001)")


    priors[[paste0("sd.", pname)]] <- ifelse(om$rel>0,
                                             paste0("sd.", pname, " ~ dunif(0, ", om$rel, ")"),
                                             paste0("sd.", pname, " ~ dunif(", om$rel, ", 0)"))
    priors[[paste0("sd.", toupper(pname))]] <- ifelse(om$rel>0,
                                                      paste0("sd.", toupper(pname), " ~ dunif(0, ", om$rel, ")"),
                                                      paste0("sd.", toupper(pname), " ~ dunif(", om$rel, ", 0)"))

    # For priors with truncated normal
    if (pname %in% c("ed50", "onset") | (pname %in% "rate" & "ditp" %in% fun$name)) {
      priors[[pname]] <- paste(priors[[pname]], "T(0,)")
      priors[[toupper(pname)]] <- paste(priors[[toupper(pname)]], "T(0,)")
    }

    # For UME models
    if (UME==TRUE) {
      priors[[pname]] <- gsub("\\[k\\]", "[k,c]", priors[[pname]])
      priors[[toupper(pname)]] <- gsub("\\[k\\]", "[k,c]", priors[[toupper(pname)]])
    }
  }

  vars <- colnames(regress.mat)
  for (i in seq_along(vars)) {
    if (any(c("common", "random") %in% regress.effect)) {
      priors[[paste0("B.", vars[i])]] <- paste0("B.", vars[i], " ~ dnorm(0,0.0001)")
    } else {
      priors[[paste0("B.", vars[i])]] <- paste0("B.", vars[i], "[k] ~ dnorm(0,0.0001)")
    }

    priors[[paste0("sd.B.", vars[i])]] <- paste0("sd.B.", vars[i], " ~ dunif(0,", om$rel, ")")
  }

  # For NMA models
  priors[["d"]] <- "d[k] ~ dnorm(0,0.0001)"
  if (UME==TRUE) {
    priors[["d"]] <- "d[k,c] ~ dnorm(0,0.0001)"
  }

  return(priors)
}
