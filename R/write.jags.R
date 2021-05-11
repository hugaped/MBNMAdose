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
                        cor=TRUE, cor.prior="wishart",
                        omega=NULL,
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
              omega=omega,
              class.effect=class.effect, UME=UME)

  model <- write.model(class.effect=class.effect, UME=UME)

  # Add dose-response function
  dosefun <- write.dose.fun(fun=fun, UME=UME)
  model <- model.insert(model, pos=which(names(model)=="arm"), x=dosefun[[1]])

  # if (length(dosefun)==2) {  # For models with multiple DR functions
  #   model <- model.insert(model, pos=which(names(model)=="study"), x=dosefun[[2]])
  # }


  # Add likelihood
  model <- write.likelihood(model, likelihood=likelihood, link=link)

  # Add treatment delta effects
  model <- write.delta(model, method=method)

  # Add treatment beta effects
  model <- write.beta(model, fun=fun, method=method, class.effect=class.effect, UME=UME)

  # Add correlation between dose-response parameters
  model <- write.cor(model, fun=fun, method=method, cor=cor, omega=omega,
                     class.effect=class.effect, UME=UME)

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
                        UME=FALSE,
                        cor.prior="wishart",
                        omega=NULL,
                        user.fun=NULL,
                        class.effect=list()) {


  # Run argument checks
  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertClass(fun, "dosefun", null.ok=FALSE, add=argcheck)
  checkmate::assertChoice(method, choices=c("common", "random"), null.ok=FALSE, add=argcheck)
  checkmate::assertLogical(UME, null.ok=FALSE, add=argcheck)
  if (method=="random") {
    checkmate::assertChoice(cor.prior, choices=c("wishart", "rho"))
  }
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

    if (any(c("rcs", "ns", "bs") %in% fun$name)) {
      warning("Class effects applied to spline function parameters may produce\n#
              non-interpretable results since knot locations will differ between agents")
    }


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
      warning("Class effects cannot be modelled with correlation between time-course relative effects - 'omega' will be ignored")
    }

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






#' Write the basic JAGS model code for MBNMA to which other lines of model
#' code can be added
#'
#' @return A character object of JAGS model code
#' @noRd
#' @examples
#' model <- write.model()
#' cat(model)
write.model <- function(UME=FALSE, class.effect=list()) {

  model <- c(
    start= "model{ 			# Begin Model Code",
    study="for(i in 1:NS){ # Run through all NS trials",
    "mu[i] ~ dnorm(0,0.001)",
    arm="for (k in 1:narm[i]){ # Run through all arms within a study",
    "}",
    "resstudydev[i] <- sum(resdev[i, 1:narm[i]])",
    te="for(k in 2:narm[i]){ # Treatment effects",
    "}",
    "}",
    trt.prior="for (k in 2:Nagent){ # Priors on relative treatment effects",
    "}"
    )

  if (length(class.effect)>0) {
    model <- append(model,  c(class.prior="for (k in 2:Nclass){ # Priors on relative class effects",
                              "}"))
  }
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
write.likelihood <- function(model, likelihood="binomial", link=NULL) {

  argcheck <- checkmate::makeAssertCollection()
  checkmate::assertChoice(likelihood, choices=c("binomial", "normal", "poisson"), add=argcheck)
  checkmate::assertChoice(link, choices=c("logit", "identity", "cloglog", "probit", "log", "smd"), add=argcheck)
  checkmate::reportAssertions(argcheck)

  if (likelihood=="binomial") {
    like <- "r[i,k] ~ dbin(psi[i,k], N[i,k])"
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
           "theta[i,k] <- mu[i] + delta[i,k]")
  if (!link %in% c("identity", "smd")) {
    glm <- gsub("(psi\\[i,k\\])(.+)", paste0(link, "(\\1)\\2"), glm)
  }

  if (link=="smd") {
    if (likelihood!="normal") {
      stop("link='smd' can only be used with likelihood='normal'")
    }
    glm[1] <- paste0(glm[1], " * pool.sd[i]")

    # Add SMD components
    smd.sub <- c(
      "sd[i,k] <- se[i,k] * pow(N[i,k],0.5)",
      "nvar[i,k] <- (N[i,k]-1) * pow(sd[i,k],2)"
    )
    model <- model.insert(model, pos=which(names(model)=="arm"), x=smd.sub)

    pool.sd <- c(
      "df[i] <- sum(N[i,1:narm[i]]) - narm[i]",
      "pool.var[i] <- sum(nvar[i,1:narm[i]])/df[i]",
      "pool.sd[i] <- pow(pool.var[i], 0.5)"
    )
    model <- model.insert(model, pos=which(names(model)=="study"), x=pool.sd)
  }

  model <- model.insert(model, pos=which(names(model)=="arm"), x=c(like, glm))


  # Add deviance contributions
  if (likelihood=="binomial") {
    resdevs <- c(
      "rhat[i,k] <- psi[i,k] * N[i,k]",
      "resdev[i,k] <- 2 * (r[i,k] * (log(r[i,k]) - log(rhat[i,k])) + (N[i,k] - r[i,k]) * (log(N[i,k] - r[i,k]) - log(N[i,k] - rhat[i,k])))"
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
write.delta <- function(model, method="common") {

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
        "sd ~ dnorm(0,0.0025) T(0,)",
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
write.beta <- function(model, fun=dloglin(), method="common",
                       class.effect=list(), UME=FALSE
) {

  if ("nonparam" %in% fun$name) {
    model <- model.insert(model, pos=which(names(model)=="start"), x="d.1[1,1] <- 0")

    if ("increasing" %in% fun$direction) {
      insert <- c("d.1[1,k] <- 0",
                  "for (c in 2:maxdose[k]) {",
                  "d.1[c,k] ~ dnorm(d.1[c-1,k],0.0001) T(d.1[c-1,k],)",
                  "}")
      model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)

    } else if ("decreasing" %in% fun$direction) {
      insert <- c("d.1[1,k] <- 0",
                  "for (c in 2:maxdose[k]) {",
                  "d.1[c,k] ~ dnorm(d.1[c-1,k],0.0001) T(,d.1[c-1,k])",
                  "}")
      model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)
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
          insert <- paste0("s.beta.", i, "[1] <- 0")
          model <- model.insert(model, pos=which(names(model)=="start"), x=insert)

          insert <- paste0("s.beta.", i, "[k] <- ", pname, "[k]")
          model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)
        }

        # If class effects are modelled for this parameter
        if (pname %in% names(class.effect)) {
          insert <- paste0(toupper(pname), "[k] ~ dnorm(0,0.001)")
          model <- model.insert(model, pos=which(names(model)=="class.prior"), x=insert)

          if (class.effect[[pname]]=="common") {
            insert <- paste0(pname, "[k]  <- ", toupper(pname), "[class[k]]")
            model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)

          } else if (class.effect[[pname]]=="random") {
            insert <- paste0(pname, "[k]  ~ dnorm(", toupper(pname), "[class[k]], tau.", toupper(pname), ")")
            model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)

            insert <- c(paste0("sd.", toupper(pname), " ~ dnorm(0,0.0025) T(0,)"),
                        paste0("tau.", toupper(pname), " <- pow(sd.", toupper(pname), ", -2)"))
            model <- model.insert(model, pos=which(names(model)=="end"), x=insert)
          }
        } else {
          if (UME==FALSE) {
            insert <- paste0(pname, "[k] ~ dnorm(0,0.001)")
            model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)

          } else if (UME==TRUE) {
            insert <- c(paste0("s.beta.", i, "[k,c] <- ", pname, "[k,c]"),
                        paste0(pname, "[k,c] ~ dnorm(0,0.001)"))
            model <- model.insert(model, pos=which(names(model)=="ume.prior"), x=insert)

            insert <- c(paste0("s.beta.", i, "[k,k] <- 0"),
                        paste0(pname, "[k,k] <- 0"))
            model <- model.insert(model, pos=which(names(model)=="ume.prior.ref"), x=insert)
          }
        }
      } else if (fun$apool[i] %in% c("common", "random")) {
        insert <- paste0(pname, " ~ dnorm(0,0.001)")
        model <- model.insert(model, pos=which(names(model)=="end"), x=insert)

        if (UME==FALSE) {
          insert <- paste0("s.beta.", i, "[1] <- 0")
          model <- model.insert(model, pos=which(names(model)=="start"), x=insert)

          if (fun$apool[i] %in% "random") {

            insert <- paste0("s.beta.", i, "[k] ~ dnorm(", pname, ", tau.", pname, ")")
            model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)

            insert <- c(paste0("sd.", pname, " ~ dnorm(0,0.0025) T(0,)"),
                        paste0("tau.", pname, " <- pow(sd.", pname, ", -2)"))
            model <- model.insert(model, pos=which(names(model)=="end"), x=insert)

          } else if (fun$apool[i] %in% "common") {
            insert <- paste0("s.beta.", i, "[k] <- ", pname)
            model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)
          }

        } else if (UME==TRUE) {
          insert <- c(paste0("s.beta.", i, "[k,k] <- 0"))
          model <- model.insert(model, pos=which(names(model)=="ume.prior.ref"), x=insert)

          if (fun$apool[i] %in% "random") {

            insert <- paste0("s.beta.", i, "[k,c] ~ dnorm(", pname, ", tau.", pname, ")")
            model <- model.insert(model, pos=which(names(model)=="ume.prior"), x=insert)

            insert <- c(paste0("sd.", pname, " ~ dnorm(0,0.0025) T(0,)"),
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
          model <- model.insert(model, pos=which(names(model)=="trt.prior"), x=insert)

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






#' Adds correlation between dose-response relative effects
#'
#' This uses a Wishart prior as default for modelling the correlation
#'
#' @inheritParams mbnma.run
#' @inheritParams write.beta
#' @inheritParams mbnma.write
#' @noRd
write.cor <- function(model, fun=dloglin(), cor=TRUE, omega=NULL,
                      method="random", class.effect=list(), UME=FALSE) {

  if (length(class.effect)>0 & cor==TRUE) {
    warning("Class effects cannot be modelled with correlation between time-course relative effects - correlation will be ignored")
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
      model <- write.cov.mat(model, sufparams=corparams, omega=omega, UME=UME)
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
#' @noRd
write.cov.mat <- function(model, sufparams,
                          omega=NULL, UME=FALSE) {

  if (UME==FALSE) {
    priorloc <- "trt.prior"
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

  jagswish <- c(paste0("for (r in 1:", mat.size, ") {"),
                "d.prior[r] <- 0",
                "}",
                "",
                paste0("inv.R ~ dwish(omega[,], ", length(sufparams), ")")
  )
  model <- model.insert(model, pos=which(names(model)=="end"), x=jagswish)

  return(model)
}







#' Write JAGS code for mbnma.nodesplit
#'
#' @inheritParams nma.run
#' @noRd
add.nodesplit <- function(model) {

  # If method=="common"
  if (any(grepl("delta\\[i\\,k\\] <-", model))) {

    match <- grep("^delta\\[i\\,k\\] <- DR", model)
    model[match] <- "delta[i,k] <- md[i,k]"
    model <- model.insert(model, pos=match, x="md[i,k] <- ifelse(split.ind[i,k]==1, direct, DR[i,k])")

    # Else if method=="random
  } else if (grepl("delta\\[i\\,k\\] ~", model)) {

    match <- grep("^md\\[i\\,k\\] <- DR", model)
    model[match] <- "md[i,k] <- ifelse(split.ind[i,k]==1, direct, DR[i,k] + sw[i,k])"
  }

  # Add prior for direct
  model <- model.insert(model, pos=which(names(model)=="end"), "direct ~ dnorm(0,0.0001)")

  return(model)
}




#' Removes empty loops from JAGS code
#'
#' @noRd
remove.loops <- function(model) {

  segs <- c("trt.prior")

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
write.nma <- function(method="common", likelihood="binomial", link="logit",
                      UME=FALSE) {
  model <- c(start="model{ 			# Begin Model Code",
             study="for(i in 1:NS){ # Run through all NS trials",
             "mu[i] ~ dnorm(0,0.001)",
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

    insert <- c("tau <- pow(sd,-2)",
                "sd ~ dnorm(0,0.0025) T(0,)")
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
                  "d[k] ~ dnorm(0,0.0001)",
                  "}")

  } else if (UME==TRUE) {

    model <- gsub("d\\[treatment\\[i,k\\]\\] (\\+|-) d\\[treatment\\[i,1\\]\\]",
                  "d[treatment[i,k],treatment[i,1]]",
                  model
    )

    te.prior <- c("for (k in 1:NT) { d[k,k] <- 0 }",
                  "for (c in 1:(NT-1)) {",
                  "for (k in (c+1):NT) {",
                  "d[k,c] ~ dnorm(0,0.0001)",
                  "}",
                  "}")
  }

  model <- model.insert(model, pos=which(names(model)=="end"), x=te.prior)

return(model)
}
