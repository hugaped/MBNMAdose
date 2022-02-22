testthat::context("Testing full set of functions")

# Includes tests for mbnma.run()

# Across a range of datasets and dose-response functions:
# Tests running
# Tests default plots of running (including fitplot and devplot)
# Tests default ranking (including cumrank)
# Tests default prediction
# Occasionally tests get.relative


# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)

# Datasets with no placebo
network <- mbnma.network(psoriasis90)
psoriasis90.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

network <- mbnma.network(ssri)
ssri.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

alldfs <- list(triptans, psoriasis90.noplac, osteopain, psoriasis75, ssri, ssri.noplac, gout)
datanams <- c("triptans", "psoriasis90.noplac", "osteopain", "psoriasis75", "ssri", "ssri.noplac", "gout")



for (dat in seq_along(alldfs)) {

  datanam <- datanams[dat]
  dataset <- alldfs[[dat]]

  print(datanam)

  ### Datasets ####
  network <- mbnma.network(dataset)


  # Make class data
  df <- dataset

  if ("class" %in% names(dataset)) {
    netclass <- mbnma.network(df)
  }


  test_that(paste("mbnma.run functions correctly for:", datanam), {
    skip_on_appveyor()
    skip_on_ci()
    skip_on_cran()

    n.iter=500
    pd <- "pv"

    # Single parameter DR functions
    result <- mbnma.run(network, fun=dpoly(degree=1), method="common",
                        pd="plugin", n.iter=n.iter)
    expect_equal(class(result), c("mbnma", "rjags"))
    expect_equal("beta.1" %in% result$parameters.to.save, TRUE)
    expect_equal(result$model.arg$pd, "plugin")
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)


    result <- mbnma.run(network, fun=dexp(), method="random",
                        pd="pd.kl", n.iter=n.iter)
    expect_equal(class(result), c("mbnma", "rjags"))
    expect_equal("sd" %in% result$parameters.to.save, TRUE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(devplot(result), NA)
    expect_error(fitplot(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)


    if ("class" %in% names(dataset)) {
      result <- expect_warning(mbnma.run(netclass, fun=dexp(), method="common",
                          pd="popt", class.effect = list(emax="random"), n.iter=n.iter))
      expect_equal(class(result), c("mbnma", "rjags"))
      expect_equal("EMAX" %in% result$parameters.to.save, TRUE)
      expect_equal("sd.EMAX" %in% result$parameters.to.save, TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(devplot(result), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result))
      expect_error(suppressWarnings(summary(result)), NA)
      expect_error(get.relative(result), NA)

    }



    # Two parameter DR functions
    result <- mbnma.run(network, fun=demax(), method="common",
                        n.iter=n.iter, pd=pd)
    expect_equal(all(c("emax", "ed50") %in% result$parameters.to.save), TRUE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(devplot(result), NA)
    expect_error(fitplot(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)


    result <- mbnma.run(network, fun=demax(ed50="common"), method="random",
                        n.iter=n.iter, pd=pd)
    expect_equal("sd" %in% result$parameters.to.save, TRUE)
    expect_equal("ed50" %in% rownames(result$BUGSoutput$summary), TRUE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)


    if ("class" %in% names(dataset)) {
      expect_error(mbnma.run(netclass, fun=demax(), method="random",
                             class.effect=list(fakeparam="common")), "corresponding to dose-response parameters")

      expect_warning(mbnma.run(netclass, fun=demax(), method="random",
                               class.effect=list(ed50="common"), n.iter=n.iter, pd=pd))
      result <- suppressWarnings(mbnma.run(netclass, fun=demax(), method="random",
                          class.effect=list(ed50="common"), n.iter=n.iter, pd=pd))
      expect_equal(all(c("emax", "ED50", "ed50", "sd") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result))
      expect_error(suppressWarnings(summary(result)), NA)

    }


    result <- mbnma.run(network, fun=demax(ed50="random"), method="common",
                        n.iter=n.iter, pd=pd)
    expect_equal(all(c("emax", "ed50", "sd.ed50") %in% result$parameters.to.save), TRUE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)


    # Three parameter DR function
    result <- mbnma.run(network, fun=demax(emax="rel", ed50="random", hill="common"),
                        method="random", n.iter=n.iter, pd=pd, priors = list(hill="dunif(-3,3)"))
    expect_equal(all(c("emax", "sd", "ed50", "sd.ed50", "hill") %in% result$parameters.to.save), TRUE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)

    result <- mbnma.run(network, fun=dspline(type="ns", knots=3, beta.1="rel", beta.2="random", beta.3="common"),
                        method="random", n.iter=n.iter, pd=pd)
    expect_equal(all(c("beta.1", "beta.2", "sd.beta.2", "sd", "beta.3") %in% result$parameters.to.save), TRUE)
    expect_equal(any(grepl("spline", result$model.arg$jagscode)), TRUE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(devplot(result), NA)
    expect_error(fitplot(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)

    result <- mbnma.run(network, fun=demax(ed50="random", hill=0.2),
                        parameters.to.save = "psi",
                        method="random", n.iter=n.iter, pd=pd)
    expect_equal("psi" %in% result$parameters.to.save, TRUE)
    expect_equal("d.1" %in% result$parameters.to.save, FALSE)
    expect_error(plot(result))
    expect_error(rank(result))
    expect_error(devplot(result), NA)
    expect_error(fitplot(result), NA)
    expect_error(predict(result))
    expect_error(summary(result))


    # Splines and polynomials
    result <- mbnma.run(network, fun=dspline(type="bs", knots=2, beta.1="common", beta.2 = "rel", beta.3="random"),
                        n.iter=n.iter, pd=pd)
    expect_equal(all(c("beta.1", "beta.2", "sd.beta.3", "beta.3") %in% result$parameters.to.save), TRUE)
    expect_equal(all(c("sd") %in% result$parameters.to.save), FALSE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)

    result <- mbnma.run(network, fun=dspline(type="ns", knots=c(0.2,0.5), beta.1="common", beta.2 = "rel", beta.3="random"),
                        n.iter=n.iter, pd=pd)
    expect_equal(all(c("beta.1", "beta.2") %in% result$parameters.to.save), TRUE)
    expect_equal(all(c("sd", "beta.3") %in% result$parameters.to.save), FALSE)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(devplot(result), NA)
    expect_error(fitplot(result), NA)
    expect_error(predict(result), NA)
    expect_error(suppressWarnings(summary(result)), NA)

    if (!datanam %in% "osteopain") {
      result <- mbnma.run(network, fun=dpoly(degree=3, beta.1="common", beta.2 = "rel", beta.3="random"),
                          method="random", n.iter=n.iter, pd=pd)
      expect_equal(all(c("beta.1", "beta.2", "sd.beta.3", "beta.3", "sd") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
    }

    # Test Omega, cor
    expect_equal("inv.R" %in% names(result$model.arg$priors), FALSE)
    expect_equal(result$model.arg$omega, NULL)

    result <- mbnma.run(network, fun=dpoly(degree=3, beta.1="rel", beta.2 = "rel", beta.3="random"),
                        method="random", n.iter=n.iter, pd=pd)
    expect_equal("inv.R" %in% names(result$model.arg$priors), TRUE)
    expect_equal(ncol(result$model.arg$omega), 2)
    expect_error(get.relative(result), NA)


    omega <- matrix(c(10,2,2,10), byrow = TRUE, ncol=2)
    result <- mbnma.run(network, fun=demax(), omega=omega, n.iter=n.iter, pd=pd)
    expect_equal(ncol(result$model.arg$omega), 2)
    expect_equal(result$model.arg$omega[1,1], omega[1,1])

    omega <- matrix(c(10,4,2,10), byrow = TRUE, ncol=2)
    expect_error(mbnma.run(network, fun=demax(), omega=omega, n.iter=n.iter, pd=pd), "omega must be")

    result <- mbnma.run(network, fun=demax(), cor=FALSE, n.iter=n.iter, pd=pd)
    expect_equal("inv.R" %in% names(result$model.arg$priors), FALSE)


    # Test UME
    result <- mbnma.run(network, fun=dpoly(degree=3, beta.1="random", beta.2 = "rel", beta.3="rel"),
                        method="random", n.iter=n.iter, pd=pd, UME=TRUE)
    expect_equal(all(c("beta.1", "beta.2", "sd.beta.1", "beta.3", "sd") %in% result$parameters.to.save), TRUE)
    expect_equal(any(grepl("beta\\.2\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
    expect_equal(any(grepl("beta\\.3\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
    expect_equal(any(grepl("beta\\.1\\[2,1\\]", rownames(result$BUGSoutput$summary))), FALSE)

    result <- mbnma.run(network, fun=dspline(knots=c(0.1, 0.3), type="ls"),
                        n.iter=n.iter, pd=pd, UME=TRUE)
    expect_equal(any(grepl("beta\\.2\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
    expect_equal(any(grepl("beta\\.3\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
    expect_equal(any(grepl("beta\\.1\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
    expect_error(get.relative(result), "cannot be used on")

    expect_error(plot(result), "model results cannot be plotted")
    expect_error(rank(result), "does not work for UME")
    expect_error(devplot(result), NA)
    expect_error(fitplot(result), NA)
    expect_error(predict(result), "does not work with UME")
    expect_error(suppressWarnings(summary(result)), NA)


    # Link functions
    if (datanam %in% c("triptans", "psoriasis90", "ssri")) {
      expect_error(mbnma.run(network, fun=demax(), link="probit", n.iter=n.iter, pd=pd), NA)
      result <- mbnma.run(network, fun=dfpoly(degree=2), link="cloglog", n.iter=n.iter, pd=pd)
      expect_equal(result$model.arg$link, "cloglog")

      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(devplot(result), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
    } else if (datanam %in% c("osteopain", "gout")) {

      if (datanam %in% "gout") {
        expect_error(mbnma.run(network, fun=dspline(knots=2, type="bs"), link="smd", n.iter=n.iter, pd=pd), "must be included in data\\.ab")
      } else if (datanam %in% "osteopain") {
        expect_error(mbnma.run(network, fun=dspline(knots=2, type="bs"), link="smd", n.iter=n.iter, pd=pd), NA)
      }

      result <- suppressWarnings(mbnma.run(network, fun=dexp(), link="log", n.iter=n.iter, pd=pd))
      expect_equal(result$model.arg$link, "log")

      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(devplot(result), NA)
      if (!datanam %in% "gout") {
        expect_error(fitplot(result), NA)
      }
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
    }


    # Add tests for non-parametric
    if (grepl("noplac", datanam)) {
      expect_error(mbnma.run(network, fun=dnonparam()), "must currently include placebo")
    } else {
      result <- mbnma.run(network, fun=dnonparam(direction = "decreasing"), method="random", n.iter=n.iter, pd=pd)
      expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
      expect_equal("sd" %in% result$parameters.to.save, TRUE)
      expect_error(summary(result))

      expect_error(result <- mbnma.run(network, fun=dpoly(), method="fixed", n.iter=n.iter, pd=pd), "Must be element of set")
    }



    # Changing priors
    result <- mbnma.run(network, fun=demax(), method="random",
                        n.iter=n.iter, pd=pd)
    prior <- list(sd="dunif(0,5)", inv.R="dwish(omega[,],5)")
    runprior <- mbnma.run(network, fun=demax(), method="random",
                          n.iter=n.iter, pd=pd, priors = prior)
    expect_equal(runprior$model.arg$priors$sd, prior$sd)
    expect_equal(runprior$model.arg$priors$inv.R, prior$inv.R)
    expect_equal(result$model.arg$priors$inv.R!=runprior$model.arg$priors$inv.R, TRUE)





    # Multiple dose-response functions
    mult <- dmulti(
      c(rep(list(dexp()),2),
      rep(list(dpoly(degree=2)),1),
      rep(list(demax()),length(network$agents)-3)
      ))

    multifun <- mbnma.run(network, fun=mult,
                          n.iter=n.iter, pd=pd)
    expect_equal(length(multifun$model.arg$fun$name)>1, TRUE)
    expect_equal(all(c("beta.1", "beta.2", "emax", "ed50", "emax.4") %in% multifun$parameters.to.save), TRUE)
    expect_error(get.relative(multifun), NA)

    # Check devdev
    expect_error(devdev(result, multifun),NA)


    mult <- dmulti(
      c(rep(list(dpoly(degree=1)),2),
        rep(list(dspline(knots = 2, type="ns", beta.1=0.2)),1),
        rep(list(dfpoly(degree=2)),length(network$agents)-3)
      ))

    multifun <- mbnma.run(network, fun=mult,
                          method="random", n.iter=n.iter, pd=pd)
    expect_equal(all(c(paste0("beta.", c(1:6)), "power.1", "power.2", "sd") %in% multifun$parameters.to.save), TRUE)
    expect_equal(length(multifun$model.arg$fun$name)>1, TRUE)

    expect_error(plot(multifun), NA)
    expect_error(rank(get.relative(multifun)), NA)
    expect_error(devplot(multifun), NA)
    expect_error(fitplot(multifun), NA)
    expect_error(predict(multifun), NA)
    expect_error(suppressWarnings(summary(multifun)), NA)


    mult <- dmulti(c(list(dloglin()),
                     list(dspline("bs", knots=2)),
                     list(dspline("ns", knots=0.5)),
                     rep(list(dloglin()), length(network$agents)-3)
    ))
    multifun <- mbnma.run(network, fun=mult, n.iter=n.iter)
    expect_equal(all(c("rate", paste0("beta.", c(1,2,3,5,6))) %in% multifun$parameters.to.save), TRUE)
    expect_error(plot(multifun), NA)
    expect_error(rank(get.relative(multifun)), NA)
    expect_error(devplot(multifun), NA)
    expect_error(fitplot(multifun), NA)
    expect_error(predict(multifun), NA)
    expect_error(suppressWarnings(summary(multifun)), NA)


    if ("class" %in% names(dataset)) {
      expect_error(mbnma.run(netclass, fun=mult,
                             class.effect = list(beta.2="common")), "single dose-response function")
    }
  })
}
