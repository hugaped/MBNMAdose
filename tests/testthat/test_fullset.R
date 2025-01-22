testthat::context("Testing full set of functions")

# Includes tests for mbnma.run()

# Across a range of datasets and dose-response functions:
# Tests running
# Tests default plots of running (including fitplot and devplot)
# Tests default ranking (including cumrank)
# Tests default prediction
# Occasionally tests get.relative
# Includes tests for meta-regression on SSRI dataset


# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)

# Datasets with no placebo
network <- mbnma.network(psoriasis90)
psoriasis90.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

network <- mbnma.network(ssri)
ssri.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

# Make SSRI a class dataset
ssri <- ssri %>% dplyr::mutate(class= dplyr::case_when(grepl("pram", agent) ~ "Pram",
                                              grepl("etine", agent) ~ "Etine",
                                              grepl("aline", agent) ~ "Aline",
                                              TRUE ~ "Placebo"))

alldfs <- list(triptans, psoriasis90.noplac, osteopain, psoriasis75, ssri, ssri.noplac, gout)
datanams <- c("triptans", "psoriasis90.noplac", "osteopain", "psoriasis75", "ssri", "ssri.noplac", "gout")



for (dat in seq_along(alldfs)) {

  datanam <- datanams[dat]
  dataset <- alldfs[[dat]]

  test_that(paste("Testing full set of functions for:", datanam), {

    # Add sdscale to osteopain
    if ("osteopain" %in% datanam[dat]) {
      dataset$standsd <- 2
      dataset$standsd[dataset$studyID %in% c(1:5)] <- 1.5
    }

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
      pD <- FALSE

      if ("stansd" %in% names(network$data.ab)) {
        sdscale <- TRUE
      } else {
        sdscale <- FALSE
      }

      # NMA
      nma <- suppressWarnings(nma.run(network, method="random", pD=pD, n.iter=500, sdscale=sdscale))
      expect_error(plot(nma), NA)

      # Nonparm MBNMA
      if (grepl("noplac", datanam)) {
        expect_error(mbnma.run(network, fun=dnonparam(), method="common",
                               pD=pD, n.iter=n.iter, sdscale=sdscale), "must currently include placebo")
      } else {
        mbnma <- mbnma.run(network, fun=dnonparam(), method="common",
                           pD=pD, n.iter=n.iter, sdscale=sdscale)
        expect_error(plot(mbnma), NA)
        expect_error(summary(mbnma), "Cannot use")
        expect_error(predict(mbnma), "does not work with non-parametric")
        expect_error(rank(mbnma), "cannot currently be performed for non-parametric")
        expect_error(get.relative(mbnma), "cannot be used on non-parametric")
      }

      # Single parameter DR functions
      result <- mbnma.run(network, fun=dpoly(degree=1), method="common",
                          pD=TRUE, n.iter=n.iter, sdscale=sdscale)
      expect_equal(class(result), c("mbnma", "rjags"))
      expect_equal("beta.1" %in% result$parameters.to.save, TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)
      expect_identical(sort(network$studyID), sort(result$model.arg$jagsdata$studyID))


      result <- mbnma.run(network, fun=dexp(), method="random",
                          pD=TRUE, n.iter=n.iter, sdscale=sdscale)
      expect_equal(class(result), c("mbnma", "rjags"))
      expect_equal("sd" %in% result$parameters.to.save, TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(devplot(result), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      if ("class" %in% names(dataset)) {
        result <- mbnma.run(netclass, fun=dexp(), method="common",
                                           pD=TRUE, class.effect = list(emax="random"), n.iter=n.iter,
                                           sdscale=sdscale)
        expect_equal(class(result), c("mbnma", "rjags"))
        expect_equal("EMAX" %in% result$parameters.to.save, TRUE)
        expect_equal("sd.EMAX" %in% result$parameters.to.save, TRUE)
        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(devplot(result), NA)
        expect_error(fitplot(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)
        expect_error(get.relative(result), NA)

        expect_warning(mbnma.run(netclass, fun=dexp(p.expon=TRUE), method="common", cor=TRUE,
                                           pD=TRUE, class.effect = list(emax="random"), n.iter=n.iter,
                                           sdscale=sdscale), "Class effects cannot be modelled with correlation")

      }



      # Two parameter DR functions
      result <- mbnma.run(network, fun=demax(), method="common",
                          n.iter=n.iter, pD=pD, sdscale=sdscale)
      expect_equal(all(c("emax", "ed50") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(devplot(result), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      result <- mbnma.run(network, fun=demax(ed50="common"), method="random",
                          n.iter=n.iter, pD=pD, sdscale=sdscale)
      expect_equal("sd" %in% result$parameters.to.save, TRUE)
      expect_equal("ed50" %in% rownames(result$BUGSoutput$summary), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      if ("class" %in% names(dataset)) {
        expect_error(mbnma.run(netclass, fun=demax(p.expon=TRUE), method="random", cor=TRUE,
                               class.effect=list(fakeparam="common"), sdscale=sdscale), "corresponding to dose-response parameters")

        expect_warning(mbnma.run(netclass, fun=demax(p.expon=TRUE), method="random", cor=TRUE,
                                 class.effect=list(ed50="common"), n.iter=n.iter, pD=pD, sdscale=sdscale))
        result <- suppressWarnings(mbnma.run(netclass, fun=demax(), method="random",
                                             class.effect=list(ed50="common"), n.iter=n.iter, pD=pD, sdscale=sdscale))
        expect_equal(all(c("emax", "ED50", "ed50", "sd") %in% result$parameters.to.save), TRUE)
        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)

      }


      result <- mbnma.run(network, fun=demax(ed50="random"), method="common",
                          n.iter=n.iter, pD=pD, sdscale=sdscale)
      expect_equal(all(c("emax", "ed50", "sd.ed50") %in% result$parameters.to.save), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      # Three parameter DR function
      result <- tryCatch(mbnma.run(network, fun=demax(emax="rel", ed50="random", hill="common"),
                                   method="random", n.iter=n.iter, pD=pD, priors = list(hill="dunif(0.1,3)"),
                                   sdscale=sdscale), error=function(e){})

      if (!is.null(result)) {
        expect_equal(all(c("emax", "sd", "ed50", "sd.ed50", "hill") %in% result$parameters.to.save), TRUE)
        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)
      }

      result <- mbnma.run(network, fun=dspline(type="ns", knots=3, beta.1="rel", beta.2="random", beta.3="common"),
                          method="random", n.iter=n.iter, pD=pD, sdscale=sdscale)
      expect_equal(all(c("beta.1", "beta.2", "sd.beta.2", "sd", "beta.3") %in% result$parameters.to.save), TRUE)
      expect_equal(any(grepl("spline", result$model.arg$jagscode)), TRUE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(devplot(result), NA)
      expect_error(fitplot(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)


      result <- tryCatch(mbnma.run(network, fun=demax(ed50="random", hill=1.2),
                          parameters.to.save = "psi", sdscale=sdscale,
                          method="random", n.iter=n.iter, pD=pD), error=function(e){})

      if (!is.null(result)) {
        expect_equal("psi" %in% result$parameters.to.save, TRUE)
        expect_equal("d.1" %in% result$parameters.to.save, FALSE)
        expect_error(plot(result))
        expect_error(rank(result))
        # expect_error(devplot(result), NA)
        # expect_error(fitplot(result), NA)
        expect_error(predict(result))
        expect_error(summary(result))
      }



      # Splines and polynomials
      result <- mbnma.run(network, fun=dspline(type="bs", knots=2, beta.1="common", beta.2 = "rel", beta.3="random"),
                          n.iter=n.iter, pD=pD, sdscale=sdscale)
      expect_equal(all(c("beta.1", "beta.2", "sd.beta.3", "beta.3") %in% result$parameters.to.save), TRUE)
      expect_equal(all(c("sd") %in% result$parameters.to.save), FALSE)
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), NA)
      expect_error(suppressWarnings(summary(result)), NA)

      result <- mbnma.run(network, fun=dspline(type="ns", knots=c(0.2,0.5), beta.1="common", beta.2 = "rel", beta.3="random"),
                          n.iter=n.iter, pD=pD, sdscale=sdscale)
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
                            method="random", n.iter=n.iter, pD=pD, sdscale=sdscale)
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
                          method="random", n.iter=n.iter, pD=pD, sdscale=sdscale, cor = TRUE)
      expect_equal("inv.R" %in% names(result$model.arg$priors), TRUE)
      expect_equal(ncol(result$model.arg$omega), 2)
      expect_error(get.relative(result), NA)


      omega <- matrix(c(10,2,2,10), byrow = TRUE, ncol=2)
      result <- mbnma.run(network, fun=demax(p.expon=TRUE), omega=omega, n.iter=n.iter, pD=pD, sdscale=sdscale, cor=TRUE)
      expect_equal(ncol(result$model.arg$omega), 2)
      expect_equal(result$model.arg$omega[1,1], omega[1,1])

      expect_error(mbnma.run(network, fun=demax(), omega=omega, n.iter=n.iter, pD=pD, sdscale=sdscale, cor=TRUE),
                   "cannot be modelled with truncated parameters")

      omega <- matrix(c(10,4,2,10), byrow = TRUE, ncol=2)
      expect_error(mbnma.run(network, fun=demax(), omega=omega, n.iter=n.iter, pD=pD, cor=TRUE,
                             sdscale=sdscale), "omega must be")

      result <- mbnma.run(network, fun=demax(), cor=FALSE, n.iter=n.iter, pD=pD, sdscale=sdscale)
      expect_equal("inv.R" %in% names(result$model.arg$priors), FALSE)


      # Test UME
      result <- mbnma.run(network, fun=dpoly(degree=3, beta.1="random", beta.2 = "rel", beta.3="rel"),
                          method="random", n.iter=n.iter, pD=pD, UME=TRUE, sdscale=sdscale)
      expect_equal(all(c("beta.1", "beta.2", "sd.beta.1", "beta.3", "sd") %in% result$parameters.to.save), TRUE)
      expect_equal(any(grepl("beta\\.2\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
      expect_equal(any(grepl("beta\\.3\\[2,1\\]", rownames(result$BUGSoutput$summary))), TRUE)
      expect_equal(any(grepl("beta\\.1\\[2,1\\]", rownames(result$BUGSoutput$summary))), FALSE)

      result <- mbnma.run(network, fun=dspline(knots=c(0.1, 0.3), type="ls"),
                          n.iter=n.iter, pD=pD, UME=TRUE, sdscale=sdscale)
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
        expect_error(mbnma.run(network, fun=demax(), link="probit", n.iter=n.iter, pD=pD,
                               sdscale=sdscale), NA)
        result <- mbnma.run(network, fun=dfpoly(degree=2), link="cloglog", n.iter=n.iter, pD=pD,
                            sdscale=sdscale)
        expect_equal(result$model.arg$link, "cloglog")

        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(devplot(result), NA)
        expect_error(fitplot(result), NA)
        expect_error(predict(result), NA)
        expect_error(suppressWarnings(summary(result)), NA)
      } else if (datanam %in% c("osteopain", "gout")) {

        if (datanam %in% "gout") {
          expect_error(mbnma.run(network, fun=dspline(knots=2, type="bs"), link="smd",
                                 n.iter=n.iter, pD=pD, sdscale=sdscale), "must be included in data\\.ab")
        } else if (datanam %in% "osteopain") {
          expect_error(mbnma.run(network, fun=dspline(knots=2, type="bs"), link="smd",
                                 n.iter=n.iter, pD=pD, sdscale=sdscale), NA)
        }

        result <- suppressWarnings(mbnma.run(network, fun=dexp(), link="log", n.iter=n.iter, pD=pD,
                                             sdscale=sdscale))
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
        result <- mbnma.run(network, fun=dnonparam(direction = "decreasing"), method="random",
                            n.iter=n.iter, pD=pD, sdscale=sdscale)
        expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
        expect_equal("sd" %in% result$parameters.to.save, TRUE)
        expect_error(summary(result))

        expect_error(result <- mbnma.run(network, fun=dpoly(), method="fixed", n.iter=n.iter, pD=pD), "Must be element of set")
      }



      # Changing priors
      result <- mbnma.run(network, fun=demax(p.expon=TRUE), method="random", cor=TRUE,
                          n.iter=n.iter, pD=pD, sdscale=sdscale)
      prior <- list(sd="dunif(0,5)", inv.R="dwish(omega[,],5)")
      runprior <- mbnma.run(network, fun=demax(p.expon=TRUE), method="random", cor=TRUE,
                            n.iter=n.iter, pD=pD, priors = prior, sdscale=sdscale)
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
                            n.iter=n.iter, pD=pD, sdscale=sdscale)
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
                            method="random", n.iter=n.iter, pD=pD, sdscale=sdscale)
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
      multifun <- mbnma.run(network, fun=mult, n.iter=n.iter, sdscale=sdscale)
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

    if (datanam=="ssri") {

      test_that(paste("regression functions run correctly for:", datanam), {

        skip_on_appveyor()
        skip_on_ci()
        skip_on_cran()

        n.iter=500
        pD <- FALSE

        # Edit dataset
        ssri.reg <- ssri %>%
          dplyr::mutate(
            class=dplyr::case_when(agent %in% c("citalopram", "escitalopram") ~ "alopram",
                                   agent %in% c("fluoxetine", "paroxetine") ~ "xetine",
                                   TRUE ~ "Sertraline")
          )

        ssri.reg <- ssri.reg %>%
          dplyr::mutate(x.weeks = weeks - mean(weeks, na.rm=TRUE)) %>% # For a continuous covariate
          dplyr::mutate(r.weeks=factor(weeks, levels=c(8,4,5,6,9,10))) # For a categorical covariate


        # Create network object
        ssrinet <- mbnma.network(ssri.reg)

        ssrimod.c <- mbnma.run(ssrinet, fun=dfpoly(degree=2),
                               regress=~r.weeks, regress.effect = "common",
                               n.iter=n.iter)

        # Regress for continuous weeks
        # Random effect modification across all agents vs Placebo
        ssrimod.r <- mbnma.run(ssrinet, fun=dfpoly(degree=2),
                               regress=~x.weeks, regress.effect = "random",
                               n.iter=n.iter)

        # Regress for continuous weeks
        # Separate effect modification for each agent vs Placebo
        ssrimod.a <- mbnma.run(ssrinet, fun=dfpoly(degree=2),
                               regress=~x.weeks, regress.effect = "agent",
                               n.iter=n.iter)

        # Separate effect modification for each agent vs Placebo
        ssrimod.c1 <- mbnma.run(ssrinet,
                                fun=dmulti(list(
                                  dfpoly(degree=2),
                                  dfpoly(degree=2),
                                  dfpoly(degree=1),
                                  dfpoly(degree=2),
                                  dfpoly(degree=2),
                                  dfpoly(degree=2)
                                )),
                                regress=~r.weeks, regress.effect = "class",
                                n.iter=n.iter)

        ssrimod.c2 <- mbnma.run(ssrinet,
                                dfpoly(degree=2),
                                class.effect = list(beta.1="random"),
                                regress = ~r.weeks,
                                regress.effect = "common",
                                n.iter=n.iter)

        modlist <- list(ssrimod.c, ssrimod.r, ssrimod.a, ssrimod.c1, ssrimod.c2)
        binlist <- list(ssrimod.c, ssrimod.c1, ssrimod.c2)
        contlist <- list(ssrimod.r, ssrimod.a)

        for (mod in seq_along(modlist)) {

          # Predict
          pred <- predict(modlist[[mod]])

          if (all.vars(modlist[[mod]]$model.arg$regress)=="r.weeks") {

            regvec <- sample(c(1,0,0,0,0), size=5)
            names(regvec) <- c("r.weeks10", "r.weeks4", "r.weeks5", "r.weeks6", "r.weeks9")

            predreg <- predict(modlist[[mod]], regress.vals = regvec)

            expect_error(predict(modlist[[mod]], regress.vals=c("pop"=1)), "must contain a single named regressor value for each covariate")

          } else if (all.vars(modlist[[mod]]$model.arg$regress)=="x.weeks") {
            predreg <- predict(modlist[[mod]], regress=c("x.weeks"=runif(1,0,10)))

          }

          # Plot predict
          expect_error(plot(pred), NA)
          expect_error(plot(predreg), NA)

          expect_warning(plot(pred, overlay.split = TRUE), "mismatches in results")

          # Print predict
          expect_error(print(pred), NA)
          expect_error(print(predreg), NA)

          # Summary predict
          expect_error(summary(pred), NA)
          expect_error(summary(predreg), NA)

          # Rank predict
          ranks <- rank(pred)
          ranksreg <- rank(predreg)

          expect_error(plot(ranks), NA)
          expect_error(plot(ranksreg), NA)

          expect_error(cumrank(ranks), NA)
          expect_error(cumrank(ranksreg), NA)

          expect_error(print(ranks), NA)
          expect_error(print(ranksreg), NA)

          # Rank MBNMA
          if (length(modlist[[mod]]$model.arg$fun$name)==1) {
            if (length(modlist[[mod]]$model.arg$class.effect)>0) {
              ranks <- rank(modlist[[mod]], level="class")
            } else {
              ranks <- rank(modlist[[mod]], level="agent")
            }

            expect_error(plot(ranks), NA)
            expect_error(cumrank(ranks), NA)
            expect_error(print(ranks), NA)
          }

          # Get relative
          if (all.vars(modlist[[mod]]$model.arg$regress)=="r.weeks") {

            regvec <- sample(c(1,0,0,0,0), size=5)
            names(regvec) <- c("r.weeks10", "r.weeks4", "r.weeks5", "r.weeks6", "r.weeks9")

            rels <- get.relative(lower.diag = modlist[[mod]],
                                 upper.diag = binlist[[sample(length(binlist), 1)]],
                                 regress.vals=regvec,
                                 lim="pred")

          } else if (all.vars(modlist[[mod]]$model.arg$regress)=="x.weeks") {

            regvec <- sample(c(1,0,0,0,0), size=5)
            names(regvec) <- c("r.weeks10", "r.weeks4", "r.weeks5", "r.weeks6", "r.weeks9")

            rels <- get.relative(lower.diag = modlist[[mod]],
                                 upper.diag = contlist[[sample(length(contlist), 1)]],
                                 regress.vals=c("x.weeks"=runif(1,0,10)),
                                 lim="pred")

            expect_error(get.relative(lower.diag = modlist[[mod]],
                                      upper.diag = binlist[[sample(length(binlist), 1)]],
                                      regress.vals=regvec,
                                      lim="pred"), "must contain a single named regressor")
          }

          expect_error(print(rels), NA)


        }


      })
    }

  })
}
