testthat::context("Testing run.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
#datanam <- "HF2PPITT"
#dataset <- HF2PPITT

### Datasets ####
network <- mbnma.network(dataset)



# Make class data
df <- dataset
df1 <- dataset

if ("class" %in% names(dataset)) {
  netclass <- mbnma.network(df)
}

# Make data with no placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)


test_that(paste("mbnma.run functions correctly for:", datanam), {
  n.iter=500

  # Single parameter DR functions
  result <- mbnma.run(network, fun="linear", beta.1="rel", method="common",
                      pd="plugin", n.iter=n.iter)
  expect_equal(class(result), c("mbnma", "rjags"))
  expect_equal("d.1" %in% result$parameters.to.save, TRUE)
  expect_equal(result$BUGSoutput$pD<0, TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
                      pd="pd.kl", n.iter=n.iter)
  expect_equal(class(result), c("mbnma", "rjags"))
  expect_equal("sd" %in% result$parameters.to.save, TRUE)
  expect_error(summary(result), NA)

  if ("class" %in% names(dataset)) {
    result <- mbnma.run(netclass, fun="exponential", beta.1="rel", method="common",
                        pd="popt", class.effect = list(beta.1="random"), n.iter=n.iter)
    expect_equal(class(result), c("mbnma", "rjags"))
    expect_equal("D.1" %in% result$parameters.to.save, TRUE)
    expect_equal("sd.D.1" %in% result$parameters.to.save, TRUE)
    expect_error(summary(result), NA)

    result <- mbnma.run(network, fun="nonparam.up", method="common", n.iter=n.iter)
    expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
    expect_error(summary(result))
  }


  result <- mbnma.run(network, fun="nonparam.down", method="random", n.iter=n.iter)
  expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
  expect_equal("sd" %in% result$parameters.to.save, TRUE)
  expect_error(summary(result))

  expect_error(result <- mbnma.run(network, fun="linear", beta.1="rel", method="fixed", n.iter=n.iter))



  # Two parameter DR functions
  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="rel", method="common",
                      n.iter=n.iter)
  expect_equal(all(c("d.1", "d.2") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="rel", method="random",
                      n.iter=n.iter)
  expect_equal("sd" %in% result$parameters.to.save, TRUE)
  expect_error(summary(result), NA)

  if ("class" %in% names(dataset)) {
    expect_warning(mbnma.run(netclass, fun="emax", beta.1="rel", beta.2="rel", method="random",
                             class.effect=list(beta.2="common"), n.iter=n.iter))
    expect_error(summary(result), NA)

    result <- mbnma.run(netclass, fun="emax", beta.1="rel", beta.2="rel", method="common",
                        class.effect=list(beta.2="common"), n.iter=n.iter)
    expect_equal(all(c("d.1", "D.2") %in% result$parameters.to.save), TRUE)
    expect_error(summary(result), NA)
  }


  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="random", method="common",
                      n.iter=n.iter)
  expect_equal(all(c("d.1", "beta.2", "sd.2") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)


  # Three parameter DR function
  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="random", beta.3="common",
                      method="random", n.iter=n.iter)
  expect_equal(all(c("d.1", "sd", "beta.2", "sd.2") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.run(network, fun="rcs", knots=4, beta.1="rel", beta.2="random", beta.3="common",
                      method="random", n.iter=n.iter)
  expect_equal(all(c("d.1", "sd", "beta.2", "sd.2") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)
  expect_equal(grepl("spline", result$model.arg$jagscode), TRUE)

  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="random", beta.3="common",
                      parameters.to.save = "psi",
                      method="random", n.iter=n.iter)
  expect_equal("psi" %in% result$parameters.to.save, TRUE)
  expect_equal("d.1" %in% result$parameters.to.save, FALSE)
  expect_error(summary(result))

  expect_error(mbnma.run(net.noplac, fun="nonparam.up"))


  # Changing priors
  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="rel", method="random",
                      n.iter=n.iter)
  prior <- list(sd="dunif(0,5)", inv.R="dwish(Omega[,],5)")
  runprior <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="rel", method="random",
                        n.iter=n.iter, priors = prior)
  expect_equal(runprior$model.arg$priors$sd, prior$sd)
  expect_equal(runprior$model.arg$priors$inv.R, prior$inv.R)
  expect_equal(result$model.arg$priors$inv.R!=runprior$model.arg$priors$inv.R, TRUE)


  # Multiple dose-response functions
  multifun <- mbnma.run(network, fun=c(rep("exponential", 2), rep("linear",1), rep("emax",length(network$agents)-3)),
                        n.iter=n.iter)
  expect_equal(length(multifun$model.arg$fun)>1, TRUE)
  expect_equal(all(c("d.1", "d.2", "d.3", "d.4") %in% multifun$parameters.to.save), TRUE)

  multifun <- mbnma.run(network, fun=c(rep("exponential", 2), rep("rcs", 1), rep("linear", length(network$agents)-3)),
                        method="random", knots=3, n.iter=n.iter)
  expect_equal(all(c("d.1", "d.2") %in% multifun$parameters.to.save), TRUE)
  expect_equal(all(c("d.3", "d.4") %in% multifun$parameters.to.save), TRUE)

  if ("class" %in% names(dataset)) {
    expect_error(mbnma.run(netclass, fun=c(rep("exponential", 3), rep("linear",1)),
                           class.effect = list(beta.2="common")), "single dose-response function")
  }


})




test_that(paste("mbnma.run wrappers function correctly for:", datanam), {
  n.iter=500

  # Single parameter DR functions
  result <- mbnma.linear(network, slope="random", n.iter=n.iter)
  expect_equal(all(c("beta.slope", "sd.slope") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.exponential(network, lambda="rel", method="common", n.iter=n.iter)
  expect_equal(all(c("d.lambda") %in% result$parameters.to.save), TRUE)
  expect_equal(all(c("sd") %in% result$parameters.to.save), FALSE)
  expect_error(summary(result), NA)

  if ("class" %in% names(dataset)) {
    # Two parameter DR functions
    result <- mbnma.emax(netclass, emax="rel", ed50="rel", method="common",
                         class.effect=list(emax="common"), n.iter=n.iter, cor = FALSE)
    expect_equal(all(c("D.emax", "d.ed50") %in% result$parameters.to.save), TRUE)
    expect_equal(all(c("d.emax") %in% result$parameters.to.save), FALSE)
    expect_error(summary(result), NA)

    # Three parameter DR functions
    if (datanam!="osteopain_2wkabs") {
      result <- mbnma.emax.hill(netclass, emax="rel", ed50="rel", hill="common",
                                method="random", n.iter=n.iter)
      expect_equal(all(c("d.emax", "d.ed50", "d.hill", "sd") %in% result$parameters.to.save), TRUE)
      expect_error(summary(result), NA)
    }
  }

})



test_that(paste("check.likelink function correctly for:", datanam), {

  if (all(c("y", "se") %in% names(dataset))) {
    expect_silent(check.likelink(df, likelihood = "normal", link="identity"))
    expect_silent(check.likelink(df, likelihood = "normal", link="logit"))

    # Expect error due to misspecified df
    expect_error(check.likelink(df, likelihood = "binomial", link="identity"))
    expect_error(check.likelink(df, likelihood = "poisson", link="identity"))

    # Expect errror due to misspecified arguments
    expect_error(check.likelink(df, likelihood = "normal", link="badger"))
    expect_error(check.likelink(df, likelihood = "test", link="identity"))

  } else if (all(c("r", "N") %in% names(dataset))) {
    expect_silent(check.likelink(df, likelihood = "binomial", link="identity"))
    expect_silent(check.likelink(df, likelihood = "binomial", link="logit"))

    # Expect error due to misspecified df
    expect_error(check.likelink(df, likelihood = "normal", link="identity"))
    expect_error(check.likelink(df, likelihood = "poisson", link="identity"))

    # Expect errror due to misspecified arguments
    expect_error(check.likelink(df, likelihood = "binomial", link="badger"))
    expect_error(check.likelink(df, likelihood = "test", link="logit"))
  }

})




test_that(paste("nma.run function correctly for:", datanam), {
  n.iter <- 500

  # expect_warning(nma.run(network, method="random", n.iter=100, warn.rhat = TRUE))

  expect_warning(nma.run(network, method="common", n.iter=n.iter, warn.rhat = FALSE), NA)

  result <- nma.run(network, method="random", n.iter=n.iter, warn.rhat = FALSE)
  expect_equal(names(result), c("jagsresult", "trt.labs"))
  expect_equal(all(c("d", "sd") %in% result$jagsresult$parameters.to.save), TRUE)

  result <- nma.run(network, method="random", n.iter=n.iter, warn.rhat = FALSE,
                    UME=TRUE)
  expect_equal("d[1,1]" %in% rownames(result$jagsresult$BUGSoutput$summary), TRUE)


  # Creating a broken network
  df.num <- mbnma.network(df1)$data.ab

  sepcomp <- mbnma.comparisons(df.num)[nrow(mbnma.comparisons(df.num)),]
  keep <- df.num$studyID[df.num$treatment %in% c(sepcomp$t1, sepcomp$t2)]
  df.num <- df.num[!(df.num$studyID %in% keep & !df.num$treatment  %in% c(sepcomp$t1, sepcomp$t2)),]

  df.num <- df.num %>% group_by(studyID) %>% mutate(narm=n())
  df.num <- df.num[df.num$narm>1,]

  fullrow <- nrow(df.num)
  network.disc <- mbnma.network(df.num)

  result.1 <- nma.run(network.disc, method="random", n.iter=n.iter, warn.rhat = FALSE,
                    UME=TRUE, drop.discon = TRUE)
  result.2 <- nma.run(network.disc, method="random", n.iter=n.iter, warn.rhat = FALSE,
                      UME=TRUE, drop.discon = FALSE)
  result.3 <- nma.run(network.disc, method="random", n.iter=n.iter, warn.rhat = FALSE,
                      UME=TRUE, drop.discon = TRUE)
  expect_equal(length(result.1$trt.labs)!=length(result.2$trt.labs), TRUE)
  expect_equal(length(result.1$trt.labs)==length(result.3$trt.labs), TRUE)
})



test_that(paste("pDcalc functions correctly for:", datanam), {
  n.iter=1000

  if (all(c("y", "se") %in% names(dataset))) {
    likelihood <- "normal"
    link <- "identity"
  } else if (all(c("r", "N") %in% names(dataset))) {
    likelihood <- "binomial"
    link <- "logit"


    # For binomial likelihood
    result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
                        parameters.to.save = c("psi", "resdev"),
                        n.iter=n.iter)

    jagsdata <- getjagsdata(network$data.ab, likelihood = likelihood, link=link)

    obs1 <- jagsdata$r
    obs2 <- jagsdata$N

    pd <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                 theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                 likelihood=likelihood, type="dose")
    expect_equal(length(pd),1)
    expect_equal(class(pd),"numeric")

    pd <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=5,
                 theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                 likelihood=likelihood, type="dose")
    expect_equal(length(pd),1)
    expect_equal(class(pd),"numeric")

    pd <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=5,
                 theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                 likelihood="poisson", type="dose")

    expect_error(pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                        theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                        likelihood="poisson", type="time"))

    expect_error(pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                        theta.result=NULL, resdev.result=result$BUGSoutput$mean$resdev,
                        likelihood="poisson", type="dose"))

    expect_error(pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=NULL,
                        theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                        likelihood="poisson", type="dose"))
  }

})





test_that(paste("mbnma.update function correctly for:", datanam), {

  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="rel", method="common",
                      n.iter=500)

  expect_error(mbnma.update(result, param="test", n.iter=100))

  update <- mbnma.update(result, param="resdev", n.iter=100)
  expect_equal(names(update), c("study", "arm", "mean", "facet", "fupdose", "groupvar"))

  update <- mbnma.update(result, param="theta", n.iter=100)
  expect_equal(names(update), c("study", "arm", "mean", "facet", "fupdose", "groupvar"))

})
