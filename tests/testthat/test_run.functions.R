testthat::context("Testing run.functions")

network <- mbnma.network(HF2PPITT)

# Make class data
df <- HF2PPITT
df$class <- ifelse(df$agent=="placebo", "placebo", "active")
netclass <- mbnma.network(df)

# Make data with no placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)


test_that("mbnma.run functions correctly", {
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

  result <- mbnma.run(netclass, fun="exponential", beta.1="rel", method="common",
                      pd="popt", class.effect = list(beta.1="random"), n.iter=n.iter)
  expect_equal(class(result), c("mbnma", "rjags"))
  expect_equal("D.1" %in% result$parameters.to.save, TRUE)
  expect_equal("sd.D.1" %in% result$parameters.to.save, TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.run(network, fun="nonparam.up", method="common", n.iter=n.iter)
  expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
  expect_error(summary(result))

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

  expect_warning(mbnma.run(netclass, fun="emax", beta.1="rel", beta.2="rel", method="random",
                      class.effect=list(beta.2="common"), n.iter=n.iter))
  expect_error(summary(result), NA)

  result <- mbnma.run(netclass, fun="emax", beta.1="rel", beta.2="rel", method="common",
                      class.effect=list(beta.2="common"), n.iter=n.iter)
  expect_equal(all(c("d.1", "D.2") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="random", method="common",
                      n.iter=n.iter)
  expect_equal(all(c("d.1", "beta.2", "sd.2") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)


  # Three parameter DR function
  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="random", beta.3="common",
                      method="random", n.iter=n.iter)
  expect_equal(all(c("d.1", "sd", "beta.2", "sd.2", "beta.3") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)


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

})




test_that("mbnma.run wrappers function correctly", {
  n.iter=500

  # Single parameter DR functions
  result <- mbnma.linear(network, slope="random", n.iter=n.iter)
  expect_equal(all(c("beta.slope", "sd.slope") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)

  result <- mbnma.exponential(network, lambda="rel", method="common", n.iter=n.iter)
  expect_equal(all(c("d.lambda") %in% result$parameters.to.save), TRUE)
  expect_equal(all(c("sd") %in% result$parameters.to.save), FALSE)
  expect_error(summary(result), NA)

  # Two parameter DR functions
  result <- mbnma.emax(netclass, emax="rel", ed50="rel", method="common",
                       class.effect=list(emax="common"), n.iter=n.iter)
  expect_equal(all(c("D.emax", "d.ed50") %in% result$parameters.to.save), TRUE)
  expect_equal(all(c("d.emax") %in% result$parameters.to.save), FALSE)
  expect_error(summary(result), NA)

  # Three parameter DR functions
  result <- mbnma.emax.hill(netclass, emax="rel", ed50="rel", hill="rel",
                            method="random", n.iter=n.iter)
  expect_equal(all(c("d.emax", "d.ed50", "d.hill", "sd") %in% result$parameters.to.save), TRUE)
  expect_error(summary(result), NA)

})



test_that("check.likelink function correctly", {

  expect_silent(check.likelink(df, likelihood = "binomial", link="identity"))
  expect_silent(check.likelink(df, likelihood = "binomial", link="logit"))

  # Expect error due to misspecified df
  expect_error(check.likelink(df, likelihood = "normal", link="identity"))
  expect_error(check.likelink(df, likelihood = "poisson", link="identity"))

  # Expect errror due to misspecified arguments
  expect_error(check.likelink(df, likelihood = "binomial", link="badger"))
  expect_error(check.likelink(df, likelihood = "test", link="logit"))

})




test_that("nma.run function correctly", {
  n.iter <- 500

  expect_warning(nma.run(network, method="common", n.iter=n.iter, warn.rhat = TRUE))

  expect_warning(nma.run(network, method="common", n.iter=n.iter, warn.rhat = FALSE), NA)

  result <- nma.run(network, method="random", n.iter=n.iter, warn.rhat = FALSE)
  expect_equal(names(result), c("jagsresult", "trt.labs"))
  expect_equal(all(c("d", "sd") %in% result$jagsresult$parameters.to.save), TRUE)

  result <- nma.run(network, method="random", n.iter=n.iter, warn.rhat = FALSE,
                    UME=TRUE)
  expect_equal("d[1,1]" %in% rownames(result$jagsresult$BUGSoutput$summary), TRUE)


  # Create broken network to test drop.discon
  df.num <- mbnma.network(df)$data.ab
  df.num$dose[df.num$studyID==3 & df.num$agent==1] <- 1
  df.num$agent[df.num$studyID==3 & df.num$agent==1] <- 5
  df.num <- df.num[!(df.num$studyID %in% c(3,11,14,16,21,29,31,37,39,40,43,44,51,63,70)),]

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



test_that("pDcalc functions correctly", {
  n.iter=1000

  # For binomial likelihood
  likelihood <- "binomial"
  result <- mbnma.run(network, fun="exponential", beta.1="rel", method="random",
                      parameters.to.save = c("psi", "resdev"),
                      n.iter=n.iter)

  jagsdata <- getjagsdata(network$data.ab)

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

})





test_that("mbnma.update function correctly", {

  result <- mbnma.run(network, fun="emax", beta.1="rel", beta.2="rel", method="common",
                      n.iter=500)

  expect_error(mbnma.update(result, param="test"))

  update <- mbnma.update(result, param="resdev")
  expect_equal(names(update), c("study", "arm", "mean", "facet", "fupdose", "groupvar"))

  update <- mbnma.update(result, param="theta")
  expect_equal(names(update), c("study", "arm", "mean", "facet", "fupdose", "groupvar"))

})
