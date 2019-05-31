testthat::context("Testing run.functions")

network <- MBNMA.network(HF2PPITT)

# Make class data
df <- HF2PPITT
df$class <- ifelse(df$agent=="placebo", "placebo", "active")
netclass <- MBNMA.network(df)


test_that("MBNMA.run functions correctly", {
  n.iter=500

  # Single parameter DR functions
  result <- MBNMA.run(network, fun="linear", beta.1="rel", method="common",
                      pd="plugin", n.iter=n.iter)
  expect_equal(class(result), c("MBNMA", "rjags"))
  expect_equal("d.1" %in% result$parameters.to.save, TRUE)
  expect_equal(result$BUGSoutput$pD<0, TRUE)

  result <- MBNMA.run(network, fun="exponential", beta.1="rel", method="random",
                      pd="pd.kl", n.iter=n.iter)
  expect_equal(class(result), c("MBNMA", "rjags"))
  expect_equal("sd" %in% result$parameters.to.save, TRUE)

  result <- MBNMA.run(netclass, fun="exponential", beta.1="rel", method="common",
                      pd="popt", class.effect = list(beta.1="random"), n.iter=n.iter)
  expect_equal(class(result), c("MBNMA", "rjags"))
  expect_equal("D.1" %in% result$parameters.to.save, TRUE)
  expect_equal("sd.D.1" %in% result$parameters.to.save, TRUE)

  result <- MBNMA.run(network, fun="nonparam.up", method="common", n.iter=n.iter)
  expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)

  result <- MBNMA.run(network, fun="nonparam.down", method="random", n.iter=n.iter)
  expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
  expect_equal("sd" %in% result$parameters.to.save, TRUE)

  expect_error(result <- MBNMA.run(network, fun="linear", beta.1="rel", method="fixed", n.iter=n.iter))



  # Two parameter DR functions
  result <- MBNMA.run(network, fun="emax", beta.1="rel", beta.2="rel", method="common",
                      n.iter=n.iter)
  expect_equal(all(c("d.1", "d.2") %in% result$parameters.to.save), TRUE)

  result <- MBNMA.run(network, fun="emax", beta.1="rel", beta.2="rel", method="random",
                      n.iter=n.iter)
  expect_equal("sd" %in% result$parameters.to.save, TRUE)

  expect_warning(MBNMA.run(netclass, fun="emax", beta.1="rel", beta.2="rel", method="random",
                      class.effect=list(beta.2="common"), n.iter=n.iter))

  result <- MBNMA.run(netclass, fun="emax", beta.1="rel", beta.2="rel", method="common",
                      class.effect=list(beta.2="common"), n.iter=n.iter)
  expect_equal(all(c("d.1", "D.2") %in% result$parameters.to.save), TRUE)

  result <- MBNMA.run(network, fun="emax", beta.1="rel", beta.2="random", method="common",
                      n.iter=n.iter)
  expect_equal(all(c("d.1", "beta.2", "sd.2") %in% result$parameters.to.save), TRUE)


  # Three parameter DR function
  result <- MBNMA.run(network, fun="emax", beta.1="rel", beta.2="random", beta.3="common",
                      method="random", n.iter=n.iter)
  expect_equal(all(c("d.1", "sd", "beta.2", "sd.2", "beta.3") %in% result$parameters.to.save), TRUE)


  result <- MBNMA.run(network, fun="emax", beta.1="rel", beta.2="random", beta.3="common",
                      parameters.to.save = "psi",
                      method="random", n.iter=n.iter)
  expect_equal("psi" %in% result$parameters.to.save, TRUE)
  expect_equal("d.1" %in% result$parameters.to.save, FALSE)

})




test_that("MBNMA.run wrappers function correctly", {
  n.iter=500

  # Single parameter DR functions
  result <- MBNMA.linear(network, slope="random", n.iter=n.iter)
  expect_equal(all(c("beta.slope", "sd.slope") %in% result$parameters.to.save), TRUE)

  result <- MBNMA.exponential(network, lambda="rel", method="common", n.iter=n.iter)
  expect_equal(all(c("d.lambda") %in% result$parameters.to.save), TRUE)
  expect_equal(all(c("sd") %in% result$parameters.to.save), FALSE)

  # Two parameter DR functions
  result <- MBNMA.emax(netclass, emax="rel", ed50="rel", method="common",
                       class.effect=list(emax="common"), n.iter=n.iter)
  expect_equal(all(c("D.emax", "d.ed50") %in% result$parameters.to.save), TRUE)
  expect_equal(all(c("d.emax") %in% result$parameters.to.save), FALSE)

  # Three parameter DR functions
  result <- MBNMA.emax.hill(netclass, emax="rel", ed50="rel", hill="rel",
                            method="random", n.iter=n.iter)
  expect_equal(all(c("d.emax", "d.ed50", "d.hill", "sd") %in% result$parameters.to.save), TRUE)

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




test_that("NMA.run function correctly", {
  n.iter <- 500

  expect_warning(NMA.run(network, method="common", n.iter=n.iter, warn.rhat = TRUE))

  expect_warning(NMA.run(network, method="common", n.iter=n.iter, warn.rhat = FALSE), NA)

  result <- NMA.run(network, method="random", n.iter=n.iter, warn.rhat = FALSE)
  expect_equal(names(result), c("jagsresult", "trt.labs"))
  expect_equal(all(c("d", "sd") %in% result$jagsresult$parameters.to.save), TRUE)

  result <- NMA.run(network, method="random", n.iter=n.iter, warn.rhat = FALSE,
                    UME=TRUE)
  expect_equal("d[1,1]" %in% rownames(result$jagsresult$BUGSoutput$summary), TRUE)


  # Create broken network to test drop.discon
  df.num <- MBNMA.network(df)$data.ab
  df.num$dose[df.num$studyID==3 & df.num$agent==1] <- 1
  df.num$agent[df.num$studyID==3 & df.num$agent==1] <- 5
  df.num <- df.num[!(df.num$studyID %in% c(3,11,14,16,21,29,31,37,39,40,43,44,51,63,70)),]

  fullrow <- nrow(df.num)
  network.disc <- MBNMA.network(df.num)

  result.1 <- NMA.run(network.disc, method="random", n.iter=n.iter, warn.rhat = FALSE,
                    UME=TRUE, drop.discon = TRUE)
  result.2 <- NMA.run(network.disc, method="random", n.iter=n.iter, warn.rhat = FALSE,
                      UME=TRUE, drop.discon = FALSE)
  result.3 <- NMA.run(network.disc, method="random", n.iter=n.iter, warn.rhat = FALSE,
                      UME=TRUE, drop.discon = TRUE)
  expect_equal(length(result.1$trt.labs)!=length(result.2$trt.labs), TRUE)
  expect_equal(length(result.1$trt.labs)==length(result.3$trt.labs), TRUE)
})



test_that("pDcalc functions correctly", {

  # For binomial likelihood
  likelihood <- "binomial"
  result <- MBNMA.run(network, fun="exponential", beta.1="rel", method="random",
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
