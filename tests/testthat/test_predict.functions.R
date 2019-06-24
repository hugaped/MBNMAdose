testthat::context("Testing predict.functions")

### Datasets ####
network <- MBNMA.network(HF2PPITT)
netgout <- MBNMA.network(GoutSUA_2wkCFB)
netalog <- MBNMA.network(alog_pcfb)
netpain <- MBNMA.network(osteopain_2wkabs)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- MBNMA.network(noplac.df)


#### Models ####

linear <- MBNMA.run(netpain, fun="linear")

emax <- MBNMA.emax(network, emax="rel", ed50="rel", method="random")
emax.gout <- MBNMA.emax(netgout, emax="rel", ed50="rel", method="random")

emax.class <- MBNMA.emax(netclass, emax="rel", ed50="random", method="common",
                         class.effect=list(emax="random"))

nonparam <- MBNMA.run(network, fun="nonparam.up")

emax.noplac <- MBNMA.emax(net.noplac, emax="rel", ed50="rel", method="random")



testthat::test_that("ref.synth functions correctly", {
  ref.df <- network$data.ab[network$data.ab$agent==1,]

  result <- ref.synth(ref.df, mbnma=emax, synth="fixed")
  expect_equal(length(result), 1)
  expect_equal(names(result)=="m.mu", TRUE)
  expect_equal(nrow(result$m.mu), emax$BUGSoutput$n.sims)

  ref.df <- netalog$data.ab[netalog$data.ab$agent==2,]
  expect_error(ref.synth(ref.df, mbnma=linear, synth="fixed"))

  ref.df <- netalog$data.ab[netalog$data.ab$agent==2 & netalog$data.ab$dose==25,]
  expect_error(ref.synth(ref.df, mbnma=linear, synth="fixed"), NA)

  expect_error(ref.synth(ref.df, mbnma=emax.noplac, synth="arndom"))

  ref.df <- network$data.ab[network$data.ab$agent==1,]
  result <- ref.synth(ref.df, mbnma=emax.noplac, synth="random", n.iter=1500)
  expect_identical(names(result), c("m.mu", "sd.mu"))
  expect_equal(nrow(result$sd.mu), 1500)

})



testthat::test_that("rescale.link functions correctly", {

  x <- c(-5,1,0)

  y <- rescale.link(x, direction="link", link="identity")
  expect_identical(x,y)

  y <- rescale.link(x, direction="natural", link="identity")
  expect_identical(x,y)

  expect_silent(rescale.link(x, direction="natural", link="logit"))
  expect_warning(rescale.link(x, direction="link", link="logit"))

  expect_warning(rescale.link(x, direction="link", link="probit"))
  expect_silent(rescale.link(x, direction="natural", link="probit"))

})



testthat::test_that("predict.MBNMA functions correctly", {
  #ref.df <- network$data.ab[network$data.ab$agent==1,]

  # Estimating E0
  ref.df <- netalog$data.ab[netalog$data.ab$agent==1,]
  pred <- predict(linear, E0 = ref.df)
  expect_identical(names(pred), c("predicts", "likelihood", "link"))
  expect_equal(linear$model.arg$likelihood, pred$likelihood)
  expect_equal(linear$model.arg$link, pred$link)
  expect_identical(names(pred$predicts), linear$agents)
  expect_silent(as.numeric(names(pred$predicts[[4]])))
  expect_equal(class(pred$predicts[[4]][[4]]), "matrix")
  expect_equal(nrow(pred$predicts[[4]][[4]]), linear$BUGSoutput$n.sims)
  expect_equal(all(pred$predicts[[2]][[2]][1] < 0), TRUE)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  # Stochastic E0 values
  expect_silent(predict(linear, E0 = "rnorm(n, 0.5,0.01)"))
  expect_silent(predict(linear, E0 = "rbeta(n, shape1=1, shape2=5)"))
  expect_error(predict(linear, E0 = "badgers(n, shape1=1, shape2=5)"))
  expect_error(predict(linear, E0 = "rbeta(badgers, shape1=1, shape2=5)"))

  # Determinsitic E0 values
  expect_silent(predict(linear, E0 = 0.01))
  expect_silent(predict(linear, E0 = 0.99))
  expect_warning(predict(emax, E0 = 1.5))

  # Changing n.doses
  pred <- predict(linear, E0=0.5, n.doses = 10)
  expect_equal(length(pred$predicts[[2]]), 10)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  # Changing max.doses
  max.doses <- list()
  for (i in seq_along(network$agents)) {
    max.doses[[length(max.doses)+1]] <- 1
  }
  pred <- predict(emax, E0=0.1, max.doses = max.doses)
  expect_identical(names(pred$predicts), emax$agents)
  expect_equal(all(pred$predicts[[2]][[2]][1] > 0), TRUE)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  names(max.doses) <- emax$agents
  expect_silent(predict(emax, E0=0.1, max.doses = max.doses))

  max.doses[[9]] <- 1
  expect_error(predict(linear, E0=0.1, max.doses = max.doses))

  max.doses <- list("eletriptan"=3, "rizatriptan"=2)
  pred <- predict(emax, E0=0.1, max.doses = max.doses, n.doses = 10)
  expect_equal(length(pred$predicts), length(max.doses))
  expect_equal(names(pred$predicts), names(max.doses))
  expect_equal(names(pred$predicts$eletriptan)[10], "3")
  expect_equal(names(pred$predicts$rizatriptan)[10], "2")
  expect_equal(all(pred$predicts[[1]][[1]][1] > 0), TRUE)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  max.doses <- list("badger"=3, "test"=1)
  expect_error(predict(linear, E0=0.1, max.doses = max.doses))

  max.doses <- list("eletriptan"="badger", "rizatriptan"=2)
  expect_error(predict(linear, E0=0.1, max.doses = max.doses, n.doses = 10))


  # Changing exact.doses
  doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  pred <- predict(emax, E0=0.1, exact.doses = doses)
  expect_identical(as.numeric(names(pred$predicts[[1]])), doses[[1]])
  expect_identical(as.numeric(names(pred$predicts[[2]])), doses[[2]])
  expect_equal(all(pred$predicts[[2]][[2]][1] > 0), TRUE)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  doses <- list(c(0,1,2,3), c(0.5,1,2))
  expect_error(predict(linear, E0=0.1, exact.doses = doses))

  dose <- c(0,0.5,1,2,4)
  doses <- list()
  for (i in seq_along(network$agents)) {
    doses[[length(doses)+1]] <- dose
  }
  expect_silent(predict(emax, E0=0.1, exact.doses = doses))

  doses <- list("eletriptan"=c("I","am","a","test"), "rizatriptan"=c(0.5,1,2))
  expect_error(predict(emax, E0=0.1, exact.doses = doses))

  doses <- list("badger"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  expect_error(predict(emax, E0=0.1, exact.doses = doses))


  #### Repeat sections with emax.noplac ####
  max.doses <- list()
  for (i in seq_along(net.noplac$agents)) {
    max.doses[[length(max.doses)+1]] <- 1
  }
  pred <- predict(emax.noplac, E0=0.1, max.doses = max.doses)
  expect_identical(names(pred$predicts), emax.noplac$agents)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

})




