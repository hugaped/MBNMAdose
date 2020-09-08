testthat::context("Testing predict.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
#datanam <- "HF2PPITT"
#dataset <- HF2PPITT

### Datasets ####
network <- mbnma.network(dataset)


### Datasets ####
# network <- mbnma.network(HF2PPITT)
# netgout <- mbnma.network(GoutSUA_2wkCFB)
# netalog <- mbnma.network(alog_pcfb)
# netpain <- mbnma.network(osteopain_2wkabs)
# netclass <- mbnma.network(osteopain_2wkabs)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)


#### Models ####

linear <- mbnma.run(network, fun="linear", n.iter=1000)

emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random", n.iter=1000)
#emax.gout <- mbnma.emax(netgout, emax="rel", ed50="rel", method="random", n.iter=1000)

if ("class" %in% names(dataset)) {
  emax.class <- suppressWarnings(mbnma.emax(netclass, emax="rel", ed50="random", method="common",
                                            class.effect=list(emax="random"), n.iter=1000))
}

nonparam <- mbnma.run(network, fun="nonparam.up", n.iter=1000)

emax.noplac <- mbnma.emax(net.noplac, emax="rel", ed50="rel", method="random", n.iter=1000)



testthat::test_that(paste0("ref.synth functions correctly for: ", datanam), {
  ref.df <- network$data.ab[network$data.ab$agent==1,]

  result <- ref.synth(ref.df, mbnma=emax, synth="fixed")
  expect_equal(length(result), 1)
  expect_equal(names(result)=="m.mu", TRUE)
  expect_equal(nrow(result$m.mu), emax$BUGSoutput$n.sims)

  ref.df <- network$data.ab[network$data.ab$agent==2,]
  if (!length(unique(ref.df$studyID)) == nrow(ref.df)) {
    expect_error(ref.synth(ref.df, mbnma=linear, synth="fixed"), "contain >1 arm")
  }

  ref.df <- network$data.ab[network$data.ab$agent==network$data.ab$agent[10] & network$data.ab$dose==network$data.ab$dose[10],]
  expect_error(ref.synth(ref.df, mbnma=linear, synth="fixed"), NA)

  expect_error(ref.synth(ref.df, mbnma=emax.noplac, synth="arndom"))

  ref.df <- network$data.ab[network$data.ab$agent==1,]
  expect_error(ref.synth(ref.df, mbnma=emax.noplac, synth="random", n.iter=500, n.burnin=5000))
  result <- ref.synth(ref.df, mbnma=emax.noplac, synth="random", n.iter=1500, n.burnin = 500)
  expect_identical(names(result), c("m.mu", "sd.mu"))
  #expect_equal(nrow(result$sd.mu), 600)

})



testthat::test_that(paste0("rescale.link functions correctly for: ", datanam), {

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



testthat::test_that(paste0("predict.mbnma functions correctly for: ", datanam), {
  #ref.df <- network$data.ab[network$data.ab$agent==1,]

  # Estimating E0
  ref.df <- network$data.ab[network$data.ab$agent==1,]
  pred <- predict(linear, E0 = ref.df)
  expect_identical(names(pred), c("predicts", "likelihood", "link", "network"))
  expect_equal(linear$model.arg$likelihood, pred$likelihood)
  expect_equal(linear$model.arg$link, pred$link)
  expect_identical(names(pred$predicts), linear$network$agents)
  expect_silent(as.numeric(names(pred$predicts[[4]])))
  expect_equal("matrix" %in% class(pred$predicts[[4]][[4]]), TRUE)
  expect_equal(nrow(pred$predicts[[4]][[4]]), linear$BUGSoutput$n.sims)
  #expect_equal(all(pred$predicts[[2]][[2]][1] < 0), TRUE)
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

  if (all(c("y", "se") %in% names(dataset))) {
    expect_warning(predict(emax, E0 = 1.5), NA)
  } else if (all(c("r", "N") %in% names(dataset))) {
    expect_warning(predict(emax, E0 = 1.5))
  }


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
  expect_identical(names(pred$predicts), emax$network$agents)
  expect_equal(all(pred$predicts[[2]][[2]][1] > 0), TRUE)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  names(max.doses) <- emax$agents
  expect_silent(predict(emax, E0=0.1, max.doses = max.doses))

  max.doses[[length(max.doses)+1]] <- 1
  expect_error(predict(emax, E0=0.1, max.doses = max.doses), "A greater number of agents")


  max.doses <- list()
  max.doses[[network$agents[2]]] <- 3
  max.doses[[network$agents[5]]] <- 2
  pred <- predict(emax, E0=0.1, max.doses = max.doses, n.doses = 10)
  expect_equal(length(pred$predicts), length(max.doses))
  expect_equal(names(pred$predicts), names(max.doses))
  expect_equal(names(pred$predicts[[network$agents[2]]])[10], "3")
  expect_equal(names(pred$predicts[[network$agents[5]]])[10], "2")
  expect_equal(all(pred$predicts[[1]][[1]][1] > 0), TRUE)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")

  max.doses <- list("badger"=3, "test"=1)
  expect_error(predict(linear, E0=0.1, max.doses = max.doses))

  max.doses <- list("eletriptan"="badger", "rizatriptan"=2)
  expect_error(predict(linear, E0=0.1, max.doses = max.doses, n.doses = 10))


  # Changing exact.doses
  doses <- list()
  doses[[network$agents[3]]] <- c(0,1,2,3)
  doses[[network$agents[4]]] <- c(0.5,1,2)
  pred <- predict(emax, E0=0.1, exact.doses = doses)
  expect_identical(as.numeric(names(pred$predicts[[1]])), doses[[1]])
  expect_identical(as.numeric(names(pred$predicts[[2]])), doses[[2]])

  if (all(c("r", "N") %in% names(dataset))) {
    expect_equal(all(pred$predicts[[2]][[2]][1] > 0), TRUE)
  }

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

  doses <- list()
  doses[[network$agents[2]]] <- c("I","am","a","test")
  doses[[network$agents[4]]] <- c(0.5,1,2)
  expect_error(predict(emax, E0=0.1, exact.doses = doses))

  doses <- list("badger"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  expect_error(predict(emax, E0=0.1, exact.doses = doses))


  #### Repeat sections with emax.noplac ####
  max.doses <- list()
  for (i in seq_along(net.noplac$agents)) {
    max.doses[[length(max.doses)+1]] <- 1
  }
  pred <- predict(emax.noplac, E0=0.1, max.doses = max.doses)
  expect_identical(names(pred$predicts), emax.noplac$network$agents)
  expect_error(print(pred), NA)
  expect_equal(class(summary(pred)), "data.frame")


  # Multiple dose-response functions
  multifun <- mbnma.run(network, fun=c(rep("exponential", 3), rep("emax",length(network$agents)-3)))
  expect_silent(predict(multifun, E0=0.2))

  doses <- list()
  doses[[network$agents[2]]] <- c(0,2,0.1)
  doses[[network$agents[4]]] <- c(0,2,0.1)
  pred <- predict(multifun, E0=0.2, exact.doses = doses)
  expect_identical(names(pred$predicts), c(network$agents[2], network$agents[4]))

})




