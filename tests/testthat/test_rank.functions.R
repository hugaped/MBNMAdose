testthat::context("Testing rank.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
#datanam <- "psoriasis"
#dataset <- psoriasis

network <- mbnma.network(dataset)

# Make class data
if ("class" %in% names(dataset)) {
  netclass <- mbnma.network(df)

  emax.class <- suppressWarnings(mbnma.emax(netclass, emax="rel", ed50="random", method="common",
                                            class.effect=list(emax="random"), n.iter=1000))
}

# Make data with no placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)

# Models
linear.run <- mbnma.run(network, fun="linear", n.iter=1000)

exponential <- mbnma.exponential(network, lambda="rel", method="common", n.iter=1000)

emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random", n.iter=1000)

nonparam <- mbnma.run(network, fun="nonparam.up", n.iter=1000)

emax.noplac <- mbnma.emax(net.noplac, emax="rel", ed50="rel", method="random", n.iter=1000)




testthat::test_that(paste0("rank.mbnma functions correctly for: ", datanam), {

  rank <- rank.mbnma(linear.run)
  expect_equal(names(rank), "d.1")
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "direction"))
  expect_equal(class(rank[[1]]$summary), "data.frame")
  expect_equal("matrix" %in% class(rank[[1]]$rank.matrix), TRUE)
  expect_equal("matrix" %in% class(rank[[1]]$prob.matrix), TRUE)
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")


  rank <- rank.mbnma(emax)
  expect_equal(sort(names(rank)), sort(c("d.emax", "d.ed50")))
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "direction"))
  expect_equal(class(rank[[2]]$summary), "data.frame")
  expect_equal("matrix" %in% class(rank[[1]]$rank.matrix), TRUE)
  expect_equal("matrix" %in% class(rank[[2]]$prob.matrix), TRUE)
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  expect_error(rank(emax, params=c("badger", "d.ed50")))

  # Checking direction=1 and direction=-1 are opposites
  rank.down <- rank(emax, direction=-1)
  expect_equal(dplyr::arrange(rank.down$d.emax$summary, '50%')$rank.param[1] %in%
                 dplyr::arrange(rank$d.emax$summary, '50%')$rank.param[nrow(rank$d.emax$summary)-1:nrow(rank$d.emax$summary)],
               TRUE)
  expect_error(print(rank.down), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  to.ranks <- c(2,5)
  rank <- rank(exponential, to.rank = to.ranks)
  expect_equal(ncol(rank$d.lambda$rank.matrix), length(to.ranks))
  expect_warning(rank.mbnma(exponential, to.rank = c(1,3,4)))
  expect_silent(rank.mbnma(exponential, to.rank = c(network$agents[2], network$agents[3])))

  # Test classes
  if ("class" %in% names(dataset)) {
    expect_error(rank.mbnma(emax, level="class"))
    expect_error(rank.mbnma(emax.class, level="agent"))
    rank <- rank.mbnma(emax.class, level="class")
    expect_equal(ncol(rank$D.emax$rank.matrix), length(unique(dataset$class[dataset$dose>0])))
    expect_error(print(rank), NA)
    expect_equal(class(summary(rank)[[1]]), "data.frame")
  }


  expect_error(rank.mbnma(nonparam))


  # Test params
  rank <- rank.mbnma(emax)
  expect_equal(sort(names(rank)), sort(c("d.emax", "d.ed50")))
  rank <- rank.mbnma(emax, params="d.ed50")
  expect_equal(names(rank), c("d.ed50"))
  expect_error(rank.mbnma(emax, params="test"))
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  # With no placebo data
  rank <- rank.mbnma(emax.noplac)
  expect_equal(sort(names(rank)), sort(c("d.emax", "d.ed50")))
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "direction"))
  expect_equal(class(rank[[2]]$summary), "data.frame")
  expect_equal("matrix" %in% class(rank[[1]]$rank.matrix), TRUE)
  expect_equal("matrix" %in% class(rank[[2]]$prob.matrix), TRUE)
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  # With multiple-dose response functions
  multifun <- mbnma.run(network, fun=c(rep("exponential", 3), rep("emax",length(network$agents)-3)))
  expect_error(rank(multifun), "multiple dose-response")

})






testthat::test_that(paste0("rank.mbnma.predict functions correctly for: ", datanam), {

  pred <- predict(linear.run, E0 = 0.5)
  rank <- rank.mbnma.predict(pred)
  expect_equal(names(rank), "Predictions")
  expect_equal(names(rank$Predictions), c("summary", "prob.matrix", "rank.matrix", "direction"))
  expect_equal(class(rank$Predictions$summary), "data.frame")
  expect_equal("matrix" %in% class(rank$Predictions$rank.matrix), TRUE)
  expect_equal("matrix" %in% class(rank$Predictions$prob.matrix), TRUE)


  #doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  doses <- list()
  doses[[network$agents[2]]] <- c(0,1,2,3)
  doses[[network$agents[4]]] <- c(0.5,1,2)
  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
                  exact.doses=doses)
  rank <- rank.mbnma.predict(pred)
  expect_equal(names(rank), "Predictions")
  expect_equal(names(rank$Predictions), c("summary", "prob.matrix", "rank.matrix", "direction"))
  expect_equal(class(rank$Predictions$summary), "data.frame")
  expect_equal("matrix" %in% class(rank$Predictions$rank.matrix), TRUE)
  expect_equal("matrix" %in% class(rank$Predictions$prob.matrix), TRUE)

  expect_equal(nrow(rank$Predictions$summary), length(unlist(doses)))

  # Test direction
  rank.up <- rank.mbnma.predict(pred, direction=-1)
  rank.down <- rank.mbnma.predict(pred, direction=1)
  expect_equal(rank.down$Predictions$summary$rank.param[rank.down$Predictions$summary$`50%`==min(rank.down$Predictions$summary$`50%`)],
               rank.up$Predictions$summary$rank.param[rank.up$Predictions$summary$`50%`==max(rank.up$Predictions$summary$`50%`)]
               )

  # Test rank.doses
  doses <- list()
  doses[[network$agents[2]]] <- c(0,1,2,3)
  doses[[network$agents[4]]] <- c(0.5,1,2)
  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
                  exact.doses=doses)

  doses[[network$agents[2]]] <- 2
  doses[[network$agents[4]]] <- 2
  rank <- rank.mbnma.predict(pred, rank.doses = doses)
  expect_equal(nrow(rank$Predictions$summary), 2)

  expect_error(rank.mbnma.predict(pred, rank.doses = list("badger"=2, "rizatriptan"=2)), "badger")

  doses[[network$agents[2]]] <- c(2, 50, 100)
  doses[[network$agents[4]]] <- 2
  expect_error(rank.mbnma.predict(pred, rank.doses = doses), "cannot be included in ranking: 50\\, 100")

})
