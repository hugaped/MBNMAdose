testthat::context("Testing rank.functions")

network <- MBNMA.network(HF2PPITT)

# Make class data
df <- HF2PPITT
df$class <- ifelse(df$agent=="placebo", "placebo", "active")
df$class <- ifelse(df$agent=="eletriptan", "active2", df$class)
netclass <- MBNMA.network(df)

# Make data with no placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- MBNMA.network(noplac.df)

# Models
linear.run <- MBNMA.run(network, fun="linear")

exponential <- MBNMA.exponential(network, lambda="rel", method="common")

emax <- MBNMA.emax(network, emax="rel", ed50="rel", method="random")

emax.class <- MBNMA.emax(netclass, emax="rel", ed50="random", method="common",
                         class.effect=list(emax="random"))

nonparam <- MBNMA.run(network, fun="nonparam.up")

emax.noplac <- MBNMA.emax(net.noplac, emax="rel", ed50="rel", method="random")

testthat::test_that("rank.MBNMA functions correctly", {

  rank <- rank.MBNMA(linear.run)
  expect_equal(names(rank), "d.1")
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix"))
  expect_equal(class(rank[[1]]$summary), "data.frame")
  expect_equal(class(rank[[1]]$rank.matrix), "matrix")
  expect_equal(class(rank[[1]]$prob.matrix), "matrix")
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")


  rank <- rank.MBNMA(emax)
  expect_equal(names(rank), c("d.emax", "d.ed50"))
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix"))
  expect_equal(class(rank[[2]]$summary), "data.frame")
  expect_equal(class(rank[[1]]$rank.matrix), "matrix")
  expect_equal(class(rank[[2]]$prob.matrix), "matrix")
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  expect_error(rank(emax, params=c("badger", "d.ed50")))

  # Checking direction=1 and direction=-1 are opposites
  rank.down <- rank.MBNMA(emax, direction=-1)
  expect_equal(rank.down$d.emax$summary$rank.param[rank.down$d.emax$summary$`50%`==1]==
                 rank$d.emax$summary$rank.param[rank$d.emax$summary$`50%`==7],
               TRUE)
  expect_error(print(rank.down), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  to.ranks <- c(2,5,6)
  rank <- rank.MBNMA(exponential, to.rank = to.ranks)
  expect_equal(ncol(rank$d.lambda$rank.matrix), length(to.ranks))
  expect_error(rank.MBNMA(exponential, to.rank = c(1,5,6)))
  expect_error(rank.MBNMA(exponential, to.rank = c("eletriptan", "sumatriptan")))

  # Test classes
  expect_error(rank.MBNMA(emax, level="class"))
  expect_error(rank.MBNMA(emax.class, level="agent"))
  rank <- rank.MBNMA(emax.class, level="class")
  expect_equal(ncol(rank$D.emax$rank.matrix), 2)
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  expect_error(rank.MBNMA(nonparam))


  # Test params
  rank <- rank.MBNMA(emax)
  expect_equal(names(rank), c("d.emax", "d.ed50"))
  rank <- rank.MBNMA(emax, params="d.ed50")
  expect_equal(names(rank), c("d.ed50"))
  expect_error(rank.MBNMA(emax, params="test"))
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

  # With no placebo data
  rank <- rank.MBNMA(emax.noplac)
  expect_equal(names(rank), c("d.emax", "d.ed50"))
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix"))
  expect_equal(class(rank[[2]]$summary), "data.frame")
  expect_equal(class(rank[[1]]$rank.matrix), "matrix")
  expect_equal(class(rank[[2]]$prob.matrix), "matrix")
  expect_error(print(rank), NA)
  expect_equal(class(summary(rank)[[1]]), "data.frame")

})






testthat::test_that("rank.MBNMA.predict functions correctly", {

  pred <- predict(linear.run, E0 = 0.5)
  rank <- rank.MBNMA.predict(pred)
  expect_equal(names(rank), "Predictions")
  expect_equal(names(rank$Predictions), c("summary", "prob.matrix", "rank.matrix", "direction"))
  expect_equal(class(rank$Predictions$summary), "data.frame")
  expect_equal(class(rank$Predictions$rank.matrix), "matrix")
  expect_equal(class(rank$Predictions$prob.matrix), "matrix")


  doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
                  exact.doses=doses)
  rank <- rank.MBNMA.predict(pred)
  expect_equal(names(rank), "Predictions")
  expect_equal(names(rank$Predictions), c("summary", "prob.matrix", "rank.matrix", "direction"))
  expect_equal(class(rank$Predictions$summary), "data.frame")
  expect_equal(class(rank$Predictions$rank.matrix), "matrix")
  expect_equal(class(rank$Predictions$prob.matrix), "matrix")

  expect_equal(nrow(rank$Predictions$summary), length(unlist(doses)))

  # Test direction
  rank.up <- rank.MBNMA.predict(pred, direction=-1)
  rank.down <- rank.MBNMA.predict(pred, direction=1)
  expect_equal(rank.down$Predictions$summary$rank.param[rank.down$Predictions$summary$`50%`==min(rank.down$Predictions$summary$`50%`)],
               rank.up$Predictions$summary$rank.param[rank.up$Predictions$summary$`50%`==max(rank.up$Predictions$summary$`50%`)]
               )

  # Test rank.doses
  doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
                  exact.doses=doses)
  rank <- rank.MBNMA.predict(pred, rank.doses = list("eletriptan"=2, "rizatriptan"=2))
  expect_equal(nrow(rank$Predictions$summary), 2)

  expect_error(rank.MBNMA.predict(pred, rank.doses = list("badger"=2, "rizatriptan"=2)), "badger")

  expect_error(rank.MBNMA.predict(pred, rank.doses = list("eletriptan"=c(2, 50, 100), "rizatriptan"=2)), "cannot be included in ranking: 50\\, 100")

})
