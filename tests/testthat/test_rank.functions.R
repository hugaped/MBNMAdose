testthat::context("Testing rank.functions")

network <- MBNMA.network(HF2PPITT)

# Make class data
df <- HF2PPITT
df$class <- ifelse(df$agent=="placebo", "placebo", "active")
df$class <- ifelse(df$agent=="eletriptan", "active2", df$class)
netclass <- MBNMA.network(df)

# Models
linear.run <- MBNMA.run(network, fun="linear")

exponential <- MBNMA.exponential(network, lambda="rel", method="common")

emax <- MBNMA.emax(network, emax="rel", ed50="rel", method="random")

emax.class <- MBNMA.emax(netclass, emax="rel", ed50="random", method="common",
                         class.effect=list(emax="random"))

nonparam <- MBNMA.run(network, fun="nonparam.up")



testthat::test_that("rank.MBNMA functions correctly", {

  rank <- rank.MBNMA(linear.run)
  expect_equal(names(rank), "d.1")
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix"))
  expect_equal(class(rank[[1]]$summary), "data.frame")
  expect_equal(class(rank[[1]]$rank.matrix), "matrix")
  expect_equal(class(rank[[1]]$prob.matrix), "matrix")


  rank <- rank.MBNMA(emax)
  expect_equal(names(rank), c("d.emax", "d.ed50"))
  expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix"))
  expect_equal(class(rank[[2]]$summary), "data.frame")
  expect_equal(class(rank[[1]]$rank.matrix), "matrix")
  expect_equal(class(rank[[2]]$prob.matrix), "matrix")

  # Checking direction=1 and direction=-1 are opposites
  rank.down <- rank.MBNMA(emax, direction=-1)
  expect_equal(rank.down$d.emax$summary$rank.param[rank.down$d.emax$summary$`50%`==1]==
                 rank$d.emax$summary$rank.param[rank$d.emax$summary$`50%`==7],
               TRUE)

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

  expect_error(rank.MBNMA(nonparam))


  # Test params
  rank <- rank.MBNMA(emax)
  expect_equal(names(rank), c("d.emax", "d.ed50"))
  rank <- rank.MBNMA(emax, params="d.ed50")
  expect_equal(names(rank), c("d.ed50"))
  expect_error(rank.MBNMA(emax, params="test"))

})
