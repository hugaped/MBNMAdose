testthat::context("Testing rank.functions")
network <- MBNMA.network(HF2PPITT)








# NEEDS TO BE CHANGED FOR MBNMAdose
testthat::test_that("rank.MBNMA functions correctly", {
  model.list <- list(emax1, emax2, quad)
  treats.list <- list(c(1,5,8,15), c(2,6,9,16:20), c(1:4), c(5:9))
  i <- 1

  rank <- rank.MBNMA(exponential, params="d.lambda",
                     direction=-1, treats=treats.list[[i]])

  testthat::expect_equal(names(rank), "d.lambda")
  testthat::expect_equal(names(rank$d.lambda), c("summary", "prob.matrix", "rank.matrix"))
  testthat::expect_equal(nrow(rank$d.lambda[["summary"]]), length(treats.list[[i]]))
  testthat::expect_equal(nrow(rank$d.lambda[["prob.matrix"]]), ncol(rank$d.lambda[["prob.matrix"]]))
  testthat::expect_equal(nrow(rank$d.lambda[["prob.matrix"]]), length(treats.list[[i]]))
  testthat::expect_equal(nrow(rank$d.lambda[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)
  testthat::expect_equal(colnames(rank$d.lambda[["rank.matrix"]]), as.character(treats.list[[i]]))



  rank <- rank.MBNMA(emax1, params=c("d.emax", "d.et50"),
                     direction=-1, treats=treats.list[[i]])

  testthat::expect_equal(sort(names(rank)), sort(c("d.emax", "d.et50")))
  testthat::expect_equal(names(rank$d.et50), c("summary", "prob.matrix", "rank.matrix"))
  testthat::expect_equal(nrow(rank$d.et50[["summary"]]), length(treats.list[[i]]))
  testthat::expect_equal(nrow(rank$d.et50[["prob.matrix"]]), ncol(rank$d.et50[["prob.matrix"]]))
  testthat::expect_equal(nrow(rank$d.et50[["prob.matrix"]]), length(treats.list[[i]]))
  testthat::expect_equal(nrow(rank$d.et50[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)
  testthat::expect_equal(colnames(rank$d.et50[["rank.matrix"]]), as.character(treats.list[[i]]))

  i <- 2
  rank <- rank.MBNMA(quad, params=c("beta.2", "d.1"),
                     direction=-1, treats=treats.list[[i]])

  testthat::expect_equal(sort(names(rank)), sort(c("beta.2", "d.1")))
  testthat::expect_equal(names(rank$d.1), c("summary", "prob.matrix", "rank.matrix"))
  testthat::expect_equal(nrow(rank$d.1[["summary"]]), length(treats.list[[i]]))
  testthat::expect_equal(nrow(rank$d.1[["prob.matrix"]]), ncol(rank$d.1[["prob.matrix"]]))
  testthat::expect_equal(nrow(rank$beta.2[["prob.matrix"]]), length(treats.list[[i]]))
  testthat::expect_equal(nrow(rank$beta.2[["rank.matrix"]]), model.list[[i]]$BUGSoutput$n.sims)
  testthat::expect_equal(colnames(rank$beta.2[["rank.matrix"]]), as.character(treats.list[[i]]))


  testthat::expect_error(rank.MBNMA(emax1, params=c("beta.1", "beta.2"),
                                    direction=-1, treats=treats.list[[i]]))
})
