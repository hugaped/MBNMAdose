testthat::context("Testing inconsistency.functions")

### Datasets ####
network <- mbnma.network(HF2PPITT)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)


testthat::test_that("test.inconsistency.loops", {
  expect_equal(nrow(inconsistency.loops(network$data.ab)), 4)
  expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), 7) # more loops since ref treatment has changed

  incon <- inconsistency.loops(network$data.ab)
  expect_identical(names(incon), c("t1", "t2", "path"))
})




testthat::test_that("test.nma.nodesplit", {

  split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
                           method="common", n.iter=1000)
  expect_equal(nrow(inconsistency.loops(network$data.ab)), length(split))
  expect_equal(class(split), "nma.nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "nma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "direct.model", "indirect.model", "nma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  split <- nma.nodesplit(net.noplac, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000)
  expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), length(split))
  expect_equal(class(split), "nma.nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "nma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "direct.model", "indirect.model", "nma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  # Test drop.discon
  split <- nma.nodesplit(net.noplac, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000, drop.discon = FALSE)
  expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), length(split))
  expect_equal(class(split), "nma.nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "nma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "direct.model", "indirect.model", "nma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  # Test comparisons
  split <- nma.nodesplit(net.noplac, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000, drop.discon = FALSE,
                           comparisons = rbind(c(9,10), c(7,11)))
  expect_equal(2, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000, drop.discon = FALSE,
                           comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5")))
  expect_equal(1, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  expect_error(nma.nodesplit(network, likelihood = "binomial", link="logit",
                               method="random", n.iter=1000, drop.discon = FALSE,
                               comparisons = rbind(c("badger","rizatriptan_0.5"))),
               "Treatment names given")

  expect_error(nma.nodesplit(network, likelihood = "binomial", link="logit",
                               method="random", n.iter=1000, drop.discon = FALSE,
                               comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                   c("zolmitriptan_4", "eletriptan_1"),
                                                   c("naratriptan_2", "Placebo_0"))))

})
