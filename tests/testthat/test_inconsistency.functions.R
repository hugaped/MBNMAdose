testthat::context("Testing inconsistency.functions")

### Datasets ####
network <- mbnma.network(HF2PPITT)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)


testthat::test_that("test.inconsistency.loops", {
  expect_equal(nrow(inconsistency.loops(network$data.ab)), 4)
  expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), 8) # more loops since ref treatment has changed

  incon <- inconsistency.loops(network$data.ab)
  expect_identical(names(incon), c("t1", "t2", "path"))
})




testthat::test_that("test.nma.nodesplit", {

  split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
                           method="common", n.iter=1000)
  expect_equal(nrow(inconsistency.loops(network$data.ab)), length(split))
  expect_equal(class(split), "nodesplit")
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
                           method="random", n.iter=1000, drop.discon = TRUE)
  expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), length(split))
  expect_equal(class(split), "nodesplit")
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
                           method="random", n.iter=1000, drop.discon = TRUE,
                           comparisons = rbind(c(17,20), c(6,2)))
  expect_equal(2, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000, drop.discon = TRUE,
                           comparisons = rbind(c("sumatriptan_1", "almotriptan_1")))
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





testthat::test_that("test.mbnma.nodesplit", {

  split <- mbnma.nodesplit(network, likelihood = "binomial", link="logit",
                         method="common", n.iter=1000)
  expect_equal(nrow(inconsistency.loops(network$data.ab)), length(split))
  expect_equal(class(split), "nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "nma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "split.model", "nma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  split <- mbnma.nodesplit(net.noplac, likelihood = "binomial", link="logit",
                         method="random", n.iter=1000, drop.discon = TRUE)
  expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), length(split))
  expect_equal(class(split), "nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "nma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "split.model", "nma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  # Test comparisons
  split <- mbnma.nodesplit(net.noplac, likelihood = "binomial", link="logit",
                         method="random", n.iter=1000, drop.discon = TRUE,
                         comparisons = rbind(c(17,20), c(6,2)))
  expect_equal(2, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  split <- mbnma.nodesplit(network, likelihood = "binomial", link="logit",
                         method="random", n.iter=1000, drop.discon = TRUE,
                         comparisons = rbind(c("sumatriptan_1", "almotriptan_1")))
  expect_equal(1, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  expect_error(mbnma.nodesplit(network, likelihood = "binomial", link="logit",
                             method="random", n.iter=1000, drop.discon = FALSE,
                             comparisons = rbind(c("badger","rizatriptan_0.5"))),
               "Treatment names given")

  expect_error(mbnma.nodesplit(network, likelihood = "binomial", link="logit",
                             method="random", n.iter=1000, drop.discon = FALSE,
                             comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                 c("zolmitriptan_4", "eletriptan_1"),
                                                 c("naratriptan_2", "Placebo_0"))))

})
