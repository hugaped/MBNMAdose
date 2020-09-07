testthat::context("Testing inconsistency.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
datanam <- "ssri"
dataset <- ssri

### Datasets ####
network <- mbnma.network(dataset)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)


testthat::test_that(paste0("test.inconsistency.loops for: ", datanam), {

  expect_error(inconsistency.loops(network$data.ab, incldr = TRUE), NA)
  comps <- inconsistency.loops(network$data.ab, incldr = TRUE)

  expect_equal(any(grepl("drparams", comps$path)), TRUE)

  compsnodr <- inconsistency.loops(network$data.ab, incldr = FALSE)
  expect_equal(nrow(compsnodr)<nrow(comps), TRUE)

  if (datanam=="HF2PPITT") {
    expect_equal(nrow(inconsistency.loops(network$data.ab)), 4)
    expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), 8) # more loops since ref treatment has changed
  }

  incon <- inconsistency.loops(network$data.ab)
  expect_identical(names(incon), c("t1", "t2", "path"))
})




testthat::test_that(paste0("test.nma.nodesplit for: ", datanam), {

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
  expect_equal(is.numeric(split[[1]]$indirect), TRUE)
  expect_equal(is.numeric(split[[1]]$direct), TRUE)
  expect_equal(length(split[[1]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[1]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  if (datanam=="HF2PPITT") {
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
  }


  # Test comparisons
  split <- nma.nodesplit(net.noplac, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000, drop.discon = TRUE,
                           comparisons = rbind(c(17,20), c(6,2)))
  expect_equal(2, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  comps <- inconsistency.loops(network$data.ab, incldr = TRUE)
  compsi <- c(comps$t1[1], comps$t2[1])
  split <- nma.nodesplit(network, likelihood = "binomial", link="logit",
                           method="random", n.iter=1000, drop.discon = TRUE,
                           comparisons = rbind(compsi))
  expect_equal(1, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  expect_error(nma.nodesplit(network, likelihood = "binomial", link="logit",
                               method="random", n.iter=1000, drop.discon = FALSE,
                               comparisons = rbind(c("badger","rizatriptan_0.5"))),
               "Treatment names given")

  if (datanam=="HF2PPITT") {
    expect_error(nma.nodesplit(network, likelihood = "binomial", link="logit",
                               method="random", n.iter=1000, drop.discon = FALSE,
                               comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                   c("zolmitriptan_4", "eletriptan_1"),
                                                   c("naratriptan_2", "Placebo_0"))))

  }

})





testthat::test_that("test.mbnma.nodesplit", {

  comps <- inconsistency.loops(network$data.ab, incldr = TRUE)

  split <- mbnma.nodesplit(network, fun="rcs", knots=3, likelihood = "binomial", link="logit",
                         method="common", n.iter=1000)
  expect_equal(nrow(comps), length(split))
  expect_equal(nrow(inconsistency.loops(network$data.ab, incldr = FALSE))==length(split), FALSE)
  expect_equal(class(split), "nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "mbnma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "split.model", "mbnma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  comps.noplac <- inconsistency.loops(net.noplac$data.ab, incldr = TRUE)
  split <- mbnma.nodesplit(net.noplac, fun="exponential", likelihood = "binomial", link="logit",
                         method="random", n.iter=1000)
  expect_equal(nrow(comps.noplac), length(split))
  expect_equal(class(split), "nodesplit")
  expect_identical(names(split[[1]]), c("comparison",
                                        "direct", "indirect", "mbnma",
                                        "overlap matrix", "p.values", "quantiles",
                                        "forest.plot", "density.plot",
                                        "split.model", "mbnma.model"))
  expect_equal(is.numeric(split[[1]]$p.values), TRUE)
  expect_equal(is.numeric(split[[2]]$indirect), TRUE)
  expect_equal(is.numeric(split[[3]]$direct), TRUE)
  expect_equal(length(split[[4]]$comparison), 2)
  expect_identical(class(split[[1]]$forest.plot), c("gg", "ggplot"))
  expect_identical(class(split[[2]]$density.plot), c("gg", "ggplot"))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")


  # Test comparisons
  split <- mbnma.nodesplit(net.noplac, fun="user", user.fun=~beta.1 * dose + beta.2 * (dose^2),
                           likelihood = "binomial", link="logit",
                         method="random", n.iter=1000,
                         comparisons = rbind(comps.noplac[1,], comps.noplac[3,]))
  expect_equal(2, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  comps <- inconsistency.loops(network$data.ab, incldr = TRUE)
  split <- mbnma.nodesplit(network, fun="emax", likelihood = "binomial", link="logit",
                         method="random", n.iter=1000,
                         comparisons = rbind(c(network$treatment[comps$t1[2]], network$treatment[comps$t2[2]])))
  expect_equal(1, length(split))
  expect_error(print(split), NA)
  expect_equal(class(summary(split)), "data.frame")

  expect_error(mbnma.nodesplit(network, fum="linear", likelihood = "binomial", link="logit",
                             method="random", n.iter=1000,
                             comparisons = rbind(c("badger","rizatriptan_0.5"))),
               "Treatment names given")

  expect_error(mbnma.nodesplit(network, fun="emax", likelihood = "binomial", link="logit",
                             method="random", n.iter=1000,
                             comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                 c("zolmitriptan_4", "eletriptan_1"),
                                                 c("naratriptan_2", "Placebo_0"))),
               "Treatment names given")

})
