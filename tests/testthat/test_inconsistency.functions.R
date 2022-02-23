testthat::context("Testing inconsistency.functions")

alldats <- c("triptans", "psoriasis75", "ssri", "osteopain", "gout")


for (dat in seq_along(alldats)) {

  # Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
  datanam <- alldats[dat]
  dataset <- get(alldats[dat])

  ### Datasets ####
  network <- mbnma.network(dataset)

  # Generate data without placebo
  noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
  net.noplac <- mbnma.network(noplac.df)


  testthat::test_that(paste0("test.inconsistency.loops for: ", datanam), {

    expect_error(inconsistency.loops(network$data.ab, incldr = TRUE), NA)
    comps <- inconsistency.loops(network$data.ab, incldr = TRUE)

    expect_equal(any(grepl("drparams", comps$path)), TRUE)

    if (!datanam %in% c("osteopain", "gout")) {
      compsnodr <- inconsistency.loops(network$data.ab, incldr = FALSE)
      expect_equal(nrow(compsnodr)<nrow(comps), TRUE)

      incon <- inconsistency.loops(network$data.ab)
      expect_identical(names(incon), c("t1", "t2", "path"))
    }

    if (datanam=="triptans") {
      expect_equal(nrow(inconsistency.loops(network$data.ab)), 4)
      expect_equal(nrow(inconsistency.loops(net.noplac$data.ab)), 8) # more loops since ref treatment has changed
    }

  })




  testthat::test_that(paste0("test.nma.nodesplit for: ", datanam), {

    if (all(c("y", "se") %in% names(dataset))) {
      like <- "normal"
      link <- "identity"
    } else if (all(c("r", "n") %in% names(dataset))) {
      like <- "binomial"
      link <- "logit"
    }

    if (!datanam %in% c("osteopain", "gout")) {
      split <- nma.nodesplit(network, likelihood = like, link=link,
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


      # Test comparisons
      comps <- inconsistency.loops(network$data.ab, incldr = FALSE)
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

    }


    if (datanam=="triptans") {
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


    if (datanam=="triptans") {
      expect_error(nma.nodesplit(network, likelihood = "binomial", link="logit",
                                 method="random", n.iter=1000, drop.discon = FALSE,
                                 comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                     c("zolmitriptan_4", "eletriptan_1"),
                                                     c("naratriptan_2", "Placebo_0"))))

    }

    # Prevent skip
    expect_equal(5,5)

  })





  testthat::test_that("test.mbnma.nodesplit", {

    if (all(c("y", "se") %in% names(dataset))) {
      like <- "normal"
      link <- "identity"
    } else if (all(c("r", "n") %in% names(dataset))) {
      like <- "binomial"
      link <- "logit"
    }

    comps <- inconsistency.loops(network$data.ab, incldr = TRUE)

    split <- mbnma.nodesplit(network, fun=dspline(type="ns", knots=2), likelihood = like, link=link,
                             method="common", n.iter=1000)
    expect_equal(nrow(comps), length(split))
    if (!datanam %in% c("osteopain", "gout")) {
      expect_equal(nrow(inconsistency.loops(network$data.ab, incldr = FALSE))==length(split), FALSE)
    }
    expect_equal(class(split), "nodesplit")
    expect_identical(names(split[[1]]), c("comparison",
                                          "direct", "indirect", "mbnma",
                                          "overlap matrix", "p.values", "quantiles",
                                          #"forest.plot", "density.plot",
                                          "split.model", "mbnma.model"))
    expect_equal(is.numeric(split[[1]]$p.values), TRUE)
    expect_equal(is.numeric(split[[2]]$indirect), TRUE)
    expect_equal(is.numeric(split[[3]]$direct), TRUE)
    expect_equal(length(split[[4]]$comparison), 2)
    expect_error(print(split), NA)
    expect_equal(class(summary(split)), "data.frame")


    comps.noplac <- inconsistency.loops(net.noplac$data.ab, incldr = TRUE)
    split <- mbnma.nodesplit(net.noplac, fun=dexp(), likelihood = like, link=link,
                             method="random", n.iter=1000)
    expect_equal(nrow(comps.noplac), length(split))
    expect_equal(class(split), "nodesplit")
    expect_identical(names(split[[1]]), c("comparison",
                                          "direct", "indirect", "mbnma",
                                          "overlap matrix", "p.values", "quantiles",
                                          #"forest.plot", "density.plot",
                                          "split.model", "mbnma.model"))
    expect_equal(is.numeric(split[[1]]$p.values), TRUE)
    expect_equal(is.numeric(split[[2]]$indirect), TRUE)
    expect_equal(is.numeric(split[[3]]$direct), TRUE)
    expect_equal(length(split[[4]]$comparison), 2)
    expect_error(print(split), NA)
    expect_equal(class(summary(split)), "data.frame")


    # Test comparisons
    split <- mbnma.nodesplit(net.noplac, fun=duser(fun= ~beta.1 * dose + beta.2 * (dose^2)),
                             likelihood = like, link=link,
                             method="random", n.iter=1000,
                             comparisons = rbind(comps.noplac[1,], comps.noplac[3,]))
    expect_equal(2, length(split))
    expect_error(print(split), NA)
    expect_equal(class(summary(split)), "data.frame")

    mult <- dmulti(c(list(dloglin()),
                     list(dspline("bs", knots=2)),
                     list(dspline("ns", knots=0.5)),
                     rep(list(dloglin()), length(network$agents)-3)
    ))

    # mult <- dmulti(
    #   c(rep(list(dexp()),2),
    #     rep(list(dpoly(degree=2)),1),
    #     rep(list(demax()),length(network$agents)-3)
    #   ))
    comps <- inconsistency.loops(network$data.ab, incldr = TRUE)
    split <- mbnma.nodesplit(network, fun=mult,
                             method="random", n.iter=1000,
                             comparisons = rbind(c(network$treatment[comps$t1[2]], network$treatment[comps$t2[2]])))
    expect_equal(1, length(split))
    expect_error(print(split), NA)
    expect_equal(class(summary(split)), "data.frame")

    expect_error(mbnma.nodesplit(network, fun="linear",
                                 method="random", n.iter=1000,
                                 comparisons = rbind(c("badger","rizatriptan_0.5"))),
                 "Treatment names given")

    if (datanam=="triptans") {
      expect_error(mbnma.nodesplit(network, fun="emax",
                                   method="random", n.iter=1000,
                                   comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                       c("zolmitriptan_4", "eletriptan_1"),
                                                       c("naratriptan_2", "Placebo_0"))))

    } else {
      expect_error(mbnma.nodesplit(network, fun="emax",
                                   method="random", n.iter=1000,
                                   comparisons = rbind(c("sumatriptan_0.5","rizatriptan_0.5"),
                                                       c("zolmitriptan_4", "eletriptan_1"),
                                                       c("naratriptan_2", "Placebo_0"))),
                   "Treatment names given")
    }



  })

}

