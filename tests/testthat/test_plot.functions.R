testthat::context("Testing plot.functions")

test_that("plot functions correctly", {

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  pd <- "pv"
  n.iter <- 1000

  # Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)

  # Datasets with no placebo
  network <- mbnma.network(psoriasis90)
  psoriasis90.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

  network <- mbnma.network(ssri)
  ssri.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

  alldfs <- list(triptans, psoriasis90.noplac, osteopain, gout, psoriasis75, ssri, ssri.noplac)
  datanams <- c("triptans", "psoriasis90.noplac", "osteopain", "gout", "psoriasis75", "ssri", "ssri.noplac")

  for (dat in seq_along(alldfs)) {

      datanam <- datanams[dat]

      network <- mbnma.network(alldfs[[dat]])

      # Models
      linear <- mbnma.run(network, fun=dpoly(), n.iter=n.iter, pd=pd)

      emax <- mbnma.run(network, fun=demax(emax="rel", ed50="rel"), method="random", n.iter=n.iter, pd=pd)

      if (!grepl("noplac", datanam)) {
        nonparam <- mbnma.run(network, fun=dnonparam(direction = "increasing"), n.iter=n.iter, pd=pd)
      }

      resdev <- mbnma.run(network, fun=dpoly(), parameters.to.save = "resdev", n.iter=n.iter, pd=pd)

      ns <- mbnma.run(network, fun=dspline(knots=c(0.5)), method="random", n.iter=n.iter, pd=pd)

      mult <- dmulti(c(list(dloglin()),
                       list(dspline("bs", knots=2)),
                       list(dspline("ns", knots=0.5)),
                       rep(list(dloglin()), length(network$agents)-3)
      ))
      multifun <- mbnma.run(network, fun=mult, n.iter=n.iter, pd=pd)

      modellist <- NULL
      modellist <- list(linear, emax, ns, multifun)

      if ("class" %in% names(alldfs[[dat]])) {
        emax.class <- suppressWarnings(mbnma.emax(network, emax="rel", ed50="random", method="common",
                                                  class.effect=list(emax="random"), n.iter=1000))

        emax.class2 <- suppressWarnings(mbnma.emax(network, emax="rel", ed50="rel", method="common",
                                                   class.effect=list(emax="random"), n.iter=1000))

        modellist[[length(modellist)+1]] <- emax.class
      }




      ###################################################
      ################## Run Tests ######################
      ###################################################

      test_that(paste0("plot.mbnma.network functions correctly for:", datanam), {

        if (!(grepl("noplac", datanam) | ("gout" %in% datanam))) {
          expect_silent(plot(network, layout = igraph::as_star(),
                             edge.scale=1, label.distance=0))

          expect_silent(plot(network, layout = igraph::with_fr(),
                             edge.scale=1, label.distance=0))

          expect_silent(plot(network, layout = igraph::in_circle(),
                             edge.scale=0.5, label.distance=10))

          expect_silent(plot(network, layout = igraph::in_circle(),
                             level="agent", remove.loops = TRUE))
        } else {
          expect_warning(plot(network, layout = igraph::as_star(),
                             edge.scale=1, label.distance=0), "not connected")

          expect_warning(plot(network, layout = igraph::with_fr(),
                             edge.scale=1, label.distance=0), "not connected")

          expect_warning(plot(network, layout = igraph::in_circle(),
                             level="agent", remove.loops = TRUE), "not connected")
        }


        g1 <- suppressWarnings(plot(network, level="treatment"))
        g2 <- suppressWarnings(plot(network, level="agent"))

        expect_silent(plot(g1))
        expect_silent(plot(g2))

        expect_equal(length(igraph::V(g1))==length(igraph::V(g2)), FALSE)


        if (grepl("noplac", datanam) | "gout" %in% datanam) {
          expect_warning(plot(network, layout=igraph::as_star(),
                              level="agent"))

          expect_warning(plot(network, layout=igraph::with_fr(),
                              level="agent", doselink = 10))

          expect_message(suppressWarnings(plot(network, layout=igraph::with_fr(),
                              level="agent", doselink = 10)), "degrees of freedom")

        } else {
          g1 <- plot(network,
                     level="treatment", v.color = "agent")

          expect_equal("Placebo_0" %in% names(igraph::V(g1)), TRUE)
          expect_equal(length(network$treatments), length(igraph::V(g1)))
          expect_equal(length(unique(igraph::V(g1)$color)), length(network$agents))

          expect_error(plot(network, layout=igraph::in_circle(),
                            level="class"))
        }


      })



      testthat::test_that(paste0("plot.mbnma functions correctly for: ", datanam), {
        for (i in seq_along(modellist)) {
          mbnma <- modellist[[i]]
          expect_silent(plot(mbnma))
        }

        if (!grepl("noplac", datanam)) {
          expect_equal("ggplot" %in% class(plot(nonparam)), TRUE)
        }

        # Test number of panels is equal to number of rel effect parameters
        g <- plot(emax)
        expect_equal(length(unique(g$data$doseparam)), 2)

        expect_error(plot(multifun), NA)
        expect_error(plot(ns), NA)

        # params argument
        expect_error(plot(emax, params="rabbit"))
        g <- plot(emax, params = "emax")
        expect_equal(length(unique(g$data$doseparam)), 1)

        # No relative effects saved
        expect_error(plot(resdev), "can be identified from the model")

        if ("class" %in% datanam) {
          g <- plot(emax.class)
          expect_equal(length(unique(g$data$doseparam)), 1)

          # Class labs
          expect_silent(
            plot(emax.class2, agent.labs = netclass$agents, class.labs=netclass$classes))
        }

      })




      testthat::test_that(paste0("plot.mbnma.predict functions correctly for: ", datanam), {
        pred <- predict(linear, E0 = 0.5)
        expect_silent(plot(pred))

        pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)")
        expect_silent(plot(pred))

        pred <- predict(ns, E0 = "rbeta(n, shape1=1, shape2=5)")
        expect_silent(plot(pred))

        pred <- predict(multifun, E0 = 0.5)
        expect_silent(plot(pred))

        # Test disp.obs

        if (!grepl("noplac", datanam)) {
          expect_message(plot(pred, disp.obs = TRUE))
        } else {
          expect_silent(plot(pred, disp.obs = TRUE))

          doses <- list()
          doses[[network$agents[2]]] <- c(0,1,2,3)
          doses[[network$agents[5]]] <- c(0.5,1,2)
          pred <- predict(emax, E0=0.1, exact.doses = doses)
          expect_message(plot(pred, disp.obs=TRUE), "placebo arms")
        }

        pred <- predict(emax, E0 = 0.5)

        # Test agent.labs
        doses <- list()
        doses[[network$agents[2]]] <- c(0,1,2,3)
        doses[[network$agents[5]]] <- c(0.5,1,2)
        pred <- predict(multifun, E0=0.1, exact.doses = doses)
        g <- plot(pred, agent.labs = c("Badger", "Ferret"))
        expect_identical(levels(g$data$agent), c("Placebo", "Badger", "Ferret"))

        expect_error(plot(pred, agent.labs = c("Badger", "Ferret", "Whippet")))


        # Test overlay.split
        if (!grepl("noplac", datanam)) {
          pred <- predict(linear, E0 = 0.5)
          expect_output(plot(pred, overlay.split = TRUE))

          pred <- predict(emax, E0 = 0.5)
          expect_output(plot(pred, overlay.split = TRUE))

          doses <- list()
          doses[[network$agents[2]]] <- c(0,1,2,3)
          doses[[network$agents[5]]] <- c(0.5,1,2)
          pred <- predict(ns, E0=0.1, exact.doses = doses)
          expect_output(suppressWarnings(plot(pred, overlay.split = TRUE)))

          doses[[network$agents[2]]] <- c(1,2,3)
          doses[[network$agents[5]]] <- c(0.5,1,2)
          pred <- predict(multifun, E0=0.1, exact.doses = doses)
          expect_output(plot(pred, overlay.split = TRUE))


          # Test method="common"
          pred <- predict(ns, E0 = 0.5)
          expect_output(plot(pred, overlay.split = TRUE, method="random"),
                        "SD")

          pred <- predict(emax, E0 = 0.5)
          expect_output(plot(pred, overlay.split = TRUE, method="random"),
                        "SD")

        } else {
          pred <- predict(linear, E0 = 0.5)
          expect_error(plot(pred, overlay.split = TRUE), "Placebo required")

        }

        # Test scales
        pred <- predict(multifun, E0 = "rbeta(n, shape1=1, shape2=5)")
        expect_silent(plot(pred, scales="fixed"))
        expect_error(plot(pred, scales="badger"))

      })



      testthat::test_that(paste0("devplot functions correctly for: ", datanam), {
        expect_message(devplot(emax, dev.type="resdev", plot.type = "scatter", n.iter=100))

        if ("class" %in% datanam) {
          expect_message(devplot(emax.class, dev.type="resdev", plot.type = "box", n.iter=100))
        }

        expect_silent(devplot(resdev, dev.type="resdev", n.iter=100))

        expect_error(devplot(emax, dev.type="dev", n.iter=100))

      })



      testthat::test_that(paste0("fitplot functions correctly for: ", datanam), {

        expect_message(fitplot(emax, disp.obs = TRUE, n.iter=100))

        if ("class" %in% datanam) {
          expect_message(fitplot(emax.class, disp.obs=FALSE, n.iter=100))
        }

        theta.run <- mbnma.run(network, fun="linear", parameters.to.save = "theta", n.iter=1000)

        if (!grepl("noplac", datanam)) {
          expect_silent(fitplot(theta.run, n.iter=100))
        }

      })


      testthat::test_that(paste0("plot.mbnma.rank functions correctly for: ", datanam), {
        rank <- rank.mbnma(emax)
        g <- plot(rank)
        expect_equal(length(g), 2)

        if (grepl("noplac", datanam)) {
          expect_error(plot(rank, treat.labs = network$agents), NA)
        } else {
          expect_error(plot(rank, treat.labs = network$agents), "same length as the number of ranked")
        }

        if ("class" %in% datanam) {
          rank <- rank.mbnma(emax.class2)
          expect_silent(plot(rank))
        }

        rank <- rank.mbnma(ns)
        g <- plot(rank, params="beta.2")
        expect_equal(length(g), 1)

        rank <- rank(get.relative(multifun))
        expect_error(plot(rank), NA)

      })


      testthat::test_that(paste0("cumrank functions correctly for: ", datanam), {
        rank <- rank.mbnma(emax)
        g <- cumrank(rank)
        expect_equal(names(g), c("cumplot", "sucra"))

        expect_silent(cumrank(rank, params="emax", sucra=FALSE))
        expect_error(cumrank(rank, params="badger", sucra=FALSE))

      })

    }
})


