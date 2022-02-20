testthat::context("Testing plot.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
#datanam <- "psoriasis"
#dataset <- psoriasis

network <- mbnma.network(dataset)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- mbnma.network(noplac.df)

netlist <- list(network, net.noplac)


# Models
linear <- mbnma.run(mbnma.network(dataset), fun="linear", n.iter=1000)

emax <- mbnma.emax(network, emax="rel", ed50="rel", method="random", n.iter=1000)

nonparam <- mbnma.run(network, fun="nonparam.up", n.iter=1000)

emax.noplac <- mbnma.emax(net.noplac, emax="rel", ed50="rel", method="random", n.iter=1000)

resdev <- mbnma.linear(network, parameters.to.save = "resdev", n.iter=1000)

modellist <- NULL
modellist <- list(linear, emax, emax.noplac)

if ("class" %in% names(dataset)) {
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

  if (datanam!="gout") {
    expect_silent(plot(network, layout = igraph::as_star(),
                       edge.scale=1, label.distance=0))

    expect_silent(plot(network, layout = igraph::with_fr(),
                       edge.scale=1, label.distance=0))

    expect_silent(plot(network, layout = igraph::in_circle(),
                       edge.scale=0.5, label.distance=10))

    expect_silent(plot(network, layout = igraph::in_circle(),
                       level="agent", remove.loops = TRUE))
  }

  g1 <- plot(network, level="treatment")
  g2 <- plot(network, level="agent")

  expect_silent(plot(g1))
  expect_silent(plot(g2))

  expect_equal(length(igraph::V(g1))==length(igraph::V(g2)), FALSE)


  if (datanam %in% c("HF2PPITT", "psoriasis", "ssri")) {
    expect_warning(plot(net.noplac, layout=igraph::as_star(),
                        level="agent"))

    expect_warning(plot(net.noplac, layout=igraph::with_fr(),
                        level="agent", doselink = 10))

  }

  g1 <- plot(network,
            level="treatment", v.color = "agent")
  g2 <- plot(net.noplac,
             level="treatment", v.color="agent", doselink = 1)

  expect_equal("Placebo_0" %in% names(igraph::V(g1)), TRUE)
  expect_equal(length(network$treatments), length(igraph::V(g1)))
  expect_equal(length(unique(igraph::V(g1)$color)), length(network$agents))

  #expect_equal("Placebo" %in% names(igraph::V(g2)), TRUE)
  expect_equal(length(net.noplac$treatments), length(igraph::V(g2))-1)
  expect_equal(length(unique(igraph::V(g2)$color)), length(net.noplac$agents)+1)
  expect_equal(length(unique(igraph::E(g2)$color)), 2)

  expect_error(plot(network, layout=igraph::in_circle(),
                              level="class"))

})



testthat::test_that(paste0("plot.mbnma functions correctly for: ", datanam), {
  for (i in seq_along(modellist)) {
    mbnma <- modellist[[i]]
    expect_silent(plot(mbnma))
  }
  expect_equal("ggplot" %in% class(plot(nonparam)), TRUE)

  # Test number of panels is equal to number of rel effect parameters
  g <- plot(emax)
  expect_equal(length(unique(g$data$doseparam)), 2)

  # params argument
  expect_error(plot(emax, params="rabbit"))
  g <- plot(emax, params = "d.emax")
  expect_equal(length(unique(g$data$doseparam)), 1)

  # Agent labs
  expect_silent(plot(emax, agent.labs = network$agents))
  expect_error(plot(emax, agent.labs = network$agents[-3]))

  # No relative effects saved
  expect_error(plot(resdev))

  if ("class" %in% names(dataset)) {
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

  # Test disp.obs
  # expect_error(plot(pred, disp.obs = TRUE))
  expect_message(plot(pred, disp.obs = TRUE))
  # expect_error(plot(pred, disp.obs = TRUE))

  pred <- predict(emax, E0 = 0.5)
  #expect_error(plot(pred, disp.obs = TRUE))

  doses <- list()
  doses[[network$agents[2]]] <- c(0,1,2,3)
  doses[[network$agents[5]]] <- c(0.5,1,2)
  pred <- predict(emax, E0=0.1, exact.doses = doses)
  expect_silent(plot(pred, disp.obs=TRUE))

  pred <- predict(emax.noplac, E0 = 0.5)
  expect_silent(plot(pred, disp.obs = TRUE))


  # Test agent.labs
  doses <- list()
  doses[[network$agents[2]]] <- c(0,1,2,3)
  doses[[network$agents[5]]] <- c(0.5,1,2)
  pred <- predict(emax, E0=0.1, exact.doses = doses)
  g <- plot(pred, agent.labs = c("Badger", "Ferret"))
  expect_identical(levels(g$data$agent), c("Badger", "Ferret"))

  expect_error(plot(pred, agent.labs = c("Badger", "Ferret", "Whippet")))


  # Test overlay.split
  pred <- predict(linear, E0 = 0.5)
  #expect_error(plot(pred, overlay.split = TRUE))
  expect_output(plot(pred, overlay.split = TRUE))

  pred <- predict(emax, E0 = 0.5)
  expect_output(plot(pred, overlay.split = TRUE))

  doses <- list()
  doses[[network$agents[2]]] <- c(0,1,2,3)
  doses[[network$agents[5]]] <- c(0.5,1,2)
  pred <- predict(emax, E0=0.1, exact.doses = doses)
  expect_output(suppressWarnings(plot(pred, overlay.split = TRUE)))

  doses[[network$agents[2]]] <- c(1,2,3)
  doses[[network$agents[5]]] <- c(0.5,1,2)
  pred <- predict(emax, E0=0.1, exact.doses = doses)
  expect_error(plot(pred, overlay.split = TRUE), "at least one agent")

  pred <- predict(emax.noplac, E0 = 0.5)
  expect_error(plot(pred, overlay.split = TRUE))


  # Test method="common"
  pred <- predict(linear, E0 = 0.5)
  expect_output(plot(pred, overlay.split = TRUE, method="random"),
                "SD")

  pred <- predict(emax, E0 = 0.5)
  expect_output(plot(pred, overlay.split = TRUE, method="random"),
                "SD")


  # Test scales
  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)")
  expect_silent(plot(pred, scales="fixed"))
  expect_error(plot(pred, scales="badger"))

})



testthat::test_that(paste0("devplot functions correctly for: ", datanam), {
  expect_message(devplot(emax, dev.type="resdev", plot.type = "scatter", n.iter=100))

  if ("class" %in% names(dataset)) {
    expect_message(devplot(emax.class, dev.type="resdev", plot.type = "box", n.iter=100))
  }

  expect_message(devplot(emax.noplac, dev.type="resdev", plot.type = "box", n.iter=100))

  expect_message(devplot(emax.noplac, dev.type="resdev", facet = FALSE, n.iter=100))

  expect_silent(devplot(resdev, dev.type="resdev", n.iter=100))

  expect_error(devplot(emax, dev.type="dev", n.iter=100))

})



testthat::test_that(paste0("fitplot functions correctly for: ", datanam), {

  expect_message(fitplot(emax, disp.obs = TRUE, n.iter=100))

  if ("class" %in% names(dataset)) {
    expect_message(fitplot(emax.class, disp.obs=FALSE, n.iter=100))
  }

  theta.run <- mbnma.run(network, fun="linear", parameters.to.save = "theta", n.iter=1000)
  expect_silent(fitplot(theta.run, n.iter=100))

})


testthat::test_that(paste0("plot.mbnma.rank functions correctly for: ", datanam), {
  rank <- rank.mbnma(emax)
  g <- plot(rank)
  expect_equal(length(g), 2)

  if ("class" %in% names(dataset)) {
    rank <- rank.mbnma(emax.class2)
    expect_silent(plot(rank))
  }

  rank <- rank.mbnma(emax.noplac)
  g <- plot(rank, params="d.emax")
  expect_equal(length(g), 1)

  expect_silent(plot(rank, treat.labs = net.noplac$agents))
  expect_error(plot(rank, treat.labs = network$agents))

  rank <- rank.mbnma(emax.noplac, to.rank=c(2,4,5))
  expect_silent(plot(rank, treat.labs = network$agents[c(2,4,5)]))
  expect_error(plot(rank, treat.labs = net.noplac$agents), "same length as the number of ranked")

})


testthat::test_that(paste0("cumrank functions correctly for: ", datanam), {
  rank <- rank.mbnma(emax)
  g <- cumrank(rank)
  expect_equal(names(g), c("cumplot", "sucra"))

  expect_silent(cumrank(rank, params="d.emax", sucra=FALSE))
  expect_error(cumrank(rank, params="badger", sucra=FALSE))

})
