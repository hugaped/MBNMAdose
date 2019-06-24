testthat::context("Testing plot.functions")

network <- MBNMA.network(HF2PPITT)
netgout <- MBNMA.network(GoutSUA_2wkCFB)
netalog <- MBNMA.network(alog_pcfb)

datalist <- list(HF2PPITT, GoutSUA_2wkCFB, alog_pcfb)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- MBNMA.network(noplac.df)

netlist <- list(network, net.noplac)


# Models
linear <- MBNMA.run(MBNMA.network(alog_pcfb), fun="linear")

emax <- MBNMA.emax(netgout, emax="rel", ed50="rel", method="random")
emax.tript <- MBNMA.emax(network, emax="rel", ed50="rel", method="random")

emax.class <- MBNMA.emax(netclass, emax="rel", ed50="random", method="common",
                         class.effect=list(emax="random"))

emax.class2 <- MBNMA.emax(netclass, emax="rel", ed50="rel", method="common",
                         class.effect=list(emax="random"))

nonparam <- MBNMA.run(network, fun="nonparam.up")

emax.noplac <- MBNMA.emax(net.noplac, emax="rel", ed50="rel", method="random")

resdev <- MBNMA.linear(network, parameters.to.save = "resdev")

modellist <- list(linear, emax, emax.class, emax.noplac)

###################################################
################## Run Tests ######################
###################################################

test_that("plot.MBNMA.network functions correctly", {

  expect_silent(plot(network, layout_in_circle = TRUE,
                               edge.scale=1, label.distance=0))

  expect_silent(plot(network, layout_in_circle = FALSE,
                               edge.scale=1, label.distance=0))

  expect_silent(plot(network, layout_in_circle = FALSE,
                               edge.scale=0.5, label.distance=10))

  g1 <- plot(network, layout_in_circle = TRUE,
             level="treatment")
  g2 <- plot(network, layout_in_circle = TRUE,
             level="agent")

  expect_silent(plot(g1))
  expect_silent(plot(g2))

  expect_equal(length(V(g1))==length(V(g2)), FALSE)

  expect_silent(plot(network, layout_in_circle = TRUE,
                     level="agent", remove.loops = TRUE))

  expect_message(plot(net.noplac, layout_in_circle = TRUE,
                     level="agent"))

  expect_warning(plot(net.noplac, layout_in_circle = TRUE,
                      level="agent", doseparam = 5))

  g1 <- plot(network, layout_in_circle = TRUE,
            level="treatment", v.color = "agent")
  g2 <- plot(net.noplac, layout_in_circle = TRUE,
             level="treatment", v.color="agent")

  expect_equal("Placebo_0" %in% names(V(g1)), TRUE)
  expect_equal(length(network$treatments), length(V(g1)))
  expect_equal(length(unique(V(g1)$color)), length(network$agents))

  expect_equal("Placebo" %in% names(V(g2)), TRUE)
  expect_equal(length(net.noplac$treatments), length(V(g2))-1)
  expect_equal(length(unique(V(g2)$color)), length(net.noplac$agents)+1)
  expect_equal(length(unique(E(g2)$color)), 2)

  expect_error(plot(network, layout_in_circle = TRUE,
                              level="class"))

})



testthat::test_that("plot.MBNMA functions correctly", {
  for (i in seq_along(modellist)) {
    expect_silent(plot(modellist[[i]]))
  }
  expect_equal("ggplot" %in% class(plot(nonparam)), TRUE)

  # Test number of panels is equal to number of rel effect parameters
  g <- plot(emax)
  expect_equal(length(unique(g$data[[names(g$facet$params$facets)]])), 2)

  g <- plot(emax.class)
  expect_equal(length(unique(g$data[[names(g$facet$params$facets)]])), 1)

  # params argument
  expect_error(plot(emax, params="rabbit"))
  g <- plot(emax, params = "d.emax")
  expect_equal(length(unique(g$data[[names(g$facet$params$facets)]])), 1)

  # Agent labs
  expect_silent(plot(emax, agent.labs = netgout$agents))
  expect_error(plot(emax, agent.labs = netgout$agents[-3]))

  # Class labs
  expect_silent(
    plot(emax.class2, agent.labs = network$agents, class.labs=netclass$classes))

  # No relative effects saved
  expect_error(plot(resdev))

})




testthat::test_that("plot.MBNMA.predict functions correctly", {
  pred <- predict(linear, E0 = 0.5)
  expect_silent(plot(pred))

  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)")
  expect_silent(plot(pred))

  # Test disp.obs
  expect_error(plot(pred, disp.obs = TRUE))
  expect_message(plot(pred, disp.obs = TRUE, network=netgout))
  expect_error(plot(pred, disp.obs = TRUE, network=net.noplac))

  pred <- predict(emax, E0 = 0.5)
  expect_error(plot(pred, disp.obs = TRUE, network=net.noplac))

  doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  pred <- predict(emax.tript, E0=0.1, exact.doses = doses)
  expect_silent(plot(pred, disp.obs=TRUE, network=network))

  pred <- predict(emax.noplac, E0 = 0.5)
  expect_silent(plot(pred, disp.obs = TRUE, network=net.noplac))


  # Test agent.labs
  doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
  pred <- predict(emax.tript, E0=0.1, exact.doses = doses)
  g <- plot(pred, agent.labs = c("Badger", "Ferret"))
  expect_identical(levels(g$data$agent), c("Badger", "Ferret"))

  expect_error(plot(pred, agent.labs = c("Badger", "Ferret", "Whippet")))


  # Test overlay.split
  pred <- predict(linear, E0 = 0.5)
  expect_error(plot(pred, overlay.split = TRUE))
  expect_output(plot(pred, overlay.split = TRUE, network=netalog))

  pred <- predict(emax, E0 = 0.5)
  expect_output(plot(pred, overlay.split = TRUE, network=netgout))

  doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0,0.5,1,2))
  pred <- predict(emax.tript, E0=0.1, exact.doses = doses)
  expect_output(plot(pred, overlay.split = TRUE, network=network))

  doses <- list("eletriptan"=c(1,2,3), "rizatriptan"=c(0.5,1,2))
  pred <- predict(emax.tript, E0=0.1, exact.doses = doses)
  expect_error(plot(pred, overlay.split = TRUE, network=network))

  pred <- predict(emax.noplac, E0 = 0.5)
  expect_error(plot(pred, overlay.split = TRUE, network=net.noplac))


  # Test method="common"
  pred <- predict(linear, E0 = 0.5)
  expect_output(plot(pred, overlay.split = TRUE, network=netalog, method="random"),
                "SD")

  pred <- predict(emax, E0 = 0.5)
  expect_output(plot(pred, overlay.split = TRUE, network=netgout, method="random"),
                "SD")


  # Test scales
  pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)")
  expect_silent(plot(pred, scales="fixed"))
  expect_error(plot(pred, scales="badger"))

})



testthat::test_that("devplot functions correctly", {
  expect_message(devplot(emax, dev.type="resdev", plot.type = "scatter"))

  expect_message(devplot(emax.class, dev.type="resdev", plot.type = "box"))

  expect_message(devplot(emax.noplac, dev.type="resdev", plot.type = "box"))

  expect_message(devplot(emax.noplac, dev.type="resdev", facet = FALSE))

  expect_silent(devplot(resdev, dev.type="resdev"))

  expect_error(devplot(emax, dev.type="dev"))

})



testthat::test_that("fitplot functions correctly", {

  expect_message(fitplot(emax, disp.obs = TRUE))

  expect_message(fitplot(emax.class, disp.obs=FALSE))

  theta.run <- MBNMA.run(network, fun="linear", parameters.to.save = "theta")
  expect_silent(fitplot(theta.run))

})


testthat::test_that("plot.MBNMA.rank functions correctly", {
  rank <- rank.MBNMA(emax)
  g <- plot(rank)
  expect_equal(length(g), 2)

  rank <- rank.MBNMA(linear.run)
  expect_silent(plot(rank))

  rank <- rank.MBNMA(emax.noplac)
  g <- plot(rank, params="d.emax")
  expect_equal(length(g), 1)

  expect_silent(plot(rank, treat.labs = net.noplac$agents))
  expect_error(plot(rank, treat.labs = network$agents))

  rank <- rank.MBNMA(emax.noplac, to.rank=c(3,5,6))
  expect_silent(plot(rank, treat.labs = network$agents[c(3,5,6)]))
  expect_error(plot(rank, treat.labs = net.noplac$agents))

})
