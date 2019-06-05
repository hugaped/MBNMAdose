testthat::context("Testing plot.functions")

network <- MBNMA.network(HF2PPITT)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- MBNMA.network(noplac.df)

netlist <- list(network, net.noplac)


# Models
linear <- MBNMA.run(network, fun="linear")

emax <- MBNMA.emax(network, emax="rel", ed50="rel", method="random")

emax.class <- MBNMA.emax(netclass, emax="rel", ed50="random", method="common",
                         class.effect=list(emax="random"))

emax.class2 <- MBNMA.emax(netclass, emax="rel", ed50="rel", method="common",
                         class.effect=list(emax="random"))

nonparam <- MBNMA.run(network, fun="nonparam.up")

emax.noplac <- MBNMA.emax(net.noplac, emax="rel", ed50="rel", method="random")

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
  expect_silent(plot(emax, agent.labs = network$agents))
  expect_error(plot(emax, agent.labs = network$agents[-3]))

  # Class labs
  expect_silent(
    plot(emax.class2, agent.labs = network$agents, class.labs=netclass$classes))

})




testthat::test_that("plot.MBNMA.predict functions correctly", {


})



testthat::test_that("devplot functions correctly", {

})



testthat::test_that("fitplot functions correctly", {

})


testthat::test_that("plot.MBNMA.rank functions correctly", {

})
