testthat::context("Testing plot.functions")

network <- MBNMA.network(HF2PPITT)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- MBNMA.network(noplac.df)

netlist <- list(network, net.noplac)

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

})




testthat::test_that("plot.MBNMA.predict functions correctly", {

})



testthat::test_that("devplot functions correctly", {

})



testthat::test_that("fitplot functions correctly", {

})


testthat::test_that("plot.MBNMA.rank functions correctly", {

})
