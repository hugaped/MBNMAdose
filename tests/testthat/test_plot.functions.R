testthat::context("Testing plot.functions")

network <- MBNMA.network(HF2PPITT)

# Generate data without placebo
noplac.df <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]
net.noplac <- MBNMA.network(noplac.df)



testthat::test_that("plot.MBNMA.network functions correctly", {

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
