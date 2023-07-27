testthat::context("Testing dose.functions")

testthat::test_that("dexp functions correctly", {
  dosefun <- dexp()
  expect_equal(dosefun$nparam, 1)
  expect_equal(dosefun$apool, c("emax"="rel"))
  expect_equal(dosefun$name, "exp")
})


testthat::test_that("demax functions correctly", {
  dosefun <- demax(emax="rel", ed50="rel")
  expect_equal(dosefun$nparam, 2)

  dosefun <- demax(emax="rel", ed50="rel", hill="rel")
  expect_equal(dosefun$nparam, 3)

  expect_message(demax(emax="rel", ed50="rel", p.expon = TRUE), "ed50")
  expect_message(demax(emax="rel", ed50="rel", hill="random", p.expon = TRUE), "hill")

  dosefun <- demax(emax="random", ed50="rel", hill="common")
  expect_equal(dosefun$apool, c(emax="random", ed50="rel", hill="common"))

  dosefun <- demax(emax="random", ed50="rel", hill=0.5)
  expect_equal(dosefun$apool, c(emax="random", ed50="rel", hill="0.5"))

  expect_error(demax(emax="random", ed50="rel", hill="test"), "hill must take either")

  expect_error(demax(emax="random", ed50="common", hill="random"), "must include at least one")
})
