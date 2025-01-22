testthat::context("Testing write.functions")


testthat::test_that("mbnma.write functions correctly", {

  write <- mbnma.write(fun=dpoly(degree=1, beta.1="rel"), method="common",
                       likelihood="binomial", link="logit")
  expect_equal(any(grepl("s\\.beta\\.1", write)), TRUE)
  expect_equal(any(grepl("s\\.beta\\.2", write)), FALSE)
  expect_equal(any(grepl("sd ", write)), FALSE)
  expect_equal(any(grepl("n\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("r\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("logit", write)), TRUE)

  write <- mbnma.write(fun=demax(emax="rel", ed50="rel"), method="random",
                       likelihood="normal", link="identity")
  expect_equal(any(grepl("s\\.beta\\.1", write)), TRUE)
  expect_equal(any(grepl("s\\.beta\\.2", write)), TRUE)
  expect_equal(any(grepl("sd ", write)), TRUE)
  expect_equal(any(grepl("y\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("se\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("logit", write)), FALSE)

  write <- mbnma.write(fun=demax(emax="rel", ed50="rel"), method="random",
                       likelihood="normal", link="identity")
  expect_equal(any(grepl("s\\.beta\\.1", write)), TRUE)
  expect_equal(any(grepl("s\\.beta\\.2", write)), TRUE)
  expect_equal(any(grepl("sd ", write)), TRUE)
  #expect_equal(any(grepl("omega", write)), TRUE)
  expect_equal(any(grepl("y\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("se\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("logit", write)), FALSE)

  write <- mbnma.write(fun=demax(emax="rel", ed50="rel", hill="random"),
                       method="random",
                       likelihood="poisson", link="cloglog")
  expect_equal(any(grepl("s\\.beta\\.1", write)), TRUE)
  expect_equal(any(grepl("s\\.beta\\.2", write)), TRUE)
  expect_equal(any(grepl("beta\\.3", write)), TRUE)
  expect_equal(any(grepl("\nd\\.3", write)), FALSE)
  expect_equal(any(grepl("sd ", write)), TRUE)
  expect_equal(any(grepl("sd\\.hill", write)), TRUE)
  expect_equal(any(grepl("sd\\.emax", write)), FALSE)
  #expect_equal(any(grepl("omega", write)), TRUE)
  expect_equal(any(grepl("r\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("E\\[i,k\\]", write)), TRUE)
  expect_equal(any(grepl("cloglog", write)), TRUE)

  write <- mbnma.write(fun=demax(emax="rel", ed50="rel", hill="random"),
                       method="random",
                       likelihood="poisson", link="cloglog",
                       omega=matrix(c(10,0,0,5), nrow=2, byrow = TRUE))
  #expect_equal(any(grepl("omega", write)), TRUE)

  expect_error(suppressWarnings(mbnma.write(fun=demax(emax="rel", ed50="rel", hill="random"),
                                            method="random",
                                            likelihood="poisson", link="cloglog",
                                            omega=matrix(c(10,0,5,5), nrow=2, byrow = TRUE))), "definite")

  write <- suppressWarnings(mbnma.write(fun=dpoly(degree=1, beta.1="rel"), method="common",
                                        likelihood="binomial", link="logit",
                                        class.effect=list(beta.1="random")))
  expect_equal(any(grepl("BETA\\.1", write)), TRUE)
  expect_equal(any(grepl("sd\\.BETA\\.1", write)), TRUE)

  expect_error(mbnma.write(fun=dpoly(degree=1, beta.1="rel"), method="common",
                           likelihood="binomial", link="logit",
                           class.effect=list(beta.2="random")))

  expect_error(suppressWarnings(mbnma.write(fun=dpoly(degree=3, beta.1="rel"), method="common",
                                            likelihood="binomial", link="logit",
                                            class.effect=list(beta.2="random"))), NA)

  write <- suppressWarnings(mbnma.write(fun=dpoly(degree=1, beta.1="rel"), method="common",
                                        likelihood="binomial", link="logit",
                                        class.effect=list(beta.1="common")))
  expect_equal(any(grepl("BETA\\.1", write)), TRUE)
  expect_equal(any(grepl("sd\\.BETA\\.1", write)), FALSE)

})
