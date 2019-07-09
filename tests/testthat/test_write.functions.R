testthat::context("Testing write.functions")


testthat::test_that("mbnma.write functions correctly", {

  write <- mbnma.write(fun="linear", beta.1="rel", method="common",
                       likelihood="binomial", link="logit")
  expect_equal(grepl("s\\.beta\\.1", write), TRUE)
  expect_equal(grepl("s\\.beta\\.2", write), FALSE)
  expect_equal(grepl("sd ", write), FALSE)
  expect_equal(grepl("N\\[i,k\\]", write), TRUE)
  expect_equal(grepl("r\\[i,k\\]", write), TRUE)
  expect_equal(grepl("logit", write), TRUE)

  write <- mbnma.write(fun="emax", beta.1="rel", beta.2="rel", method="random",
                       likelihood="normal", link="identity")
  expect_equal(grepl("s\\.beta\\.1", write), TRUE)
  expect_equal(grepl("s\\.beta\\.2", write), TRUE)
  expect_equal(grepl("sd ", write), TRUE)
  expect_equal(grepl("y\\[i,k\\]", write), TRUE)
  expect_equal(grepl("se\\[i,k\\]", write), TRUE)
  expect_equal(grepl("logit", write), FALSE)

  write <- mbnma.write(fun="emax", beta.1="rel", beta.2="rel", method="random",
                       likelihood="normal", link="identity")
  expect_equal(grepl("s\\.beta\\.1", write), TRUE)
  expect_equal(grepl("s\\.beta\\.2", write), TRUE)
  expect_equal(grepl("sd ", write), TRUE)
  expect_equal(grepl("Omega", write), TRUE)
  expect_equal(grepl("y\\[i,k\\]", write), TRUE)
  expect_equal(grepl("se\\[i,k\\]", write), TRUE)
  expect_equal(grepl("logit", write), FALSE)

  write <- mbnma.write(fun="emax.hill", beta.1="rel", beta.2="rel", beta.3="random",
                       method="random",
                       likelihood="poisson", link="cloglog")
  expect_equal(grepl("s\\.beta\\.1", write), TRUE)
  expect_equal(grepl("s\\.beta\\.2", write), TRUE)
  expect_equal(grepl("beta\\.3,", write), TRUE)
  expect_equal(grepl("\nd\\.3,", write), FALSE)
  expect_equal(grepl("sd ", write), TRUE)
  expect_equal(grepl("sd\\.3", write), TRUE)
  expect_equal(grepl("sd\\.1", write), FALSE)
  expect_equal(grepl("Omega", write), TRUE)
  expect_equal(grepl("r\\[i,k\\]", write), TRUE)
  expect_equal(grepl("E\\[i,k\\]", write), TRUE)
  expect_equal(grepl("cloglog", write), TRUE)

  write <- mbnma.write(fun="emax.hill", beta.1="rel", beta.2="rel", beta.3="random",
                       method="random",
                       likelihood="poisson", link="cloglog",
                       var.scale=c(10,1))
  expect_equal(grepl("Omega\\[1,1\\] <- 10", write), TRUE)
  expect_equal(grepl("Omega\\[2,2\\] <- 1", write), TRUE)

  expect_error(mbnma.write(fun="emax.hill", beta.1="rel", beta.2="rel", beta.3="random",
                       method="random",
                       likelihood="poisson", link="cloglog",
                       var.scale=100))

  write <- mbnma.write(fun="linear", beta.1="rel", method="common",
                       likelihood="binomial", link="logit",
                       class.effect=list(beta.1="random"))
  expect_equal(grepl("D\\.1", write), TRUE)
  expect_equal(grepl("sd\\.D\\.1", write), TRUE)

  expect_error(mbnma.write(fun="linear", beta.1="rel", method="common",
                       likelihood="binomial", link="logit",
                       class.effect=list(beta.2="random")))

  write <- mbnma.write(fun="linear", beta.1="rel", method="common",
                       likelihood="binomial", link="logit",
                       class.effect=list(beta.1="common"))
  expect_equal(grepl("D\\.1", write), TRUE)
  expect_equal(grepl("sd\\.D\\.1", write), FALSE)

})
