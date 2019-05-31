testthat::context("Testing run.functions")

network <- MBNMA.network(HF2PPITT)

# Make class data
df <- HF2PPITT
df$class <- ifelse(df$agent=="placebo", "placecbo", "active")
netclass <- MBNMA.network(df)


test_that("MBNMA.run functions correctly", {
  n.iter=500

  # Single parameter DR functions
  result <- MBNMA.run(network, fun="linear", beta.1="rel", method="common",
                      pd="pv", n.iter=n.iter)
  expect_equal(class(result), c("MBNMA", "rjags"))
  expect_equal("d.1" %in% result$parameters.to.save, TRUE)

  result <- MBNMA.run(network, fun="exponential", beta.1="rel", method="random",
                      pd="pd.kl", n.iter=n.iter)
  expect_equal(class(result), c("MBNMA", "rjags"))
  expect_equal("sd" %in% result$parameters.to.save, TRUE)

  result <- MBNMA.run(netclass, fun="exponential", beta.1="rel", method="common",
                      pd="popt", class.effect = list(beta.1="random"), n.iter=n.iter)
  expect_equal(class(result), c("MBNMA", "rjags"))
  expect_equal("D.1" %in% result$parameters.to.save, TRUE)
  expect_equal("sd.D.1" %in% result$parameters.to.save, TRUE)

  result <- MBNMA.run(network, fun="nonparam.up", method="common", n.iter=n.iter)
  expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)

  result <- MBNMA.run(network, fun="nonparam.down", method="random", n.iter=n.iter)
  expect_equal("d.1[1,1]" %in% rownames(result$BUGSoutput$summary), TRUE)
  expect_equal("sd" %in% result$parameters.to.save, TRUE)



  expect_error(result <- MBNMA.run(network, fun="linear", beta.1="rel", method="fixed", n.iter=n.iter))

})
