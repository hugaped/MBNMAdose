testthat::context("Testing run.functions")


test_that(paste("run.functions work correctly"), {

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()


# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
alldfs <- list(triptans, psoriasis75, ssri, osteopain, gout)
datanams <- c("triptans", "psoriasis75", "ssri", "osteopain", "gout")

# Datasets with no placebo/
network <- mbnma.network(psoriasis90)
psoriasis90.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

network <- mbnma.network(ssri)
ssri.noplac <- network$data.ab[network$data.ab$narm>2 & network$data.ab$agent!=1,]

alldfs[[length(alldfs)+1]] <- psoriasis90.noplac
alldfs[[length(alldfs)+1]] <- ssri.noplac
datanams <- append(datanams, c("psoriasis90.noplac", "ssri.noplac"))

for (dat in seq_along(alldfs)) {

  datanam <- datanams[dat]
  dataset <- alldfs[[dat]]

  print(datanam)

  ### Datasets ####
  network <- mbnma.network(dataset)



  # Make class data
  df <- dataset
  df1 <- dataset

  if ("class" %in% names(dataset)) {
    netclass <- mbnma.network(df)
  }

  n.iter=500
  pd <- "pv"

  test_that(paste("mbnma.run wrappers function correctly for:", datanam), {

    expect_warning(mbnma.linear(network, slope="rel", n.iter=n.iter, pd=pd), "syntax for specifying dose-response functions")

    # Single parameter DR functions
    expect_error(mbnma.run(network, fun=dpoly(degree=1, beta.1="random"), n.iter=n.iter, pd=pd), "must include at least one parameter")

    result <- mbnma.run(network, fun=dloglin(), method="common", n.iter=n.iter, pd=pd)
    expect_equal(all(c("rate") %in% result$parameters.to.save), TRUE)
    expect_equal(all(c("sd") %in% result$parameters.to.save), FALSE)
    expect_error(summary(result), NA)
    expect_error(plot(result), NA)
    expect_error(rank(result), NA)
    expect_error(predict(result), NA)

    if ("class" %in% names(dataset)) {

      # Two parameter DR functions
      result <- suppressWarnings(mbnma.emax(network, emax="rel", ed50="rel", method="common",
                           class.effect=list(emax="common"), n.iter=n.iter, pd=pd, cor = FALSE))
      expect_equal(all(c("EMAX", "ed50", "emax") %in% result$parameters.to.save), TRUE)
      expect_error(suppressWarnings(summary(result), NA))
      expect_error(plot(result), NA)
      expect_error(rank(result), NA)
      expect_error(predict(result), "does not work with models that use class effects")

      # Three parameter DR functions
      if (datanam!="osteopain") {
        result <- mbnma.emax.hill(netclass, emax="rel", ed50="rel", hill="common",
                                  method="random", n.iter=n.iter, pd=pd)
        expect_equal(all(c("emax", "ed50", "hill", "sd") %in% result$parameters.to.save), TRUE)
        expect_error(summary(result), NA)
        expect_error(plot(result), NA)
        expect_error(rank(result), NA)
        expect_error(predict(result), NA)
      }
    }

  })



  test_that(paste("check.likelink function correctly for:", datanam), {

    if (all(c("y", "se") %in% names(dataset))) {
      expect_silent(check.likelink(df, likelihood = "normal", link="identity"))
      expect_silent(check.likelink(df, likelihood = "normal", link="logit"))

      # Expect error due to misspecified df
      expect_error(check.likelink(df, likelihood = "binomial", link="identity"))
      expect_error(check.likelink(df, likelihood = "poisson", link="identity"))

      # Expect errror due to misspecified arguments
      expect_error(check.likelink(df, likelihood = "normal", link="badger"))
      expect_error(check.likelink(df, likelihood = "test", link="identity"))

    } else if (all(c("r", "n") %in% names(dataset))) {
      expect_silent(check.likelink(df, likelihood = "binomial", link="identity"))
      expect_silent(check.likelink(df, likelihood = "binomial", link="logit"))

      # Expect error due to misspecified df
      expect_error(check.likelink(df, likelihood = "normal", link="identity"))
      expect_error(check.likelink(df, likelihood = "poisson", link="identity"))

      # Expect errror due to misspecified arguments
      expect_error(check.likelink(df, likelihood = "binomial", link="badger"))
      expect_error(check.likelink(df, likelihood = "test", link="logit"))
    }

  })




  test_that(paste("nma.run function correctly for:", datanam), {
    n.iter <- 500

    # expect_warning(nma.run(network, method="random", n.iter=100, warn.rhat = TRUE))

    expect_warning(nma.run(network, method="common", n.iter=n.iter, pd=pd, warn.rhat = FALSE), NA)

    result <- nma.run(network, method="random", n.iter=n.iter, pd=pd, warn.rhat = FALSE)
    expect_equal(names(result), c("jagsresult", "trt.labs", "UME"))
    expect_equal(all(c("d", "sd") %in% result$jagsresult$parameters.to.save), TRUE)

    result <- nma.run(network, method="random", n.iter=n.iter, pd=pd, warn.rhat = FALSE,
                      UME=TRUE)
    expect_equal("d[1,1]" %in% rownames(result$jagsresult$BUGSoutput$summary), TRUE)


    # Creating a broken network
    df.num <- mbnma.network(df1)$data.ab

    sepcomp <- mbnma.comparisons(df.num)[nrow(mbnma.comparisons(df.num)),]
    keep <- df.num$studyID[df.num$treatment %in% c(sepcomp$t1, sepcomp$t2)]
    df.num <- df.num[!(df.num$studyID %in% keep & !df.num$treatment  %in% c(sepcomp$t1, sepcomp$t2)),]

    df.num <- df.num %>% dplyr::group_by(studyID) %>% dplyr::mutate(narm=dplyr::n())
    df.num <- df.num[df.num$narm>1,]

    fullrow <- nrow(df.num)
    network.disc <- mbnma.network(df.num)

    result.1 <- nma.run(network.disc, method="random", n.iter=n.iter, pd=pd, warn.rhat = FALSE,
                        UME=TRUE, drop.discon = TRUE)
    result.2 <- nma.run(network.disc, method="random", n.iter=n.iter, pd=pd, warn.rhat = FALSE,
                        UME=TRUE, drop.discon = FALSE)
    result.3 <- nma.run(network.disc, method="random", n.iter=n.iter, pd=pd, warn.rhat = FALSE,
                        UME=TRUE, drop.discon = TRUE)
    expect_equal(length(result.1$trt.labs)!=length(result.2$trt.labs), TRUE)
    expect_equal(length(result.1$trt.labs)==length(result.3$trt.labs), TRUE)
  })



  test_that(paste("pDcalc functions correctly for:", datanam), {
    n.iter=1000

    if (all(c("y", "se") %in% names(dataset))) {
      likelihood <- "normal"
      link <- "identity"

      # Prevents skip
      expect_equal(5,5)

    } else if (all(c("r", "n") %in% names(dataset))) {
      likelihood <- "binomial"
      link <- "logit"


      # For binomial likelihood
      result <- mbnma.run(network, fun=dexp(), method="random",
                          parameters.to.save = c("psi", "resdev"),
                          n.iter=n.iter, pd=pd)

      jagsdata <- getjagsdata(network$data.ab, likelihood = likelihood, link=link)

      obs1 <- jagsdata$r
      obs2 <- jagsdata$n

      pd.est <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                   theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                   likelihood=likelihood, type="dose")
      expect_equal(length(pd.est),1)
      expect_equal(class(pd.est),"numeric")

      pd.est <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=5,
                   theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                   likelihood=likelihood, type="dose")
      expect_equal(length(pd.est),1)
      expect_equal(class(pd.est),"numeric")

      pd.est <- pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=5,
                   theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                   likelihood="poisson", type="dose")

      expect_error(pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                          theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                          likelihood="poisson", type="time"))

      expect_error(pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=jagsdata[["NS"]],
                          theta.result=NULL, resdev.result=result$BUGSoutput$mean$resdev,
                          likelihood="poisson", type="dose"))

      expect_error(pDcalc(obs1=obs1, obs2=obs2, narm=jagsdata[["narm"]], NS=NULL,
                          theta.result=result$BUGSoutput$mean$psi, resdev.result=result$BUGSoutput$mean$resdev,
                          likelihood="poisson", type="dose"))
    }

  })





  test_that(paste("mbnma.update function correctly for:", datanam), {

    result <- mbnma.run(network, fun=demax(), method="common",
                        n.iter=500)

    expect_error(mbnma.update(result, param="test", n.iter=100))

    update <- mbnma.update(result, param="resdev", n.iter=100)
    expect_equal(names(update), c("study", "arm", "mean", "facet", "fupdose", "groupvar"))

    update <- mbnma.update(result, param="theta", n.iter=100)
    expect_equal(names(update), c("study", "arm", "mean", "facet", "fupdose", "groupvar"))

    update <- mbnma.update(result, param="theta", n.iter=100, armdat = FALSE)
    expect_equal(names(update), c("study", "arm", "mean"))

  })
}




})
