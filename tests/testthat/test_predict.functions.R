testthat::context("Testing predict.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
alldfs <- list(triptans, psoriasis75, ssri, osteopain, gout)
datanams <- c("triptans", "psoriasis75", "ssri", "osteopain", "gout")

# Datasets with no placebo
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

  if ("class" %in% names(dataset)) {
    netclass <- mbnma.network(df)
  }


  pd <- "pv"
  n.iter <- 1000

  testthat::test_that(paste0("predict.functions works correctly for: ", datanam), {

    skip_on_appveyor()
    skip_on_ci()
    skip_on_cran()

    #### Models ####

    linear <- mbnma.run(network, fun=dpoly(), n.iter=n.iter, pd=pd)

    emax <- mbnma.run(network, demax(), method="random", n.iter=n.iter, pd=pd)

    if ("class" %in% names(dataset)) {
      emax.class <- suppressWarnings(mbnma.run(network, demax(emax="rel", ed50="random"), method="common",
                                                class.effect=list(emax="random"), n.iter=n.iter, pd=pd))
    }

    mult <- dmulti(
      c(rep(list(dexp()),2),
        rep(list(dpoly(degree=2)),1),
        rep(list(demax()),length(network$agents)-3)
      ))

    multifun <- mbnma.run(network, fun=mult,
                          n.iter=n.iter, pd=pd)




    testthat::test_that(paste0("ref.synth functions correctly for: ", datanam), {
      ref.df <- network$data.ab[network$data.ab$agent==1,]
      ref.df <- ref.df[!duplicated(ref.df$studyID),]

      result <- suppressWarnings(ref.synth(ref.df, mbnma=emax, synth="fixed"))
      expect_equal(names(result), c("jagsmod", "m.mu"))
      expect_equal(nrow(result$m.mu), emax$BUGSoutput$n.sims)

      ref.df <- network$data.ab[network$data.ab$agent==2,]
      if (!length(unique(ref.df$studyID)) == nrow(ref.df)) {
        expect_error(ref.synth(ref.df, mbnma=linear, synth="fixed"), "contain >1 arm")
      }

      ref.df <- network$data.ab[network$data.ab$agent==network$data.ab$agent[2] & network$data.ab$dose==network$data.ab$dose[2],]
      ref.df <- ref.df[!duplicated(ref.df$studyID),]
      expect_error(suppressWarnings(ref.synth(ref.df, mbnma=linear, synth="fixed")), NA)

      expect_error(ref.synth(ref.df, mbnma=emax, synth="arndom"))

      ref.df <- network$data.ab[network$data.ab$agent==1,]
      ref.df <- ref.df[!duplicated(ref.df$studyID),]
      result <- suppressWarnings(ref.synth(ref.df, mbnma=multifun, synth="random", n.iter=1500, n.burnin=500))
      expect_identical(names(result), c("jagsmod", "m.mu", "sd.mu"))

    })



    testthat::test_that(paste0("rescale.link functions correctly for: ", datanam), {

      x <- c(-5,1,0)

      y <- rescale.link(x, direction="link", link="identity")
      expect_identical(x,y)

      y <- rescale.link(x, direction="natural", link="identity")
      expect_identical(x,y)

      expect_silent(rescale.link(x, direction="natural", link="logit"))
      expect_warning(rescale.link(x, direction="link", link="logit"))

      expect_warning(rescale.link(x, direction="link", link="probit"))
      expect_silent(rescale.link(x, direction="natural", link="probit"))

    })



    testthat::test_that(paste0("predict.mbnma functions correctly for: ", datanam), {
      #ref.df <- network$data.ab[network$data.ab$agent==1,]

      # Estimating E0
      ref.df <- network$data.ab[network$data.ab$agent==1,]
      ref.df <- ref.df[!duplicated(ref.df$studyID),]
      pred <- suppressWarnings(predict(linear, E0 = ref.df))
      expect_identical(names(pred), c("predicts", "likelihood", "link", "network", "lim", "E0"))
      expect_equal(linear$model.arg$likelihood, pred$likelihood)
      expect_equal(linear$model.arg$link, pred$link)
      expect_identical(names(pred$predicts), linear$network$agents)
      expect_silent(as.numeric(names(pred$predicts[[4]])))
      expect_equal("matrix" %in% class(pred$predicts[[4]][[4]]), TRUE)
      expect_equal(nrow(pred$predicts[[4]][[4]]), linear$BUGSoutput$n.sims)
      #expect_equal(all(pred$predicts[[2]][[2]][1] < 0), TRUE)
      expect_error(print(pred), NA)
      expect_equal(class(summary(pred)), "data.frame")

      # Stochastic E0 values
      expect_silent(predict(linear, E0 = "rnorm(n, 0.5,0.01)"))
      expect_silent(predict(linear, E0 = "rbeta(n, shape1=1, shape2=5)"))
      expect_error(predict(linear, E0 = "badgers(n, shape1=1, shape2=5)"))
      expect_error(predict(linear, E0 = "rbeta(badgers, shape1=1, shape2=5)"))

      # Determinsitic E0 values
      expect_silent(predict(linear, E0 = 0.01))
      expect_silent(predict(linear, E0 = 0.99))

      if (all(c("y", "se") %in% names(dataset))) {
        expect_warning(predict(emax, E0 = 1.5), NA)
      } else if (all(c("r", "N") %in% names(dataset))) {
        expect_warning(predict(emax, E0 = 1.5))
      }


      # Changing n.doses
      pred <- predict(linear, E0=0.5, n.doses = 10)
      expect_equal(length(pred$predicts[[2]]), 10)
      expect_error(print(pred), NA)
      expect_equal(class(summary(pred)), "data.frame")


      # Changing exact.doses
      if (length(network$agents)>=4) {
        doses <- list()
        doses[[network$agents[3]]] <- c(0,1,2,3)
        doses[[network$agents[4]]] <- c(0.5,1,2)
        pred <- predict(emax, E0=0.1, exact.doses = doses)
        expect_identical(as.numeric(names(pred$predicts[[2]])), doses[[1]])
        expect_identical(as.numeric(names(pred$predicts[[3]])), doses[[2]])

        if (all(c("r", "n") %in% names(dataset))) {
          expect_equal(all(pred$predicts[[2]][[2]][1] > 0), TRUE)
        }

        expect_error(print(pred), NA)
        expect_equal(class(summary(pred)), "data.frame")

        doses <- list(c(0,1,2,3), c(0.5,1,2))
        expect_error(predict(linear, E0=0.1, exact.doses = doses))

        dose <- c(0,0.5,1,2,4)
        doses <- list()
        for (i in seq_along(network$agents)) {
          doses[[length(doses)+1]] <- dose
        }
        expect_silent(predict(emax, E0=0.1, exact.doses = doses))

        doses <- list()
        doses[[network$agents[2]]] <- c("I","am","a","test")
        doses[[network$agents[4]]] <- c(0.5,1,2)
        expect_error(predict(emax, E0=0.1, exact.doses = doses))

        doses <- list("badger"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
        expect_error(predict(emax, E0=0.1, exact.doses = doses))
      }


      # Multiple dose-response functions
      ref.df <- network$data.ab[network$data.ab$agent==1,]
      ref.df <- ref.df[!duplicated(ref.df$studyID),]
      expect_error(suppressWarnings(predict(multifun, E0=ref.df)), NA)


      if (length(network$agents)>=4) {
        doses <- list()
        doses[[network$agents[2]]] <- c(0,2,0.1,0.05, 0.3)
        doses[[network$agents[4]]] <- c(0,2,0.1,0.5,0.9,1)
        pred <- predict(multifun, E0=0.2, exact.doses = doses)
        expect_identical(names(pred$predicts), c("Placebo", network$agents[2], network$agents[4]))

      }

    })

  })
}



