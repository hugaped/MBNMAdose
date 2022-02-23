testthat::context("Testing rank.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)


test_that(paste("rank.functions work correctly"), {

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

    df <- alldfs[[dat]]
    dataset <- df
    datanam <- datanams[dat]

    network <- mbnma.network(df)

    # Make class data
    if ("class" %in% names(df)) {
      netclass <- mbnma.network(df)

      emax.class <- suppressWarnings(mbnma.run(netclass, demax(), method="random", n.iter=1000,
                              class.effect = list(ed50="random")))
    }

    # Models
    quad <- mbnma.run(network, fun=dpoly(degree=2, beta.1="rel", beta.2="random"), n.iter=1000)

    exponential <- mbnma.run(network, fun=dexp(onset="rel"), method="common", n.iter=1000)

    emax <- mbnma.run(network, demax(), method="random", n.iter=1000)

    if (!grepl("noplac", datanam)) {
      nonparam <- mbnma.run(network, fun=dnonparam(direction="increasing"), n.iter=1000)
    }

    spline <- mbnma.run(network, fun=dspline(type="bs", knots=c(0.1,0.8)), n.iter=1000)


    mult <- dmulti(c(list(dloglin()),
                     list(dspline("bs", knots=2)),
                     list(dspline("ns", knots=0.5)),
                     rep(list(dloglin()), length(network$agents)-3)
    ))
    multifun <- mbnma.run(network, fun=mult, n.iter=1000)



    testthat::test_that(paste0("rank.mbnma functions correctly for: ", datanam), {

      rank <- rank.mbnma(quad)
      expect_equal(names(rank), "beta.1")
      expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "lower_better"))
      expect_equal(class(rank[[1]]$summary), "data.frame")
      expect_equal("matrix" %in% class(rank[[1]]$rank.matrix), TRUE)
      expect_equal("matrix" %in% class(rank[[1]]$prob.matrix), TRUE)
      expect_error(print(rank), NA)
      expect_equal(class(summary(rank)[[1]]), "data.frame")


      rank <- rank.mbnma(emax)
      expect_equal(sort(names(rank)), sort(c("emax", "ed50")))
      expect_equal(names(rank[[1]]), c("summary", "prob.matrix", "rank.matrix", "cum.matrix", "lower_better"))
      expect_equal(class(rank[[2]]$summary), "data.frame")
      expect_equal("matrix" %in% class(rank[[1]]$rank.matrix), TRUE)
      expect_equal("matrix" %in% class(rank[[2]]$prob.matrix), TRUE)
      expect_error(print(rank), NA)
      expect_equal(class(summary(rank)[[1]]), "data.frame")

      expect_error(rank(emax, params=c("badger", "d.ed50")), "has not been monitored by the model")

      # Checking direction=1 and direction=-1 are opposites
      rank.down <- rank(emax, lower_better=TRUE)
      expect_equal(dplyr::arrange(rank.down$emax$summary, '50%')$rank.param[1] %in%
                     dplyr::arrange(rank$emax$summary, '50%')$rank.param[nrow(rank$emax$summary)-1:nrow(rank$emax$summary)],
                   TRUE)
      expect_error(print(rank.down), NA)
      expect_equal(class(summary(rank)[[1]]), "data.frame")

      to.ranks <- c(2,4)
      rank <- rank(exponential, to.rank = to.ranks)
      expect_equal(ncol(rank$emax$rank.matrix), length(to.ranks))

      if (grepl("noplac", datanam)) {
        expect_silent(rank.mbnma(exponential, to.rank = c(1,3,4)))
      } else {
        expect_warning(rank.mbnma(exponential, to.rank = c(1,3,4)), "Placebo \\(d\\[1\\] or D\\[1\\]\\) cannot be included in the ranking")
      }
      expect_silent(rank.mbnma(exponential, to.rank = c(network$agents[2], network$agents[3])))

      # Test classes
      if ("class" %in% names(dataset)) {
        expect_error(rank.mbnma(emax, level="class"), "classes have not been used")
        expect_error(rank.mbnma(emax.class, level="agent"), NA)

        rank <- rank.mbnma(emax.class, level="class")
        expect_equal(ncol(rank$ED50$rank.matrix), length(unique(dataset$class[dataset$dose>0])))
        expect_error(print(rank), NA)
        expect_equal(class(summary(rank)[[1]]), "data.frame")
      }


      if (!grepl("noplac", datanam)) {
        expect_error(rank.mbnma(nonparam), "Ranking cannot currently be performed")
      }


      # Test params
      rank <- rank.mbnma(emax)
      expect_equal(sort(names(rank)), sort(c("emax", "ed50")))
      rank <- rank.mbnma(emax, params="ed50")
      expect_equal(names(rank), c("ed50"))
      expect_error(rank.mbnma(emax, params="test"))
      expect_error(print(rank), NA)
      expect_equal(class(summary(rank)[[1]]), "data.frame")

      # With multiple-dose response functions

      expect_error(rank(multifun), "Ranking cannot be performed for models with agent-specific")

    })


    testthat::test_that(paste0("rank.mbnma.predict functions correctly for: ", datanam), {

      pred <- predict(quad, E0 = 0.5)
      rank <- rank.mbnma.predict(pred)
      expect_equal(names(rank), "Predictions")
      expect_equal(names(rank$Predictions), c("summary", "prob.matrix", "rank.matrix", "lower_better"))
      expect_equal(class(rank$Predictions$summary), "data.frame")
      expect_equal("matrix" %in% class(rank$Predictions$rank.matrix), TRUE)
      expect_equal("matrix" %in% class(rank$Predictions$prob.matrix), TRUE)


      #doses <- list("eletriptan"=c(0,1,2,3), "rizatriptan"=c(0.5,1,2))
      doses <- list()
      doses[[network$agents[2]]] <- c(0,1,2,3)
      doses[[network$agents[4]]] <- c(0.5,1,2)
      pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
                      exact.doses=doses)
      rank <- rank.mbnma.predict(pred)
      expect_equal(names(rank), "Predictions")
      expect_equal(names(rank$Predictions), c("summary", "prob.matrix", "rank.matrix", "lower_better"))
      expect_equal(class(rank$Predictions$summary), "data.frame")
      expect_equal("matrix" %in% class(rank$Predictions$rank.matrix), TRUE)
      expect_equal("matrix" %in% class(rank$Predictions$prob.matrix), TRUE)

      expect_equal(nrow(rank$Predictions$summary), length(unlist(doses)))

      # Test direction
      rank.up <- rank.mbnma.predict(pred, lower_better=TRUE)
      rank.down <- rank.mbnma.predict(pred, lower_better=FALSE)
      expect_equal(rank.down$Predictions$summary$rank.param[rank.down$Predictions$summary$`50%`==min(rank.down$Predictions$summary$`50%`)],
                   rank.up$Predictions$summary$rank.param[rank.up$Predictions$summary$`50%`==max(rank.up$Predictions$summary$`50%`)]
      )

      # Test rank.doses
      doses <- list()
      doses[[network$agents[2]]] <- c(0,1,2,3)
      doses[[network$agents[4]]] <- c(0.5,1,2)
      pred <- predict(emax, E0 = "rbeta(n, shape1=1, shape2=5)",
                      exact.doses=doses)

      doses[[network$agents[2]]] <- 2
      doses[[network$agents[4]]] <- 2
      rank <- rank.mbnma.predict(pred, rank.doses = doses)
      expect_equal(nrow(rank$Predictions$summary), 2)

      expect_error(rank.mbnma.predict(pred, rank.doses = list("badger"=2, "rizatriptan"=2)), "Agent badger not in `predicts`")

      doses[[network$agents[2]]] <- c(2, 50, 100)
      doses[[network$agents[4]]] <- 2
      expect_error(rank.mbnma.predict(pred, rank.doses = doses), "cannot be included in ranking: 50\\, 100")

    })
  }

})
















