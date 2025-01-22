testthat::context("Testing get.relative")


test_that(paste("get.relative functions work correctly"), {

  skip_on_appveyor()
  skip_on_ci()
  skip_on_cran()

  n.iter <- 1000
  pD <- FALSE

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

  # alldfs <- list(triptans)
  # datanams <- c("triptans")

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

    emax <- mbnma.run(network, fun=demax(), method="random", n.iter=n.iter, pD=pD)

    emax2 <- mbnma.run(network, fun=demax(hill=0.2), method="random", n.iter=n.iter, pD=pD)

    bs <- mbnma.run(network, fun=dspline(knots=2), n.iter=n.iter, pD=pD)

    ns <- mbnma.run(network, fun=dspline(knots=c(0.5)), method="random", n.iter=n.iter, pD=pD)

    mult <- dmulti(c(list(dloglin()),
                     list(dspline("bs", knots=2)),
                     list(dspline("ns", knots=0.5)),
                     rep(list(dloglin()), length(network$agents)-3)
    ))
    multifun1 <- mbnma.run(network, fun=mult, n.iter=n.iter, pD=pD)


    mult <- dmulti(
      c(rep(list(dpoly(degree=1)),2),
        rep(list(dspline(knots = 2, type="ns", beta.1=0.2)),1),
        rep(list(dfpoly(degree=2)),length(network$agents)-3)
      ))

    multifun2 <- mbnma.run(network, fun=mult,
                           method="random", n.iter=n.iter, pD=pD)


    test_that(paste("get.relative works correctly for:", datanam), {

      expect_error(get.relative(emax, treatments=list("Placebo"=0, "Badger"=c(5,10))), "are not all named agents in")

      treatments <- list()

      temp <- get.relative(emax, treatments = treatments)
      expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      expect_equal(round(temp$mean[3,1] - temp$mean[2,1], 1), round(temp$mean[3,2], 1))
      expect_error(rank(temp), NA)

      temp <- get.relative(emax2, treatments=treatments)
      expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      expect_equal(round(temp$mean[3,1] - temp$mean[2,1], 1), round(temp$mean[3,2], 1))
      expect_error(rank(temp), NA)

      temp <- get.relative(bs, treatments=treatments)
      expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      expect_equal(round(temp$mean[3,1] - temp$mean[2,1], 1), round(temp$mean[3,2], 1))
      expect_error(rank(temp), NA)

      temp <- get.relative(ns, treatments=treatments)
      expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      expect_equal(round(temp$mean[3,1] - temp$mean[2,1], 1), round(temp$mean[3,2], 1))
      expect_error(rank(temp), NA)

      temp <- get.relative(multifun1, treatments=treatments)
      expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      if (!grepl("noplac", datanam)) {
        expect_equal(round(temp$mean[3,1] - temp$mean[2,1], 1), round(temp$mean[3,2], 1))
      }
      expect_error(rank(temp), NA)

      temp <- get.relative(multifun2, treatments=treatments)
      expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      if (!grepl("noplac", datanam)) {
        expect_equal(round(temp$mean[3,1] - temp$mean[2,1], 1), round(temp$mean[3,2], 1))
      }
      expect_error(rank(temp), NA)

      if (datanam %in% "triptans") {
        temp <- get.relative(ns, treatments=list(Placebo=0, eletriptan=1))
        expect_equal(anyNA(temp$relarray[2,1,]), FALSE)

        temp <- get.relative(ns, treatments=list(Placebo=0, zolmitriptan=10))
        expect_equal(anyNA(temp$relarray[2,1,]), FALSE)

        temp <- get.relative(multifun2, treatments=list(zolmitriptan=1, eletriptan=1))
        expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      }

      # Check prediction intervals
      # For common effects model
      temp <- get.relative(bs)
      temp2 <- get.relative(bs, lim="pred")
      expect_equal(temp$se[2], temp2$se[2])

      # For random effects model
      temp <- get.relative(ns)
      temp2 <- get.relative(ns, lim="pred")
      expect_gte(temp2$se[2], temp$se[2])

      if (datanam %in% "osteopain") {
        temp <- get.relative(bs, treatments=list("Celebrex"=c(0,100,250,400,500)))
        expect_equal(anyNA(temp$relarray[2,1,]), FALSE)
      }

      # Datasets with logit link
      if (datanam %in% c("triptans", "psoriasis90.noplac", "psoriasis75", "ssri")) {

        # Check eform
        temp <- get.relative(emax, treatments=treatments, eform=TRUE)
        expect_equal(all(temp$relarray>0, na.rm=TRUE), TRUE)

      }
    })

  }

})
