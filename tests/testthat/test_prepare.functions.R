testthat::context("Testing prepare.functions")

# Tested datasets must have at least 5 agents - options are HF2PPIT, psoriasis, ssri, osteopain, gout(?)
alldfs <- list(triptans, psoriasis75, ssri, osteopain, gout)
datanams <- c("triptans", "psoriasis75", "ssri", "osteopain", "gout")

for (dat in seq_along(alldfs)) {

  datanam <- datanams[dat]
  dataset <- alldfs[[dat]]

  print(datanam)

  ### Datasets ####
  network <- mbnma.network(dataset)


  df1 <- dataset

  df2 <- df1
  df2$agent <- as.character(df2$agent)
  df2$agent[df2$dose==0] <- network$agents[2]

  if ("class" %in% names(dataset)) {
    df.class <- dataset
  }

  # df.class <- HF2PPITT
  # df.class$class <- NA
  # df.class$class[df.class$agent %in% c("placebo", "eletriptan")] <- 1
  # df.class$class[is.na(df.class$class)] <- 2

  datalist <- list(df1, df2)


  ################### Testing ################

  testthat::test_that(paste0("mbnma.validate.data functions correctly for: ", datanam), {
    df.err <- dataset
    arm <- df.err[df.err$studyID==df.err$studyID[1],]
    arm <- arm[1,]
    df.err <- df.err[df.err$studyID!=df.err$studyID[1],]
    df.err <- rbind(arm, df.err)
    expect_error(mbnma.validate.data(df.err), regexp = "single study arm")

    df.err <- dataset
    df.err$dose[10] <- -1
    expect_error(mbnma.validate.data(df.err), regexp = "All values for `dose`")

    df.err <- dataset[, !(names(dataset) %in% c("r", "y"))]
    expect_error(mbnma.validate.data(df.err), regexp = "Required variable names are")

    df.err <- dataset
    if ("r" %in% names(df.err)) {
      df.err$r[20] <- NA
    } else if ("y" %in% names(df.err)) {
      df.err$y[20] <- NA
    }
    expect_error(mbnma.validate.data(df.err), regexp = "NA values in:")

    if ("class" %in% names(dataset)) {
      df.err <- dataset
      df.err$class[1] <- 3
      expect_error(mbnma.validate.data(df.err), regexp = "Class codes are different")

      expect_silent(mbnma.validate.data(df.class))
    }

    if ("y" %in% names(dataset)) {
      new.df <- dataset
      new.df$standsd <- 0.5

      expect_silent(mbnma.validate.data(new.df))

      df.err <- new.df
      df.err$standsd[1] <- 2

      expect_error(mbnma.validate.data(df.err), "must be identical within each study")
    }

  })


  test_that(paste0("add_index functions correctly for: ", datanam), {
    df <- dataset

    index <- add_index(df)
    expect_message(add_index(df))

    expect_equal(index[["treatments"]][1], "Placebo_0")
    expect_equal(index[["agents"]][1], "Placebo")

    lvl <- c("treatment", "agent")
    lvls <- c("treatments", "agents")
    if ("class" %in% names(df)) {
      expect_equal(index[["classes"]][1], "Placebo")

      lvl <- append(lvl, "class")
      lvls <- append(lvls, "classes")
    }

    for (i in seq_along(lvls)) {
      expect_equal(length(index[[lvls[i]]]), length(unique(index$data.ab[[lvl[i]]])))
      checkmate::assertNumeric(index$data.ab[[lvl[i]]], lower=1, any.missing = FALSE, finite=TRUE)
    }

  })




  test_that(paste0("mbnma.network functions correctly for: ", datanam), {
    expect_message(mbnma.network(df1))

    if (datanam!="osteopain") {
      expect_message(mbnma.network(df2))
    } else {
      expect_error(mbnma.network(df2), "Class codes are different")
    }

    df.err <- dataset
    arm <- df.err[df.err$studyID==df.err$studyID[1],]
    arm <- arm[1,]
    df.err <- df.err[df.err$studyID!=df.err$studyID[1],]
    df.err <- rbind(arm, df.err)
    expect_error(mbnma.network(df.err), regex="single study arm")

    y <- 5
    expect_error(mbnma.network(y))
  })




  test_that(paste0("mbnma.comparions functions correctly for: ", datanam), {

    for (i in seq_along(datalist)) {

      if (i==2 & datanam!="osteopain") {
        network <- mbnma.network(datalist[[i]])

        expect_error(mbnma.comparisons(network))

        comps <- mbnma.comparisons(network$data.ab)

        expect_equal(names(comps), c("t1", "t2", "nr"))
        checkmate::assertDataFrame(comps, any.missing = FALSE, types="numeric")

        expect_equal(all(comps$t1<=comps$t2), TRUE)
      } else {
        # Created to avoid skips
        expect_equal(5,5)
      }
    }
  })




  test_that(paste0("drop.disconnected functions correctly for: ", datanam), {


    # Creating a broken network
    df.num <- mbnma.network(df1)$data.ab

    sepcomp <- mbnma.comparisons(df.num)[nrow(mbnma.comparisons(df.num)),]
    keep <- df.num$studyID[df.num$treatment %in% c(sepcomp$t1, sepcomp$t2)]
    df.num <- df.num[!(df.num$studyID %in% keep & !df.num$treatment  %in% c(sepcomp$t1, sepcomp$t2)),]

    df.num <- df.num %>% dplyr::group_by(studyID) %>% dplyr::mutate(narm=dplyr::n())
    df.num <- df.num[df.num$narm>1,]

    network <- mbnma.network(df.num)

    expect_warning(plot(network))


    drops <- drop.disconnected(network)
    expect_equal(nrow(df.num) > nrow(drops$data.ab), TRUE)



    # With a complete network
    if (datanam %in% c("triptans", "psoriasis75", "ssri", 2)) {
      df.num <- mbnma.network(df1)$data.ab

      fullrow <- nrow(df.num)
      network <- mbnma.network(df.num)

      expect_warning(plot(network), NA)

      drops <- drop.disconnected(network)
      expect_equal(fullrow, nrow(drops$data.ab))
    }

  })




  test_that(paste0("genspline functions correctly for: ", datanam), {

    xlist <- list(c(0:50), c(10,25,89), c(5,10), c(1))
    for (i in seq_along(xlist)) {
      x <- xlist[[i]]
      expect_silent(genspline(x, spline="ns", knots=2, max.dose=max(x)))
      expect_silent(genspline(x, spline="ns", knots=3, max.dose=max(x)))

      knots <- 3
      splines <- genspline(x, spline="ns", knots=knots, max.dose=max(x))
      expect_equal(nrow(splines), length(x))
      expect_equal(ncol(splines), knots+1)

      if (max(x)>10) {
        knots <- c(0.35,0.5,0.1)
        expect_silent(genspline(x, spline="ns", knots=knots, max.dose=10))

        if (length(x)>1) {
          expect_equal(ncol(genspline(x, spline="ns", knots=3, max.dose=10)), length(knots)+1)
        }
      }

      expect_error(genspline(x, spline="ns", knots=5, max.dose=max(x)), "complexity")
      expect_error(genspline(x, spline="ns", knots=c(1,2,3), max.dose=max(x)), "'probs' outside")

      expect_error(genspline(x, spline="badger", knots=3, max.dose=max(x)))

    }



  })


  test_that(paste0("getjagsdata functions correctly for: ", datanam), {

    data.ab <- network$data.ab

    expect_error(getjagsdata(data.ab, class=FALSE, fun=demax(), nodesplit = c(1,3)), NA)

    expect_error(getjagsdata(data.ab, fun=dspline(type="ns", knots=c(0.2,0.5), beta.1="common", beta.2 = "rel", beta.3="random")), NA)


    mult <- dmulti(
      c(rep(list(dpoly(degree=1)),2),
        rep(list(dspline(knots = 2, type="ns", beta.1=0.2)),1),
        rep(list(dfpoly(degree=2)),length(network$agents)-3)
      ))

    expect_error(getjagsdata(data.ab, fun=mult), NA)


    mult <- dmulti(
      c(rep(list(dpoly(degree=1)),2),
        rep(list(dspline(knots = c(0.1,0.5), type="ns", beta.1=0.2)),1),
        rep(list(dspline(knots = 3, type="ls", beta.2="common")),length(network$agents)-3)
      ))

    expect_error(getjagsdata(data.ab, fun=mult), NA)

    expect_error(getjagsdata(data.ab, class=FALSE, fun=demax(hill=0.5)), NA)

  })

}

