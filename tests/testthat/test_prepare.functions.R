testthat::context("Testing prepare.functions")
df1 <- HF2PPITT

df2 <- df1
df2$agent <- as.character(df2$agent)
df2$agent[df2$dose==0] <- "almotriptan"

df.class <- HF2PPITT
df.class <- NA
df.class$class[df.class$agent %in% c("placebo", "eletriptan")] <- 1
df.class$class[is.na(df.class$class)] <- 2

datalist <- list(HF2PPITT, osteopain_2wkabs, alog_pcfb)


################### Testing ################

testthat::test_that("mbnma.validate.data functions correctly", {
  df.err <- HF2PPITT[-1,]
  expect_error(mbnma.validate.data(df.err), regexp = "single study arm")

  df.err <- HF2PPITT
  df.err$dose[10] <- -1
  expect_error(mbnma.validate.data(df.err), regexp = "All values for `dose`")

  df.err <- HF2PPITT[, names(HF2PPITT)!="r"]
  expect_error(mbnma.validate.data(df.err), regexp = "Required variable names are")

  df.err <- HF2PPITT
  df.err$r[20] <- NA
  expect_error(mbnma.validate.data(df.err), regexp = "NA values in:")

  df.err <- df.class
  df.err$class[1] <- 2
  expect_error(mbnma.validate.data(df.err), regexp = "Class codes are different")

  expect_silent(mbnma.validate.data(df.class))
})


test_that("add_index functions correctly", {

  temp.df <- df.class
  temp.df$agent <- as.character(temp.df$agent)
  df.list <- list(df.class, temp.df)

  for (df in seq_along(df.list)) {
    index <- add_index(df.list[[df]])
    expect_message(add_index(df.list[[df]]))

    expect_equal(index[["treatments"]][1], "Placebo_0")
    expect_equal(index[["agents"]][1], "Placebo")
    expect_equal(index[["classes"]][1], "Placebo")

    lvls <- c("treatments", "agents", "classes")
    lvl <- c("treatment", "agent", "class")
    for (i in seq_along(lvls)) {
      expect_equal(length(index[[lvls[i]]]), length(unique(index$data.ab[[lvl[i]]])))
      checkmate::assertNumeric(index$data.ab[[lvl[i]]], lower=1, any.missing = FALSE, finite=TRUE)
    }
  }

})




test_that("mbnma.network functions correctly", {
  expect_message(mbnma.network(df1))

  expect_message(mbnma.network(df2))

  expect_message(mbnma.network(df.class))

  df.err <- HF2PPITT[-1,]
  expect_error(mbnma.network(df.err), regex="single study arm")

  y <- 5
  expect_error(mbnma.network(y))

  expect_message(mbnma.network(alog_pcfb))
  expect_message(mbnma.network(osteopain_2wkabs))
  expect_message(mbnma.network(GoutSUA_2wkCFB))
})




test_that("mbnma.comparions functions correctly", {

  for (i in seq_along(datalist)) {
    network <- mbnma.network(datalist[[i]])

    expect_error(mbnma.comparisons(network))

    comps <- mbnma.comparisons(network$data.ab)

    expect_equal(names(comps), c("t1", "t2", "nr"))
    checkmate::assertDataFrame(comps, any.missing = FALSE, types="numeric")

    expect_equal(all(comps$t1<=comps$t2), TRUE)
  }
})




test_that("drop.disconnected functions correctly", {


  # Creating a broken network
  df.num <- mbnma.network(df1)$data.ab
  df.num$dose[df.num$studyID==3 & df.num$agent==1] <- 1
  df.num$agent[df.num$studyID==3 & df.num$agent==1] <- 5
  df.num <- df.num[!(df.num$studyID %in% c(3,11,14,16,21,29,31,37,39,40,43,44,51,63,70)),]

  fullrow <- nrow(df.num)
  network <- mbnma.network(df.num)

  expect_warning(plot(network))


  drops <- drop.disconnected(network)
  expect_equal(fullrow, nrow(drops$data.ab)+2)



  # With a complete network
  df.num <- mbnma.network(df1)$data.ab

  fullrow <- nrow(df.num)
  network <- mbnma.network(df.num)

  expect_warning(plot(network), NA)

  drops <- drop.disconnected(network)
  expect_equal(fullrow, nrow(drops$data.ab))

})
