library(checkmate)
library(testthat)
library(MBNMAdose)
library(igraph)
library(dplyr)

testthat::test_check("MBNMAdose")

psoriasis.r <- psoriasis
psoriasis.r$r <- psoriasis.r$r75

datalist <- list("HF2PPITT"=HF2PPITT, "psoriasis"=psoriasis.r, "ssri"=ssri,
                 "osteopain_2wkabs"=osteopain_2wkabs, "GoutSUA_2wkCFB"=GoutSUA_2wkCFB)



datalist <- list("psoriasis"=psoriasis.r)
for (z in seq_along(datalist)) {
  print(paste("Testing dataset:", names(datalist)[z]))

  datanam <- names(datalist)[z]
  dataset <- datalist[[z]]

  # testsel <- test_env()
  # with(testsel, currentdat <- list(
  #      datanam = names(datalist)[z],
  #      dataset = datalist[[z]]
  # )
  #      )

  # testthat::test_dir("tests/testthat", filter="inconsistency", env=testsel)
  testthat::test_dir("tests/testthat")
}

