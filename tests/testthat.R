library(checkmate)
library(testthat)
library(MBNMAdose)
library(igraph)
library(dplyr)

psoriasis.r <- psoriasis
psoriasis.r$r <- psoriasis.r$r75

datalist <- list("HF2PPITT"=HF2PPITT,
                 "psoriasis"=psoriasis.r,
                 "ssri"=ssri,
                 "osteopain_2wkabs"=osteopain_2wkabs,
                 "GoutSUA_2wkCFB"=GoutSUA_2wkCFB
                 )



datalist <- list("osteopain_2wkabs"=osteopain_2wkabs)
for (dataset.number in seq_along(datalist)) {
  print(paste("Testing dataset:", names(datalist)[dataset.number]))

  datanam <- names(datalist)[dataset.number]
  dataset <- datalist[[dataset.number]]

  testthat::test_dir("tests/testthat")
}

