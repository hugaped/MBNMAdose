library(checkmate)
library(testthat)
library(MBNMAdose)
library(igraph)
library(dplyr)
library(zoo)


datalist <- list("triptans"=triptans,
                 "psoriasis75"=psoriasis75,
                 "ssri"=ssri,
                 "osteopain"=osteopain,
                 "gout"=gout
                 )

