## ----setup, include = FALSE---------------------------------------------------
library(MBNMAdose)
#devtools::load_all()
library(rmarkdown)
library(knitr)
library(dplyr)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 5,
  include=TRUE,
  tidy.opts=list(width.cutoff=80),
  tidy=TRUE
)

## ---- echo=FALSE--------------------------------------------------------------
kable(head(triptans), digits=2) 

