## Resubmission
This is a resubmission. In this version I have:

* Changed the license to be standard GPL-3
* Added files to .Rbuildignore so that only standard things are in the check directory
* Added explanation of additional NOTE in R CMD check results (below)


## Test environments

* ubuntu 16.04 (on travis-ci), R 3.2.4
* local Windows, R 3.6.1 (devel and release)
* Fedora Linux, R-devel, clang, gfortran (on rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on rhub)
* Cannot test on macOS (on rhub or travis-ci) since requires system installation of JAGS (not available on CRAN)


## R CMD check results

There were 2 NOTEs:

* New submission

* Overall checktime 30 min > 10 min
  - Many functions within the package use an MCMC sampler (JAGS) which is slow to compile and run. I have set examples not to run for the purposes of checking (though I have also tested that they do run successfully), but for the tests and the vignette building I cannot reduce the time taken for these models to run so the checktime will be slower than for many other packages.


There were no ERRORs or WARNINGs

First submission of the package

R CMD check succeeded


## Downstream dependencies

There are no downstream dependencies (yet!)
