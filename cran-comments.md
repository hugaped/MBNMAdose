## Resubmission
This is a resubmission. In this version I have:

* Changed running of models in tests to ensure total build time is <10 min
* Reduced the number of models evaluated in vignette for the purposes of CRAN checks
* Added a URL to DESCRIPTION for the academic paper from which the methods have been developed
* Added explanation of additional NOTE in R CMD check results (below)


## Test environments

* ubuntu 16.04 (on travis-ci), R 3.2.4
* local Windows, R 3.6.1 (devel and release)
* Fedora Linux, R-devel, clang, gfortran (on rhub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (on rhub)
* Cannot test on macOS (on rhub or travis-ci) since requires system installation of JAGS (not available on CRAN)


## R CMD check results

There was 1 NOTEs:

* New submission


There were no ERRORs or WARNINGs

First submission of the package

R CMD check succeeded


## Downstream dependencies

There are no downstream dependencies (yet!)
