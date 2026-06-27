## Resubmission

This is a submission of an update to MBNMAdose (version 0.5.1).


## R CMD check results

There were no ERRORs, WARNINGs or NOTEs.


## Test environments

* Local Windows 11, R 4.6.0 (release) -- `R CMD check --as-cran`
* win-builder, R-devel
* win-builder, R-release
* GitHub Actions:
  - macOS-latest (release)
  - Windows Server (release)
  - Ubuntu-latest (devel, release, oldrel-1)
* R-hub: linux, windows, macos

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Notes for CRAN

* This package requires the external library JAGS (>= 4.3.0,
  https://mcmc-jags.sourceforge.net/) via the rjags / R2jags packages.
  JAGS 4.3.1 was used for the checks above. Examples and tests that fit
  JAGS models are wrapped so that they are skipped where JAGS is not
  available.
