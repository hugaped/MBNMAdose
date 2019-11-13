# MBNMAdose 0.2.5

## Documentation changes

- Added Lujin Li to list of authors as a reviewer thanks to his considerable help in identifying issues in the previous version.

## Additions

- `mbnma.network` objects returned from `plot.mbnma.network` now have specific igraph attributes assigned to them, which can be easily changed by the user.

## Bug fixes

### Major
- Exponential function models were not working previously but the dose-response function has been rewritten so that it runs the model correctly.
- DIC reported correctly in output when using plugin (`pd="plugin"`), or Kullback-Leibler divergence (`pd="pd.kl"`) 


# MBNMAdose 0.2.4

## Documentation changes

- Update to README to ensure package workflow image works correctly on CRAN.
- Added JAGS as a System Requirement (JAGS >= 4.3.0) to the DESCRIPTION

## Bug fixes

### Major
- Fixed incorrect ordering of treatment codes in mbnma.network, which also led to problems in subsequent commands/plots

### Minor
- Fixed minor with numeric vs character coding of treatments in arguments for ranking functions
- Fixed issue with `nma.nodesplit()` that prevented the model running if disconnected treatments were included in the analysis (`drop.discon=FALSE`)


# MBNMAdose 0.2.3

## Documentation changes

- Made corrections to some arguments specified in documentation
- Fixed incorrect vignette reference
- Allowed more examples in the vignette to run so that plots can be created which better illustrates how MBNMA.run() functions are used.


# MBNMAdose 0.2.2

## First release of package

Welcome to MBNMAdose. Ready for release into the world. I hope it can be of service to you! For time-course MBNMA, also check out the sister package, MBNMAtime.
