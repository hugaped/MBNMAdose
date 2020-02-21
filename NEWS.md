# MBNMAdose 0.2.7



# MBNMAdose 0.2.6

## Additions/Changes

- Ensured non-parametric functions are properly monotonic by setting default initial values. Each agent now includes an index dose level of 1, which corresponds to the reference treatment effect (placebo)



# MBNMAdose 0.2.5

## Additions/Changes

- Added Lujin Li to list of authors as a reviewer thanks to his considerable help in identifying issues in the previous version.
- Each agent can be assigned a different dose-response function (by assigning a vector of functions to `fun` in `mbnma.run()`) so that multiple functions can be modelled simultaneously. Some downstream package functions still may not yet work with these models though.
- `mbnma.network` objects returned from `plot.mbnma.network` now have specific igraph attributes assigned to them, which can be easily changed by the user.
- `user.fun` now takes a formula as an argument (for example `~ (beta.1 * dose) + (beta.2 * dose^2)`) rather than a string.
- `plot.mbnma.network()` now uses a `layout` argument that takes an igraph layout function instead of `layout_in_circle` (which was a logical argument). This allows any igraph layout to be plotted rather than just a circle (e.g. `igraph::as_star()`)


## Bug fixes

- Changed `if {class(x)=="matrix"}` statements to `if {is.matrix(x)}` to address R development changes


### Major
- Exponential function models were not working previously but the dose-response function has been rewritten so that it runs the model correctly.
- DIC reported correctly in output when using plugin (`pd="plugin"`), or Kullback-Leibler divergence (`pd="pd.kl"`)
- Using the argument `parallel=TRUE` in `mbnma.run()` (or wrapper functions) now properly runs JAGS in parallel on multiple cores.

### Minor
- Downstream mbnma-related objects now store `mbnma.network` in their output rather than just treatment and agent names.


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
