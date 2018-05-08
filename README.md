# TSexploreR

Package containing miscellaneous tools for Time Series (TS) analysis. A main emphasis is given on peak and oscillation analysis. It contains convenient wrappers for existing implementation along with re-implementation of existing techniques and completely novel methods.

The package is still under development and was mainly developed on Ubuntu 16.04LTS, R 3.4.4 and Python 3.6.4.

## Installation

Installation can be performed directly from github.

```R
library(devtools)
install_github("majpark21/TSexploreR")
```

Note that some functions rely on Python, specifically numpy and scipy. Therefore make sure that these modules are part of your python installation. TSexploreR uses reticulate package to call Python functions: https://rstudio.github.io/reticulate/. If you encounter python-related errors, please follow the documentation to make sure that you have correctly configured the reticulate package.

## Overview of the contents of the package

### Peaks and oscillations analysis
* Peak extraction: 

Can be performed using `extract_peak()`. It is a wrapper for 5 different methods using rolling windows approaches and wavelet tranform.

* Single-Peak features:

Simple geometric features such as height, width, growth rate... can be extracted for single peaks with: `FeatAllFeat()`. All individual single-peak features can be obtain by the call of the corresponding `Feat*()` function.


### Oscillation analysis
* Trends and season decomposition

`classical.decomposition()` provides a decomposition between trends and season using an additive model: $y(t) = trend(t) + seasonal(t) + remainder(t)$. The decomposition can be performed robustly, such that the seasonal component is estimated by using the median (as opposed to the mean) of the detrended TS over seasons.

* Synchrony of a pair of TS

`synchrony.measures()` returns 4 measures of synchrony of oscillating trajectories: 3 types of correlations (Pearson, Kendall, Spearman) and "overlap.clip", a novel approach. The "overlap.clip" method relies on the prior transformation of the oscillating signals into binary "up-down" representations through clipping. The 2 transformed signals are then overlapped and the proportion of points were they are identical is reported as a measure of synchrony.

* Simulation of noisy populations of oscillating signals

Some convenience functions for generating sinusoids and exponential decays are regrouped in the `sim_*()` functions. Each function generates a population with a certain type of noise between individual TS: phase-shift, Gaussian white noise, response lag...



* Noise in a set of time series
* Pairwise Separability of an ensemble of time series sets
* Utils
