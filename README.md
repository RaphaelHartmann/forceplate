
# forceplate <a><img src="man/figures/forceplate_logo.png" width="200" align="right" /> </a>

[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/forceplate)](https://cran.r-project.org/package=forceplate)
[![License](https://img.shields.io/badge/license-GPL(>=2)-0A7FC4.svg)](http://www.gnu.org/licenses/gpl-2.0.html)

Process raw force-plate time series data (txt-files) by segmenting them into trials and, if needed, calculating (user-defined) descriptive statistics of variables for user-defined time bins (relative to trigger onsets) for each trial. When segmenting, the data a baseline correction, a filter, and a data imputation can be applied if needed. Experimental data can also be processed and combined with the segmented force-plate data. This procedure is suggested by Johannsen et al. (2023) and some of the options (e.g., choice of low-pass filter) are also suggested by Winter (2009).

## Installation
``` r
# The easiest and most stable way is to use:
install.packages("forceplate")

# Alternatively, you might install it from GitHub via devtools:
install.packages("devtools")
devtools::install_github("RaphaelHartmann/forceplate"))
```

## References
Johannsen, L., Stephan, D. N., Straub, E., Döhring, F., Kiesel, A., Koch, I., & Müller, H. (2023). Assessing the influence of cognitive response conflict on balance control: An event-related approach using response-aligned force-plate time series data. *Psychological Research, 87*, 2297–2315.

Winter, D. A. (2009). *Biomechanics and Motor Control of Human Movement*. John Wiley & Sons.
