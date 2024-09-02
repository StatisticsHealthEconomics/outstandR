# outstandR: *Out*come regression *stand*ardisation
<!-- <img align="right" src="mime.png" width="100"> -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/StatisticsHealthEconomics/outstandR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/StatisticsHealthEconomics/outstandR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

> Model-based Standardisation with G-estimation

## Overview

`outstandR` is an R package designed to facilitate **outcome regression standardisation** using model-based approaches, particularly focusing on G-estimation. The package provides tools to apply standardisation techniques for indirect treatment comparisons, especially in scenarios with limited individual patient data.

Key features include:
- **G-computation methods**
- **Multiple imputation marginalisation**
- **MAIC, STC** for comparison

This R package contains code originally written for the papers:

> Antonio Remiro-Azócar, Anna Heath, Gianluca Baio. _Parametric G-computation for Compatible Indirect Treatment Comparisons with Limited Individual Patient Data._
> Res Synth Methods. 2022;1–31.

and

> Remiro-Azócar, A., Heath, A., & Baio, G. (2023). _Model-based standardization using multiple imputation. BMC Medical Research Methodology_, 1–15. https://doi.org/10.1186/s12874-024-02157-x

## Installation
Install the [development version from GitHub](https://github.com/StatisticsHealthEconomics/) using `remotes`:

```r
remotes::install_github("StatisticsHealthEconomics/outstandR")
```

## License
This package is licensed under the GPLv3. For more information, see [LICENSE](https://www.gnu.org/licenses/gpl-3.0).

## Contributing
We welcome contributions! Please submit contributions through `Pull Requests`, following the [contributing guidelines](https://github.com/n8thangreen/BCEA/blob/dev/CONTRIBUTING.md).
To report issues and/or seek support, please file a new ticket in the
[issue](https://github.com/StatisticsHealthEconomics/outstandR/issues) tracker.

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/n8thangreen/BCEA/blob/dev/CONDUCT.md). By participating in this project you agree to abide by its terms.
