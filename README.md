
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Covariance Estimation for Random Surfaces beyond Separability

<!-- badges: start -->
<!-- badges: end -->

The purpose of this repo is to allow for reproducibility of the
simulation studies and applications in papers

-   Random Surface Covariance Estimation by Shifted Partial Tracing
    [\[arXiv:1912.12870v2\]](https://arxiv.org/abs/1912.12870) (aka
    separable-plus-banded model),
-   Principal Separable Component Analysis via the Partial Inner Product
    [\[arXiv:2007.12175v1\]](https://arxiv.org/abs/2007.12175) (aka
    separable component decomposition),
-   Inference and Computation for Sparsely Sampled Random Surfaces
    [\[arXiv:2103.10077v1\]](https://arxiv.org/abs/2103.10077) (aka
    sparsely observed separable model),

which constitute the doctoral thesis of the author.

While the separable-plus-banded model and the separable component
decomposition form the basis for the `surfcov` package available
[here](https://github.com/TMasak/surfcov), the package was created only
later. The simulation studies in the paper are instead produced by the
scripts included here, which contain all the functions. As such, this
repo provides the code for reproducing the results in the aforementioned
papers. A total majority of these were computed on a
[cluster](https://www.epfl.ch/research/facilities/scitas/). The
remainder of this page shows examples how to reproduce (chunks of) the
results on a single machine, running R version 4.0.4 (2021-02-15).

There may be errors in the code (especially in parts not explicitly
covered below). If a mistake is found, please contact the author at
tomas\[tecka\]masak\[zavinac\]epfl\[tecka\]ch.

## Separable-plus-banded Model

In the respective folder,
[library.R](separable_plus_banded_model/library.R) contains functions
implementing all the algorithms from the
[paper](https://arxiv.org/abs/1912.12870).

Script [demo](separable_plus_banded_model/demo.R) shows how the
functions can be used in form of a sample simulation run.

### Simulation Study

Most of the results in Section 5.1 of the
[paper](https://arxiv.org/abs/1912.12870) corresponding to the scripts
in the [separate folder](separable_plus_banded_model/scripts) were
computed on a cluster. The resulting data are stored in the [data
folder](separable_plus_banded_model/data).

To reproduce, for example, the data underlying the left panel of Figure
2 on a single machine, open
[script\_legendre\_ar\_buf.R](separable_plus_banded_model/scripts/script_legendre_ar_buf.R),
comment out the following two lines at the bottom of the script

``` r
library(parallel)
data <- mclapply(1:25,simulate,mc.cores=25)
```

and uncomment

``` r
# data <- lapply(1:25,simulate)
```

The script can then be run to reproduce
[legendre\_ar\_buf.RData](separable_plus_banded_model/data/legendre_ar_buf.RData).

The resulting data are in turn used by
[results.R](separable_plus_banded_model/results.R) to create the left
panel of Figure 2, apart from the bias curve. All the bias curves are
computed separately (and locally) using
[script\_norms.R](separable_plus_banded_model/script_norms.R)

TODO: We have apparently used the triangular kernel instead of the
Epanechnikov kernel in the respective simulations. Of course, this does
not affect the results very much (differences should not be visible by
eye), but the next version of the paper should mention the triangular
kernel instead of Epanechnikov.

### Mortality Data Analysis

The raw data were downloaded from the [Human Mortality
Database](https://www.mortality.org/) on 12/4/209/ The .txt files were
saved in a folder, labeled “death\_rates”. Then, the script
[mortality\_data\_read.R](separable_plus_banded_model/mortality_data_read.R)
was used to produce `mortality_data.RData`. As requested by the Human
Mortality Database, we do not distribute the data set, hence users have
to go through this process themselves.

After `mortality_data.RData` is created,
[mortality\_analysis.R](separable_plus_banded_model/mortality_analysis.R)
can be used to reproduce the graphical outputs of Section 5.2 of the
paper. The empirical bootstrap tests are conducted in a separate script:
[testing.R](separable_plus_banded_model/testing.R).

## Separable Component Decomposition

All the results with the exception of runtimes
([script\_fourier\_runtimes.R](separable_component_decomposition/script_fourier_times.R))
