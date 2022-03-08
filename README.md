Using Helicopter Counts to Estimate the Abundance of Himalayan Tahr in
New Zealand’s Southern Alps
================

## Overview

This repository contains data and code from:

Ramsey, D.S.L., Forsyth, D.M., Perry, M., Thomas, P., McKay, M., and
Wright, E. (2022). “Using Helicopter Counts to Estimate the Abundance of
Himalayan Tahr in New Zealand’s Southern Alps” *Journal of Wildlife
Management*

### File descriptions:

-   `r/tahr_models.r` reads data and fits dynamic N-mixture models to
    helicopter counts from 117 plots using `nimble`.
-   `r/misc_functions.r` contains various functions require by the main
    script.
-   `data/tahr_data.rds` helicopter count data for 117 plots by 3
    occasions including metadata.
-   `data/plots.rds` sf object of plot boundaries.
-   `data/MU.rds` sf object of tahr management unit boundaries.
-   `data/PCL.rds` sf object of public conservaton land boundaries.
-   `data/nzsouth.rds` sf object of the south island of New Zealand.

## Prerequisites

The script require packages `tidyverse`, `nimble`, `sf`, `nngeo`,
`MCMCviz`, `bayesplot`, `gridExtra`,`progress`,`ggspatial` and
`ggvoronoi`.
