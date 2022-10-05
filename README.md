# bpa-jfdi-p2-publication

This repository contains scripts and documentation associated with the publication:

> Hernandez, K.M., Bramlett, K., Agius, P., ..., and Leiman, L.C. 2022. Contrived materials and a dataset for the evaluation of liquid biopsy tests: A Blood Profiling Atlas in Cancer (BLOODPAC) community study

## The datasets

Results from the manuscript can be reproduced from a single CSV file that is registered to the [BLOODPAC Data Commons](https://data.bloodpac.org/).
The GUID for this file is `dg.5B0D/b47b8599-1488-49d1-9047-d5b3cc661806`. (_TODO: Add in discovery page permalink and instructions_).

## Analysis

Most figures and tables from the manuscript can be reproduced by running the R-script provided in this repository.

The main script is `bpa-jfdi-p2-results.R`. To use this script, you need to have the dataset downloaded
and make sure the required R packages are installed.

The R package requirements include:

* `tidyverse`
* `ggpubr`
* `combat`

R session information:

```
R version 4.1.0 (2021-05-18)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 12.6

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base

other attached packages:
 [1] combinat_0.0-8  ggpubr_0.4.0    forcats_0.5.1   stringr_1.4.0   dplyr_1.0.7
 [6] purrr_0.3.4     readr_2.0.1     tidyr_1.1.3     tibble_3.1.4    ggplot2_3.3.5
[11] tidyverse_1.3.1

loaded via a namespace (and not attached):
 [1] tidyselect_1.1.1  haven_2.4.3       carData_3.0-4     colorspace_2.0-2
 [5] vctrs_0.3.8       generics_0.1.0    utf8_1.2.2        rlang_0.4.11
 [9] pillar_1.6.2      foreign_0.8-81    glue_1.4.2        withr_2.4.2
[13] DBI_1.1.1         dbplyr_2.1.1      modelr_0.1.8      readxl_1.3.1
[17] lifecycle_1.0.0   munsell_0.5.0     ggsignif_0.6.2    gtable_0.3.0
[21] cellranger_1.1.0  zip_2.2.0         rvest_1.0.1       rio_0.5.27
[25] tzdb_0.1.2        curl_4.3.2        fansi_0.5.0       broom_0.7.9
[29] Rcpp_1.0.7        scales_1.1.1      backports_1.2.1   jsonlite_1.7.2
[33] abind_1.4-5       fs_1.5.0          hms_1.1.0         openxlsx_4.2.4
[37] stringi_1.7.4     rstatix_0.7.0     grid_4.1.0        cli_3.0.1
[41] tools_4.1.0       magrittr_2.0.1    crayon_1.4.1      car_3.0-11
[45] pkgconfig_2.0.3   ellipsis_0.3.2    data.table_1.14.0 xml2_1.3.2
[49] reprex_2.0.1      lubridate_1.7.10  assertthat_0.2.1  httr_1.4.2
[53] rstudioapi_0.13   R6_2.5.1          compiler_4.1.0
```

After making sure you have the required files and packages, you can generate the figures by running:

```
Rscript bpa-jfdi-p2-results.R </path/to/dataset.csv>
```

The script will create a directory `figures/` where all the figures will be saved as PNGs.
