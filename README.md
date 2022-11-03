> This repository houses the code and data for the article: _Thermal performance
> of the electron transport system Complex III in seven Alabama fishes_ by L. Horne, D. DeVries, R. Wright, E. Irwin, B. Staton, H. Abdelrahman, and J. Stoeckel accepted for publication in the _Journal of Experimental Zoology: Series A_

[![ArticleDOI](https://img.shields.io/badge/Article%20DOI-10.1002%2Fjez.2267-blue)](https://doi.org/10.1002/jez.2667)

[![GitHub Repo Archive DOI](https://img.shields.io/badge/GitHub%20Repo%20Archive%20DOI-10.5281%2Fzenodo.7187078-blue)](https://doi.org/10.5281/zenodo.7187078)

## Repository Structure

*  `analysis.R` contains the R code to fit all GAMs, bootstrap them, obtain critical temperature values from the bootstrapped curves, and perform the ANOVAs comparing temperature values across species, all as described in the article. This file also contains code to reproduce all plots found in the manuscript and prints or exports files containing the numerical values contained in the results text.
*  `functions.R` contains several custom R functions that perform the GAM fitting, bootstrapping, and critical temperature value derivation. Is `source()`-ed by `analysis.R` prior to performing these tasks.
*  `ETS-vs-temp-data.csv` contains the raw data for the analysis. Each row is one individual of one species observed at one temperature and contains four variables:  
    -  **Fish**: an identifier for individual nested within species
    -  **Temp**: the temperature (&#176;C) of the observation
    -  **ETS**: the activity measurement of the observation (mL O<sub>2</sub> $\cdot$ gWW<sup>-1</sup> $\cdot$ hr<sup>-1</sup>)
    -  **Species**: an identifier for unique species; the species codes are as follows:
    <br>
    
    | Scientific Name           | Common Name            | Code in Data |
    | ------------------------- | ---------------------- | ------------ |
    | *Campostoma oligolepis*   | Largescale Stoneroller | `stone`      |
    | *Notropis baileyi*        | Rough Shiner           | `rough`      |
    | *Cyprinella gibbsi*       | Tallapoosa Shiner      | `poosa`      |
    | *Lepomis macrochirus*     | Bluegill Sunfish       | `bluegill`   |
    | *Cyprinella venusta*      | Blacktail Shiner       | `black`      |
    | *Semotilus atromaculatus* | Creek Chub             | `creek`      |
    | *Cottus carolinae*        | Banded Sculpin         | `sculpin`    |
    

## Model Output

This repository does not track the output files, since they can be reproduced from source at any time with the above files.

To replicate the analysis, open the file `ETS-vs-temp-ms-analysis.Rproj` file in RStudio and then execute the entire `analysis.R` script. You will have a new `output` subdirectory created that will contain all output generated for the article. An addition PDF file will be created (`output/boot-gams.pdf`) that displays the fitted and bootstrap GAMs for each individual, as well as the derived critical temperature values. 



## Session Info

For reproducibility purposes, all analyses were conducted using this configuration:

```
- Session info --------------------------------------------------------------------------
 setting  value                       
 version  R version 4.0.2 (2020-06-22)
 os       Windows 10 x64              
 system   x86_64, mingw32             
 ui       RStudio                     
 language (EN)                        
 collate  English_United States.1252  
 ctype    English_United States.1252  
 tz       America/Los_Angeles         
 date     2022-10-06                  

- Packages ------------------------------------------------------------------------------
 package      * version  date       lib source        
 cli            3.3.0    2022-04-25 [1] CRAN (R 4.0.2)
 coda           0.19-3   2019-07-05 [1] CRAN (R 4.0.2)
 codetools      0.2-16   2018-12-24 [1] CRAN (R 4.0.2)
 colorspace     2.0-3    2022-02-21 [1] CRAN (R 4.0.5)
 emmeans        1.7.2    2022-01-04 [1] CRAN (R 4.0.5)
 estimability   1.3      2018-02-11 [1] CRAN (R 4.0.3)
 farver         2.1.0    2021-02-28 [1] CRAN (R 4.0.5)
 latex2exp      0.4.0    2015-11-30 [1] CRAN (R 4.0.2)
 lattice        0.20-41  2020-04-02 [1] CRAN (R 4.0.2)
 lifecycle      1.0.1    2021-09-24 [1] CRAN (R 4.0.5)
 magrittr       2.0.3    2022-03-30 [1] CRAN (R 4.0.5)
 MASS           7.3-51.6 2020-04-26 [1] CRAN (R 4.0.2)
 Matrix         1.2-18   2019-11-27 [1] CRAN (R 4.0.2)
 mgcv         * 1.8-31   2019-11-09 [1] CRAN (R 4.0.2)
 multcomp       1.4-18   2022-01-04 [1] CRAN (R 4.0.5)
 multcompView   0.1-8    2019-12-19 [1] CRAN (R 4.0.3)
 munsell        0.5.0    2018-06-12 [1] CRAN (R 4.0.2)
 mvtnorm        1.1-1    2020-06-09 [1] CRAN (R 4.0.0)
 nlme         * 3.1-148  2020-05-24 [1] CRAN (R 4.0.2)
 R6             2.5.1    2021-08-19 [1] CRAN (R 4.0.5)
 rlang          1.0.2    2022-03-04 [1] CRAN (R 4.0.5)
 sandwich       3.0-1    2021-05-18 [1] CRAN (R 4.0.5)
 scales         1.2.0    2022-04-13 [1] CRAN (R 4.0.5)
 sessioninfo    1.1.1    2018-11-05 [1] CRAN (R 4.0.2)
 stringi        1.7.6    2021-11-29 [1] CRAN (R 4.0.5)
 stringr        1.4.0    2019-02-10 [1] CRAN (R 4.0.2)
 survival       3.1-12   2020-04-10 [1] CRAN (R 4.0.2)
 TH.data        1.1-0    2021-09-27 [1] CRAN (R 4.0.5)
 withr          2.2.0    2020-04-20 [1] CRAN (R 4.0.2)
 xtable         1.8-4    2019-04-21 [1] CRAN (R 4.0.2)
 zoo            1.8-8    2020-05-02 [1] CRAN (R 4.0.2)

[1] C:/Users/bstaton/Documents/R/R-4.0.2/library
```
