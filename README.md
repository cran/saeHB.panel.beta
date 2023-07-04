
<!-- README.md is generated from README.Rmd. Please edit that file -->

# saeHB.panel.beta

Several functions are provided for small area estimation at the area
level using the hierarchical bayesian (HB) method with panel data under
beta distribution for variable interest. This package also provides a
dataset produced by data generation. The ‘rjags’ package is employed to
obtain parameter estimates. Model-based estimators involve the HB
estimators, which include the mean and the variation of the mean. For
the reference, see Rao and Molina (2015, <ISBN:978-1-118-73578-7>).

## Author

Dian Rahmawati Salis, Azka Ubaidillah

## Maintaner

Dian Rahmawati Salis <dianrahmawatisalis03@gmail.com>

## Function

- `RaoYuAr1.beta()` This function gives estimation of y using
  Hierarchical Bayesian Rao Yu Model under Beta distribution
- `Panel.beta()` This function gives estimation of y using Hierarchical
  Bayesian Rao Yu Model under Beta distribution when rho = 0

## Installation

You can install the development version of saeHB.panel.beta from GitHub
with:

``` r
# install.packages("devtools")
devtools::install_github("DianRahmawatiSalis/saeHB.panel.beta")
#> Downloading GitHub repo DianRahmawatiSalis/saeHB.panel.beta@HEAD
#> rlang  (1.0.6 -> 1.1.1) [CRAN]
#> cli    (3.6.0 -> 3.6.1) [CRAN]
#> vctrs  (0.5.2 -> 0.6.3) [CRAN]
#> tibble (3.1.8 -> 3.2.1) [CRAN]
#> rjags  (4-13  -> 4-14 ) [CRAN]
#> dplyr  (1.1.0 -> 1.1.2) [CRAN]
#> Installing 6 packages: rlang, cli, vctrs, tibble, rjags, dplyr
#> Installing packages into 'C:/Users/LENOVO/AppData/Local/R/win-library/4.2'
#> (as 'lib' is unspecified)
#> package 'rlang' successfully unpacked and MD5 sums checked
#> Warning: cannot remove prior installation of package 'rlang'
#> Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\00LOCK\rlang\libs\x64\rlang.dll
#> to C:\Users\LENOVO\AppData\Local\R\win-library\4.2\rlang\libs\x64\rlang.dll:
#> Permission denied
#> Warning: restored 'rlang'
#> package 'cli' successfully unpacked and MD5 sums checked
#> Warning: cannot remove prior installation of package 'cli'
#> Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\00LOCK\cli\libs\x64\cli.dll to
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\cli\libs\x64\cli.dll: Permission
#> denied
#> Warning: restored 'cli'
#> package 'vctrs' successfully unpacked and MD5 sums checked
#> Warning: cannot remove prior installation of package 'vctrs'
#> Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\00LOCK\vctrs\libs\x64\vctrs.dll
#> to C:\Users\LENOVO\AppData\Local\R\win-library\4.2\vctrs\libs\x64\vctrs.dll:
#> Permission denied
#> Warning: restored 'vctrs'
#> package 'tibble' successfully unpacked and MD5 sums checked
#> Warning: cannot remove prior installation of package 'tibble'
#> Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\00LOCK\tibble\libs\x64\tibble.dll
#> to C:\Users\LENOVO\AppData\Local\R\win-library\4.2\tibble\libs\x64\tibble.dll:
#> Permission denied
#> Warning: restored 'tibble'
#> package 'rjags' successfully unpacked and MD5 sums checked
#> Warning: cannot remove prior installation of package 'rjags'
#> Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\00LOCK\rjags\libs\x64\rjags.dll
#> to C:\Users\LENOVO\AppData\Local\R\win-library\4.2\rjags\libs\x64\rjags.dll:
#> Permission denied
#> Warning: restored 'rjags'
#> package 'dplyr' successfully unpacked and MD5 sums checked
#> Warning: cannot remove prior installation of package 'dplyr'
#> Warning in file.copy(savedcopy, lib, recursive = TRUE): problem copying
#> C:\Users\LENOVO\AppData\Local\R\win-library\4.2\00LOCK\dplyr\libs\x64\dplyr.dll
#> to C:\Users\LENOVO\AppData\Local\R\win-library\4.2\dplyr\libs\x64\dplyr.dll:
#> Permission denied
#> Warning: restored 'dplyr'
#> 
#> The downloaded binary packages are in
#>  C:\Users\LENOVO\AppData\Local\Temp\Rtmp6VfJ5T\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\LENOVO\AppData\Local\Temp\Rtmp6VfJ5T\remotes1b885e6f22c1\DianRahmawatiSalis-saeHB.panel.beta-b359c9b/DESCRIPTION' ...     checking for file 'C:\Users\LENOVO\AppData\Local\Temp\Rtmp6VfJ5T\remotes1b885e6f22c1\DianRahmawatiSalis-saeHB.panel.beta-b359c9b/DESCRIPTION' ...   ✔  checking for file 'C:\Users\LENOVO\AppData\Local\Temp\Rtmp6VfJ5T\remotes1b885e6f22c1\DianRahmawatiSalis-saeHB.panel.beta-b359c9b/DESCRIPTION' (4.7s)
#>       ─  preparing 'saeHB.panel.beta': (1.5s)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#>       ─  checking for LF line-endings in source and make files and shell scripts (999ms)
#>               checking for empty or unneeded directories  ─  checking for empty or unneeded directories
#>       ─  building 'saeHB.panel.beta_0.1.2.tar.gz'
#>      
#> 
#> Installing package into 'C:/Users/LENOVO/AppData/Local/R/win-library/4.2'
#> (as 'lib' is unspecified)
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(saeHB.panel.beta)
data("dataPanelbeta")
dataPanelbeta <- dataPanelbeta[1:25,] #for the example only use part of the dataset
formula <- ydi~xdi1+xdi2 
area <- max(dataPanelbeta[,2])
period <- max(dataPanelbeta[,3])
result<-Panel.beta(formula,area=area, period=period ,iter.mcmc = 10000,thin=5,burn.in = 1000,data=dataPanelbeta)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 25
#>    Unobserved stochastic nodes: 62
#>    Total graph size: 359
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 25
#>    Unobserved stochastic nodes: 62
#>    Total graph size: 359
#> 
#> Initializing model
#> 
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 25
#>    Unobserved stochastic nodes: 62
#>    Total graph size: 359
#> 
#> Initializing model
```

Extract area mean estimation

``` r
result$Est
#>              MEAN         SD      2.5%       25%       50%       75%     97.5%
#> mu[1,1] 0.9717543 0.02049596 0.9190452 0.9647537 0.9771319 0.9855882 0.9944996
#> mu[2,1] 0.9510178 0.03288811 0.8726664 0.9375707 0.9584745 0.9742376 0.9896039
#> mu[3,1] 0.9423424 0.04229746 0.8342821 0.9257270 0.9531654 0.9709367 0.9884046
#> mu[4,1] 0.9685811 0.02368207 0.9106474 0.9610285 0.9744508 0.9837632 0.9934777
#> mu[5,1] 0.9404871 0.04833029 0.8108703 0.9242932 0.9548771 0.9729015 0.9889299
#> mu[1,2] 0.9712634 0.02151995 0.9176578 0.9634467 0.9766530 0.9855381 0.9940899
#> mu[2,2] 0.9623792 0.02698444 0.8903765 0.9515899 0.9698480 0.9802553 0.9930563
#> mu[3,2] 0.9199737 0.05673965 0.7725415 0.8974167 0.9355057 0.9591863 0.9829062
#> mu[4,2] 0.9785373 0.01764545 0.9313410 0.9731022 0.9830868 0.9898914 0.9964443
#> mu[5,2] 0.9380432 0.04478322 0.8197203 0.9205959 0.9494996 0.9685290 0.9868380
#> mu[1,3] 0.9706437 0.02230411 0.9062298 0.9628493 0.9761609 0.9859619 0.9954731
#> mu[2,3] 0.8655381 0.07748549 0.6717318 0.8303173 0.8821522 0.9202044 0.9650687
#> mu[3,3] 0.9514076 0.03298669 0.8613422 0.9374960 0.9599521 0.9751765 0.9904608
#> mu[4,3] 0.9581918 0.02852203 0.8840658 0.9465256 0.9652434 0.9774888 0.9910127
#> mu[5,3] 0.9168573 0.05704092 0.7595892 0.8949962 0.9311539 0.9560012 0.9824765
#> mu[1,4] 0.9551522 0.02996238 0.8751646 0.9425921 0.9623560 0.9763032 0.9909293
#> mu[2,4] 0.9342440 0.04342608 0.8190859 0.9170173 0.9451454 0.9635272 0.9851920
#> mu[3,4] 0.9334419 0.04256771 0.8268696 0.9137840 0.9447703 0.9633846 0.9844740
#> mu[4,4] 0.9757125 0.02018988 0.9194181 0.9701157 0.9809731 0.9882532 0.9956180
#> mu[5,4] 0.8548727 0.09644667 0.5998419 0.8110210 0.8829172 0.9249182 0.9678862
#> mu[1,5] 0.9682527 0.02270194 0.9080987 0.9594752 0.9742140 0.9836890 0.9938743
#> mu[2,5] 0.8867797 0.07030962 0.7065786 0.8544728 0.9037534 0.9359809 0.9757281
#> mu[3,5] 0.9575245 0.02978026 0.8859338 0.9445934 0.9654878 0.9779015 0.9917703
#> mu[4,5] 0.9300983 0.04640135 0.8169823 0.9101969 0.9413470 0.9620436 0.9852102
#> mu[5,5] 0.8637439 0.08470990 0.6407232 0.8247754 0.8853290 0.9257943 0.9666170
```

Extract coefficient estimation

``` r
result$coefficient
#>          Mean        SD      2.5%       25%      50%      75%    97.5%
#> b[0] 1.945340 0.3914832 1.1681623 1.6812227 1.951651 2.214936 2.712910
#> b[1] 1.169466 0.5381308 0.1030402 0.8199215 1.161090 1.528560 2.200846
#> b[2] 1.137731 0.4537707 0.2633356 0.8324130 1.125769 1.435226 2.020636
```

Extract area random effect variance

``` r
result$refVar
#> [1] 0.4626572
```

Extract MSE

``` r
MSE_HB<-result$Est$SD^2
summary(MSE_HB)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0003114 0.0005608 0.0010881 0.0021821 0.0023358 0.0093020
```

Extract RSE

``` r
RSE_HB<-sqrt(MSE_HB)/result$Est$MEAN*100
summary(RSE_HB)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.803   2.445   3.467   4.528   5.139  11.282
```

Extract convergence diagnostic using geweke test

``` r
result$covergence.test
#> NULL
```

## References

- Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
  York: John Wiley and Sons, Inc.
- Torabi, M., & Shokoohi, F. (2012). Likelihood inference in small area
  estimation by combining time-series and cross-sectional data. Journal
  of Multivariate Analysis, 111, 213–221.
  <https://doi.org/10.1016/j.jmva.2012.05.016>
