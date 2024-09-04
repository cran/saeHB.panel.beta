
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
#> Skipping install of 'saeHB.panel.beta' from a github remote, the SHA1 (fe67bb61) has not changed since last install.
#>   Use `force = TRUE` to force installation
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
#> mu[1,1] 0.9745402 0.02043167 0.9242318 0.9671592 0.9798158 0.9874685 0.9962959
#> mu[2,1] 0.9529850 0.03378111 0.8687035 0.9391639 0.9606148 0.9762169 0.9920275
#> mu[3,1] 0.9416476 0.04334522 0.8257088 0.9271068 0.9515223 0.9695644 0.9886169
#> mu[4,1] 0.9707650 0.02343253 0.9100078 0.9631909 0.9768535 0.9858918 0.9956301
#> mu[5,1] 0.9392371 0.05230504 0.7937674 0.9226590 0.9552521 0.9731984 0.9904665
#> mu[1,2] 0.9730519 0.02075604 0.9151144 0.9651103 0.9784012 0.9872296 0.9955669
#> mu[2,2] 0.9644632 0.02716892 0.8916524 0.9553131 0.9715651 0.9825200 0.9941719
#> mu[3,2] 0.9190929 0.05909138 0.7582093 0.8974757 0.9337880 0.9586956 0.9846804
#> mu[4,2] 0.9806928 0.01626167 0.9376644 0.9753268 0.9851461 0.9914809 0.9977165
#> mu[5,2] 0.9414686 0.04347486 0.8300576 0.9253633 0.9528456 0.9703039 0.9899214
#> mu[1,3] 0.9727516 0.02296491 0.9084448 0.9653176 0.9785098 0.9877710 0.9961441
#> mu[2,3] 0.8650717 0.08258513 0.6512432 0.8283691 0.8813968 0.9243896 0.9664624
#> mu[3,3] 0.9547238 0.03083291 0.8773217 0.9424085 0.9616815 0.9760577 0.9921776
#> mu[4,3] 0.9604722 0.02744216 0.8909641 0.9491892 0.9671605 0.9790931 0.9939556
#> mu[5,3] 0.9185326 0.05607900 0.7741958 0.8964084 0.9307007 0.9567351 0.9842448
#> mu[1,4] 0.9584364 0.03167149 0.8737031 0.9456576 0.9671949 0.9791662 0.9929711
#> mu[2,4] 0.9360060 0.04438870 0.8225671 0.9184701 0.9461601 0.9665405 0.9865605
#> mu[3,4] 0.9350573 0.04267082 0.8231444 0.9169837 0.9452816 0.9656088 0.9877537
#> mu[4,4] 0.9774635 0.01875693 0.9297391 0.9713515 0.9826057 0.9896004 0.9971276
#> mu[5,4] 0.8457488 0.10996071 0.5617364 0.8037121 0.8748879 0.9238268 0.9705472
#> mu[1,5] 0.9700137 0.02447969 0.9072987 0.9622983 0.9762001 0.9854658 0.9954380
#> mu[2,5] 0.8870202 0.07028000 0.7023668 0.8552202 0.9036043 0.9369899 0.9772847
#> mu[3,5] 0.9605815 0.03019414 0.8833168 0.9502352 0.9675795 0.9797505 0.9936401
#> mu[4,5] 0.9366477 0.04333110 0.8299269 0.9170241 0.9469230 0.9665555 0.9886723
#> mu[5,5] 0.8613721 0.08544666 0.6459957 0.8197358 0.8803638 0.9229781 0.9675223
```

Extract coefficient estimation

``` r
result$coefficient
#>          Mean        SD      2.5%       25%      50%      75%    97.5%
#> b[0] 1.932909 0.3961688 1.1738158 1.6717925 1.933390 2.195974 2.704684
#> b[1] 1.188665 0.5270641 0.1656371 0.8391137 1.179424 1.531993 2.244406
#> b[2] 1.206062 0.4730080 0.3134756 0.8882160 1.197545 1.528039 2.165126
```

Extract area random effect variance

``` r
result$refVar
#> [1] 0.5076331
```

Extract MSE

``` r
MSE_HB<-result$Est$SD^2
summary(MSE_HB)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0002644 0.0005993 0.0011412 0.0023440 0.0027358 0.0120914
```

Extract RSE

``` r
RSE_HB<-sqrt(MSE_HB)/result$Est$MEAN*100
summary(RSE_HB)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   1.658   2.524   3.545   4.626   5.569  13.002
```

Extract convergence diagnostic using geweke test

``` r
result$convergence.test
#>              b[0]      b[1]     b[2]
#> Z-score 0.7387093 0.1672595 1.760278
```

## References

- Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
  York: John Wiley and Sons, Inc.
- Torabi, M., & Shokoohi, F. (2012). Likelihood inference in small area
  estimation by combining time-series and cross-sectional data. Journal
  of Multivariate Analysis, 111, 213–221.
  <https://doi.org/10.1016/j.jmva.2012.05.016>
