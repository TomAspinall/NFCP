# NFCP
N-Factor Commodity Pricing Through Term Structure Estimation

<!-- badges: start -->
<!-- badges: end -->

Commodity pricing models are (systems of) stochastic differential equations that are utilized for the valuation and hedging of commodity contingent claims (i.e. derivative products on the commodity) and other commodity related investments. Parameters of commodity pricing models are estimated through maximum likelihood estimation, using available term structure futures data of a commodity. 'NFCP' (n-factor commodity pricing) provides a framework for the modeling, parameter estimation, probabilistic forecasting, option valuation and simulation of commodity prices through state space and Monte Carlo methods, risk-neutral valuation and Kalman filtering. 'NFCP' allows the commodity pricing model to consist of n correlated factors, with both random walk and mean-reverting elements. Commodity pricing models that capture market dynamics are of great importance to commodity market participants in order to exercise sound investment and risk-management strategies. The n-factor commodity pricing model framework was first presented in the work of Cortazar and Naranjo (2006). Examples presented in 'NFCP' replicate the two-factor crude oil commodity pricing model presented in the prolific work of Schwartz and Smith (2000) with the approximate term structure futures data applied within this study provided in the 'NFCP' package. Kalman filtering in 'NFCP' is performed using sequential processing through the 'FKF.SP' package to optimise computational efficiency. Parameter estimation of n-factor models is performed using genetic algorithms through the 'rGenoud' package to maximise the likelihood that a global maximum is reached during maximum likelihood estimation.

Primary features of 'NFCP' include:

- Parameter estimation of n-factor commodity pricing models through state space methods, Kalman filtering and maximum likelihood estimation.

- Analytic pricing of European call and put options under estimated n-factor commodity pricing models

- Probabilistic forecasting and Monte Carlo simulation of future commodity price paths.

## Installation

You can install the released version of NFCP from [CRAN](https://CRAN.R-project.org) with:

```
install.packages("NFCP")
```

And the development version from [GitHub](https://github.com/) with:

```
devtools::install_github("TomAspinall/NFCP")
```
which contains source code for the package starting with version 0.1.0.
