
#' @name NFCP
#' N-factor Commodity Term Structure Model estimation, analysis, forecasting and derivative pricing
#'
#' The \code{NFCP} package provides tools to develop Term Structure Commodity Pricing Models. These models can feature both Random Walk and Mean-Reverting Behaviour,
#' characterized as the exponential sum of N unobservable state variables. The package allows for parameter estimation of N-factor stochastic models through
#' the Kalman filter and maximum likelihood estimation. Expected Spot and Futures prices can be probabilistically forecasted through analytic expressions and European
#' Call and Put Options can be values analytically. The package also allows risk-neutral spot and futures prices to be simulated through Monte-Carlo methods.
#'
#' The primary features of \code{NFCP} are:
#'
#' -	Estimate Commodity Pricing Models through Maximum Likelihood Estimation using the Kalman Filter and Term Structure time-series data.
#' -	Evaluate the fit and robustness of Commodity Pricing models to Term Structure data.
#' -	Probabilistically forecast future spot and futures prices analytically
#' -  Value European Call and Put options under an N-factor model analytically
#' -  Simulate risk-neutral spot and futures price paths of commodities through Monte-Carlo Simulation
#'
#' @author
#' Thomas Aspinall tomaspinall2512@gmail.com (0000-0002-6968-1989)
#' Adrian Gepp (0000-0003-1666-5501)
#' Geoff Harris (0000-0003-4284-8619)
#' Simone Kelly (0000-0002-6528-8557)
#' Colette Southam (0000-0001-7263-2347)
#' Bruce Vanstone (0000-0002-3977-2468)
