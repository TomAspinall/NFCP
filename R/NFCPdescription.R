
#' @name NFCP
#' N-factor Commodity Pricing Through Term Structure Estimation
#'
#'Commodity pricing models are (systems of) stochastic differential equations that are utilized for the valuation and hedging of commodity contingent claims (i.e. derivative products on the commodity) and
#'other commodity related investments. Commodity pricing models that capture market dynamics are of great importance to commodity market participants in order to exercise sound investment and risk-management
#'strategies. Parameters of commodity pricing models are estimated through maximum likelihood estimation, using available term structure futures data of a commodity.
#''NFCP' (n-factor commodity pricing) provides a framework for the modeling, parameter estimation, probabilistic forecasting, option valuation and simulation of commodity prices through state space and
#'Monte Carlo methods, risk-neutral valuation and Kalman filtering. 'NFCP' allows the commodity pricing model to consist of n correlated factors, with both random walk and mean-reverting elements.
#'The n-factor commodity pricing model framework was first presented in the work of Cortazar and Naranjo (2006) <doi:10.1002/fut.20198>.
#'Examples presented in 'NFCP' replicate the two-factor crude oil commodity pricing model presented in the prolific work of Schwartz and Smith (2000) <doi:10.1287/mnsc.46.7.893.12034> with the approximate
#'term structure futures data applied within this study provided in the 'NFCP' package.
#'
#' The primary features of \code{NFCP} are:
#'
#' -	Estimate commodity pricing models through maximum likelihood estimation using the Kalman Filter and term structure time-series data.
#' -	Evaluate the fit and robustness of commodity pricing models to term structure data.
#' -	Probabilistically forecast future spot and futures prices analytically.
#' -  Value European call and put options under an N-factor model analytically.
#' -  Simulate risk-neutral spot and futures price paths of commodities through Monte Carlo Simulation.
#' -  Value American put options under an N-factor model numerically through the least-squares Monte Carlo simulation method.
#'
#' @author
#' Thomas Aspinall tomaspinall2512@gmail.com (0000-0002-6968-1989)
#' Adrian Gepp (0000-0003-1666-5501)
#' Geoff Harris (0000-0003-4284-8619)
#' Simone Kelly (0000-0002-6528-8557)
#' Colette Southam (0000-0001-7263-2347)
#' Bruce Vanstone (0000-0002-3977-2468)
