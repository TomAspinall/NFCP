#'N-Factor Model American Put Option Pricing
#'
#' @description Value American Put Options under the parameters of an N-factor model through the Least-Squares Monte Carlo (LSM) Simulation Method.
#' This function is a wrapper to the 'LSM.AmericanOption' function of the 'LSMRealOptions' package.
#'
#' @param X.0 Initial values of the state vector.
#' @param parameters Named vector of parameter values of a specified N-factor model. Function \code{NFCP.Parameters} is recommended.
#' @param n total number of simulated price paths
#' @param t Time to expiration of the option (in years)
#' @param dt discrete time step of simulation
#' @param K Strike price of the American put option
#' @param r Risk-free interest rate.
#' @param orthogonal The orthogonal polynomial used to approximate the continuation value of the option in the LSM simulation method.
#' Orthogonal polynomial arguments available are: "Power", "Laguerre", "Jacobi", "Legendre", "Chebyshev", "Hermite". See \code{help(LSM.AmericanOption)}
#' @param degree 	The degree of polynomials used in the least squares fit. See \code{help(LSM.AmericanOption)}
#' @param verbose \code{logical} Should additional information be output?
#' @param debugging \code{logical} Should additional simulation information be output?
#'
#'@details
#'
#'The 'American.Option.Value' function is a wrapper to the 'Spot.Price.Simulate' and 'LSM.AmericanOption' of the 'LSMRealOptions' package that
#'returns the value of American put options under a given N-factor model.
#'
#'The least-squares Monte Carlo (LSM) simulation method is an option valuation method first presented by Longstaff and Schwartz (2001) that
#'approximates the value of American options.
#'
#'Methods to solve for the value of options with early exercise opportunities include partial differential equations, lattice-based methods
#'and Monte-Carlo simulation. LSM simulation is the optimal solution method to value American options under an N-factor model due
#'to the multiple factors that can make up the spot price process and influence the option value. Whilst lattice and partial differential equation
#'approaches suffer from the 'curse of dimensionality', LSM simulation may be readily applied under multi-factor settings.
#'
#'Longstaff and Schwartz (2001) state that as the conditional expectation of the continuation value belongs to a Hilbert space,
#'it can be represented by a combination of orthogonal basis functions. Increasing the number of stochastic state variables
#'therefore increases the number of required basis functions exponentially.
#'
#'@return
#'
#'The 'American.Option.Value' function by default returns a \code{numeric} object corresponding to the calculated value of the American put option.
#'
#'When \code{verbose = T}, 6 objects are returned within a \code{list} class object. The objects returned are:
#'
#'\tabular{ll}{
#'
#' \code{Value} \tab The calculated option value. \cr
#' \code{Standard Error} \tab The standard error of the calculated option value. \cr
#' \code{Expected Timing} \tab The expected time of early exercise.. \cr
#' \code{Expected Timing SE} \tab The standard error of the expected time of early exercise. \cr
#' \code{Exercise Probability} \tab The probability of early exercise of the option being exercised. \cr
#' \code{Cumulative Exercise Probability} \tab \code{vector}. The cumulative probability of option exercise at each discrete observation point \cr
#' }
#'
#'When \code{debugging = T}, an additional 2 objects are returned within the \code{list} class object. These are the results output by both the 'Spot.Price.Simulate' and
#''LSM.AmericanOption' of the 'LSMRealOptions' package respectively. The objects returned are:
#'
#' \tabular{ll}{
#' \code{State_Variables} \tab A matrix of simulated state variables for each factor is returned when \code{verbose = T}. The number of factors returned corresponds to the number of factors in the specified N-factor model. \cr
#' \code{Prices} \tab A matrix of simulated price paths. Each column represents one simulated price path and each row represents one simulated observation. \cr
#' }
#'
#'@references
#'
#'Longstaff, F.A., and E.S. Schwartz. 2001. Valuing American Options by Simulation: A Simple Least-Squares Approach. The Review of Financial Studies. 14:113-147.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Aspinall, T., A. Gepp, G. Harris, S. Kelly, C. Southam, and B. Vanstone, (2021). LSMRealOptions: Value
#'American and Real Options Through LSM Simulation. R package version 0.1.0.
#'
#'@examples
#'
#' ##Example 1 - An American put option on a stock following 'GBM'
#' ##growing at the risk-free rate:
#' American.Option.Value(X.0 = log(36),
#'     parameters = c(mu_star = (0.06 - (1/2) * 0.2^2), sigma_1 = 0.2),
#'     n = 1e2,
#'     t = 1,
#'     dt = 1/50,
#'     K = 40,
#'     r = 0.05,
#'     verbose = FALSE,
#'     orthogonal = "Laguerre",
#'     degree = 3)
#'
#'##Example 2 - An American put option under a two-factor crude oil model:
#'
#'##Step 1 - Obtain current (i.e. most recent) state vector by filtering the
#'##two-factor oil model:
#'Schwartz.Smith.Oil <- NFCP.Kalman.filter(parameter.values = SS.Oil$Two.Factor,
#'                                         parameters = names(SS.Oil$Two.Factor),
#'                                         log.futures = log(SS.Oil$Stitched.Futures),
#'                                         dt = SS.Oil$dt,
#'                                         TTM = SS.Oil$Stitched.TTM,
#'                                         verbose = TRUE)
#'
#'##Step 2 - Calculate 'put' option price:
#'American.Option.Value(X.0 = Schwartz.Smith.Oil$X.t,
#'                      parameters = SS.Oil$Two.Factor,
#'                      n = 1e2,
#'                      t = 1,
#'                      dt = 1/12,
#'                      K = 20,
#'                      r = 0.05,
#'                      verbose = FALSE,
#'                      orthogonal = "Power",
#'                      degree = 2)
#' @export
American.Option.Value <- function(X.0, parameters, n, t, dt, K, r, orthogonal = "Power", degree = 2, verbose = FALSE, debugging = FALSE){

  ## Step 1 - Simulate the N-factor model:
  Simulated.States <- Spot.Price.Simulate(X.0 = X.0, parameters = parameters, t = t, dt = dt, n = n, verbose = TRUE)

  ## Step 2 - Apply the 'LSM.AmericanOption' function from the 'LSMRealOptions' to value the American put option:
  output <- LSMRealOptions::LSM.AmericanOption(state.variables = Simulated.States$State_Variables,
                                               payoff = Simulated.States$Prices,
                                               K = K,
                                               dt = dt,
                                               rf = r,
                                               orthogonal = orthogonal,
                                               degree = degree,
                                               verbose = verbose)
  if(debugging){
    return(c(output, Simulated.States))
  } else {
    return(output)
  }
}
