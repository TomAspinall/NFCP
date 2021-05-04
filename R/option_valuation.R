#'N-factor model American put option pricing
#'
#' @description Value American put options under the parameters of an N-factor model through the Least-Squares Monte Carlo (LSM) Simulation Method.
#' This function is a wrapper to the 'LSM_American_option' function of the 'LSMRealOptions' package.
#'
#' @param x_0 Initial values of the state vector.
#' @param parameters Named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#' @param N_simulations total number of simulated price paths
#' @param option_maturity Time to expiration of the option (in years)
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
#'The 'American_option_value' function is a wrapper to the 'spot_price_simulate' and 'LSM_American_option' of the 'LSMRealOptions' package that
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
#'The 'American_option_value' function by default returns a \code{numeric} object corresponding to the calculated value of the American put option.
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
#'# Example 1 - An American put option on a stock following 'GBM'
#'# growing at the risk-free rate:
#'American_option_value(x_0 = log(36),
#'                      parameters = c(mu_rn = (0.06 - (1/2) * 0.2^2), sigma_1 = 0.2),
#'                      N_simulations = 1e2,
#'                      option_maturity = 1,
#'                      dt = 1/50,
#'                      K = 40,
#'                      r = 0.05,
#'                      verbose = FALSE,
#'                      orthogonal = "Laguerre",
#'                      degree = 3)
#'
#'# Example 2 - An American put option under a two-factor crude oil model:
#'
#'## Step 1 - Obtain current (i.e. most recent) state vector by filtering the
#'## two-factor oil model:
#'Schwartz_Smith_oil <- NFCP_Kalman_filter(parameter_values = SS_oil$two_factor,
#'                                         parameter_names = names(SS_oil$two_factor),
#'                                         log_futures = log(SS_oil$stitched_futures),
#'                                         dt = SS.Oil$dt,
#'                                         futures_TTM = SS_oil$stitched_TTM,
#'                                         verbose = TRUE)
#'
#'##Step 2 - Calculate 'put' option price:
#'American_option_value(x_0 = Schwartz_Smith_oil$x_t,
#'                      parameters = SS_oil$two_factor,
#'                      N_simulations = 1e2,
#'                      option_maturity = 1,
#'                      dt = 1/12,
#'                      K = 20,
#'                      r = 0.05,
#'                      verbose = FALSE,
#'                      orthogonal = "Power",
#'                      degree = 2)
#'@export
American_option_value <- function(x_0, parameters, N_simulations, option_maturity, dt, K, r, orthogonal = "Power", degree = 2, verbose = FALSE, debugging = FALSE){

  if("mu_star" %in% names(parameters)){
    warning("'mu_star' is deprecated. please rename this parameter to 'mu_rn'")
    parameters["mu_rn"] <- parameters["mu_star"]
    parameters <- parameters[!names(parameters) %in% "mu_star"]
  }

  ## Step 1 - Simulate the N-factor model:
  simulated_states <- spot_price_simulate(x_0 = x_0, parameters = parameters, t = option_maturity, dt = dt, N_simulations = N_simulations, verbose = TRUE)

  ## Step 2 - Apply the 'LSM.AmericanOption' function from the 'LSMRealOptions' to value the American put option:
  output <- suppressWarnings( LSMRealOptions::LSM.AmericanOption(state.variables = simulated_states$state_variables,
                                               payoff = simulated_states$spot_prices,
                                               K = K,
                                               dt = dt,
                                               rf = r,
                                               orthogonal = orthogonal,
                                               degree = degree,
                                               verbose = verbose))
  if(debugging){
    return(c(output, simulated_states))
  } else {
    return(output)
  }
}





#'N-factor model European option pricing
#'@description Value European Option Put and Calls under the parameters of an N-factor model.
#'
#'@param x_0 Initial values of the state vector.
#'@param parameters Named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param futures_maturity Time, in years, when the underlying futures contract matures.
#'@param option_maturity Time, in years,  when the European option expires.
#'@param K Strike price of the European Option
#'@param r Risk-free interest rate.
#'@param call \code{logical} is the European option a call or put option?
#'@param verbose \code{logical}. Should additional information be output? see \bold{details}
#'
#'@details
#'
#'\loadmathjax
#'
#'The \code{European_option_value} function calculates analytic expressions of the value of European call and put options on futures contracts within the N-factor model. Under the assumption that future futures prices
#'are log-normally distributed under the risk-neutral process, there exist analytic expressions of the value of European call and put options on futures contracts. The following analytic expression follows from that presented by Schwartz and Smith (2000) extended to the N-factor framework. The value of a European option on a futures contract
#'is given by calculating its expected future value using the risk-neutral process and subsequently discounting at the risk-free rate.
#'
#'One can verify that under the risk-neutral process, the expected futures price at time \mjeqn{t}{t} is:
#'
#'\mjdeqn{E^*[F_{T,t}] = exp(\sum_{i=1}^Ne^{-\kappa_iT}x_{i}(0) + \mu^*t +  A(T-t) +
#'\frac{1}{2}(\sigma_1^2t+\sum_{i.j\neq1}{e^{-\left(\kappa_i+\kappa_j\right)\left(T-t\right)}\sigma_i\sigma_j\rho_{i,j}.\frac{1-e^{-\left(\kappa_i+\kappa_j\right)t}}{\kappa_i+\kappa_j}})) \equiv F_{T,0}}{
#'E^*[F[T,t]] = exp(sum_{i=1}^N (e^(-kappa[i]*T) * x[i](0) + mu^* * t +  A(T-t) +  1/2 (sigma[1]^2 + sum_[i.j != 1] (e^(- (kappa[i] + kappa[j])(T-t)) * sigma[i] * sigma[j] * rho[i,j] * (1 - e^(-(kappa[i] * kappa[j])t))
#'/(kappa[i] + kappa[j]))) equiv F[T,0]}
#'
#'This follows from the derivation provided within the vignette of the NFCP package as well as the details of the \code{futures_price_forecast} package.
#'The equality of expected futures price at time \mjeqn{t}{t} being  equal to the time-\mjeqn{t}{t} current futures price \mjeqn{F_{T,0}}{F[T,0]} is proven by Futures prices
#'being given by expected spot prices under the risk-neutral process
#'\mjeqn{(F_{T,t}=E_t^\ast\left[S_T\right])}{F[T,t] = E[t]^*(S[T])} and the law of iterated expectations \mjeqn{\left(E^\ast\left[E_t^\ast\left[S_T\right]\right]=E^\ast\left[S_T\right]\right)}{E^*(E[t]^*(S[T])) = E^*(S[T])}
#'
#'Because future futures prices are log-normally distributed under the risk-neutral process, we can write a closed-form expression for valuing European put and call options on these futures. When \mjeqn{T=0}{T=0}
#'these are European options on the spot price of the commodity itself. The value of a European call option on a futures contract maturing at time
#'\mjeqn{T}{T}, with strike price \mjeqn{K}{K}, and with time \mjeqn{t}{t} until the option expires, is:
#'
#'\mjdeqn{e^{-rt}E^\ast\left[\max{\left(F_{T,t}-K,0\right)}\right]}{e^(-rt) * E^*(max(F[T,t] - K, 0))}
#'\mjdeqn{= e^{-rt}( F_{T,0}N(d) - KN(d-\sigma_\phi(t,T)))}{e^(-rt) * (F[T,0] * N(d) - K * N(d - sigma[phi](t,T)))}
#'
#'Where:
#'\mjdeqn{d=\frac{\ln(F/K)}{\sigma_\phi(t,T)}+\frac{1}{2}\sigma_\phi\left(t,T\right)}{d = ln(F/K) / sigma[phi](t,T) + 1/2 sigma[phi](t,T)}
#'
#'and:
#'
#'\mjdeqn{\sigma_\phi\left(t,T\right) = \sqrt(\sigma_1^2t+\sum_{i.j\neq1}{e^{-\left(\kappa_i+\kappa_j\right)\left(T-t\right)}\sigma_i\sigma_j\rho_{i,j}.
#'\frac{1-e^{-\left(\kappa_i+\kappa_j\right)t}}{\kappa_i+\kappa_j}})}{
#'sigma[phi](t,T) = sqrt( sigma[1]^2 + sum_[i.j != 1]( e^(-(kappa[i] + kappa[j])(T-t)) sigma[i] sigma[j] rho[i,j] * (1 - e^(-(kappa[i] + kappa[j])t))/(kappa[i] + kappa[j])))}
#'
#'Parameter \mjeqn{ N(d) }{N(d)} indicates cumulative probabilities for the standard normal distribution (i.e. \mjeqn{P(Z<d)}{P(Z<d)}).
#'
#'Similarly, the value of a European put with the same parameters is given by:
#'
#'\mjdeqn{e^{-rt} E^*[max(K-F_{T,t},0)]}{e^(-rt) E^*(max(K - F[T,t],0))}
#'
#'\mjdeqn{=e^{-rt}\left(-F_{T,0}N\left(-d\right)+KN\left(\sigma_\phi\left(t,T\right)-d\right)\right)}{e^(-rt) * (- F[T,0] * N(-d) + K * N(sigma[phi](t,T) - d))}
#'
#'
#'The presented option valuation formulas are analogous to the Black-Scholes formulas for valuing European options on stocks that do not pay dividends
#'
#'Under this terminology, the stock price corresponds to the present value of the futures commitment \mjeqn{(e^{-rt}F_{T,0})}{e^(-rt) F[T,0]}  and the equivalent annualized volatility would be \mjeqn{\sigma_\phi(t,T)/\sqrt t}{sigma[phi](t,T) / sqrt(t)}
#'
#'When \code{verbose = T}, the \code{European_option_value} function numerically calculates the sensitivity of option prices to the underlying parameters specified within the N-factor model, as well as some of the most common
#'"Greeks" related to European put and call option pricing. All gradients are calculated numerically by calling the \code{grad} function from the \code{numDeriv} package.
#'
#'@return
#'The \code{European_option_value} function returns a numeric value corresponding to the present value of an option when \code{verbose = F}.
#'When \code{verbose = T}, \code{European_option_value} returns a list with three objects:
#'
#'\tabular{ll}{
#'
#' \code{option value} \tab Present value of the option. \cr
#'
#'\code{annualized volatility} \tab Annualized volatility of the option. \cr
#'
#'\code{parameter sensitivity} \tab Sensitivity of the option value to each parameter of the N-factor model. \cr
#'
#'\code{greeks} \tab Sensitivity of the option value to different option parameters. \cr
#'
#' }
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'##Example 1 - A European 'put' option on a stock following 'GBM'
#'##growing at the risk-free rate:
#'
#'### Risk-free rate:
#'rf <- 0.05
#'### Stock price:
#'S_0 <- 20
#'### Stock volatility:
#'S_sigma <- 0.2
#'### Option maturity:
#'Tt <- 1
#'### Exercise price:
#'K <- 20
#'### Calculate 'put' option price:
#'European_option_value(x_0 = log(S_0), parameters = c(mu_rn = rf, sigma_1 = S_sigma),
#'                      futures_maturity = Tt, option_maturity = 1,
#'                      K = K, r = rf, call = FALSE, verbose = TRUE)
#'
#'##Example 2 - A European call option under a two-factor crude oil model:
#'
#'##Step 1 - Obtain current (i.e. most recent) state vector by filtering the
#'##two-factor oil model:
#'Schwartz_Smith_oil <- NFCP_Kalman_filter(parameter_values = SS_oil$two_factor,
#'                                         parameter_names = names(SS_oil$two_factor),
#'                                         log_futures = log(SS_oil$stitched_futures),
#'                                         dt = SS_oil$dt,
#'                                         futures_TTM = SS_oil$stitched_TTM,
#'                                         verbose = TRUE)
#'
#'##Step 2 - Calculate 'call' option price:
#'European_option_value(x_0 = Schwartz_Smith_oil$x_t,
#'                      parameters = SS_oil$two_factor,
#'                      futures_maturity = 1,
#'                      option_maturity = 1,
#'                      K = 20,
#'                      r = 0.05,
#'                      call = TRUE,
#'                      verbose = FALSE)
#'@export
European_option_value <- function(x_0, parameters, futures_maturity, option_maturity, K, r, call, verbose = FALSE){

  if("mu_star" %in% names(parameters)){
    warning("'mu_star' is deprecated. please rename this parameter to 'mu_rn'")
    parameters["mu_rn"] <- parameters["mu_star"]
    parameters <- parameters[!names(parameters) %in% "mu_star"]
  }


  ##Remove unnecessary parameters:
  parameters <- parameters[!names(parameters) %in% c("mu")]
  parameters <- parameters[!grepl("ME", names(parameters))]

  x_0 <- c(x_0, use.names = F)
  names(x_0) <- paste0("x_", 1:length(x_0), "_0")

  ###Inputs:
  parameters <- c(futures_maturity = futures_maturity, option_maturity = option_maturity, K = K, r = r, x_0, parameters)

  if(futures_maturity < option_maturity) stop("futures_maturity must be greater than option_maturity!")

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
  GBM <- !("kappa_1" %in% names(parameters))
  if(GBM) parameters["kappa_1"] <- 0

  ##Calculate Initial Values:
  ###Calculate expected Futures price:
  F_Tt <- NFCP::futures_price_forecast(x_0 = x_0, parameters = parameters, t = 0, futures_TTM = futures_maturity)

  ### Variance Calculation:
  covariance <- NFCP::cov_func(parameters, option_maturity)
  sigma_Tt <- 0
  for(i in 1:N_factors){
    for(j in 1:N_factors){
      sigma_Tt <- sigma_Tt + covariance[i,j] * exp(- (parameters[paste0("kappa_", i)] + parameters[paste0("kappa_", j)]) * (futures_maturity - option_maturity))
    }
  }
  sd <- c(sqrt(sigma_Tt), use.names = FALSE)
  if(GBM) parameters <- parameters[!names(parameters) %in% "kappa_1"]

  parameters <- c(parameters, F_Tt = F_Tt, sd = sd)

  ###Define the European option value function:
  European_option_calc <- function(X = 0, parameters, gradit = 0, call){

    if(gradit>0 && gradit<=length(parameters)) parameters[gradit] <- X
    N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
    GBM <- !("kappa_1" %in% names(parameters))
    if(GBM) parameters["kappa_1"] <- 0


    x_0 <- parameters[grepl("x_", names(parameters))]
    ## Current expected futures price:
    F_Tt <- NFCP::futures_price_forecast(x_0 = x_0, parameters = parameters, t = 0, futures_TTM = parameters["futures_maturity"])

    ### Underlying Volatility:
    covariance <- NFCP::cov_func(parameters, parameters["option_maturity"])
    sigma_Tt <- 0
    for(i in 1:N_factors){
      for(j in 1:N_factors){
        sigma_Tt <- sigma_Tt + covariance[i,j] * exp(- (parameters[paste0("kappa_", i)] + parameters[paste0("kappa_", j)]) * (parameters["futures_maturity"] - parameters["option_maturity"]))
      }
    }

    parameters["F_Tt"] <- F_Tt
    parameters["sd"] <- c(sqrt(sigma_Tt), use.names = FALSE)

    # Are we evaluating the sensitivity to the futures / underlying volatility?
    if(length(parameters) - gradit < 3) parameters[gradit] <- X

    d1 <-  (log(parameters["F_Tt"] / parameters["K"]) + (0.5 * parameters["sd"]^2)) / parameters["sd"]
    d2 <- d1 - parameters["sd"]


    if(call){
      value <- as.numeric(exp(-parameters["r"]*(parameters["option_maturity"])) * (parameters["F_Tt"] * stats::pnorm(d1) - parameters["K"] * stats::pnorm(d2)))
    } else {
      value <- as.numeric(exp(-parameters["r"]*(parameters["option_maturity"])) * (parameters["K"] * stats::pnorm(-d2) - parameters["F_Tt"] * stats::pnorm(-d1)))
    }

    return(value)
  }
  Value <- European_option_calc(parameters = parameters,call = call)
  if(!verbose){
    return(Value)
  } else {

    ###Calculate Gradients:
    TheGreeks <- rep(0, length(parameters))
    names(TheGreeks) <- names(parameters)

    for(gradit in 1:length(parameters)){
      if(futures_maturity - option_maturity < 0.02){
        if(names(parameters)[gradit] %in% c("option_maturity", "futures_maturity")){
          if(names(parameters)[gradit] %in% "option_maturity")   TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call, side = -1)
          if(names(parameters)[gradit] %in% "futures_maturity") TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call, side = 1)
        }
        else {
          TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call)
        }
      }
      else {
        TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call)
      }
    }
    names(TheGreeks)[(length(TheGreeks)-1):length(TheGreeks)] <- c("underlying futures price", "underlying volatility")


    final_output <- list(
      "option value" = Value,
      "annualized volatility" = sd/sqrt(option_maturity),
      "parameter sensitivity" = TheGreeks[!names(TheGreeks) %in% c("underlying volatility", "underlying futures price", "option_maturity", "futures_maturity", "K", "r")],
      "greeks" = TheGreeks[c("underlying volatility", "underlying futures price", "option_maturity", "futures_maturity", "K", "r")]
    )
    return(final_output)
  }
}



