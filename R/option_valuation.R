#'N-factor model American options on futures contracts valuation
#'
#'@description Value American options on futures contracts under the parameters of an N-factor model
#'
#'@param x_0 \code{vector}. Initial values of the state variables, where the length must correspond to the number of factors specified in the parameters.
#'@param parameters \code{vector}. A named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param futures_maturity \code{numeric}. Time, in years, when the underlying futures contract matures.
#'@param option_maturity \code{numeric}. Time, in years,  when the American option expires.
#'@param K \code{numeric}. Strike price of the American Option.
#'@param r \code{numeric}. Annualized risk-free interest rate.
#'@param call \code{logical}. Is the American option a call or put option?
#'@param N_simulations \code{numeric}. Total number of simulated price paths.
#'@param dt \code{numeric}. Discrete time step, in years, of the Monte Carlo simulation.
#'@param orthogonal \code{character}. The orthogonal polynomial used to approximate the continuation value of the option in the LSM simulation method.
#' Orthogonal polynomial arguments available are: "Power", "Laguerre", "Jacobi", "Legendre", "Chebyshev", "Hermite".
#'@param degree \code{numeric}. The degree of polynomials used in the least squares fit.
#'@param verbose \code{logical}. Should additional option value information be output? see \bold{details}.
#'@param debugging \code{logical} Should the simulated state variables and futures prices be output?
#'
#'@details
#'
#'The \code{American_option_value} function calculates numerically the value of American options on futures contracts within the N-factor model. An American option on a commodity futures contract gives the holder
#'the right, but not the obligation, to buy (call) or sell (put) the underlying asset at any time before option maturity. If the American option is exercised, the option devolves into buying or selling of the underlying
#'futures asset at the exercise price.
#'
#'The 'American_option_value' function uses Monte Carlo simulation and the Least-Squares Monte Carlo (LSM) simulation approach to numerically calculate the value of American options on futures contracts under the N-factor model. LSM simulation is a method
#'that values options with early exercise opportunities, first presented by Longstaff and Schwartz (2001). LSM simulation uses discrete time steps to approximate the value of the American option and thus technically values Bermudan-style options,
#'converging to American option values as the size of the time step approaches zero. For more information on LSM simulation, see help('LSM_American_option') of the 'LSMRealOption' package or Longstaff and Schwartz (2001).
#'
#'For a provided N-factor model,the 'American_option_value' function simulates state variables under the N-factor framework through the 'spot_price_simulate' function, developing expected futures prices from
#'these simulated state variables. The function then uses
#' the 'LSM_American_option' of the 'LSMRealOption' package to calculate the value of the American option with early exercise opportunities.
#'
#'The number of simulations has a large influence on the standard error and accuracy of calculated option values at the cost of computational expense. Large numbers of simulations are suggested to converge upon appropriate values.
#'
#'Orthogonal polynomials are used in the LSM simulation method to approximate the value of continuing to hold the American option. In general, increasing the degree of orthogonal polynomials used should increase the accuracy of results, at the cost of increased computational expense.
#'
#'
#'@return
#'
#'The 'American_option_value' function by default returns a \code{numeric} object corresponding to the calculated value of the American option.
#'
#'When \code{verbose = T}, 6 objects related to the American option value are returned within a \code{list} class object. The objects returned are:
#'
#'\tabular{ll}{
#'
#' \code{Value} \tab The calculated option value. \cr
#' \code{Standard Error} \tab The standard error of the calculated option value. \cr
#' \code{Expected Timing} \tab The expected time of early exercise. \cr
#' \code{Expected Timing SE} \tab The standard error of the expected time of early exercise. \cr
#' \code{Exercise Probability} \tab The probability of early exercise of the option being exercised. \cr
#' \code{Cumulative Exercise Probability} \tab \code{vector}. The cumulative probability of option exercise at each discrete observation point. \cr
#' }
#'
#'When \code{debugging = T}, an additional 2 simulation objects are returned within the \code{list} class object. These objects can have high dimensions and thus memory usage, so caution should be applied. The objects returned are:
#'
#' \tabular{ll}{
#' \code{Simulated State Variables} \tab An array of simulated state variables. The three dimensions of the array correspond to a discrete time observation, simulated price path, and factor of the N-factor model, respectively. \cr
#' \code{Simulated Futures Prices} \tab A matrix of simulated futures contract price paths. Each row represents one simulated discrete time observation and each column represents one simulated price path \cr
#' }
#'
#'@references
#'
#'Longstaff, F.A., and E.S. Schwartz, (2001). Valuing American Options by Simulation: A Simple Least-Squares Approach. \emph{The Review of Financial Studies.}, 14:113-147.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Aspinall, T., A. Gepp, G. Harris, S. Kelly, C. Southam, and B. Vanstone, (2021). LSMRealOptions: Value
#'American and Real Options Through LSM Simulation. R package version 0.1.1.
#'
#'@examples
#'
#'# Example 1 - An American put option on a futures contract following 'GBM'
#'American_option_value(x_0 = log(36),
#'                      parameters = c(mu_rn = 0.06, sigma_1 = 0.2),
#'                      N_simulations = 1e2,
#'                      futures_maturity = 1,
#'                      option_maturity = 1,
#'                      dt = 1/50,
#'                      K = 40,
#'                      r = 0.06,
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
#'                                         dt = SS_oil$dt,
#'                                         futures_TTM = SS_oil$stitched_TTM,
#'                                         verbose = TRUE)
#'
#'##Step 2 - Calculate 'put' option price:
#'American_option_value(x_0 = Schwartz_Smith_oil$x_t,
#'                      parameters = SS_oil$two_factor,
#'                      futures_maturity = 2,
#'                      option_maturity = 1,
#'                      K = 20,
#'                      r = 0.05,
#'                      call = FALSE,
#'                      N_simulations = 1e2,
#'                      dt = 1/12,
#'                      verbose = TRUE,
#'                      orthogonal = "Power",
#'                      degree = 2)
#'@export
American_option_value <- function(x_0, parameters, futures_maturity, option_maturity, K, r, call = FALSE, N_simulations, dt, orthogonal = "Power", degree = 2, verbose = FALSE, debugging = FALSE){

  ## How many factors are specified?
  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("sigma_",1:length(parameters))],
                                                                                                FUN = is.numeric) & !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))
  ## Futures maturity:
  futures_TTM <- seq(futures_maturity, futures_maturity - option_maturity, -dt)

  ## How many seasonality factors are specified?
  N_season <- length(grep("season", names(parameters)))/2
  ## Incorporate deterministic seasonality (if there is any):
  seasonality <- 0
  if(N_season > 0){
    if(is.na(x_0["seasonal_trend"])) warning("Deterministic seasonality considered, however a 'seasonal_trend' element not specified within 'x_0'. Starting point of 0 along the seasonal curve is assumed.")
    seas_trend <- ifelse(is.na(x_0["seasonal_trend"]), 0, x_0["seasonal_trend"])
    seasonal_offset <- futures_TTM - floor(futures_TTM) + seas_trend
    for(i in 1:N_season) seasonality <- seasonality + parameters[paste0("season_", i, "_1")] * cos(2 * i * pi * seasonal_offset) + parameters[paste0("season_", i, "_2")] * sin(2 * i * pi * seasonal_offset)
  }
  # if(N_season > 0) for(i in 1:N_season) seasonality <- seasonality + parameters[paste0("season_", i, "_1")] * cos(2 * i * pi * futures_TTM) + parameters[paste0("season_", i, "_2")] * sin(2 * i * pi * futures_TTM)

  ## From spot to expected futures:
  futures_adjustment <- seasonality + A_T(parameters, futures_TTM)

  ## Step 1 - Simulate the N-factor model:
  suppressWarnings(simulated_states <- spot_price_simulate(x_0 = x_0, parameters = parameters, t = option_maturity, dt = dt, N_simulations = N_simulations, verbose = TRUE)$state_variables)

  # Step 2 - Adjust from simulated state variables to expected futures prices:
  GBM <- "mu_rn" %in% names(parameters)
  if(GBM) parameters[c("kappa_1", "E")] <- 0
  simulated_futures <- simulated_states[,,1] * exp(-parameters["kappa_1"] * futures_TTM) + futures_adjustment
  if(N_factors > 1) for(i in 2:N_factors) simulated_futures <- simulated_futures + simulated_states[,,i] * exp(-parameters[paste0("kappa_", i)] * futures_TTM)
  simulated_futures <- exp(parameters["E"] + simulated_futures)

  ## Step 2 - Apply the 'LSM_american_option' function from the 'LSMRealOptions' to value the American option:
  output <- LSMRealOptions::LSM_american_option(state_variables = simulated_futures,
                                                payoff = simulated_futures,
                                                K = K,
                                                dt = dt,
                                                rf = r,
                                                call = call,
                                                orthogonal = orthogonal,
                                                degree = degree,
                                                verbose = verbose)
  if(debugging){
    return(c(output, list("Simulated State Variables" = simulated_states, "Simulated Futures Prices" = simulated_futures)))
  } else {
    return(output)
  }
}


#'N-factor model European options on futures contracts valuation
#'@description Value European Option Put and Calls under the parameters of an N-factor model.
#'
#'@param x_0 \code{vector}. Initial values of the state variables, where the length must correspond to the number of factors specified in the parameters.
#'@param parameters \code{vector}. A named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param futures_maturity \code{numeric}. Time, in years, when the underlying futures contract matures.
#'@param option_maturity \code{numeric}. Time, in years,  when the American option expires.
#'@param K \code{numeric}. Strike price of the American Option.
#'@param r \code{numeric}. Annualized risk-free interest rate.
#'@param call \code{logical}. Is the American option a call or put option?
#'@param verbose \code{logical}. Should additional option value information be output? see \bold{details}.
#'
#'@details
#'
#'\loadmathjax
#'
#'The \code{European_option_value} function calculates analytic expressions of the value of European call and put options on futures contracts within the N-factor model. A European option on a commodity futures contract gives the holder
#'the right, but not the obligation, to buy (call) or sell (put) the underlying asset at option maturity. If the European option is exercised, the option devolves into buying or selling of the underlying futures asset.
#'
#'State variables (i.e., the states of the factors of an N-factor model) are generally unobservable. Filtering the commodity pricing model using term structure data will provide the most recent optimal estimates of state variables, which can then be used to forecast and value European options.
#'
#'Under the assumption that future futures prices
#'are log-normally distributed under the risk-neutral process, there exist analytic expressions of the value of European call and put options on futures contracts. The value of a European option on a futures contract
#'is given by calculating the current expected futures price and the average instantaneous variance of the futures return innovations over the life of the option.
#'
#'Consider a European option with strike price \mjeqn{K}{K} and a risk-free interest rate of \mjeqn{r_f}{r_f}. The option maturity is at time \mjeqn{T_0}{T_0} and futures maturity at time \mjeqn{T_1}{T_1}.
#'The particular model features a state vector of length \mjeqn{N}{N} (i.e., N-factors) \mjeqn{x(t)}{x(t)}
#'
#'The value of a European call option would thus be:
#'
#'\mjdeqn{e^{-r T_0} . E^*[max(F(x(T_0),T_0,T_1) - K, 0)]}{e^(-r T_0) . E^*[max(F(x(T_0),T_0,T_1) - K, 0)]}
#'
#'The analytic solution to call and put options are given by:
#'
#'Call options:
#'\mjdeqn{e^{-r T_0}(F(x(0), 0, T_1) . N(d_1) - K.N(d_2))}{e^(-r T_0) (F(x(0), 0, T_1) . N(d_1) - K.N(d_2))}
#'
#'Put options:
#'\mjdeqn{e^{-r T_0}(K.N(-d_2) - F(x(0), 0, T_1) . N(-d_1))}{e^(-r T_0) (K.N(-d_2) - F(x(0), 0, T_1) . N(-d_1))}
#'
#'Where:
#'
#'Where:
#'\mjdeqn{d_1 = \frac{\ln(F/K) + \frac{1}{2} v^2}{v}}{d_1 = (ln(F/K) + 1/2 v^2) / v}
#'
#'\mjdeqn{d_2 = d_1 - v}{d_2 = d_1 - v}
#'
#'Parameter \mjeqn{ N(d) }{N(d)} indicates cumulative probabilities for the standard normal distribution (i.e. \mjeqn{P(Z<d)}{P(Z<d)}).
#'
#'Finally, parameter \mjeqn{v}{v}, the annualized option volatility, is given by:
#'
#'\mjdeqn{Var^*[\ln(F(x(T_0), T_0, T_1))] \equiv v^2 = \sum_{i.j=1} e^{(-\kappa_i + \kappa_j)(T_1 - T_0)}Cov^*(x_i(T_0), x_j(T_0))}{Var^*[ln(F(x(T_0), T_0, T_0))] = v^2 = sum_(i.j=1) e^((-kappa[i] + kappa[j])(T_1 - T_0)) Cov(x_i(T_0), x_j(T_0))}
#'
#'The annualized option volatility approaches \mjeqn{\sigma_1^2 T_0}{sigma[1]^2 T_0} as both \mjeqn{T_0}{T_0} and \mjeqn{T_1}{T_1} increase, as most uncertainty about spot prices at futures contract maturity and option expiration
#'are a result of uncertainty about spot prices, rather than the cost of carry (Schwartz and Smith, 2000).
#'
#'The presented option valuation formulas are analogous to the Black-Scholes formulas for valuing European options on stocks that do not pay dividends
#'
#'When \code{verbose = T}, the \code{European_option_value} function numerically calculates the sensitivity of option prices to underlying option and model parameters. Gradients are calculated numerically through the
#' \code{grad} function of the \code{numDeriv} package.
#'
#'@return
#'The \code{European_option_value} function returns a numeric value corresponding to the present value of an option when \code{verbose = F}.
#'When \code{verbose = T}, \code{European_option_value} returns a list with three objects:
#'
#'\tabular{ll}{
#'
#'\code{option value} \tab Present value of the option. \cr
#'
#'\code{annualized volatility} \tab Annualized volatility of the option. \cr
#'
#'\code{parameter sensitivity} \tab Sensitivity of option value to model parameters. \cr
#'
#'\code{greeks} \tab Sensitivity of option value to option parameters. \cr
#'
#' }
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Paul Gilbert and Ravi Varadhan (2016). numDeriv: Accurate Numerical Derivatives. R package version 2016.8-1. https://CRAN.R-project.org/package=numDeriv
#'
#'@examples
#'##Example 1 - A European 'put' option on a futures contract following 'GBM'
#'
#'European_option_value(x_0 = log(20), parameters = c(mu_rn = 0.06, sigma_1 = 0.2),
#'                      futures_maturity = 1, option_maturity = 1,
#'                      K = 20, r = 0.06, call = FALSE, verbose = TRUE)
#'
#'##Example 2 - A European put option under a two-factor crude oil model:
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
#'                      futures_maturity = 2,
#'                      option_maturity = 1,
#'                      K = 20,
#'                      r = 0.05,
#'                      call = FALSE,
#'                      verbose = FALSE)
#'@export
European_option_value <- function(x_0, parameters, futures_maturity, option_maturity, K, r, call = FALSE, verbose = FALSE){

  if(is.null(names(parameters))) stop("parameters must be a named vector. NFCP_parameters function is suggested")

  ##Remove unnecessary parameters:
  parameters <- parameters[!names(parameters) %in% c("mu")]
  parameters <- parameters[!grepl("ME", names(parameters))]

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
  N_season <- length(grep("season", names(parameters)))/2

  seas_trend <- ifelse(is.na(x_0["seasonal_trend"]), 0, x_0["seasonal_trend"])

  x_0 <- c(x_0, use.names = F)
  names(x_0) <- paste0("x_", 1:length(x_0), "_0")
  if(N_season > 0) x_0 <- c(x_0, seasonal_trend = c(seas_trend, use.names = FALSE))

  ###Inputs:
  parameters <- c(futures_maturity = futures_maturity, option_maturity = option_maturity, K = K, r = r, x_0, parameters)

  if(futures_maturity < option_maturity) stop("futures_maturity must be greater than option_maturity!")

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

    seas_trend <- ifelse(is.na(parameters["seasonal_trend"]), 0, parameters["seasonal_trend"])

    x_0 <- parameters[grepl("x_", names(parameters))]
    if(N_season > 0 ) x_0 <- c(x_0, seasonal_trend = c(seas_trend, use.names = FALSE))


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

    params <- parameters[!names(parameters) %in% "seasonal_trend"]

    ###Calculate Gradients:
    TheGreeks <- rep(0, length(params))
    names(TheGreeks) <- names(params)

    for(gradit in 1:length(params)){

      if(futures_maturity - option_maturity < 0.02){
        if(names(params)[gradit] %in% c("option_maturity", "futures_maturity")){
          if(names(params)[gradit] %in% "option_maturity")   TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = params[gradit], gradit = gradit, call = call, side = -1)
          if(names(params)[gradit] %in% "futures_maturity") TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = params[gradit], gradit = gradit, call = call, side = 1)
        }
        else {
          TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = params[gradit], gradit = gradit, call = call)
        }
      }
      else {
        TheGreeks[gradit] <- numDeriv::grad(func = European_option_calc, parameters = parameters, x = params[gradit], gradit = gradit, call = call)
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



