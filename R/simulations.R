#' Forecast N-Factor Model Spot Prices
#' @description Analytically forecast expected spot prices following the "true" process of a given n-factor stochastic model
#'
#'@param x_0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param t a vector of discrete time points to forecast
#'@param percentiles Optional. A vector of percentiles to include probabilistic forecasting intervals.
#'
#'@details
#'Future expected spot prices under the N-factor model can be forecasted through the analytic expression of expected future prices under the "true" N-factor process.
#'
#'Given that the log of the spot price is equal to the sum of the state variables (equation 1), the spot price is log-normally distributed with the expected prices given by:
#'
#'\loadmathjax
#'\mjdeqn{E[S_t] = exp(E[ln(S_t)] + \frac{1}{2}Var[ln(S_t)])}{exp(E[ln(S[t])] + 1/2 Var[ln(S[t])])}
#'Where:
#'\mjdeqn{E[ln(S_t)] = \sum_{i=1}^Ne^{-(\kappa_it)}x_i(0) + \mu t}{E[ln(S[t])] = sum_{i=1}^N (e^(-(kappa[i] t)) x[i,0] + mu * t)}
#'
#'Where \mjeqn{\kappa_i = 0}{kappa[i] = 0} when \code{GBM=T} and \mjeqn{\mu = 0}{mu = 0} when \code{GBM = F}
#'
#'\mjdeqn{Var[ln(S_t)] =  \sigma_1^2t + \sum_{i.j\neq1}\sigma_i\sigma_j\rho_{i,j}\frac{1-e^{-(\kappa_i+\kappa_j)t}}{\kappa_i+\kappa_j}}{
#'Var[ln(S[t])] = sigma[1]^2 * t + sum_{i.j != 1} (sigma[i] sigma[j] rho[i,j] (1 - e^(-(kappa[i] + kappa[j])t)) / (kappa[i] + kappa[j]) )}
#'
#'and thus:
#'
#'\mjdeqn{E[S_t] = exp(\sum_{i=1}^N e^{-\kappa_it}x_i(0) + (\mu + \frac{1}{2}\sigma_1^2)t + \frac{1}{2}\sum_{i.j\neq1} \sigma_i\sigma_j\rho_{i,j}\frac{1-e^{-(\kappa_i+\kappa_j)t}}{\kappa_i+\kappa_j})}{
#'E[S[t]] = exp( sum_{i=1}^N e^(-kappa[i] t) x[i,0] + (mu + 1/2 sigma[1]^2)t + 1/2 (sum_{i.j != 1}( sigma[i] sigma[j] rho[i,j] (1 - e^(-(kappa[i] + kappa[j])t)) / (kappa[i] + kappa[j]))) )}
#'
#'Under the assumption that the first factor follows a Brownian Motion, in the long-run expected spot prices grow over time at a constant rate of \mjeqn{\mu + \frac{1}{2}\sigma_1^2}{mu + 1/2 sigma[1]} as the \mjeqn{e^{-\kappa_it}}{e^(-kappa[i] * t)} and \mjeqn{e^{-(\kappa_i + \kappa_j)t}}{e^(-(kappa[i] + kappa[j]))} terms approach zero.
#'
#'An important consideration when forecasting spot prices using parameters estimated through maximum likelihood estimation is that the parameter estimation process takes the assumption of risk-neutrality and thus the true process growth rate \mjeqn{\mu}{mu} is not estimated with a high level of precision. This can be shown from the higher standard error for \mjeqn{\mu}{mu} than other estimated parameters, such as the risk-neutral growth rate \mjeqn{\mu^*}{mu^*}. See Schwartz and Smith (2000) for more details.
#'
#'@return \code{spot_price_forecast} returns a vector of expected future spot prices under a given N-factor model at specified discrete future time points. When \code{percentiles} are specified, the function returns a matrix with the corresponding confidence bands in each column of the matrix.
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'# Forecast the Schwartz and Smith (2000) two-factor oil model:
#'
#'##Step 1 - Kalman filter of the two-factor oil model:
#'SS_2F_filtered <- NFCP_Kalman_filter(SS_oil$two_factor,
#'                                     names(SS_oil$two_factor),
#'                                     log(SS_oil$stitched_futures),
#'                                     SS_oil$dt,
#'                                     SS_oil$stitched_TTM,
#'                                     verbose = TRUE)
#'
#'##Step 2 - Probabilistic forecast of N-factor stochastic differential equation (SDE):
#'spot_price_forecast(x_0 = SS_2F_filtered$x_t,
#'                    parameters = SS_oil$two_factor,
#'                    t = seq(0,9,1/12),
#'                    percentiles = c(0.1, 0.9))
#'
#'@export
spot_price_forecast <- function(x_0, parameters, t, percentiles = NULL){

  ###This is the forecast of the spot price (ie. not the risk-neutral version!)

  if(is.null(names(parameters))) stop("parameters must be a named vector")

  if("mu_star" %in% names(parameters)){
    warning("'mu_star' is deprecated. please rename this parameter to 'mu_rn'")
    parameters["mu_rn"] <- parameters["mu_star"]
    parameters <- parameters[!names(parameters) %in% "mu_star"]
  }

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))

  GBM <- "mu" %in% names(parameters)

  if(GBM){
    parameters["kappa_1"] <- 0
    parameters["E"] <- 0
  } else {
    parameters["mu"] <- 0
  }

  ###Expected Spot Price Forecasting:
  ##Save the expected future spot prices:

  spot_price_forecast <- parameters["E"] + parameters["mu"] * t
  for(i in 1:N_factors) spot_price_forecast <- spot_price_forecast + x_0[i] * exp(-parameters[paste0("kappa_",i)] * (t))

  ###Instantaneous Volatility:
  var.spot_prices <- rep(0, length(t))
  for(i in 1:length(t)) var.spot_prices[i] <- sum(cov_func(parameters, t[i]))

  ###Add to Spot Price Forecast:
  spot_price_forecast <- spot_price_forecast + 0.5 * var.spot_prices

  if(!is.null(percentiles)){

    if(min(percentiles)<0 || max(percentiles) > 1) stop("percentiles are limited between zero and one")

    percentiles <- c(0.5, percentiles)[order(c(0.5, percentiles))]

    output <- matrix(spot_price_forecast, nrow = length(t), ncol = length(percentiles))

    for(i in 1:length(percentiles)) output[,i] <-  output[,i] + sqrt(var.spot_prices) * stats::qnorm(percentiles[i])

    spot_price_forecast <- output
    colnames(spot_price_forecast) <- percentiles
  }

  saved_forecasts <- as.matrix(exp(spot_price_forecast))
  rownames(saved_forecasts) <- t

  return(saved_forecasts)
}

#'Forecast N-factor Model Futures Prices
#'@description Analytically forecast future expected Futures prices under the risk-neutral version of a specified N-factor model.
#'@param x_0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param t a numeric specifying the time point at which to forecast futures prices
#'@param futures_TTM a vector specifying the time to maturity of futures contracts to value.
#'@param percentiles Optional. A vector of percentiles to include probabilistic forecasting intervals.
#'
#'
#'
#'@details
#'\loadmathjax
#'Under the assumption or risk-neutrality, futures prices are equal to the expected future spot price. Additionally, under deterministic interest rates, forward prices are equal to futures prices. Let \mjeqn{F_{T,t}}{F[T,t]}
#'denote the market price of a futures contract at time \mjeqn{t}{t} with time \mjeqn{T}{T} until maturity. let * denote the risk-neutral expectation and variance of futures prices.
#'The following equations assume that the first factor follows a Brownian Motion.
#'
#'\mjdeqn{E^*[ln(F_{T,t})] =\sum_{i=1}^Ne^{-\kappa_iT}x_{i}(0) + \mu^*t +  A(T-t)}{E^*[ln(F[T,t])] = sum_{i=1}^N (e^(-kappa[i] T) x[i,0] + mu * t + A(T-t))}
#'
#'Where:
#'\mjdeqn{A(T-t) = \mu^*(T-t)-\sum_{i=1}^N - \frac{1-e^{-\kappa_i (T-t)}\lambda_i}{\kappa_i}+\frac{1}{2}(\sigma_1^2(T-t) + \sum_{i.j\neq 1} \sigma_i \sigma_j \rho_{i,j} \frac{1-e^{-(\kappa_i+\kappa_j)(T-t)}}{\kappa_i+\kappa_j})}{
#'A(T-t) = mu^* (T-t) - sum_{i=1}^N ( - (1 - e^(-kappa[i] (T-t))lambda[i]) / kappa[i]) + 1/2 sigma[1]^2 (T-t) + sum_{i.j != 1} sigma[i] sigma[j] rho[i,j] (1 - e^(-(kappa[i] + kappa[j])(T-t))) / (kappa[i] + kappa[j])}
#'The variance is given by:
#'\mjdeqn{Var^*[ln(F_{T,t})]= \sigma_1^2t + \sum_{i.j\neq1} e^{-(\kappa_i + \kappa_j)(T-t)}\sigma_i\sigma_j\rho_{i,j}\frac{1-e^{-(\kappa_i+\kappa_j)t}}{\kappa_i+\kappa_j}}{
#'Var^*[ln(F[T,t])] = sigma[1]^2 * t + sum_{i.j != 1} e^(-(kappa[i] + kappa[j])(T-t)) sigma[i] sigma[j] rho[i,j] (1 - e^(-(kappa[i] + kappa[j])t))/(kappa[i] + kappa[j])}
#'
#'@return
#'\code{futures_price_forecast} returns a vector of expected Futures prices under a given N-factor model with specified time to maturities at time \mjeqn{t}{t}.  When \code{percentiles} are specified, the function returns a matrix with the corresponding confidence bands in each column of the matrix.
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'
#'# Forecast futures prices of the Schwartz and Smith (2000) two-factor oil model:
#'## Step 1 - Run the Kalman filter for the two-factor oil model:
#'SS_2F_filtered <- NFCP_Kalman_filter(parameter_values = SS_oil$two_factor,
#'                                     parameter_names = names(SS_oil$two_factor),
#'                                     log_futures = log(SS_oil$stitched_futures),
#'                                     dt = SS_oil$dt,
#'                                     futures_TTM = SS_oil$stitched_TTM,
#'                                     verbose = TRUE)
#'
#'## Step 2 - Probabilistic forecast of the risk-neutral two-factor
#'## stochastic differential equation (SDE):
#'futures_price_forecast(x_0 = SS_2F_filtered$x_t,
#'                       parameters = SS_oil$two_factor,
#'                       t = 0,
#'                       futures_TTM = seq(0,9,1/12),
#'                       percentiles = c(0.1, 0.9))
#'@export
futures_price_forecast <- function(x_0, parameters, t = 0, futures_TTM = 1:10, percentiles = NULL){

  if(any(futures_TTM<t)) stop("futures_TTM must be greater than t")

  ###Expected Futures Prices according to the risk-neutral equation

  if(is.null(names(parameters))) stop("parameters must be a named vector")

  if("mu_star" %in% names(parameters)){
    warning("'mu_star' is deprecated. please rename this parameter to 'mu_rn'")
    parameters["mu_rn"] <- parameters["mu_star"]
    parameters <- parameters[!names(parameters) %in% "mu_star"]
  }

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))

  GBM <- "mu_rn" %in% names(parameters)

  if(GBM){
    parameters["kappa_1"] <- 0
    parameters["E"] <- 0
  } else {
    parameters["mu_rn"] <- 0
  }

  ###Expected Futures Price Forecasting:
  ##Save the expected Futures prices:

  ###Equilibrium / Risk-neutral growth up to time t:
  futures_price_forecast <- parameters["E"] + parameters["mu_rn"]*t + rep(0, length(futures_TTM))
  ###Factors expected movement from time zero to time t:
  for(i in 1:N_factors) futures_price_forecast <- futures_price_forecast + x_0[i] * exp(-parameters[paste0("kappa_",i)] * (futures_TTM))

  ###Take risk-premiums and volatility into account:
  futures_price_forecast <- futures_price_forecast + A_T(parameters, futures_TTM-t)

  ###Instantaneous Volatility:
  var.futures_prices <- rep(0, length(futures_TTM))
  for(i in 1:length(futures_TTM)) var.futures_prices[i] <- sum(cov_func(parameters, futures_TTM[i]))

  if(!is.null(percentiles)){

    if(min(percentiles)<0 || max(percentiles) > 1) stop("percentiles are limited between zero and one")

    percentiles <- c(0.5, percentiles)[order(c(0.5, percentiles))]

    output <- matrix(futures_price_forecast, nrow = length(futures_TTM), ncol = length(percentiles))

    for(i in 1:length(percentiles)){
      output[,i] <-  output[,i] + sqrt(var.futures_prices) * stats::qnorm(percentiles[i])
    }
    futures_price_forecast <- output
    colnames(futures_price_forecast) <- percentiles
  }

  saved_forecasts <- as.matrix(exp(futures_price_forecast))
  rownames(saved_forecasts) <- futures_TTM

  return(saved_forecasts)
}


#'Simulate N-Factor Model Spot Prices
#'
#'@description Simulate risk-neutral price paths of an an N-factor commodity pricing model through Monte Carlo Simulation.
#'
#'@param x_0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param t the number of years to simulate
#'@param dt discrete time step of simulation
#'@param N_simulations total number of simulations
#'@param antithetic \code{logical}. Should antithetic price paths be simulated?
#'@param verbose \code{logical}. Should simulated state variables be output? see \bold{returns}
#'
#'@details
#'\loadmathjax
#'The \code{spot_price_simulate} function is able to quickly simulate a large number of risk-neutral price paths of a commodity following the N-factor model.
#'Simulating risk-neutral price paths of a commodity under an N-factor model through Monte Carlo simulations allows for the
#'valuation of commodity related investments and derivatives, such as American Options and Real Options through dynamic programming methods.
#'The \code{spot_price_simulate} function quickly and efficiently simulates an N-factor model over a specified number of years, simulating antithetic price paths as a simple variance reduction technique.
#'The \code{spot_price_simulate} function uses the \code{mvrnorm} function from the \code{MASS} package to draw from a multivariate normal distribution for the simulation shocks.
#'
#'The N-factor model stochastic differential equation is given by:
#'
#'Brownian Motion processes (ie. factor one when \code{GBM = T}) are simulated using the following solution:
#'
#' \mjdeqn{x_{1,t+1} = x_{1,t} + \mu^*\Delta t + \sigma_1 \Delta t Z_{t+1}}{x[1,t+1] = x[1,t] + mu^* * Delta t + sigma[1] * Delta t * Z[t+1]}
#'
#'Where \mjeqn{\Delta t}{Delta t} is the discrete time step, \mjeqn{\mu^*}{mu^*} is the risk-neutral growth rate and \mjeqn{\sigma_1}{sigma[1]} is the instantaneous volatility. \mjeqn{Z_t}{Z[t]} represents the
#'independent standard normal at time \mjeqn{t}{t}.
#'
#'Ornstein-Uhlenbeck Processes are simulated using the following solution:
#'
#'\mjdeqn{x_{i,t} = x_{i,0}e^{-\kappa_it}-\frac{\lambda_i}{\kappa_i}(1-e^{-\kappa_it})+\int_0^t\sigma_ie^{\kappa_is}dW_s}{x[i,t] = x[i,0] * e^(-kappa[i] * t) - lambda[i]/kappa[i] * (1 - e^(-kappa[i] * t)) + int_0^t (sigma[i] *
#'e^(kappa[i] * s) dW[s])}
#'
#'Where a numerical solution is obtained by numerically discretising and approximating the integral term using the Euler-Maruyama integration scheme:
#'\mjdeqn{\int_0^t\sigma_ie^{\kappa_is}dW_s = \sum_{j=0}^t \sigma_ie^{\kappa_ij}dW_s}{int_0^t ( sigma[i] e^(kappa[i] * s) dw[s])}
#'
#'@return
#'\code{spot_price_simulate} returns a list when \code{verbose = T} and a matrix of simulated price paths when \code{verbose = F}. The returned objects in the list are:
#'
#'\tabular{ll}{
#'
#' \code{State_Variables} \tab A matrix of simulated state variables for each factor is returned when \code{verbose = T}. The number of factors returned corresponds to the number of factors in the specified N-factor model. \cr
#'
#' \code{Prices} \tab A matrix of simulated price paths. Each column represents one simulated price path and each row represents one simulated observation. \cr
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
#'
#'
#'# Example 1
#'###Simulate a geometric Brownian motion (GBM) process:
#'## starting price of 20, with a growth of 5% p.a. and
#'## volatility of 20% p.a.
#'simulated_spot_prices <- spot_price_simulate(
#'  x_0 = log(20),
#'  parameters = c(mu_rn = (0.05 - (1/2) * 0.2^2), sigma_1 = 0.2),
#'  t = 1,
#'  dt = 1/12,
#'  N_simulations = 1e3)
#'
#'# Example 2
#'###Simulate future spot price paths under Risk-Neutrality and under the
#'###Schwartz - Smith two factor model:
#'
#'##Step 1 - Run the Kalman Filter for the Two-Factor Oil Model:
#'SS_2F_filtered <- NFCP_Kalman_filter(parameter_values = SS_oil$two_factor,
#'                                     parameter_names = names(SS_oil$two_factor),
#'                                     log_futures = log(SS_oil$stitched_futures),
#'                                     dt = SS_oil$dt,
#'                                     futures_TTM = SS_oil$stitched_TTM,
#'                                     verbose = TRUE)
#'
#'#Step 2 - Simulate spot prices:
#'##100 antithetic simulations of one year of monthly observations
#'simulated_spot_prices <- spot_price_simulate(
#'  x_0 = SS_2F_filtered$x_t,
#'  parameters = SS_oil$two_factor,
#'  t = 1,
#'  dt = 1/12,
#'  N_simulations = 1e3,
#'  antithetic = TRUE,
#'  verbose = TRUE)
#'@export
spot_price_simulate <- function(x_0, parameters, t = 1, dt = 1, N_simulations = 2, antithetic = TRUE, verbose = FALSE){

  if(length(t) > 1) stop("'t' must be length 1!")
  if(length(dt) > 1) stop("'dt' must be length 1!")
  if(length(N_simulations) > 1) stop("'N_simulations' must be length 1!")


  if("mu_star" %in% names(parameters)){
    warning("'mu_star' is deprecated. please rename this parameter to 'mu_rn'")
    parameters["mu_rn"] <- parameters["mu_star"]
    parameters <- parameters[!names(parameters) %in% "mu_star"]
  }

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) &
                           sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.numeric) &
                           !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))

  if("mu_rn" %in% names(parameters)) parameters["E"] <- 0

  model_type <- rep("", N_factors)
  for(i in 1:N_factors)  model_type[i] <- ifelse(paste0("kappa_",i) %in% names(parameters), "MR", "GBM")


  #Calculations:
  n <- round(N_simulations)
  N_sim <- ifelse(n %% 2 == 0 && antithetic, n, n + 1)
  nloops <- ifelse(antithetic, N_sim/2, N_sim)
  time_periods <- seq(0, t, dt)
  nsteps <- length(time_periods) - 1
  ##Correlations:
  covariance <- diag(N_factors)
  for(i in 1:N_factors) for(j in i:N_factors) if(i != j) covariance[i,j] = covariance[j,i] = parameters[paste("rho",min(i,j),max(i,j),sep="_")]

  #Correlated Brownian Process (ie. standard normal):
  shocks <- MASS::mvrnorm(n = (nsteps) * nloops, mu = rep(0, N_factors), Sigma = covariance) * sqrt(dt)

  ##Instantiate save state matrix of simulations:
  state_matrix <- array(dim = c(nsteps+1, N_sim, N_factors))

  prices <- matrix(0, nrow = nsteps+1, ncol = N_sim)
  if(verbose) output <- list()
  X <- matrix(0, nrow = nsteps+1, ncol = N_sim)
  ##Simulate the Factors:
  for(i in 1:N_factors){

    model <- model_type[i]
    x_0.i <- x_0[i]
    state_matrix[1,,i] <- x_0.i

    ###Mean-Reverting/Ornstein-Uhlenbeck:
    if(model == "MR"){

      MR_sigma  <- parameters[paste0("sigma_", i)]
      MR_kappa  <- parameters[paste0("kappa_", i)]
      MR_lambda <- parameters[paste0("lambda_", i)]

      ## Method One - Not working at the moment:

      # mu <- - MR_lambda / MR_kappa
      # chi_x <- exp(- MR_kappa * time_periods[-1])
      # ##Drift:
      # Drift <- x_0.i * chi_x + mu * (1 - chi_x)
      # #Shock:
      # ###The MR Shock is used in accordance to a cum sum:
      # ##This is the interior of the integral / summation:
      # shock_val <-  shocks[,i] * exp(MR_kappa * time_periods[-1])
      #
      # ##The cumsum is the integral / summation, the second multiplication is the shock part before it.
      # MR_shock <- MR_sigma * apply(matrix(shock_val, nrow = nsteps), MARGIN = 2, cumsum) * chi_x
      #
      # #State Matrix:
      # state_matrix[2:(nsteps+1),seq(1, N_sim, ifelse(antithetic,2,1)),i] <- Drift + MR_shock
      # #Antithetic Values:
      # if(antithetic) state_matrix[2:(nsteps+1),seq(2, N_sim, 2), i] <- Drift - MR_shock

      # Method Two:
      mu <- - MR_lambda / MR_kappa
      shock_val <-  matrix(shocks[,i], nrow = nsteps) * MR_sigma

      for(loop in 2:(nsteps+1)){
        state_matrix[loop, seq(1, N_sim, ifelse(antithetic,2,1)), i] <- state_matrix[loop-1, seq(1, N_sim, ifelse(antithetic,2,1)), i] +
          MR_kappa * (mu - state_matrix[loop-1, seq(1, N_sim, ifelse(antithetic,2,1)), i]) * dt +
          shock_val[loop-1,]

        if(antithetic) state_matrix[loop, seq(2, N_sim, ifelse(antithetic,2,1)), i] <- state_matrix[loop-1, seq(2, N_sim, ifelse(antithetic,2,1)), i] +
          MR_kappa * (mu - state_matrix[loop-1, seq(2, N_sim, ifelse(antithetic,2,1)), i]) * dt -
          shock_val[loop-1,]
      }

    }
    ###Brownian Motion:
    if(model == "GBM"){

      GBM_sigma   <- parameters[paste0("sigma_", i)]
      GBM_lambda  <- parameters[paste0("lambda_", i)]
      GBM_mu_rn <- parameters[paste0("mu_rn")]

      ##Drift:
      drift <- GBM_mu_rn * dt
      ##Shock:
      shock <- matrix(shocks[,i], nrow = nsteps) * GBM_sigma

      #Values:
      state_matrix[2:(nsteps+1), seq(1, N_sim, ifelse(antithetic,2,1)),i] <- x_0.i + apply(drift + shock, MARGIN = 2, cumsum)
      #Antithetic Values:
      if(antithetic) state_matrix[2:(nsteps+1), seq(2, N_sim, 2),i] <- x_0.i + apply(drift - shock, MARGIN = 2, cumsum)


    }
    #update log(Prices)
    prices <- prices + state_matrix[,,i]
  }
  if(N_sim != n) prices <- prices[,-ncol(prices)]

  ##Return Simulated Values:

  if(!verbose) return(exp(parameters["E"] + prices))
   else return(list(state_variables = state_matrix, spot_prices = exp(parameters["E"] + prices)))

  }


#'Simulate N-Factor Model Futures Prices
#'@description Simulate Futures price data with dynamics that follow the parameters of an N-factor model through Monte Carlo simulation.
#'
#'@param x_0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param dt discrete time step of simulation
#'@param N_obs The number of observations to simulate
#'@param futures_TTM A vector or matrix of the time to maturity of futures contracts to simulate. See \bold{details}
#'@param ME_TTM vector of maturity groupings to consider for simulated futures prices. The length of \code{ME_TTM} must be equal to the number of 'ME' parameter values.
#'@param verbose \code{logical}. Should the simulated state variables and associated prices be output?
#'
#'@details
#'\loadmathjax
#'The \code{futures_price_simulate} function simulates futures price data using the Kalman Filter algorithm, drawing from a normal
#'distribution for the shocks in the transition and measurement equations at each discrete time step. At each discrete time point,
#'an observation of the state vector is generated through the transition equation, drawing from a normal distribution with a covariance equal to \mjeqn{Q_t}{Q[t]}.
#'Following this, simulated futures prices are generated through the measurement equation, drawing from a normal distribution with covariance matrix equal to \mjeqn{H}{H}.
#'
#'Input \code{futures_TTM} can be either a matrix specifying the constant time to maturity of futures contracts to simulate, or it can be a matrix where \code{nrow(futures_TTM) == N_obs} for the time-varying time to maturity of the futures contracts to simulate. This allows for the simulation of both aggregate stitched data and individual futures contracts.
#'
#'@return
#'\code{futures_price_simulate} returns a list with three objects when \code{verbose = T} and a matrix of simulated futures prices when \code{verbose = F}. The list objects returned are:
#'
#'#'\tabular{ll}{
#'
#' \code{state_vector} \tab  A \code{matrix} of Simulated state variables at each discrete time point. The columns represent each factor of the N-factor model and the rows represent
#'the simulated values at each discrete simulated time point. \cr
#'
#' \code{futures_prices} \tab  A \code{matrix} of Simulated futures prices, with each column representing a simulated futures contract. \cr
#'
#' \code{spot_prices} \tab A vector of simulated spot prices \cr
#' }
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'
#'
#'##Example 1 - Simulate Crude Oil Stitched futures prices
#'##under a Two-Factor model, assuming a constant time to maturity:
#'
#'simulated_futures <- futures_price_simulate(x_0 = c(log(SS_oil$spot[1,1]), 0),
#'                                            parameters = SS_oil$two_factor,
#'                                            dt = SS_oil$dt,
#'                                            N_obs = nrow(SS_oil$stitched_futures),
#'                                            futures_TTM = SS_oil$stitched_TTM)
#'
#'##Example 2 - Simulate Crude Oil Contract Prices under a Two-Factor model,
#'##using a rolling-window of measurement errors:
#'
#'simulated_futures_prices <- futures_price_simulate(x_0 = c(log(SS_oil$spot[1,1]), 0),
#'                                                   parameters = SS_2F,
#'                                                   dt = SS_oil$dt,
#'                                                   N_obs = nrow(SS_oil$contracts),
#'                                                   futures_TTM = SS_oil$contract_maturities,
#'                                                   ME_TTM = c(1/4, 1/2, 1, 2, 5))
#'@export
futures_price_simulate <- function(x_0, parameters, dt, N_obs, futures_TTM, ME_TTM = NULL, verbose = TRUE){

  if(is.null(names(parameters))) stop("parameters must be a named vector. NFCP_parameters function is suggested")
  # transposeme = FALSE


  ## Temporary
  if("mu_star" %in% names(parameters)){
    warning("'mu_star' is deprecated. please rename this parameter to 'mu_rn'")
    parameters["mu_rn"] <- parameters["mu_star"]
    parameters <- parameters[!names(parameters) %in% "mu_star"]
  }

  ##If it's a constant futures_TTM, develop the maturity matrix:
  futures_TTM <- as.matrix(futures_TTM)
  ##If a maturity matrix has been provided:
  if(all(dim(futures_TTM) > 1)){
    maturity_matrix <- futures_TTM
    if(dim(maturity_matrix)[1]!=N_obs) stop("futures_TTM does not match number of observations")
  }

    GBM <- "mu" %in% names(parameters)

    if(GBM){
      parameters["kappa_1"] <- 0
      parameters["E"] <- 0
    } else {
      parameters["mu"] <- 0
    }

    ###The highest sigma input would be our number of factors
    N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.numeric) & !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))
    N_contracts <- ncol(maturity_matrix)
    ## How many measurement errors are specified?
    N_ME <- max(which(paste0("ME_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("ME_",1:length(parameters))], FUN = is.numeric) & !sapply(parameters[paste0("ME_",1:length(parameters))], FUN = is.na)))

    ## Have enough ME maturity terms been specified?
    if(N_ME > 1 && N_ME < N_contracts){
      if(is.null(ME_TTM)) stop("Multiple measurement error (ME) terms have been specified but the maturity terms of the measurement error (ME_TTM) have not.")
      if(length(ME_TTM) != N_ME) stop("Number of measurement error (ME) terms specified does not match the length of argument 'ME_TTM'")
    }

    #Initialize State Variable
    if(length(x_0) != N_factors) stop("Initial state vector length not equal to N factors")

    # Step one - evaluate the measurement error matrix:
    homogeneous_H <- TRUE

    #Measurement Errors matrix:
    H <- matrix(rep(0, N_contracts))
    ## Case 1 - Only one ME has been specified:
    if(N_ME == 1) H[1:N_contracts] <- parameters["ME_1"]^2
    ## Case 2 - an ME for every contract has been specified:
    if(N_ME == N_contracts) H[1:N_contracts] <- sapply(1:N_ME, FUN = function(x) parameters[paste0("ME_",x)]^2)
    ## Case 3 - Multiple ME have been specified, corresponding to particular maturity groupings:
    if(!is.null(ME_TTM) && N_ME > 1 && N_ME < N_contracts){

      # H is now time inhomogenous:
      homogeneous_H <- FALSE

      H <- matrix(NA, nrow = N_obs, ncol = N_contracts)

      for(loop in length(ME_TTM):1) H[futures_TTM < ME_TTM[loop]] = parameters[paste0("ME_",loop)]^2

      H <- t(H)
    }
    H[H < 1.01e-10] <- 0
    H[is.na(H)] <- 0

    X_t <- matrix(x_0)

    ##Discrete Time Steps:

    #Insantiate State Vectors ie. the outputs:
    X <- matrix(0  , nrow = N_obs, ncol = N_factors)
    Y <- matrix(0  , nrow = N_obs, ncol = N_contracts)

    d <- A_T(parameters, maturity_matrix)

    Z <- array(NA, dim = c(N_obs, N_contracts, N_factors))
    Z[,,1:N_factors] <- sapply(1:N_factors, FUN = function(X) exp(- parameters[paste0("kappa_", X)] * maturity_matrix))

    #Multivariate Normal observation - Transition Equation:
    omega <- MASS::mvrnorm(n = N_obs, mu = rep(0, N_factors), Sigma = cov_func(parameters, dt))


    for(t in 1:N_obs){

      # H_t <- diag(H[,ifelse(homogeneous_H, 1, t)])

      #Multivariate Normal observation - Measurement Equation:
      # V <- MASS::mvrnorm(n = N_obs, mu = rep(0, N_contracts), Sigma = diag(H[,ifelse(homogeneous_H, 1, t)]))
      v_t <- MASS::mvrnorm(n = 1, mu = rep(0, N_contracts), Sigma = diag(H[,ifelse(homogeneous_H, 1, t)]))

      #Measurement Equation:
      #y = d + Z * x_t + v_t
      #d
      d_t <- matrix(d[t,])

      #Z:
      Z_t <- Z[t,,]

      #Simulated Measurement Error:
      # v_t <- V[t,]

      Y_t <- parameters["E"] + d_t + Z_t %*% X_t + v_t

      ##Record Observations:
      X[t,] <- X_t
      Y[t,] <- Y_t

      if(t != N_obs){

        #Omega - Transition Equation:
        omega_t <- omega[t,]
        #Transition Equation:
        #x_t = c + G * x_tm1 + omega
        c_t <- matrix(c(parameters["mu"] * dt, rep(0, N_factors-1)))
        #G:
        G_t <- diag(N_factors)
        diag(G_t) <- sapply(1:N_factors, FUN = function(X) exp(- parameters[paste0("kappa_", X)] * dt))

        X_t <- c_t + G_t %*% X_t + omega_t

      }
    }
    X_output <- X
    Y_output <- exp(Y)
    spot_prices <- as.matrix(exp(rowSums(parameters["E"] + X)))
    rownames(Y_output) <- rownames(spot_prices) <- rownames(X_output) <- rownames(maturity_matrix)

    colnames(spot_prices) <- "Spot"

    if(any(dim(futures_TTM)==1)) colnames(Y_output) <- as.character(round(futures_TTM,4))
    colnames(X_output) <- paste("Factor", 1:N_factors)

    if(verbose){
      return(list(state_vector = X_output, futures_prices = Y_output, spot_prices = spot_prices))
    } else {
      return(Y_output)
    }

  }
