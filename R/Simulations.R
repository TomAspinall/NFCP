#' Forecast N-Factor Model Spot Prices
#' @description Analytically forecast expected spot prices following the "true" process of a given n-factor stochastic model
#'
#'@param X.0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param t a vector of discrete time points to forecast
#'@param Percentiles Optional. A vector of percentiles to include probabilistic forecasting intervals.
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
#'@return \code{Spot.Price.Forecast} returns a vector of expected future spot prices under a given N-factor model at specified discrete future time points. When \code{percentiles} are specified, the function returns a matrix with the corresponding confidence bands in each column of the matrix.
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'# Forecast the Schwartz and Smith (2000) two-factor oil model:
#'
#'##Step 1 - Run the Kalman filter for the two-factor oil model:
#'Schwartz.Smith.Oil <- NFCP.Kalman.filter(SS.Oil$Two.Factor,
#'                                       names(SS.Oil$Two.Factor),
#'                                       log(SS.Oil$Stitched.Futures),
#'                                       SS.Oil$dt,
#'                                       SS.Oil$Stitched.TTM,
#'                                       verbose = TRUE)
#'
#'##Step 2 - Probabilistic forecast of n-factor stochastic differential equation (SDE):
#'E.Spot <- Spot.Price.Forecast(X.0 = Schwartz.Smith.Oil$X.t,
#'                                  parameters = SS.Oil$Two.Factor,
#'                                  t = seq(0,9,1/12),
#'                                  Percentiles = c(0.1, 0.9))
#'@export
Spot.Price.Forecast <- function(X.0, parameters, t, Percentiles = NULL){

  ###This is the forecast of the spot price (ie. not the risk-neutral version!)

  if(is.null(names(parameters))) stop("parameters must be a named vector")

  N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))

  GBM <- "mu" %in% names(parameters)

  if(GBM){
    parameters["kappa_1"] <- 0
    parameters["E"] <- 0
  } else {
    parameters["mu"] <- 0
  }

  ###Expected Spot Price Forecasting:
  ##Save the expected future spot prices:

  Spot.Price.Forecast <- parameters["E"] + parameters["mu"] * t
  for(i in 1:N.factors) Spot.Price.Forecast <- Spot.Price.Forecast + X.0[i] * exp(-parameters[paste0("kappa_",i)] * (t))

  ###Instantaneous Volatility:
  var.Spot <- rep(0, length(t))
  for(i in 1:length(t)) var.Spot[i] <- sum(cov_func(parameters, t[i]))

  ###Add to Spot Price Forecast:
  Spot.Price.Forecast <- Spot.Price.Forecast + 0.5 * var.Spot

  if(!is.null(Percentiles)){

    if(min(Percentiles)<0 || max(Percentiles) > 1) stop("Percentiles are limited between zero and one")

    Percentiles <- c(0.5, Percentiles)[order(c(0.5, Percentiles))]

    output <- matrix(Spot.Price.Forecast, nrow = length(t), ncol = length(Percentiles))

    for(i in 1:length(Percentiles)) output[,i] <-  output[,i] + sqrt(var.Spot) * stats::qnorm(Percentiles[i])

    Spot.Price.Forecast <- output
    colnames(Spot.Price.Forecast) <- Percentiles
  }

  saved_forecasts <- as.matrix(exp(Spot.Price.Forecast))
  rownames(saved_forecasts) <- t

  return(saved_forecasts)
}

#'Forecast N-factor Model Futures Prices
#'@description Analytically forecast future expected Futures prices under the risk-neutral version of a specified N-factor model.
#'@param X.0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param t a numeric specifying the time point at which to forecast futures prices
#'@param TTM a vector specifying the time to maturity of futures contracts to value.
#'@param Percentiles Optional. A vector of percentiles to include probabilistic forecasting intervals.
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
#'\code{Futures.Price.Forecast} returns a vector of expected Futures prices under a given N-factor model with specified time to maturities at time \mjeqn{t}{t}.  When \code{percentiles} are specified, the function returns a matrix with the corresponding confidence bands in each column of the matrix.
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'# Forecast futures prices of the Schwartz and Smith (2000) two-factor oil model:
#'## Step 1 - Run the Kalman filter for the two-factor oil model:
#'Schwartz.Smith.Oil = NFCP.Kalman.filter(parameter.values = SS.Oil$Two.Factor,
#'                                       parameters = names(SS.Oil$Two.Factor),
#'                                       log.futures = log(SS.Oil$Stitched.Futures),
#'                                       dt = SS.Oil$dt,
#'                                       TTM = SS.Oil$Stitched.TTM,
#'                                       verbose = TRUE)
#'
#'## Step 2 - Probabilistic forecast of the risk-neutral two-factor
#'## stochastic differential equation (SDE):
#'E.Futures = Futures.Price.Forecast(X.0 = Schwartz.Smith.Oil$X.t,
#'                                        parameters = SS.Oil$Two.Factor,
#'                                        t = 0,
#'                                        TTM = seq(0,9,1/12),
#'                                        Percentiles = c(0.1, 0.9))
#'@export
Futures.Price.Forecast <- function(X.0, parameters, t = 0, TTM = 1:10, Percentiles = NULL){

  if(any(TTM<t)) stop("TTM must be greater than t")

  ###Expected Futures Prices according to the risk-neutral equation

  if(is.null(names(parameters))) stop("parameters must be a named vector")

  N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))

  GBM <- "mu_star" %in% names(parameters)

  if(GBM){
    parameters["kappa_1"] <- 0
    parameters["E"] <- 0
  } else {
    parameters["mu_star"] <- 0
  }

  ###Expected Futures Price Forecasting:
  ##Save the expected Futures prices:

  ###Equilibrium / Risk-neutral growth up to time t:
  Futures.Price.Forecast <- parameters["E"] + parameters["mu_star"]*t + rep(0, length(TTM))
  ###Factors expected movement from time zero to time t:
  for(i in 1:N.factors) Futures.Price.Forecast <- Futures.Price.Forecast + X.0[i] * exp(-parameters[paste0("kappa_",i)] * (TTM))

  ###Take risk-premiums and volatility into account:
  Futures.Price.Forecast <- Futures.Price.Forecast + A_T(parameters, TTM-t)

  ###Instantaneous Volatility:
  var.Futures <- rep(0, length(TTM))
  for(i in 1:length(TTM)) var.Futures[i] <- sum(cov_func(parameters, TTM[i]))

  if(!is.null(Percentiles)){

    if(min(Percentiles)<0 || max(Percentiles) > 1) stop("Percentiles are limited between zero and one")

    Percentiles <- c(0.5, Percentiles)[order(c(0.5, Percentiles))]

    output <- matrix(Futures.Price.Forecast, nrow = length(TTM), ncol = length(Percentiles))

    for(i in 1:length(Percentiles)){
      output[,i] <-  output[,i] + sqrt(var.Futures) * stats::qnorm(Percentiles[i])
    }
    Futures.Price.Forecast <- output
    colnames(Futures.Price.Forecast) <- Percentiles
  }

  saved_forecasts <- as.matrix(exp(Futures.Price.Forecast))
  rownames(saved_forecasts) <- TTM

  return(saved_forecasts)
}


#'Simulate N-Factor Model Spot Prices
#'
#'@description Simulate risk-neutral price paths of an an N-factor commodity pricing model through Monte Carlo Simulation.
#'
#'@param X.0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param t the number of years to simulate
#'@param dt discrete time step of simulation
#'@param n total number of simulations
#'@param antithetic \code{logical}. Should antithetic price paths be simulated?
#'@param verbose \code{logical}. Should simulated state variables be output? see \bold{returns}
#'
#'@details
#'\loadmathjax
#'The \code{Spot.Price.Simulate} function is able to quickly simulate a large number of risk-neutral price paths of a commodity following the N-factor model.
#'Simulating risk-neutral price paths of a commodity under an N-factor model through Monte Carlo simulations allows for the
#'valuation of commodity related investments and derivatives, such as American Options and Real Options through dynamic programming methods.
#'The \code{Spot.Price.Simulate} function quickly and efficiently simulates an N-factor model over a specified number of years, simulating antithetic price paths as a simple variance reduction technique.
#'The \code{Spot.Price.Simulate} function uses the \code{mvrnorm} function from the \code{MASS} package to draw from a multivariate normal distribution for the simulation shocks.
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
#'\code{Spot.Price.Simulate} returns a list when \code{verbose = T} and a matrix of simulated price paths when \code{verbose = F}. The returned objects in the list are:
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
#'# Example 1
#'###Simulate a Geometric Brownian Motion (GBM) process:
#'## Starting price of 20, with a growth of 5% p.a. and
#'## volatility of 20% p.a.
#'Simulated.Spot.Prices <- Spot.Price.Simulate(
#'  X.0 = log(20),
#'  parameters = c(mu_star = (0.05 - (1/2) * 0.2^2), sigma_1 = 0.2),
#'  t = 1,
#'  dt = 1/12,
#'  n = 1e3)
#'
#'# Example 2
#'###Simulate future spot price paths under Risk-Neutrality and under the
#'###Schwartz - Smith two factor model:
#'
#'##Step 1 - Run the Kalman Filter for the Two-Factor Oil Model:
#'Schwartz.Smith.Oil <- NFCP.Kalman.filter(parameter.values = SS.Oil$Two.Factor,
#'                                       parameters = names(SS.Oil$Two.Factor),
#'                                       log.futures = log(SS.Oil$Stitched.Futures),
#'                                       dt = SS.Oil$dt,
#'                                       TTM = SS.Oil$Stitched.TTM,
#'                                       verbose = TRUE)
#'
#'#Step 2 - Simulate spot prices:
#'##100 antithetic simulations of one year of monthly observations
#'Simulated.Spot.Prices <- Spot.Price.Simulate(
#'  X.0 = Schwartz.Smith.Oil$X.t,
#'  parameters = SS.Oil$Two.Factor,
#'  t = 1,
#'  dt = 1/12,
#'  n = 1e3,
#'  antithetic = TRUE,
#'  verbose = TRUE)
#'
#'@export
Spot.Price.Simulate <- function(X.0, parameters, t = 1, dt = 1, n = 2, antithetic = TRUE, verbose = FALSE){

  if(length(t) > 1) stop("'t' must be length 1!")
  if(length(dt) > 1) stop("'dt' must be length 1!")
  if(length(n) > 1) stop("'n' must be length 1!")

  N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.numeric) & !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))

  if("mu_star" %in% names(parameters)) parameters["E"] <- 0

  model.type <- rep("", N.factors)
  for(i in 1:N.factors)  model.type[i] <- ifelse(paste0("kappa_",i) %in% names(parameters), "MR", "GBM")


  #Calculations:
  n <- round(n)
  N_sim <- ifelse(n %% 2 == 0 && antithetic, n, n + 1)
  nloops <- ifelse(antithetic, N_sim/2, N_sim)
  time_periods <- seq(0, t, dt)
  nsteps <- length(time_periods) - 1
  ##Correlations:
  covariance <- diag(N.factors)
  for(i in 1:N.factors) for(j in i:N.factors) if(i != j) covariance[i,j] = covariance[j,i] = parameters[paste("rho",min(i,j),max(i,j),sep="_")]

  #Correlated Brownian Process (ie. standard normal):
  shocks <- MASS::mvrnorm(n = (nsteps) * nloops, mu = rep(0, N.factors), Sigma = covariance) * sqrt(dt)

  ##Instantiate save state matrix of simulations:
  State.Matrix <- array(dim = c(nsteps+1, N_sim, N.factors))

  prices <- matrix(0, nrow = nsteps+1, ncol = N_sim)
  if(verbose) output <- list()
  X <- matrix(0, nrow = nsteps+1, ncol = N_sim)
  ##Simulate the Factors:
  for(i in 1:N.factors){

    model <- model.type[i]
    X.0.i <- X.0[i]
    State.Matrix[1,,i] <- X.0.i

    ###Mean-Reverting/Ornstein-Uhlenbeck:
    if(model == "MR"){

      MR_sigma  <- parameters[paste0("sigma_", i)]
      MR_kappa  <- parameters[paste0("kappa_", i)]
      MR_lambda <- parameters[paste0("lambda_", i)]

      mu <- - MR_lambda / MR_kappa
      chi_x <- exp(- MR_kappa * time_periods[-1])
      ##Drift:
      Drift <- X.0.i * chi_x + mu * (1 - chi_x)
      #Shock:
      ###The MR Shock is used in accordance to a cum sum:
      ##This is the interior of the integral / summation:
      shock_val <-  shocks[,i] * exp(MR_kappa * time_periods[-1])

      ##The cumsum is the integral / summation, the second multiplication is the shock part before it.
      MR_shock <- MR_sigma * apply(matrix(shock_val, nrow = nsteps), MARGIN = 2, cumsum) * chi_x

      #State Matrix:
      State.Matrix[2:(nsteps+1),seq(1, N_sim, ifelse(antithetic,2,1)),i] <- Drift + MR_shock
      #Antithetic Values:
      if(antithetic) State.Matrix[2:(nsteps+1),seq(2, N_sim, 2), i] <- Drift - MR_shock

    }
    ###Brownian Motion:
    if(model == "GBM"){

      GBM_sigma   <- parameters[paste0("sigma_", i)]
      GBM_lambda  <- parameters[paste0("lambda_", i)]
      GBM_mu_star <- parameters[paste0("mu_star")]

      ##Drift:
      drift <- GBM_mu_star * dt
      ##Shock:
      shock <- matrix(shocks[,i], nrow = nsteps) * GBM_sigma

      #Values:
      State.Matrix[2:(nsteps+1), seq(1, N_sim, ifelse(antithetic,2,1)),i] <- X.0.i + apply(drift + shock, MARGIN = 2, cumsum)
      #Antithetic Values:
      if(antithetic) State.Matrix[2:(nsteps+1), seq(2, N_sim, 2),i] <- X.0.i + apply(drift - shock, MARGIN = 2, cumsum)


    }
    #update log(Prices)
    prices <- prices + State.Matrix[,,i]
    ##If an uneven number was specified in antithetic:
    # if(verbose) if(N_sim != n) output[[paste("Factor",i)]] = X[,-ncol(X)] else output[[paste("Factor",i)]] = X
  }
  if(N_sim != n) prices <- prices[,-ncol(prices)]

  ##Return Simulated Values:

  if(!verbose) return(exp(parameters["E"] + prices))
   else return(list(State_Variables = State.Matrix, Prices = exp(parameters["E"] + prices)))

  }


#'Simulate N-Factor Model Futures Prices
#'@description Simulate Futures price data with dynamics that follow the parameters of an N-factor model through Monte Carlo simulation.
#'
#'@param X.0 Initial values of the state vector.
#'@param parameters A named vector of parameter values of a specified N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param dt discrete time step of simulation
#'@param N.obs The number of observations to simulate
#'@param TTM A vector or matrix of the time to maturity of futures contracts to simulate. See \bold{details}
#'@param verbose \code{logical}. Should the simulated state variables and associated prices be output?
#'
#'@details
#'\loadmathjax
#'The \code{Futures.Price.Simulate} function simulates futures price data using the Kalman Filter algorithm, drawing from a normal
#'distribution for the shocks in the transition and measurement equations at each discrete time step. At each discrete time point,
#'an observation of the state vector is generated through the transition equation, drawing from a normal distribution with a covariance equal to \mjeqn{Q_t}{Q[t]}.
#'Following this, simulated futures prices are generated through the measurement equation, drawing from a normal distribution with covariance matrix equal to \mjeqn{H}{H}.
#'
#'Input \code{TTM} can be either a matrix specifying the constant time to maturity of futures contracts to simulate, or it can be a matrix where \code{nrow(TTM) == N.obs} for the time-varying time to maturity of the futures contracts to simulate. This allows for the simulation of both aggregate stitched data and individual futures contracts.
#'
#'@return
#'\code{Futures.Price.Simulate} returns a list with three objects when \code{verbose = T} and a matrix of simulated futures prices when \code{verbose = F}. The list objects returned are:
#'
#'#'\tabular{ll}{
#'
#' \code{State.Vector} \tab  A \code{matrix} of Simulated state variables at each discrete time point. The columns represent each factor of the N-factor model and the rows represent
#'the simulated values at each discrete simulated time point. \cr
#'
#' \code{Futures} \tab  A \code{matrix} of Simulated futures prices, with each column representing a simulated futures contract. \cr
#'
#' \code{Spot} \tab A vector of simulated spot prices \cr
#' }
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'##Example 1 - Simulate Crude Oil Stitched futures prices
#'##under a Two-Factor model, assuming a constant time to maturity:
#'
#'Simulated.Stitched.Futures <- Futures.Price.Simulate(X.0 = c(log(SS.Oil$Spot[1,1]), 0),
#'                                                    parameters = SS.Oil$Two.Factor,
#'                                                    dt = SS.Oil$dt,
#'                                                    N.obs = nrow(SS.Oil$Stitched.Futures),
#'                                                    TTM = SS.Oil$Stitched.TTM)
#'
#'##Example 2 - Simulate Crude Oil Contract Prices under a Two-Factor model
#'
#'###Assume constant white noise in parameters of 1%:
#'SS.Oil.Two.Factor <- SS.Oil$Two.Factor
#'SS.Oil.Two.Factor <- SS.Oil.Two.Factor[!grepl("contract", names(SS.Oil.Two.Factor))]
#'SS.Oil.Two.Factor["sigma.contracts"] <- 0.01
#'
#'Simulated.Contracts <- Futures.Price.Simulate(X.0 = c(log(SS.Oil$Spot[1,1]), 0),
#'                                             parameters = SS.Oil.Two.Factor,
#'                                             dt = SS.Oil$dt,
#'                                             N.obs = nrow(SS.Oil$Contracts),
#'                                             TTM = SS.Oil$Contract.Maturities)
#'
#'@export
Futures.Price.Simulate <- function(X.0, parameters, dt, N.obs, TTM, verbose = TRUE){

  if(is.null(names(parameters))) stop("parameters must be a named vector. parametersgenerate function is suggested")
  transposeme = FALSE

  ##If it's a constant TTM, develop the maturity matrix:
  TTM <- as.matrix(TTM)
  ##If a maturity matrix has been provided:
  if(all(dim(TTM) > 1)){
    maturity.matrix <- TTM
    if(!any(dim(maturity.matrix)==N.obs)) stop("TTM does not match number of observations")
    if(dim(TTM)[1] != N.obs) {
      transposeme <- TRUE
      maturity.matrix <- t(maturity.matrix)
    }} else {
      maturity.matrix <- matrix(c(TTM), nrow = N.obs, ncol = length(c(TTM)), byrow = TRUE)
    }

    GBM <- "mu" %in% names(parameters)

    if(GBM){
      parameters["kappa_1"] <- 0
      parameters["E"] <- 0
    } else {
      parameters["mu"] <- 0
    }

    ###The highest sigma input would be our number of factors
    N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.numeric) & !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))
    N.contracts <- ncol(maturity.matrix)

    ###Contract White Noise:
    H <- diag(N.contracts)

    if("sigma.contracts" %in% names(parameters)){
      diag(H) <- parameters["sigma.contracts"]^2
    } else {
      diag(H) <- sapply(1:N.contracts, FUN = function(x) ifelse(parameters[paste0("sigma.contract_",x)]^2 < 1.01e-10 , 0, parameters[paste0("sigma.contract_",x)]^2))
      if(anyNA(H)) stop("Insufficient number of individual contract errors specified")
    }

    #Initialize State Variable
    if(length(X.0) != N.factors) stop("Initial state vector length not equal to N factors")
    X_t <- matrix(X.0)

    ##Discrete Time Steps:

    #Insantiate State Vectors ie. the outputs:
    X <- matrix(0  , nrow = N.obs, ncol = N.factors)
    Y <- matrix(0  , nrow = N.obs, ncol = N.contracts)

    d <- A_T(parameters, maturity.matrix)

    Z <- array(NA, dim = c(N.obs, N.contracts, N.factors))
    Z[,,1:N.factors] <- sapply(1:N.factors, FUN = function(X) exp(- parameters[paste0("kappa_", X)] * maturity.matrix))

    #Multivariate Normal observation - Measurement Equation:
    V <- MASS::mvrnorm(n = N.obs, mu = rep(0, N.contracts), Sigma = H)
    #Multivariate Normal observation - Transition Equation:
    omega <- MASS::mvrnorm(n = N.obs, mu = rep(0, N.factors), Sigma = cov_func(parameters, dt))


    for(t in 1:N.obs){

      #Measurement Equation:
      #y = d + Z * x_t + v_t
      #d
      d_t <- matrix(d[t,])

      #Z:
      Z_t <- Z[t,,]

      #Simulated Measurement Error:
      v_t <- V[t,]

      Y_t <- parameters["E"] + d_t + Z_t %*% X_t + v_t

      ##Record Observations:
      X[t,] <- X_t
      Y[t,] <- Y_t

      if(t != N.obs){

        #Omega - Transition Equation:
        omega_t <- omega[t,]
        #Transition Equation:
        #x_t = c + G * x_tm1 + omega
        c_t <- matrix(c(parameters["mu"] * dt, rep(0, N.factors-1)))
        #G:
        G_t <- diag(N.factors)
        diag(G_t) <- sapply(1:N.factors, FUN = function(X) exp(- parameters[paste0("kappa_", X)] * dt))

        X_t <- c_t + G_t %*% X_t + omega_t

      }
    }
    X_output <- X
    Y_output <- exp(Y)
    Spot <- as.matrix(exp(rowSums(parameters["E"] + X)))
    rownames(Y_output) <- rownames(Spot) <- rownames(X_output) <- rownames(maturity.matrix)

    colnames(Spot) <- "Spot"

    if(any(dim(TTM)==1)) colnames(Y_output) <- as.character(round(TTM,4))
    colnames(X_output) <- paste("Factor", 1:N.factors)

    if(verbose){
      if(transposeme)   return(list(State.Vector = t(X_output), Futures = t(Y_output), Spot = t(Spot)))
      return(list(State.Vector = X_output, Futures = Y_output, Spot = Spot))
    } else {
      if(transposeme)   return(t(Y_output))
      return(Y_output)
    }

  }
