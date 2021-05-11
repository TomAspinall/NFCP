#'Filter an N-factor commodity pricing model though the Kalman filter
#'
#'@description
#'\loadmathjax
#'Given a set of parameters of the N-factor model, filter term structure data using the Kalman filter.
#'
#'@param parameter_values \code{vector}. Numeric parameter values of an N-factor model.
#'
#'@param parameter_names \code{vector}. Parameter names, where each element of \code{parameter_names} must correspond to its respective value
#'element in object \code{parameter_values}.
#'
#'@param log_futures \code{matrix}. The natural logarithm of observed futures prices. Each row must correspond to quoted futures prices at a particular date and every column must correspond to a unique futures contract.
#'NA values are allowed.
#'
#'@param dt \code{numeric}. Constant, discrete time step of observations, in years.
#'
#'@param futures_TTM \code{vector} or \code{matrix}.  The time-to-maturity of observed futures contracts, in years, at a given observation date. This time-to-maturity can either be constant (ie. class 'vector') or variable (ie. class 'matrix') across observations.
#'The number of columns of 'futures_TTM' must be identical to the number of columns of object 'log_futures'. The number of rows of object 'futures_TTM' must be either 1 or equal to the number of rows of object 'log_futures'.
#'
#'@param ME_TTM \code{vector}. the time-to-maturity groupings to consider for observed futures prices. The length of \code{ME_TTM} must be equal to the number of 'ME' parameters specified in object 'parameter_names'. The maximum of 'ME_TTM' must be greater than the maximum value of 'futures_TTM'.
#'When the number of 'ME' parameter values is equal to one or the number of columns of object 'log_futures', this argument is ignored.
#'
#'@param verbose \code{logical}. Should additional information be output? see \bold{values}. When \code{verbose = F}, the \code{NFCP_Kalman_filter} function is significantly faster, see \bold{details}.
#'
#'@param debugging \code{logical}. Should additional filtering information be output? see \bold{values}.
#'
#'@details
#'
#'\code{NFCP_Kalman_filter} applies the Kalman filter algorithm for observable \code{log_futures} prices against the input parameters of an N-factor model
#'provided through the \code{parameter_values} and \code{parameter_names} input vectors.
#'
#'The \code{NFCP_Kalman_filter} function is
#'designed for subsequent input into optimization functions and is called within the N-factor parameter estimation function \code{NFCP_MLE}. The first input to the
#'\code{NFCP_Kalman_filter} function is a vector of parameters of an
#'N-factor model, with elements of this vector corresponding to the parameter names within the elements of input vector \code{parameter_names}.
#'When \code{logical} input \code{verbose = F}, the \code{NFCP_Kalman_filter} function calls the \code{fkf_SP} function of the \code{FKF_SP} package, which itself is a wrapper
#'of a routine of the Kalman filter written in C utilizing Sequential Processing for maximum computational efficiency (see \code{fkf_SP} for more details). When \code{verbose = T},
#'the \code{NFCP_Kalman_filter} instead applies a Kalman filter algorithm written in base \code{R} and outputs several other \code{list objects}, including filtered values and
#'measures for model fit and robustness (see \bold{Returns})
#'
#'\bold{The N-factor model}
#'The N-factor framework was first presented in the work of Cortazar and Naranjo (2006, equations 1-3).
#'It is a risk-premium class of commodity pricing model, in which futures prices are given by discounted expected future spot prices,
#'where these spot prices are discounted at a given level of risk-premium, known as the cost-of-carry.
#'
#'The N-factor framework describes the spot price process of a commodity as the correlated sum of \mjeqn{N}{N} state variables \mjeqn{x_t}{x[t]}. The 'NFCP' package also allows for a deterministic,
#'cyclical seasonal function \mjeqn{season(t)}{season(t)} to be considered.
#'
#'When \code{GBM = TRUE}:
#'\mjdeqn{log(S_{t}) = season(t) + \sum_{i=1}^N x_{i,t}}{log(S[t]) = season(t) + sum_{i=1}^n x[i,t]}
#'When \code{GBM = FALSE}:
#'\mjdeqn{log(S_{t}) = E + season(t) \sum_{i=1}^N x_{i,t}}{log(S[t]) = E + season(t) + sum_{i=1}^n x[i,t]}
#'
#'Where \code{GBM} determines whether the first factor follows a Brownian Motion or Ornstein-Uhlenbeck process to induce a unit root in the spot price process.
#'
#'When \code{GBM = TRUE}, the first factor corresponds to the spot price, and additional N-1 factors model the cost-of-carry.
#'
#'When \code{GBM = FALSE}, the commodity model assumes that there is a long-term equilibrium the commodity price will tend towards over time, with model volatility a decreasing function of time. This is not the standard approach made in the commodity pricing literature (Cortazar and Naranjo, 2006).
#'
#'State variables are thus assumed to follow the following processes:
#'
#'When \code{GBM = TRUE}:
#'\mjdeqn{dx_{1,t} = \mu^*dt + \sigma_{1} dw_{1}t}{dx[1,t] = mu^* dt + sigma[1] dw[1]t}
#'
#'When \code{GBM = FALSE}:
#'\mjdeqn{dx_{1,t} = - (\lambda_{1} + \kappa_{1}x_{1,t})dt + \sigma_{1} dw_{1}t}{dx[1,t] = - (lambda[1] + kappa[1] x[1,t]) dt + sigma[1] dw[t]t}
#'
#'And:
#'\mjdeqn{dx_{i,t} =_{i\neq 1} - (\lambda_{i} + \kappa_{i}x_{i,t})dt + \sigma_{i} dw_{i}t}{dx[i,t] =_(i != 1) - (lambda[i] + kappa[i] x[i,t]dt + sigma[i] dw[i]t)}
#'
#'where:
#'\mjdeqn{E(w_{i})E(w_{j}) = \rho_{i,j}}{E(w[i])E(w[j])}
#'
#'Additionally, the deterministic seasonal function (if specified) is given by:
#'
#'\mjdeqn{season(t) = \sum_{i=1} ( season_{i,1} cos(2i\pi) + season_{i,2} sin(2i\pi)}{season(t) = sum_{i=1} ( season_{i,1} cos(2.i.pi) + season_{i,2} sin(2.i.pi)}
#'
#'The addition of deterministic, cyclical seasonality as a function of trigonometric variables was first suggested by Hannan, Terrell, and Tuckwell (1970) and first applied to model commodities by Sørensen (2002).
#'
#'The following constant parameters are defined as:
#'
#'\code{var} \mjeqn{\mu}{mu}:  long-term growth rate of the Brownian Motion process.
#'
#'\code{var} \mjeqn{E}{E}: Constant equilibrium level.
#'
#'\code{var} \mjeqn{\mu^*=\mu-\lambda_1}{mu^* = mu-lambda[1]}: Long-term risk-neutral growth rate
#'
#'\code{var} \mjeqn{\lambda_{i}}{lambda[i]}: Risk premium of state variable \mjeqn{i}{i}.
#'
#'\code{var} \mjeqn{\kappa_{i}}{kappa[i]}: Reversion rate of state variable \mjeqn{i}{i}.
#'
#'\code{var} \mjeqn{\sigma_{i}}{sigma[i]}: Instantaneous volatility of state variable \mjeqn{i}{i}.
#'
#'\code{var} \mjeqn{\rho_{i,j} \in [-1,1]}{rho[i,j] in [-1,1]}: Instantaneous correlation between state variables \mjeqn{i}{i} and \mjeqn{j}{j}.
#'
#'Including additional factors within the spot-price process allow for additional flexibility (and possibly fit) to the term structure of a commodity.
#'The N-factor model nests simpler models within its framework, allowing for the fit of different N-factor models (applied to the same term structure data),
#'represented by the log-likelihood, to be directly compared with statistical testing possible through a chi-squared test.
#'
#'\bold{Disturbances - Measurement Error}:
#'
#'The Kalman filtering algorithm assumes a given measure of measurement error or disturbance in the measurement equation (ie. matrix \mjeqn{H}{H}). Measurement errors can be interpreted as error in the
#'model's fit to observed prices, or as errors in the reporting of prices (Schwartz and Smith, 2000). These disturbances are typically assumed independent.
#'
#'\code{var} \mjeqn{ME_i}{ME[i]} measurement error of contract \mjeqn{i}{i}.
#'
#'where the measurement error of futures contracts \mjeqn{ME_i}{ME[i]} is equal to \code{'ME_'} [i] (i.e. \code{'ME_1'}, \code{'ME_2'}, ...) specified in arguments \code{parameter_values} and \code{parameter_names}.
#'
#'There are three particular cases on how the measurement error of observations can be treated in the \code{NFCP_Kalman_filter} function:
#'
#'\bold{Case 1:} Only one ME is specified. The Kalman filter assumes that the measurement error of observations are independent and identical.
#'
#'\bold{Case 2:} One ME is specified for every observed futures contract. The Kalman filter assumes that the measurement error of observations are independent and unique.
#'
#'\bold{Case 3:} A series of ME's are specified for a given grouping of maturities of futures contracts. The Kalman filter assumes that the measurement error of observations are independent and unique to their respective time-to-maturity.
#'
#'Grouping of maturities for case 3 is specified through the \code{ME_TTM} argument. This is a vector that specifies the maximum maturity to consider for each respective ME parameter argument.
#'
#'in other words, ME_1 is considered for observations with TTM less than ME_TTM[1], ME_2 is considered for observations with TTM less than ME_TTM[2], ..., etc.
#'
#'The first case is clearly the simplest to estimate, but can be a restrictive assumption. The second case is clearly the most difficult to estimate, but can be an infeasible assumption when considering all available futures contracts that make up the term structure of a commodity.
#'
#'Case 3 thus serves to ease the restriction of case 1, and allow the user to make the modeling of measurement error as simple or complex as desired for a given set of maturities.
#'
#'\bold{Kalman Filtering}
#'
#'The following section describes the Kalman filter equations used to filter the N-factor model.
#'
#'The Kalman filter iteration is characterised by a transition and measurement equation.
#'The transition equation develops the vector of state variables between discretised time steps (whilst considering a given level of covariance between state variables over time).
#'The measurement equation relates the unobservable state vector to a vector of observable measurements (whilst also considering a
#'given level of measurement error). The typical Kalman filter algorithm is a Gaussian process state space model.
#'
#'Transition Equation:
#'\mjdeqn{\hat x_{t|t-1} = c_t + G_t \hat x_{t-1} + Q_t \eta_t }{hat(x)[t|t-1] = c[t] + G[t] * hat(x)[t-1]}
#'Measurement Equation:
#'\mjdeqn{\hat y_t = d_t + Z_t \hat x_{t|t-1} + H_t \epsilon_t}{hat(y)[t] = d[t] + Z[t] * hat(x)[t|t-1]}
#'
#'\mjdeqn{t = 1, \cdots, n }{t = 1, ..., n}
#'
#'Where \mjeqn{\eta_t}{eta[t]} and \mjeqn{\epsilon_t}{epsilon[t]} are IID \mjeqn{N(0,I(m))}{N(0,I(m))} and iid \mjeqn{N(0,I(d))}{N(0,I(d))} respectively.
#'
#'The state vector follows a normal distribution, \mjeqn{x_1 \sim N(a_1, P_1)}{x[1] ~ N(a[1], P[1])}, with \mjeqn{a_1}{a[1]} and \mjeqn{P_1}{P[1]} as the mean vector and variance matrix of
#'the initial state vector \mjeqn{x_1}{x[1]}, respectively.
#'
#'The Kalman filter can be used for parameter estimation through the maximization of the Log-Likelihood value. See \code{NFCP_MLE}.
#'
#'\bold{Filtering the N-factor model}
#'
#'let \mjeqn{m}{m} represent the number of observations at time \mjeqn{t}{t}
#'
#'let \mjeqn{n}{n} represent the number of factors in the N-factor model
#'
#'observable futures prices: \mjeqn{y_t = [ln(F(t,T_1)), ln(F(t,T_2)), \cdots, ln(F(t,T_m))]'}{y[t] = [ln(F(t,T[1])), ln(F(t,T[2])), ..., ln(F(t,T[m]))]'}
#'
#'State vector: \mjeqn{x_t=[x_1t,x_2t,\cdots,x_nt ]'}{x[t] = [x[1t], x[2t], ..., x[nt]]'}
#'
#'Measurement error: \mjeqn{diag(H) = [ME_{1}^2, ME_{2}^2, \cdots, ME_{n}^2]}{diag(H) = [ME[1]^2, ME[2]^2, ..., ME[n]^2]}
#'
#'When the number of specified ME terms is one, \mjeqn{s_1 = s_2 = \cdots = s_n = }{s[1] = s[2] = ... = s[n] = } \mjeqn{ME_1^2}{ME[1]^2}
#'
#'\code{var} \mjeqn{Z}{Z} is an \mjeqn{m \times n}{m X n} matrix, where each element \mjeqn{[i,j]}{[i,j]} is equal to:
#'
#'\mjdeqn{Z_{i,j} = e^{-\kappa_i T_j}}{Z[i,j] = e^(-kappa[i] * T[j])}
#'
#'\code{var} \mjeqn{d_t}{d[t]} is an \mjeqn{m \times 1}{m X 1} vector:
#'
#'\mjdeqn{d_t=[season(T_1) + A(T_1), season(T_2) + A(T_2), \cdots, season(T_m) + A(T_m)]'}{d[t]=[season(T[1]) + A(T[1]), season(T[2]) + A(T[2]), ..., season(T[m]) + A(T[m])]'}
#'
#'Under the assumption that Factor 1 follows a Brownian Motion, \eqn{A(T)} is given by:
#'\mjdeqn{A(T) = \mu^*T-\sum_{i=1}^N - \frac{1-e^{-\kappa_i T}\lambda_i}{\kappa_i}+\frac{1}{2}(\sigma_1^2T +
#'\sum_{i.j\neq 1} \sigma_i \sigma_j \rho_{i,j} \frac{1-e^{-(\kappa_i+\kappa_j)T}}{\kappa_i+\kappa_j})}{A(T) = mu^* * T - sum_{i=1}^N (1-e^(-kappa[i] T) lambda[i])/(kappa[i]) + 1/2 (sigma[1]^2 * T)
#' + sum_{i.j != 1} sigma[i] sigma[j] rho[i,j] (1 - e^(-(kappa[i] + kappa[j]) * T)) / (kappa[i] + kappa[j])}
#'
#'\code{var} \mjeqn{v_t}{v[t]} is a \mjeqn{n \times 1}{n X 1} vector of serially uncorrelated Guassian disturbances with
#'\mjeqn{E(V_t) = 0}{E(V[t]) = 0}  and \mjeqn{cov(v_t)=R^2}{cov(v[t])=R^2}
#'
#'Where:
#'
#'\mjeqn{diag(G_t) = [e^{-\kappa_1 \tau}, e^{-\kappa_2 \tau}, \cdots, e^{-\kappa_n \tau}]}{diag(G[t]) = [e^{-kappa[1] tau}, e^{-kappa[2] tau}, ..., e^{-kappa[n] tau}}
#'
#'
#'Where \mjeqn{ \tau =T-t}{tau = T - t}
#'
#'\code{var} \mjeqn{w_t}{w[t]} is an \mjeqn{n \times 1}{n X 1} vector of serially uncorrelated Guassian disturbances where:
#'\mjdeqn{E(w_t) = 0}{E(w[t]) = 0} and \mjeqn{cov(w_t) = Q_t}{cov(w[t]) = Q[t]}
#'
#'\code{var} \mjeqn{c_t=[\mu \Delta t,0,\cdots,0]'}{c[t] = [mu * Delta t, 0, ..., 0]'} is an \mjeqn{N \times 1}{N X 1} vector of the intercept of the transition equation.
#'
#'\code{var} \mjeqn{Q_t}{Q[t]} is equal to the covariance function, given by:
#'
#'\mjdeqn{Cov_{1,1}(x_{1,t},x_{1,t}) = \sigma_1^2t}{Cov[1,1](x[1,t],x[1,t]) = sigma[1]^2 * t}
#'\mjdeqn{Cov_{i,j}(x_{i,t},x_{j,t}) = \sigma_i\sigma_j\rho_{i,j}\frac{1-e^{-(\kappa_i+\kappa_j)t}}{\kappa_i+\kappa_j}}{Cov[i,j](x[i,t],x[j,t]) = sigma[i] sigma[j] rho[i,j] (1-e^{-(kappa[i]+kappa[j]) * t) / (kappa[i] + kappa[j])}}
#'(see also \code{cov_func})
#'
#'\bold{Penalising poorly specified models}
#'
#'The Kalman filter returns non-real log-likelihood scores when the prediction error variance matrix becomes singular or its determinant becomes negative. This generally occurs when a poorly specified parameter set is input, such as when measurement error is zero.
#'Non-real log-likelihood scores can break optimization and gradients algorithms and functions. To circumvent this, the \code{NFCP_Kalman_filter}
#'returns a heavily penalized log-likelihood score when \code{verbose = F}. Penalized log-likelihood scores are calculated by:
#'
#'\code{stats::runif(1, -2e6, -1e6)}
#'
#'\bold{Diffuse Kalman filtering}
#'
#'If the initial values of the state vector are not supplied within the \code{parameter_names} and \code{parameter_values} vectors, a 'diffuse' assumption is used within the Kalman filtering algorithm.
#'Initial states of factors that follow an Ornstein-Uhlenbeck are assumed to equal zero.
#'The initial state of the first factor, given that it follows a Brownian motion, is assumed equal to the first element of \code{log_futures}. This is an
#'assumption that the initial estimate of the spot price is equal to the closest to maturity observed futures price.
#'
#'The initial states of factors that follow an Ornstein-Uhlenbeck have a transient effect on future observations. This makes the diffuse assumption reasonable and further means that initial states cannot generally be accurately estimated.
#'
#'@return
#'\code{NFCP_Kalman_filter} returns a \code{numeric} object when \code{verbose = F}, which corresponds to the log-likelihood of observations.
#'When \code{verbose = T}, the \code{NFCP_Kalman_filter} function returns a \code{list} object of length seven with the following objects:
#'
#'\tabular{ll}{
#'
#'\code{Log-Likelihood} \tab Log-Likelihood of observations. \cr
#'
#'\code{Information Criteria} \tab \code{vector}. The Akaikie and Bayesian Information Criterion. \cr
#'
#'\code{X_t} \tab \code{vector}. The final observation of the state vector. \cr
#'
#'\code{X} \tab \code{matrix}. Optimal one-step-ahead state vector. \cr
#'
#'\code{Y} \tab \code{matrix}. Estimated futures prices. \cr
#'
#'\code{V} \tab \code{matrix}. Estimation error. \cr
#'
#'\code{Filtered Error} \tab \code{matrix}. positive mean error (high bias), negative mean error (low bias),
#'
#'mean error (bias) and root mean squared error (RMSE)
#'
#'of the filtered values to observed futures prices.  \cr
#'
#'\code{Term Structure Fit} \tab \code{matrix}. The mean error (Bias), mean absolute error, standard deviation of error
#'
#'  and root mean squared error (RMSE) of each observed futures contract. \cr
#'
#'\code{Term Structure Volatility Fit} \tab \code{matrix}. Theoretical and empirical volatility of observed futures contract returns. \cr
#' }
#'
#'When \code{debugging = T}, 9 objects are returned in addition to those returned when \code{verbose = T}:
#'
#'\tabular{ll}{
#'
#'\code{P_t} \tab \code{array}. Covariance matrix of state variables, with the third dimension indexing across time \cr
#'
#'\code{F_t} \tab \code{vector}. Prediction error variance matrix, with the third dimension indexing across time \cr
#'
#'\code{K_t} \tab \code{matrix}. Kalman Gain, with the third dimension indexing across time \cr
#'
#'\code{d} \tab \code{matrix}.  \mjeqn{d_t}{d[t]} (see \bold{details}) \cr
#'
#'\code{Z} \tab \code{matrix}.  \mjeqn{Z_t}{z[t]} (see \bold{details}) \cr
#'
#'\code{G_t} \tab \code{matrix}.  \mjeqn{G_t}{G[t]} (see \bold{details})  \cr
#'
#'\code{c_t} \tab \code{vector}.  \mjeqn{C_t}{c[t]} (see \bold{details}) \cr
#'
#'\code{Q_t} \tab \code{matrix}. \mjeqn{Q_t}{Q[t]}  (see \bold{details}) \cr
#'
#'\code{H} \tab \code{matrix}. \mjeqn{H}{H}  (see \bold{details}) \cr
#' }
#'
#'
#'
#'@references
#'
#'Hannan, E. J., et al. (1970). "The seasonal adjustment of economic time series." \emph{International economic review}, 11(1): 24-52.
#'
#'Anderson, B. D. O. and J. B. Moore, (1979). \emph{Optimal filtering} Englewood Cliffs: Prentice-Hall.
#'
#'Fahrmeir, L. and G. tutz,(1994) \emph{Multivariate Statistical Modelling Based on Generalized Linear Models.} Berlin: Springer.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Sørensen, C. (2002). "Modeling seasonality in agricultural commodity futures." \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products} 22(5): 393-426.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Durbin, J., and S. J. Koopman, (2012). \emph{Time series analysis by state space methods.} Oxford university press.
#'
#'
#'@examples
#'
#'
#'##Example 1 - complete, stitched data.
#'##Replicating the Schwartz and Smith (2000)
#'##Two-Factor commodity pricing model applied to crude oil:
#'
#'SS_stitched_filtered <- NFCP_Kalman_filter(
#'parameter_values = SS_oil$two_factor,
#'parameter_names = names(SS_oil$two_factor),
#'log_futures = log(SS_oil$stitched_futures),
#'futures_TTM = SS_oil$stitched_TTM,
#'## maturity groupings need not be considered here:
#'ME_TTM = NULL,
#'dt = SS_oil$dt,
#'verbose = FALSE)
#'
#'
#'##Example 2 - incomplete, contract data.
#'##Replicating the Schwartz and Smith (2000)
#'##Two-Factor commodity pricing model applied to all available
#'##crude oil contracts:
#'
#'SS_2F <- SS_oil$two_factor
#'##omit stitched contract white noise
#'SS_2F <- SS_2F[!grepl("ME",
#'              names(SS_2F))]
#'
#'# Evaluate two different measurement errors
#'SS_2F[c("ME_1", "ME_2")] <- c(0.01, 0.04)
#'
#'## Separate measurement error into two different maturity groupings
#'SS_ME_TTM <- c(1,3)
#'## ME_1 is applied for observed contracts with less than one year
#'## maturity, whilst ME_2 considers contracts with maturity greater
#'## than one year, and less than three years
#'
#'#Kalman filter
#'SS_contract_filtered <- NFCP_Kalman_filter(
#'parameter_values = SS_2F,
#'parameter_names = names(SS_2F),
#'## All available contracts are considered
#'log_futures = log(SS_oil$contracts),
#'## Respective 'futures_TTM' of these contracts are input:
#'futures_TTM = SS_oil$contract_maturities,
#'ME_TTM = SS_ME_TTM,
#'dt = SS_oil$dt,
#'verbose = FALSE)
#'@export
NFCP_Kalman_filter = function(parameter_values, parameter_names, log_futures, dt, futures_TTM, ME_TTM = NULL, verbose = FALSE, debugging = FALSE){

  ## Named Parameter Vector:
  params <- parameter_values
  names(params) <- parameter_names

  ## How many measurement errors (ME) are specified?
  N_ME <- max(which(paste0("ME_", 1:length(params)) %in% names(params) & sapply(params[paste0("ME_",1:length(params))], FUN = is.numeric) & !sapply(params[paste0("ME_",1:length(params))], FUN = is.na)))

  ## How many factors are specified?
  N_factors <- max(which(paste0("sigma_", 1:length(params)) %in% names(params) & sapply(params[paste0("sigma_",1:length(params))], FUN = is.numeric) & !sapply(params[paste0("sigma_",1:length(params))], FUN = is.na)))

  ## How many seasonality factors are specified?
  N_season <- length(grep("season", parameter_names))/2

  ##Standardize format:
  log_futures <- as.matrix(log_futures)
  futures_TTM <- as.matrix(futures_TTM)

  ##Dimensions of observations:
  N_obs <- nrow(log_futures)
  N_contracts <- ncol(log_futures)

  ## "Contract" data or "Aggregate" Data?
  contract_data <- all(dim(futures_TTM)>1)
  if(contract_data && !all(dim(log_futures) == dim(futures_TTM))) stop("log_futures and futures_TTM have different dimensions")
  if(!contract_data && length(futures_TTM)!=N_contracts) stop("Aggregate futures data, however ncol(log_futures) and length(futures_TTM) have different dimensions")

  ## Are we using ME maturity groupings?
  inhomogeneous_H <- N_ME > 1 && N_ME < ncol(log_futures)

  ## If so, have enough been specified?
  if(inhomogeneous_H){
    if(is.null(ME_TTM)) stop("Multiple measurement error (ME) terms have been specified but the maturity terms of the measurement error (ME_TTM) have not.")
    if(length(ME_TTM) != N_ME) stop("Number of measurement error (ME) terms specified does not match the length of argument 'ME_TTM'")
    if(max(futures_TTM, na.rm = TRUE) > max(ME_TTM, na.rm = TRUE)) stop("Maximum observed contract maturity (futures_TTM) is greater than the max specified maturity grouping for the measurement error (ME_TTM)")
  }

  ##Is factor one a GBM or MR Process?
  GBM <- "mu" %in% names(params)
  if(GBM) zerome <- c("kappa_1", "E") else zerome <- "mu"
  params[zerome] <- 0

  #x_t is our vector of state variables:
  x_t <- matrix(rep(0, N_factors))
  #This is our best estimate if we're not specifying anything:
  if(GBM) x_t[1] <- log_futures[1,1]
  #But if we are:
  if("x_0_1" %in% names(params)) x_t <- matrix(sapply(1:N_factors, FUN = function(x) if(paste0("x_0_",x) %in% names(params)) params[paste0("x_0_",x)]))


  #c:
  c_t <- c(params["mu"] * dt, rep(0, N_factors-1))

  #Variance of omega:
  # P_t <- Q_t <- cov_func(params, dt)
  Q_t <- cov_func(params, dt)
  ## Diffuse assumption:
  P_t <- diag(100, N_factors)

  #G:
  G_t <- diag(N_factors)
  diag(G_t) <- sapply(1:N_factors, FUN = function(X) exp(- params[paste0("kappa_", X)] * dt))

  ## Adjust data into required format
  futures_TTM <- t(futures_TTM)

  Zt <- array(NA, dim = c(N_contracts, N_factors, ifelse(contract_data, N_obs,1)))
  for(i in 1:N_factors) Zt[,i,] <- exp(- params[paste0("kappa_", i)] * futures_TTM)

  ## Incorporate Seasonality (if there is any):
  seasonality <- 0
  if(N_season > 0) for(i in 1:N_season) seasonality <- seasonality + params[paste0("season_", i, "_1")] * cos(2 * i * pi * futures_TTM) + params[paste0("season_", i, "_2")] * sin(2 * i * pi * futures_TTM)

  ## Place data into required format:
  dtT <- matrix(seasonality + params["E"] + A_T(params, futures_TTM), nrow = N_contracts, ncol = ifelse(contract_data,N_obs,1))

  # Measurement error - diagonal elements:
  Ht <- matrix(rep(0, N_contracts), nrow = N_contracts, ncol = ifelse(inhomogeneous_H, N_obs, 1))

  ## Case 1 - Only one ME has been specified:
  if(N_ME == 1) Ht[1:N_contracts] <- params["ME_1"]^2

  ## Case 2 - an ME for every contract has been specified:
  if(N_ME == N_contracts) Ht[1:N_contracts] <- sapply(1:N_ME, FUN = function(x) params[paste0("ME_",x)]^2)
  ## Case 3 - Multiple ME have been specified, corresponding to particular maturity groupings:
  if(inhomogeneous_H) for(loop in N_ME:1) Ht[futures_TTM < ME_TTM[loop]] = params[paste0("ME_",loop)]^2
  ## Zero out minimum values:
  Ht[Ht<1.01e-10] <- 0

  if(debugging) verbose <- TRUE
  # -------------------------------------
  ##Kalman filter in C (faster):
  if(!verbose){

      log_likelihood <- suppressWarnings(FKF.SP::fkf.SP(a0 = c(x_t),
                   P0 = P_t,
                   dt = c_t,
                   ct = dtT,
                   Tt = G_t,
                   Zt = Zt,
                   HHt = Q_t,
                   GGt = Ht,
                   yt = t(log_futures)))

    ## If The model was poorly specified, the log-likelihood returns NA. We need to return a heavily penalized score for the gradient function.
    return(ifelse(is.na(log_likelihood),stats::runif(1, -2e6, -1e6), log_likelihood))
  } else {
  # -------------------------------------
    # Kalman filter in R

    #Verbose variables:
    save_X <- matrix(0, nrow = N_obs, ncol = N_factors)
    save_X_SD <- matrix(0, nrow = N_obs, ncol = N_factors)
    save_V <- matrix(NA, nrow = N_obs, ncol = ncol(log_futures))
    save_Y <- matrix(NA, nrow = N_obs, ncol = ncol(log_futures))
    rownames(save_X) <- rownames(save_X_SD) <- rownames(save_V) <- rownames(save_Y) <- rownames(log_futures)

    #Debugging variables:
    if(debugging){
      max_obs <- max(rowSums(!is.na(log_futures)))
      save_P <- array(NA, dim = c(N_factors, N_factors, N_obs))
      save_F <- array(NA, dim = c(max_obs, max_obs, N_obs))
      save_K <- array(NA, dim = c(N_factors, max_obs, N_obs))
      LL_t <- rep(0, N_obs)
    }

    #Initialize
    log_likelihood <- 0
    I <- diag(N_factors)

    #####BEGIN Kalman filter:
    for(t in 1:N_obs){

      ##Transition Equation:
      x_tp1 <- c_t + G_t %*% x_t

      ##Covariance Transition Equation:
      P_tp1 <- G_t %*% P_t %*% t(G_t) + Q_t

      ##How many contracts are we observing this iteration?
      obs_contracts_t <- which(!is.na(log_futures[t,]))
      if(length(obs_contracts_t)>0){

        N_obs_contracts_t <- length(obs_contracts_t)

        ##Update matrices:
        d_t <- matrix(dtT[obs_contracts_t,ifelse(contract_data, t, 1)])
        Z_t <- matrix(Zt[obs_contracts_t,,ifelse(contract_data, t, 1)], ncol = N_factors)

        #Measurement Errors matrix:
        #Expectation of epsilon_t is 0
        #Variance    of epsilon_t is H
        H_t <- diag(N_obs_contracts_t)
        diag(H_t) <- Ht[obs_contracts_t,ifelse(inhomogeneous_H, t, 1)]

        #Function of covariance matrix:
        F_t  <- Z_t %*% P_tp1 %*% t(Z_t) + H_t

        det_F_t <- suppressWarnings(log(det(F_t)))

        ##Numeric Stability - Poorly Conditioned params:
        inverse_F_t <- try(solve(F_t))
        if(is.na(det_F_t))                         stop("Negative determinant in Kalman filter Covariance matrix. Parameters may be poorly specified.")
        if(any(class(inverse_F_t) == "try-error")) stop("Singular Kalman filter Covariance Matrix. Parameters may be poorly specified.")

        #Kalman Gain:
        K_t <- P_tp1 %*% t(Z_t) %*% inverse_F_t
        P_tp1 <- (I - K_t %*% Z_t) %*% P_tp1

        P_t <- P_tp1

        ##Measurement Equation:
        y_bar_t <- d_t + Z_t %*% x_tp1
        #Actual Futures Prices:
        y_t <- as.matrix(log_futures[t,obs_contracts_t])
        #Prediction Error:
        v_t <- y_t - y_bar_t
        #Correction based upon prediction error:
        x_tp1 <- x_tp1 + K_t %*% v_t
        ###Update, iteration begins anew:
        x_t <- x_tp1
        P_t <- P_tp1

        #Update Concurrent Log Likelihood of Observations:
        log_likelihood <- sum(log_likelihood, - (1/2) * sum(N_obs_contracts_t * log(2*pi), det_F_t, t(v_t) %*% inverse_F_t %*% v_t))
        #-----------------------

        #Record estimated variables
        #Updated Error terms:
        y_tt <- d_t + Z_t %*% x_t
        v_tt <- y_tt - y_t

        save_X[t,] <- x_t
        save_X_SD[t,] <- diag(P_t)
        save_Y[t,obs_contracts_t] <- y_tt
        save_V[t,obs_contracts_t] <- v_tt
        log_likelihood_result <- log_likelihood
        if(debugging){

          save_P[,,t] <- P_t
          save_F[1:N_obs_contracts_t,1:N_obs_contracts_t,t] <- F_t
          save_K[,1:N_obs_contracts_t,t] <- K_t

          LL_t[t] <- log_likelihood
        }
      }
    }

    #Final observations, for forecasting purposes:
    X.t <- c(x_t)
    names(X.t) <- paste0("x_", 1:N_factors, "_t")

    ####Term Structure Analysis:
    #Save the filtered Observations:
    Y_output <- exp(cbind(params["E"] + rowSums(save_X),save_Y))
    if(!is.null(colnames(log_futures))) colnames(Y_output) <- c("filtered Spot", colnames(log_futures))
    colnames(save_X) <- colnames(save_X_SD) <- paste("Factor", 1:N_factors)
    rownames(Y_output) <- rownames(save_X)

    ###Term Structure Fit:
    Term_Structure_Fit <- matrix(0, nrow = 4, ncol = ncol(log_futures))

    ##Mean Error:
    Term_Structure_Fit[1,] <- colMeans(save_V, na.rm = TRUE)
    ##Mean Absolute Error
    Term_Structure_Fit[2,] <- colMeans(abs(save_V), na.rm = TRUE)
    ##SD of Error:
    Term_Structure_Fit[3,] <- apply(save_V, MARGIN = 2, FUN = function(x) stats::sd(x, na.rm = TRUE))
    ##RMSE of each contract:
    Term_Structure_Fit[4,] <- sqrt(colMeans(save_V^2, na.rm = TRUE))

    rownames(Term_Structure_Fit) <- c("Mean Error", "Mean Absolute Error", "SD Error", "RMSE")
    colnames(Term_Structure_Fit) <- colnames(save_V) <- colnames(log_futures)

    ### Term structure fit to all observations:
    Filtered_Error <- c(`High Bias` = mean(save_V[save_V > 0], na.rm = TRUE), `Low Bias` =  mean(save_V[save_V < 0], na.rm = TRUE), `Bias` = mean(save_V, na.rm = TRUE),
                            `RMSE` = sqrt(mean(save_V^2, na.rm = TRUE)))

    ## Information Criteria:
    ### Number of univariate observations
    np_IC <- sum(!is.na(c(log_futures)))
    ### Akaike Information Criterion:
    AIC <- ((2 * length(parameter_values)) - 2 * log_likelihood)/np_IC
    ### Bayesian Information Criterion:
    BIC <- (length(parameter_values) * log(np_IC) - 2 * log_likelihood)/np_IC
    IC <- c("AIC" = AIC, "BIC" = BIC)

    ###Volatility TSFit:
    if(contract_data) {
      Volatility_TSFit <- TSfit_volatility(params, exp(log_futures), futures_TTM[,ncol(futures_TTM)], dt)
    } else {
      Volatility_TSFit <- TSfit_volatility(params, exp(log_futures), futures_TTM, dt) }

    ##Verbose List
    output = list("Log-Likelihood" = log_likelihood, "Information Criteria" = IC , x_t = X.t, X = save_X, Y = Y_output,
                  V = save_V, "Filtered Error" = Filtered_Error, "Term Structure Fit" = Term_Structure_Fit, "Term Structure Volatility Fit" = Volatility_TSFit)

    ##Debugging List:
    if(debugging) output <- c(output, list(LL_t = LL_t, P_t = save_P, F_t = save_F, K_t = save_K, d_t = dtT, Z_t = Zt, G_t = G_t, c_t = c_t, Q_t = Q_t, H_t = Ht))

    #Return Output value:
    return(output)
  }
  }

