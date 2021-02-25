#'N-Factor Commodity Pricing Kalman Filter
#'
#'@description
#'\loadmathjax
#'Given a set of parameters of the N-factor model, filter term structure data using the Kalman filter.
#'
#'@param parameter.values Vector of parameter values of an N-factor model. The \code{NFCP.Kalman.filter} function is designed for
#'application to \code{optim} type functions, and thus parameter values and
#'corresponding parameters different inputs within the function.
#'
#'@param parameters Vector of parameter names. Each element of \code{parameters} must correspond to its respective value
#'element in object \code{parameter.values}.
#'
#'@param log.futures Object of class \code{matrix} corresponding to the natural logarithm of observable futures prices.
#'NA's are allowed within the \code{matrix}. Every column of the matrix must correspond to a particular futures contract,
#'with each row corresponding to a quoted price on a given date.
#'
#'@param TTM Object of class 'vector' or 'matrix' that specifies the time to maturity of observed futures contracts.
#'time to maturity can be either constant (ie. class 'vector') or time homogeneous (ie. class 'matrix').
#'When the time to maturity of observed futures contracts is time homogeneous, the dimensions of
#'\code{TTM} must be identical to that of \code{log.futures}. Every element of \code{TTM}
#'corresponds to the time to maturity, in years, of a futures contract at a given observation date.
#'
#'@param dt Constant, discrete time step of observations
#'
#'@param verbose \code{logical}. Should additional information be output? see \bold{values}. When \code{verbose = F}, the \code{NFCP.Kalman.filter} function is significantly faster, see \bold{details}
#'
#'@param debugging \code{logical}. Should additional filtering information be output? see \bold{values}
#'
#'@details
#'
#'\code{NFCP.Kalman.filter} applies the Kalman filter algorithm for observable \code{log.futures} prices against the input parameters of an N-factor model
#'provided through the \code{parameter.values} and \code{parameters} input vectors.
#'
#'The \code{NFCP.Kalman.filter} function is
#'designed for subsequent input into optimization functions and is called within the N-factor parameter estimation function \code{NFCP.MLE}. The first input to the
#'\code{NFCP.Kalman.filter} function is a vector of parameters of an
#'N-factor model, with elements of this vector corresponding to the parameter names within the elements of input vector \code{parameters}.
#'When \code{logical} input \code{verbose = F}, the \code{NFCP.Kalman.filter} function calls the \code{fkf_SP} function of the \code{FKF_SP} package, which itself is a wrapper
#'of a routine of the Kalman filter written in C utilizing Sequential Processing for maximum computational efficiency (see \code{fkf_SP} for more details). When \code{verbose = T},
#'the \code{NFCP.Kalman.filter} instead applies a Kalman filter algorithm written in base \code{R} and outputs several other \code{list objects}, including filtered values and
#'measures for model fit and robustness (see \bold{Returns})
#'
#'\bold{The N-factor model}
#'The N-factor model was first presented in the work of Cortazar and Naranjo (2006, equations 1-3). The N-factor framework
#'describes the spot price process of a commodity as the correlated sum of \mjeqn{N}{N} state variables \mjeqn{x_t}{x[t]}.
#'
#'When \code{GBM = TRUE}:
#'\mjdeqn{log(S_{t}) = \sum_{i=1}^N x_{i,t}}{log(S[t]) = sum_{i=1}^n x[i,t]}
#'When \code{GBM = FALSE}:
#'\mjdeqn{log(S_{t}) = E + \sum_{i=1}^N x_{i,t}}{log(S[t]) = E + sum_{i=1}^n x[i,t]}
#'
#'Additional factors within the spot-price process are designed to result in additional flexibility, and possibly fit to the observable term structure, in
#' the spot price process of a commodity. The fit of different N-factor models, represented by the log-likelihood can be directly compared with statistical
#' testing possible through a chi-squared test.
#'
#'Flexibility in the spot price under the N-factor framework allows the first factor to follow a Brownian Motion or Ornstein-Uhlenbeck process to induce a unit root.
#'In general, an N-factor model where \code{GBM = T}
#'allows for non-reversible behaviour within the price of a commodity, whilst \code{GBM = F} assumes that there is a long-run equilibrium that
#'the commodity price will revert to in the long-term.
#'
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
#'The following constant parameters are defined as:
#'
#'\code{param} \mjeqn{\mu}{mu}:  long-term growth rate of the Brownian Motion process.
#'
#'\code{param} \mjeqn{E}{E}: Constant equilibrium level.
#'
#'\code{param} \mjeqn{\mu^*=\mu-\lambda_1}{mu^* = mu-lambda[1]}: Long-term risk-neutral growth rate
#'
#'\code{param} \mjeqn{\lambda_{i}}{lambda[i]}: Risk premium of state variable \mjeqn{i}{i}.
#'
#'\code{param} \mjeqn{\kappa_{i}}{kappa[i]}: Reversion rate of state variable \mjeqn{i}{i}.
#'
#'\code{param} \mjeqn{\sigma_{i}}{sigma[i]}: Instantaneous volatility of state variable \mjeqn{i}{i}.
#'
#'\code{param} \mjeqn{\rho_{i,j} \in [-1,1]}{rho[i,j] in [-1,1]}: Instantaneous correlation between state variables \mjeqn{i}{i} and \mjeqn{j}{j}.
#'
#'\bold{Measurement Error}:
#'
#'The Kalman filtering algorithm assumes a given measure of measurement error (ie. matrix \mjeqn{H}{H}). Measurement errors can be interpreted as error in the
#'model's fit to observed prices, or as errors in the reporting of prices (Schwartz and Smith, 2000) and are assumed independent.
#'
#'\code{var} \mjeqn{s_i}{s[i]} Observation error of contract \mjeqn{i}{i}.
#'
#'When \code{S.Constant = T}, the values of parameter \mjeqn{s_i}{s[i]} are assumed constant and equal to parameter object 'sigma.contracts'. When \code{S.Constant = F}, observation errors are assumed unique, where the error of
#'futures contracts \mjeqn{s_i}{s[i]} is equal to parameter object \code{'sigma.contract_'} [i] (i.e. \code{'sigma.contract_1'}, \code{'sigma.contract_2'}, ...).
#'
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
#'When \code{verbose = F}, the \code{NFCP.Kalman.filter} function performs Kalman filtering through the \code{fkf.SP} function of the \code{FKF.SP} package for maximum filtering efficiency, which is
#'crucial when filtering and estimating parameters of term structure data under the N-factor model, which generally has many observations, contracts and unknown parameters. When \code{verbose = T},
#'the \code{NFCP.Kalman.filter} function instead performs Kalman filtering in R, returning greater details of filtered objects (see \bold{Value})
#'
#'The Kalman filter can be used for parameter estimation through the maximization of the Log-Likelihood value. See \code{NFCP.MLE}.
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
#'Measurement error: \mjeqn{diag(H) = [s_{1}^2, s_{2}^2, \cdots, s_{n}^2]}{diag(H) = [s[1]^2, s[2]^2, ..., s[n]^2]}
#'
#'Where \mjeqn{s_i}{s[i]} == \code{"sigma.contract_"} [i] when the \code{logical} of function \code{NFCP.Parameters} \code{S.constant = F}
#'
#'When \code{S.constant = T}, \mjeqn{s_1 = s_2 = \cdots = s_n = }{s[1] = s[2] = ... = s[n] = } \code{"sigma.contracts"}
#'
#'\code{var} \mjeqn{Z}{Z} is an \mjeqn{m \times n}{m X n} matrix, where each element \mjeqn{[i,j]}{[i,j]} is equal to:
#'
#'\mjdeqn{Z_{i,j} = e^{-\kappa_i T_j}}{Z[i,j] = e^(-kappa[i] * T[j])}
#'
#'\code{var} \mjeqn{d_t}{d[t]} is an \mjeqn{m \times 1}{m X 1} vector:
#'
#'\mjdeqn{d_t=[A(T_1), A(T_2), \cdots, A(T_m)]'}{d[t]=[A(T[1]), A(T[2]), ..., A(T[m])]'}
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
#'The Kalman filter returns non-real log-likelihood scores when the function of the covariance matrix becomes singular or its determinant becomes negative.
#'This occurs when a poorly specified parameter set is input. Non-real log-likelihood scores can break optimization algorithms. To circumvent this, the \code{NFCP.Kalman.filter}
#'returns a heavily penalized log-likelihood score whilst also returning a warning. Penalized log-likelihood scores are calculated by:
#'
#'\code{stats::runif(1, -1.5e6, -1e6)}
#'
#'\bold{Diffuse Kalman filtering}
#'
#'If the initial values of the state vector are not supplied within the \code{parameters} and \code{parameter.values} vectors (ie. \code{Initial.State = F} within the
#'\code{NFCP.Parameters} function), a 'diffuse' assumption is used within the Kalman filtering algorithm. Factors that follow an Ornstein-Uhlenbeck are assumed to equal zero. When
#'\code{Estimate.Initial.State = F} and \code{GBM = T}, the initial value of the first state variable is assumed to equal the first element of \code{log.futures}. This is an
#'assumption that the initial estimate of the spot price is equal to the closest to maturity observed futures price.
#'
#'The initial covariance of the state vector for the Kalman filtering algorithm assumed to be equal to matrix \mjeqn{Q}{Q}
#'
#'Initial states of factors that follow an Ornstein-Uhlenbeck have a transient effect on future
#'observations, however the initial value of a random walk variable persists across observations and therefore influencing model fit more (see Schwartz and Smith (2000) for more details).
#'
#'
#'@return
#'\code{NFCP.Kalman.filter} returns a \code{numeric} object when \code{verbose = F}, which corresponds to the log-likelihood of observations.
#'When \code{verbose = T}, the \code{NFCP.Kalman.filter} function returns a \code{list} object of length seven with the following objects:
#'
#'\tabular{ll}{
#'
#' \code{LL} \tab Log-Likelihood of observations \cr
#'
#'\code{X.t} \tab \code{vector}. The final observation of the state vector \cr
#'
#'\code{X} \tab \code{matrix}. All observations of the state vector, after the updating equation has been applied \cr
#'
#'\code{Y} \tab \code{matrix}. Estimated futures prices at each observation \cr
#'
#'\code{V} \tab \code{matrix}. Estimation error of each futures contracts at each observation \cr
#'
#'\code{Filtered.Error} \tab \code{matrix}. positive mean error (high bias), negative mean error (low bias), mean error (bias) and root mean squared error (RMSE) of the filtered values to observed futures prices.  \cr
#'
#'\code{TSFit.Error} \tab \code{matrix}. The Mean Error (Bias), Mean Absolute Error, Standard Deviation of Error and Root Mean Squared Error (RMSE) of each
#'observed contract, matching the column names of \code{log.futures}  \cr
#'
#'\code{TSFit.Volatility} \tab \code{matrix}. The theoretical and empirical volatility of futures returns for each observed contract as returned from the \code{TSFit.Volatility} function \cr
#' }
#'
#' When \code{debugging = T}, 9 objects are returned in addition to those returned when \code{verbose = T}:
#'
#'\tabular{ll}{
#'
#' \code{P_t} \tab \code{array}. The covariance matrix at each observation point, with the third dimension indexing across time \cr
#'
#'\code{F_t} \tab \code{vector}. The function of the Kalman filter covariance matrix at each observation point, with the third dimension indexing across time \cr
#'
#'\code{K_t} \tab \code{matrix}. The Kalman Gain at each observation point, with the third dimension indexing across time \cr
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
#'Anderson, B. D. O. and J. B. Moore, (1979). \emph{Optimal filtering} Englewood Cliffs: Prentice-Hall.
#'
#'Fahrmeir, L. and G. tutz,(1994) \emph{Multivariate Statistical Modelling Based on Generalized Linear Models.} Berlin: Springer.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Durbin, J., and S. J. Koopman, (2012). \emph{Time series analysis by state space methods.} Oxford university press.
#'
#'
#'@examples
#'
#'##Example 1 - complete, stitched data.
#'##Replicating the Schwartz and Smith (2000)
#'##Two-Factor commodity pricing model applied to crude oil:
#'
#'Schwartz.Smith.Oil.Stitched <- NFCP.Kalman.filter(
#'  parameter.values = SS.Oil$Two.Factor,
#'  parameters = names(SS.Oil$Two.Factor),
#'  log.futures = log(SS.Oil$Stitched.Futures),
#'  TTM = SS.Oil$Stitched.TTM,
#'  dt = SS.Oil$dt,
#'  verbose = FALSE)
#'
#'
#'##Example 2 - incomplete, contract data.
#'##Replicating the Schwartz and Smith (2000)
#'##Two-Factor commodity pricing model applied to all available
#'##crude oil contracts:
#'
#'SS.Oil.2F <- SS.Oil$Two.Factor
#'##omit stitched contract white noise
#'SS.Oil.2F <- SS.Oil.2F[!grepl("sigma.contract",
#'                                      names(SS.Oil.2F))]
#'##Assume constant white noise in observable contracts:
#'SS.Oil.2F["sigma.contracts"] <- 0.01
#'
#'#Kalman filter
#'Schwartz.Smith.Oil.Contracts <- NFCP.Kalman.filter(
#'  parameter.values = SS.Oil.2F,
#'  parameters = names(SS.Oil.2F),
#'  ## All available contracts are considered
#'  log.futures = log(SS.Oil$Contracts),
#'  ## Respective 'TTM' of these contracts are input
#'  TTM = SS.Oil$Contract.Maturities,
#'  dt = SS.Oil$dt,
#'  verbose = FALSE)
#'@export
NFCP.Kalman.filter = function(parameter.values, parameters, log.futures, dt, TTM, verbose = FALSE, debugging = FALSE){

  if(debugging) verbose <- TRUE

  params <- parameter.values
  names(params) <- parameters

  ##Standardize format:
  log.futures <- as.matrix(log.futures)
  TTM <- as.matrix(TTM)

  ## "Contract" data or "Aggregate" Data?
  Contract.Data <- all(dim(TTM)>1)
  if(Contract.Data && !all(dim(log.futures) == dim(TTM))) stop("Observations and Maturity matrix have different dimensions")
  if(!Contract.Data && length(TTM)!=ncol(log.futures)) stop("Number of observations does not equal number of time homogeneous TTM's")

  ##GBM or MR Process?
  GBM <- "mu" %in% names(params)

  ###First factor GBM or MR?
  if(GBM){
    params["kappa_1"] <- 0
    params["E"] <- 0
  } else {
    params["mu"] <- 0
  }
  if(!"sigma.contracts" %in% names(params)) if(!paste0("sigma.contract_", ncol(log.futures)) %in% names(params)) stop("Not enough Individual Contract white noises have been specified")

  ###The highest sigma input would be our number of factors
  N.factors <- max(which(paste0("sigma_", 1:length(params)) %in% names(params) & sapply(params[paste0("sigma_",1:length(params))], FUN = is.numeric) & !sapply(params[paste0("sigma_",1:length(params))], FUN = is.na)))

  ##Dimensions of observations:
  N.obs <- nrow(log.futures)
  N.contracts <- ncol(log.futures)

  #X_t is our vector of state variables:
  X_t <- matrix(rep(0, N.factors))
  #This is our best estimate if we're not specifying anything:
  if(GBM) X_t[1] <- log.futures[1,1]
  #But if we are:
  if("X.0_1" %in% names(params)) X_t <- matrix(sapply(1:N.factors, FUN = function(x) if(paste0("X.0_",x) %in% names(params)) params[paste0("X.0_",x)]))


  #Things to calculate before the iteration:
  d <- params["E"] + A_T(params, TTM)

  ##If Aggregate (Contract) Data, matrix Z is time (in)homogenous
  if(Contract.Data){
    Z <- array(NA, dim = c(N.obs, N.contracts, N.factors))
    Z[,,1:N.factors] <- sapply(1:N.factors, FUN = function(X) exp(- params[paste0("kappa_", X)] * TTM))
  } else {
    Z <- as.matrix(sapply(1:N.factors, FUN = function(X) exp(- params[paste0("kappa_", X)] * TTM)))

    if(N.contracts == 1) Z = t(as.matrix(sapply(1:N.factors, FUN = function(X) exp(- params[paste0("kappa_", X)] * TTM))))
  }


  #G:
  G_t <- diag(N.factors)
  diag(G_t) <- sapply(1:N.factors, FUN = function(X) exp(- params[paste0("kappa_", X)] * dt))

  #c:
  c_t <- matrix(c(params["mu"] * dt, rep(0, N.factors-1)))

  #Variance of omega:
  # P_t <- Q_t <- cov_func(params, dt)
  Q_t <- cov_func(params, dt)
  P_t <- diag(100, N.factors)

  #Measurement Errors matrix:
  H <- diag(N.contracts)
  if("sigma.contracts" %in% names(params)){
    diag(H) <- params["sigma.contracts"]^2
  } else {
    if(!paste0("sigma.contract_", ncol(log.futures)) %in% names(params)) stop("Not enough Individual Contract white noises have been specified")
    diag(H) <- sapply(1:N.contracts, FUN = function(x) params[paste0("sigma.contract_",x)]^2)
  }


  ##If not verbose, run it with the fkf.SP function (faster):
  if(!verbose){

    ##These are required structures for the fkf.SP inputs:
    if(Contract.Data){
    Zt <- array(NA, dim = c(N.contracts, N.factors, N.obs))
    for(i in 1:N.factors) Zt[,i,] <- t(Z[,,i])
    d <- t(d)
    } else {
    Zt <- Z
    }
      log_likelihood <- suppressWarnings(FKF.SP::fkf.SP(a0 = c(X_t),
                   P0 = P_t,
                   dt = c_t,
                   ct = d,
                   Tt = G_t,
                   Zt = Zt,
                   HHt = Q_t,
                   GGt = as.matrix(diag(H)),
                   yt = t(log.futures)))

    ## If The model was poorly specified, the log-likelihood returns NA. We need to return a heavily penalised score for the gradient function.
    return(ifelse(is.na(log_likelihood),stats::runif(1, -1.5e6, -1e6), log_likelihood))
  } else {
    ##Kalman filter in R

    #Variables to save:
    save_X <- matrix(0, nrow = N.obs, ncol = N.factors)
    save_X_SD <- matrix(0, nrow = N.obs, ncol = N.factors)
    save_V <- matrix(NA, nrow = N.obs, ncol = ncol(log.futures))
    save_Y <- matrix(NA, nrow = N.obs, ncol = ncol(log.futures))
    rownames(save_X) <- rownames(save_X_SD) <- rownames(save_V) <- rownames(save_Y) <- rownames(log.futures)

    if(debugging){
      max.obs <- max(rowSums(!is.na(log.futures)))

      save_P <- array(NA, dim = c(N.factors, N.factors, N.obs))
      save_F <- array(NA, dim = c(max.obs, max.obs, N.obs))
      save_K <- array(NA, dim = c(N.factors, max.obs, N.obs))
      LL_t <- rep(0, N.obs)
    }

    #Initialize
    log_likelihood <- 0
    I <- diag(N.factors)
    converged <- FALSE

    if(!Contract.Data){
      d_t <- d
      Z_t <- Z
    }

    #####BEGIN Kalman filter:
    for(t in 1:N.obs){

      ##Transition Equation:
      X_tp1 <- c_t + G_t %*% X_t

      ##Covariance Transition Equation:
      P_tp1 <- G_t %*% P_t %*% t(G_t) + Q_t

      ##How many contracts are we observing this iteration?
      Observed.Contracts <- which(!is.na(log.futures[t,]))
      if(length(Observed.Contracts)>0){
        N.Observed.contracts <- length(Observed.Contracts)

      ##Time Inhomogeneous - Update:
      if(Contract.Data){
        #Measurement Errors matrix:
        #Expectation of epsilon_t is 0
        #Variance    of epsilon_t is H
        H <- diag(N.Observed.contracts)
        if("sigma.contracts" %in% names(params)){
          diag(H) <- params["sigma.contracts"]^2
        } else {
          diag(H) <- sapply(Observed.Contracts, FUN = function(x) params[paste0("sigma.contract_",x)]^2)
          ##If the model's estimated something to push it to zero, set it to zero:
          diag(H)[which(diag(H)<1.01e-10)] <- 0
        }

        #Step 1 - Calculate Required Values:

        #d
        d_t <- matrix(d[t,Observed.Contracts])

        #Z:
        Z_t <- as.matrix(Z[t,Observed.Contracts,])
        if(length(Observed.Contracts)==1) Z_t <- t(Z[t,Observed.Contracts,])
      }

      if(Contract.Data || !converged){

        #Function of covariance matrix:
        F_t  <- Z_t %*% P_tp1 %*% t(Z_t) + H

        det_F_t <- suppressWarnings(log(det(F_t)))

        ##Numeric Stability - Poorly Conditioned params:
        inverse_F_t <- try(solve(F_t))
        if(is.na(det_F_t))                         stop("Negative determinant in Kalman filter Covariance matrix. Theta may be poorly specified.")
        if(any(class(inverse_F_t) == "try-error")) stop("Singular Kalman filter Covariance Matrix. Theta may be poorly specified.")

        #Kalman Gain:
        K_t <- P_tp1 %*% t(Z_t) %*% inverse_F_t
        P_tp1 <- (I - K_t %*% Z_t) %*% P_tp1

        ###Check if the values have converged, if so, we can increase the efficiency of the algorithm:
        if(!converged && t > 3) converged <- abs(sum(P_tp1 - P_t)) < 1e-07

        P_t <- P_tp1
        convergance_time <- t

      }

        ##Measurement Equation:
        y_bar_t <- d_t + Z_t %*% X_tp1
        #Actual Futures Prices:
        y_t <- as.matrix(log.futures[t,Observed.Contracts])
        #Prediction Error:
        v_t <- y_t - y_bar_t
        #Correction based upon prediction error:
        X_tp1 <- X_tp1 + K_t %*% v_t
        ###Update, iteration begins anew:
        X_t <- X_tp1
        P_t <- P_tp1

        #Update Concurrent Log Likelihood of Observations:
        log_likelihood <- sum(log_likelihood, - (1/2) * sum(N.Observed.contracts * log(2*pi), det_F_t, t(v_t) %*% inverse_F_t %*% v_t))
        #-----------------------

        ##Verbose Saving:

        #Record our estimated variables
        #Updated Error terms:
        y_tt <- d_t + Z_t %*% X_t
        v_tt <- y_tt - y_t

        save_X[t,] <- X_t
        save_X_SD[t,] <- diag(P_t)
        save_Y[t,Observed.Contracts] <- y_tt
        save_V[t,Observed.Contracts] <- v_tt
        log_likelihood_result <- log_likelihood
        if(debugging){

          save_P[,,t] <- P_t
          save_F[1:N.Observed.contracts,1:N.Observed.contracts,t] <- F_t
          save_K[,1:N.Observed.contracts,t] <- K_t

          LL_t[t] <- log_likelihood
        }
      }
    }

    #Final observations, for forecasting purposes:
    X.t <- c(X_t)
    names(X.t) <- paste0("X.", 1:N.factors, "_t")

    ####Term Structure Analysis:
    #Save the filtered Observations:
    Y_output <- exp(cbind(params["E"] + rowSums(save_X),save_Y))
    if(!is.null(colnames(log.futures))) colnames(Y_output) <- c("filtered Spot", colnames(log.futures))
    colnames(save_X) <- colnames(save_X_SD) <- paste("Factor", 1:N.factors)
    rownames(Y_output) <- rownames(save_X)

    ###Term Structure Fit:
    Term_Structure_Fit <- matrix(0, nrow = 4, ncol = ncol(log.futures))

    ##Mean Error:
    Term_Structure_Fit[1,] <- colMeans(save_V, na.rm = TRUE)
    ##Mean Absolute Error
    Term_Structure_Fit[2,] <- colMeans(abs(save_V), na.rm = TRUE)
    ##SD of Error:
    Term_Structure_Fit[3,] <- apply(save_V, MARGIN = 2, FUN = function(x) stats::sd(x, na.rm = TRUE))
    ##RMSE of each contract:
    Term_Structure_Fit[4,] <- sqrt(colMeans(save_V^2, na.rm = TRUE))

    rownames(Term_Structure_Fit) <- c("Mean Error", "Mean Absolute Error", "SD Error", "RMSE")
    colnames(Term_Structure_Fit) <- colnames(save_V) <- colnames(log.futures)

    ### Term structure fit to all observations:
    Filtered_Error <- c(`High Bias` = mean(save_V[save_V > 0], na.rm = TRUE), `Low Bias` =  mean(save_V[save_V < 0], na.rm = TRUE), `Bias` = mean(save_V, na.rm = TRUE),
                            `RMSE` = sqrt(mean(save_V^2, na.rm = TRUE)))


    ###Volatility TSFit:
    if(Contract.Data) {
      Volatility_TSFit <- TSFit.Volatility(params, exp(log.futures), TTM[nrow(TTM),], dt)
    } else {
      Volatility_TSFit <- TSFit.Volatility(params, exp(log.futures), TTM, dt) }


    ##Verbose List
    output = list(LL = log_likelihood, X.t = X.t, X = save_X, Y = Y_output,
                  V = save_V, Filtered.Error = Filtered_Error, TSFit.Error = Term_Structure_Fit, TSFit.Volatility = Volatility_TSFit)

    ##Debugging List:
    if(debugging) output <- c(output, list(LL_t = LL_t, P_t = save_P, F_t = save_F, K_t = save_K, d = d, Z = Z, G_t = G_t, c_t = c_t, Q_t = Q_t, H = H))

    #Return Output value:
    return(output)
  }
  }

