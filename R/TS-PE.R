###This is a wrapper to use the rgenoud optimization algorithm for parameter estimation:

#' N-Factor Model Parameter Estimation through the Kalman Filter and Maximum Likelihood Estimation
#'
#'@description
#'\loadmathjax
#'the \code{NFCP.MLE} function performs parameter estimation of an n-factor model given observable term structure futures data through maximum likelihood estimation.
#' \code{NFCP.MLE} allows for missing observations as well as constant or variable time to maturity of observed futures contracts.
#'
#'@param log.futures Object of class \code{matrix} corresponding to the natural logarithm of observable futures prices.
#'NA's are allowed within the \code{matrix}. Every column of the matrix must correspond to a particular futures contract,
#'with each row corresponding to a quoted price on a given date.
#'
#'@param TTM Object of class \code{vector} or \code{matrix} that specifies the time to maturity of observed futures contracts.
#'time to maturity can be either constant (i.e. class \code{vector}) or time dependent (i.e. class \code{matrix}).
#'When the time to maturity of observed futures contracts is time dependent, the dimensions of
#'\code{TTM} must be identical to that of \code{log.futures}. Every element of \code{TTM}
#'corresponds to the time to maturity, in years, of a futures contract at a given observation date.
#'
#'@param dt Constant, discrete time step of observations
#'
#'@param N.factors \code{numeric}. Number of state variables in the spot price process.
#'
#'@param GBM \code{logical}. When \code{TRUE}, factor 1 of the model is assumed to follow a Brownian Motion, inducing a unit-root in the spot price process.
#'
#'@param S.Constant \code{logical}. When \code{TRUE}, the white noise of observable contracts are assumed identical and independent.
#'
#'@param Estimate.Initial.State  \code{logical}. When \code{TRUE}, the initial state vector is specified as unknown parameters of the commodity pricing model. When \code{FALSE},
#'a "diffuse" assumption is taken instead (see \bold{details})
#'
#'@param Richardsons.Extrapolation \code{logical}. When \code{TRUE}, the \code{grad} function from the \code{numDeriv} package is called to
#'approximate the gradient within the \code{genoud} optimization algorithm.
#'
#'@param cluster 	an optional object of the 'cluster' class returned by one of the makeCluster commands in the \code{parallel} package to allow for parameter estimation
#'to be performed across multiple cluster nodes.
#'
#'@param Domains an optional \code{matrix} of the lower and upper bounds for the parameter estimation process. The \code{NFCP.Domains} function is highly recommended.
#'When \code{Domains} is not specified, the standard bounds specified within the \code{NFCP.Domains} function are used.
#'
#'@param ... additional arguments to be passed into the \code{genoud} function. See \code{help(genoud)}
#'
#'@details
#'
#'\code{NFCP.MLE} is a wrapper function that uses the genetic algorithm optimization function \code{genoud} from the \code{rgenoud}
#'package to optimize the log-likelihood score returned from the \code{NFCP.Kalman.filter} function. When \code{Richardsons.Extrapolation = TRUE}, gradients are approximated
#'numerically within the optimization algorithm through the \code{grad} function from the \code{numDeriv} package. \code{NFCP.MLE} is designed
#'to perform parameter estimation as efficiently as possible, ensuring a global optimum is reached even with a large number of unknown parameters and state variables. Arguments
#'passed to the \code{genoud} function can greatly influence estimated parameters and must be considered when performing parameter estimation. Recommended arguments to pass
#'into the \code{genoud} function are included within the vignette of \code{NFCP}. All arguments of the \code{genoud} function may be passed through the \code{NFCP.MLE}
#'function (except for \code{gradient.check}, which is hard set to false).
#'
#'\code{NFCP.MLE} performs boundary constrained optimization of log-likelihood scores and does not allow does not allow for out-of-bounds evaluations within
#'the \code{genoud} optimization process, preventing candidates from straying beyond the bounds provided by \code{Domains}. When \code{Domains} is not specified, the default
#'bounds specified by the \code{NFCP.Domains} function are used.
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
#'the commodity price will revert to.
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
#'\code{var} \mjeqn{\mu}{mu}:  long-term growth rate of the brownian motion process.
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
#'\bold{Diffuse Kalman Filtering}
#'
#'If \code{Estimate.Initial.State = F}, a 'diffuse' assumption is used within the Kalman filtering algorithm. Factors that follow an Ornstein-Uhlenbeck are assumed to equal zero. When
#'\code{Estimate.Initial.State = F} and \code{GBM = T}, the initial value of the first state variable is assumed to equal the first element of \code{log.futures}. This is an
#'assumption that the initial estimate of the spot price is equal to the closest to maturity observed futures price.
#'
#'The initial covariance of the state vector for the Kalman Filtering algorithm assumed to be equal to matrix \mjeqn{Q}{Q}
#'
#'Initial states of factors that follow an Ornstein-Uhlenbeck process are generally not estimated with a high level of precision, due to the transient effect of the initial state vector on future
#'observations, however the initial value of a random walk variable persists across observations (see Schwartz and Smith (2000) for more details).
#'
#'@return
#'\code{NFCP.MLE} returns a \code{list} with 10 objects. 9 objects are returned when the user has specified not to calculate the hessian matrix at solution.
#'
#'\tabular{ll}{
#'
#'\code{MLE} \tab \code{numeric} The Maximum-Likelihood-Estimate of the solution \cr
#'
#'\code{Estimated.Parameters} \tab \code{vector}. The estimated parameters \cr
#'
#'\code{Standard.Errors} \tab \code{vector}. Standard error of the estimated parameters. Returned only when \code{hessian = T} is specified  \cr
#'
#'\code{X.t} \tab \code{vector}. The final observation of the state vector \cr
#'
#'\code{X} \tab \code{matrix}. All observations of the state vector, after the updating equation has been applied \cr
#'
#'\code{Y} \tab \code{matrix}. Estimated futures prices at each observation \cr
#'
#'\code{V} \tab \code{matrix}. Estimation error of each futures contracts at each observation \cr
#'
#'\code{TSFit.Error} \tab \code{matrix}. The Mean Error (Bias), Mean Absolute Error, Standard Deviation of Error and Root Mean Squared Error (RMSE) of each
#'observed contract, matching the column names of \code{log.futures}  \cr
#'
#'\code{TSFit.Volatility} \tab \code{matrix}. The theoretical and empirical volatility of futures returns for each observed contract as returned from the \code{TSFit.Volatility} function \cr
#'
#'\code{proc_time} \tab \code{list}. The real and CPU time (in seconds) the \code{NFCP.MLE} function has taken. \cr
#'
#'\code{genoud.value} \tab \code{list}. The output of the called \code{genoud} function.
#'
#' }
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Mebane, W. R., and J. S. Sekhon, (2011). Genetic Optimization Using Derivatives: The rgenoud Package for R.
#'\emph{Journal of Statistical Software}, 42(11), 1-26. URL http://www.jstatsoft.org/v42/i11/.
#'
#' @examples
#'
#'##Example 1 - Perform One Generation of Maximum Likelihood Estimation on the
#'##first 20 weekly observations of the Schwartz and Smith (2000) Crude Oil Data:
#'##Example 1 - Complete, Stitched Data:
#'Schwartz.Smith.Two.Factor.Stitched <- NFCP.MLE(
#'  ####Arguments
#'  log.futures = log(SS.Oil$Stitched.Futures)[1:5,1],
#'  dt = SS.Oil$dt,
#'  TTM= SS.Oil$Stitched.TTM[1],
#'  S.Constant = FALSE, N.factors = 1, GBM = TRUE,
#'  ####Genoud arguments:
#'  hessian = TRUE,
#'  Richardsons.Extrapolation = FALSE,
#'  pop.size = 4, optim.method = "L-BFGS-B", print.level = 0,
#'  max.generations = 0, solution.tolerance = 10)
#'
#'##Example 2 - Incomplete, Contract Data:
#'Schwartz.Smith.Two.Factor.Contract <- NFCP.MLE(
#'  ####Arguments
#'  log.futures = log(SS.Oil$Contracts)[1:20,1:5],
#'  dt = SS.Oil$dt,
#'  TTM= SS.Oil$Contract.Maturities[1:20,1:5],
#'  S.Constant = TRUE,
#'  N.factors = 1, GBM = TRUE,
#'  ####Genoud arguments:
#'  hessian = TRUE,
#'  Richardsons.Extrapolation = FALSE,
#'  pop.size = 4, optim.method = "L-BFGS-B", print.level = 0,
#'  max.generations = 0, solution.tolerance = 10)
#'
#'
#'@export
NFCP.MLE <- function(log.futures, dt, TTM, N.factors, GBM = TRUE, S.Constant = TRUE, Estimate.Initial.State = FALSE,
                                  Richardsons.Extrapolation = TRUE, cluster = FALSE, Domains = NULL, ...){

  time.0 <- proc.time()

  ##Standardize format:
  log.futures <- as.matrix(log.futures)
  TTM <- as.matrix(TTM)

  ##Dimensions of observations:
  N.obs <- nrow(log.futures)
  N.contracts <- ncol(log.futures)

  ##"Contract" data or "Aggregate" Data?
  Contract.Data <- all(dim(TTM)>1)
  if(Contract.Data && !all(dim(log.futures) == dim(TTM))) stop("Observations and Maturity matrix have different dimensions")
  if(!Contract.Data && length(TTM)!=ncol(log.futures)) stop("Number of observations does not equal number of time homogeneous TTM's")

  ##Unknown Parameters:
  parameters <- NFCP::NFCP.Parameters(N.factors, GBM, Estimate.Initial.State, S.Constant, N.contracts)

  cat("---------------------------------------------------------------- \n")
  cat("Term Structure Estimation: \n")
  cat(paste(length(parameters), "Unknown Parameters \n"))
  cat(paste(nrow(log.futures), "Observations \n"))
  cat(paste(ncol(log.futures)), "Futures Contracts \n")

  ##Gradient Function?
  if(Richardsons.Extrapolation){
    ####Richardsons Extrapolation:
    gr <- function(x,...) numDeriv::grad(func = NFCP::NFCP.Kalman.filter, x, parameters = parameters, log.futures = log.futures,
                                        dt = dt, TTM = TTM)
  } else gr <- NULL

  ##Domain Boundaries?
  if(is.null(Domains)) Domains <- NFCP::NFCP.Domains(parameters)

  ##Parallel Processing?
  if(!any(class(cluster)=="cluster" | class(cluster)=="SOCKcluster")) cluster <- F

  ##Run the Genetic Algorithm Parameter Estimation:
  NFCP.Output <- rgenoud::genoud(NFCP.Kalman.filter, nvars = length(parameters), parameters = parameters,
                               log.futures = log.futures, dt = dt, TTM = TTM,
                               max = T, gr = gr, Domains = Domains,
                               boundary.enforcement = 2, gradient.check = F, cluster = cluster, ...)

  ###Close the cluster:
  if(any(class(cluster)=="cluster" | class(cluster)=="SOCKcluster")) parallel::stopCluster(cluster)

  ####Parameter Estimation Complete:
  Estimated.Parameters <- matrix(NFCP.Output$par)
  rownames(Estimated.Parameters) <- parameters

  # Standard Errors:
  if("hessian" %in% names(NFCP.Output)){
    SE <- suppressWarnings(try(sqrt(abs(diag(solve(- NFCP.Output$hessian)))), silent = TRUE))
    if(class(SE)[1] != "try-error"){
      names(SE) <- parameters
      Estimated.Parameters <- cbind(Estimated.Parameters, SE)
    }
  }

  ### Which parameters are the Kappa?
  parameter_index <- which(grepl("kappa", parameters))

  if(length(parameter_index) > 1){

  ## Sort the estimated parameters in terms of increasing Kappa's
  Ordered <- order(Estimated.Parameters[parameter_index,1])

  ## Order the Kappas:
  Estimated.Parameters[grepl("kappa", parameters),] <- Estimated.Parameters[grepl("kappa", parameters),][Ordered,]
  ## Order the Sigma's:
  if(GBM){
  Estimated.Parameters[grepl("sigma_", parameters),][-1,] <- Estimated.Parameters[grepl("sigma_", parameters),][-1,][Ordered,]
  } else {
    Estimated.Parameters[grepl("sigma_", parameters),] <- Estimated.Parameters[grepl("sigma_", parameters),][Ordered,]
  }
  ## Order the lambda's:
  Estimated.Parameters[grepl("lambda_", parameters),] <- Estimated.Parameters[grepl("lambda_", parameters),][Ordered,]

  ## Order the Rho's:
  Rho <- Estimated.Parameters[grepl("rho_", parameters),]
  if(GBM) Ordered <- c(1, Ordered + 1)

  for(iter in 1:(N.factors-1)){
    for(iter_2 in 2:N.factors){
      if(iter != iter_2){
        Estimated.Parameters[paste("rho", iter, iter_2, sep = "_"),] <- Rho[paste("rho",
                                                                                  min(Ordered[iter], Ordered[iter_2]),
                                                                                  max(Ordered[iter], Ordered[iter_2]), sep = "_"),]
      }
    }}

  }
  ## Ordered outputs:
  if(exists("SE") & class(SE)[1] != "try-error") SE <- Estimated.Parameters[,2]
  Estimated.Parameters <- as.numeric(Estimated.Parameters[,1])
  names(Estimated.Parameters) <- parameters

  #Output List:
  NFCP.List <- list(MLE = NFCP.Output$value, Estimated.Parameters = Estimated.Parameters)
  if(exists("SE")) NFCP.List <- c(NFCP.List, list(Standard.Errors = SE))

  return(c(NFCP.List,
           NFCP::NFCP.Kalman.filter(parameter.values = Estimated.Parameters, parameters = parameters, log.futures = log.futures, TTM = TTM,
                                  dt = dt, verbose = TRUE)[-1], proc_time = list(proc.time() - time.0), list(genoud.value = NFCP.Output)))
}
