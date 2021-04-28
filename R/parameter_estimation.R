###This is a wrapper to use the rgenoud optimization algorithm for parameter estimation:

#' N-Factor Model Parameter Estimation through the Kalman Filter and Maximum Likelihood Estimation
#'
#'@description
#'\loadmathjax
#'the \code{NFCP_MLE} function performs parameter estimation of an n-factor model given observable term structure futures data through maximum likelihood estimation.
#' \code{NFCP_MLE} allows for missing observations as well as constant or variable time to maturity of observed futures contracts.
#'
#'@param log_futures Object of class \code{matrix} corresponding to the natural logarithm of observable futures prices.
#'NA's are allowed within the \code{matrix}. Every column of the matrix must correspond to a particular futures contract,
#'with each row corresponding to a quoted price on a given date.
#'
#'@param futures_TTM Object of class \code{vector} or \code{matrix} that specifies the time to maturity of observed futures contracts.
#'time to maturity can be either constant (i.e. class \code{vector}) or time dependent (i.e. class \code{matrix}).
#'When the time to maturity of observed futures contracts is time dependent, the dimensions of
#'\code{futures_TTM} must be identical to that of \code{log_futures}. Every element of \code{futures_TTM}
#'corresponds to the time to maturity, in years, of a futures contract at a given observation date.
#'
#'@param dt Constant, discrete time step of observations
#'
#'@param N_factors \code{numeric}. Number of state variables in the spot price process.
#'
#' @param N_ME \code{numeric}. The number of independent measuring errors of observable futures contracts to consider in the Kalman filter.
#'
#'@param ME_TTM vector of maturity groupings to consider for observed futures prices. The length of \code{ME_TTM} must be equal to the number of 'ME' parameter values, and the maximum of ME_TTM must be greater than the maximum observed time-to-maturiy of a futures contract. When the number of 'ME' parameter values is equal to one or
#' the total number of contracts (i.e., columns of \code{log_futures}), this argument is optional and not considered. The measurement error of an observation is highly influenced by its time-to-maturity, see \bold{details}.
#'
#'@param GBM \code{logical}. When \code{TRUE}, factor 1 of the model is assumed to follow a Brownian Motion, inducing a unit-root in the spot price process.
#'
#'@param estimate_initial_state  \code{logical}. When \code{TRUE}, the initial state vector is specified as unknown parameters of the commodity pricing model. When \code{FALSE},
#'a "diffuse" assumption is taken instead (see \bold{details})
#'
#'@param Richardsons_extrapolation \code{logical}. When \code{TRUE}, the \code{grad} function from the \code{numDeriv} package is called to
#'approximate the gradient within the \code{genoud} optimization algorithm.
#'
#'@param cluster 	an optional object of the 'cluster' class returned by one of the makeCluster commands in the \code{parallel} package to allow for parameter estimation
#'to be performed across multiple cluster nodes.
#'
#'@param Domains an optional \code{matrix} of the lower and upper bounds for the parameter estimation process. The \code{NFCP_domains} function is highly recommended.
#'When \code{Domains} is not specified, the standard bounds specified within the \code{NFCP_domains} function are used.
#'
#'@param ... additional arguments to be passed into the \code{genoud} function. See \code{help(genoud)}
#'
#'@details
#'
#'\code{NFCP_MLE} is a wrapper function that uses the genetic algorithm optimization function \code{genoud} from the \code{rgenoud}
#'package to optimize the log-likelihood score returned from the \code{NFCP_Kalman_filter} function. When \code{Richardsons_extrapolation = TRUE}, gradients are approximated
#'numerically within the optimization algorithm through the \code{grad} function from the \code{numDeriv} package. \code{NFCP_MLE} is designed
#'to perform parameter estimation as efficiently as possible, ensuring a global optimum is reached even with a large number of unknown parameters and state variables. Arguments
#'passed to the \code{genoud} function can greatly influence estimated parameters and must be considered when performing parameter estimation. Recommended arguments to pass
#'into the \code{genoud} function are included within the vignette of \code{NFCP}. All arguments of the \code{genoud} function may be passed through the \code{NFCP_MLE}
#'function (except for \code{gradient.check}, which is hard set to false).
#'
#'\code{NFCP_MLE} performs boundary constrained optimization of log-likelihood scores and does not allow does not allow for out-of-bounds evaluations within
#'the \code{genoud} optimization process, preventing candidates from straying beyond the bounds provided by \code{Domains}. When \code{Domains} is not specified, the default
#'bounds specified by the \code{NFCP_domains} function are used.
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
#'\bold{Diffuse Kalman Filtering}
#'
#'If \code{estimate_initial_state = F}, a 'diffuse' assumption is used within the Kalman filtering algorithm. Factors that follow an Ornstein-Uhlenbeck are assumed to equal zero. When
#'\code{estimate_initial_state = F} and \code{GBM = T}, the initial value of the first state variable is assumed to equal the first element of \code{log_futures}. This is an
#'assumption that the initial estimate of the spot price is equal to the closest to maturity observed futures price.
#'
#'The initial covariance of the state vector for the Kalman Filtering algorithm assumed to be equal to matrix \mjeqn{Q}{Q}
#'
#'Initial states of factors that follow an Ornstein-Uhlenbeck process are generally not estimated with a high level of precision, due to the transient effect of the initial state vector on future
#'observations, however the initial value of a random walk variable persists across observations (see Schwartz and Smith (2000) for more details).
#'
#'@return
#'\code{NFCP_MLE} returns a \code{list} with 10 objects. 9 objects are returned when the user has specified not to calculate the hessian matrix at solution.
#'
#'\tabular{ll}{
#'
#'\code{MLE} \tab \code{numeric} The Maximum-Likelihood-Estimate of the solution \cr
#'
#'\code{estimated_parameters} \tab \code{vector}. The estimated parameters \cr
#'
#'\code{standard_errors} \tab \code{vector}. Standard error of the estimated parameters. Returned only when \code{hessian = T} is specified  \cr
#'
#'\code{x_t} \tab \code{vector}. The final observation of the state vector \cr
#'
#'\code{X} \tab \code{matrix}. All observations of the state vector, after the updating equation has been applied \cr
#'
#'\code{Y} \tab \code{matrix}. Estimated futures prices at each observation \cr
#'
#'\code{V} \tab \code{matrix}. Estimation error of each futures contracts at each observation \cr
#'
#'\code{Filtered Error} \tab \code{matrix}. The Mean Error (Bias), Mean Absolute Error, Standard Deviation of Error and Root Mean Squared Error (RMSE) of each
#'observed contract, matching the column names of \code{log_futures}  \cr
#'
#'\code{Term Structure Volatility Fit} \tab \code{matrix}. The theoretical and empirical volatility of futures returns for each observed contract as returned from the \code{TSFit.Volatility} function \cr
#'
#'\code{proc_time} \tab \code{list}. The real and CPU time (in seconds) the \code{NFCP_MLE} function has taken. \cr
#'
#'\code{genoud_value} \tab \code{list}. The output of the called \code{genoud} function.
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
#'##Perform One Generation of Maximum Likelihood Estimation on the
#'##first 20 weekly observations of the Schwartz and Smith (2000) Crude Oil Data:
#'SS_2F_estimated_model <- NFCP_MLE(
#'####Arguments
#'log_futures = log(SS_oil$contracts)[1:20,1:5],
#'dt = SS_oil$dt,
#'futures_TTM= SS_oil$contract_maturities[1:20,1:5],
#'fixed_ME = TRUE,
#'N_factors = 1, GBM = TRUE,
#'####Genoud arguments:
#'hessian = TRUE,
#'Richardsons_extrapolation = FALSE,
#'pop.size = 4, optim.method = "L-BFGS-B", print.level = 0,
#'max.generations = 0, solution.tolerance = 10)
#'
#'@export
NFCP_MLE <- function(log_futures, dt, futures_TTM, N_factors, N_ME = 1, ME_TTM = NULL, GBM = TRUE, estimate_initial_state = FALSE,
                                  Richardsons_extrapolation = TRUE, cluster = FALSE, Domains = NULL, ...){

  time.0 <- proc.time()

  ##Standardize format:
  log_futures <- as.matrix(log_futures)
  futures_TTM <- as.matrix(futures_TTM)

  ## "Contract" data or "Aggregate" Data?
  contract_data <- all(dim(futures_TTM)>1)
  if(contract_data && !all(dim(log_futures) == dim(futures_TTM))) stop("log_futures and futures_TTM have different dimensions")
  if(!contract_data && length(futures_TTM)!=ncol(log_futures)) stop("Aggregate futures data, however ncol(log_futures) and length(futures_TTM) have different dimensions")

  ME_TTM_used <- FALSE
  ## Have enough ME maturity terms been specified?
  if(N_ME > 1 && N_ME < ncol(log_futures)){
    ME_TTM_used <- TRUE
    if(is.null(ME_TTM)) stop("Multiple measurement error (ME) terms have been specified but the maturity terms of the measurement error (ME_TTM) have not.")
    if(length(ME_TTM) != N_ME) stop("Number of measurement error (ME) terms specified does not match the length of argument 'ME_TTM'")
    if(max(futures_TTM) > max(ME_TTM)) stop("Maximum observed contract maturity (futures_TTM) is greater than the max specified maturity grouping for the measurement error (ME_TTM)")
  }

  ##Unknown Parameters:
  parameters <- NFCP::NFCP_parameters(N_factors, GBM, estimate_initial_state, N_ME)

  cat("----------------------------------------------------------------
Term Structure Estimation: \n")
  cat(paste(length(parameters), "Unknown Parameters \n"))
  cat(paste(nrow(log_futures), "observations \n"))
  cat(paste(ncol(log_futures)), "futures contracts \n")
  if(ME_TTM_used){
    cat(paste(length(ME_TTM), "maturity groupings used to consider the measurement error of observations \n"))
  } else {
    cat("\n")
  }


  ##Gradient Function?
  if(Richardsons_extrapolation){
    ####Richardsons Extrapolation:
    gr <- function(x,...) numDeriv::grad(func = NFCP::NFCP_Kalman_filter, x, parameter_names = parameters, log_futures = log_futures,
                                        dt = dt, futures_TTM = futures_TTM, ME_TTM = ME_TTM)
  } else gr <- NULL

  ##Domain Boundaries?
  if(is.null(Domains)) Domains <- NFCP::NFCP_domains(parameters)

  ##Parallel Processing?
  if(!any(class(cluster)=="cluster" | class(cluster)=="SOCKcluster")) cluster <- F

  ##Run the Genetic Algorithm Parameter Estimation:
  NFCP_output <- rgenoud::genoud(NFCP_Kalman_filter, nvars = length(parameters), parameter_names = parameters,
                               log_futures = log_futures, dt = dt, futures_TTM = futures_TTM, ME_TTM = ME_TTM,
                               max = T, gr = gr, Domains = Domains,
                               boundary.enforcement = 2, gradient.check = F, cluster = cluster, ...)

  ###Close the cluster:
  if(any(class(cluster)=="cluster" | class(cluster)=="SOCKcluster")) parallel::stopCluster(cluster)

  ####Parameter Estimation Complete:
  estimated_parameters <- matrix(NFCP_output$par)
  rownames(estimated_parameters) <- parameters

  # Standard Errors:
  if("hessian" %in% names(NFCP_output)){
    SE <- suppressWarnings(try(sqrt(abs(diag(solve(- NFCP_output$hessian)))), silent = TRUE))
    if(class(SE)[1] != "try-error"){
      names(SE) <- parameters
      estimated_parameters <- cbind(estimated_parameters, SE)
    }
  }

  ### Which parameters are the Kappa?
  parameter_index <- which(grepl("kappa", parameters))

  if(length(parameter_index) > 1){

  ## Sort the estimated parameters in terms of increasing Kappa's
  Ordered <- order(estimated_parameters[parameter_index,1])

  ## Order the Kappas:
  estimated_parameters[grepl("kappa", parameters),] <- estimated_parameters[grepl("kappa", parameters),][Ordered,]
  ## Order the Sigma's:
  if(GBM){
  estimated_parameters[grepl("sigma_", parameters),][-1,] <- estimated_parameters[grepl("sigma_", parameters),][-1,][Ordered,]
  } else {
    estimated_parameters[grepl("sigma_", parameters),] <- estimated_parameters[grepl("sigma_", parameters),][Ordered,]
  }
  ## Order the lambda's:
  estimated_parameters[grepl("lambda_", parameters),] <- estimated_parameters[grepl("lambda_", parameters),][Ordered,]

  ## Order the Rho's:
  Rho <- estimated_parameters[grepl("rho_", parameters),]
  if(GBM) Ordered <- c(1, Ordered + 1)

  for(iter in 1:(N_factors-1)){
    for(iter_2 in (iter+1):N_factors){
        estimated_parameters[paste("rho", iter, iter_2, sep = "_"),] <- Rho[paste("rho",
                                                                                  min(Ordered[iter], Ordered[iter_2]),
                                                                                  max(Ordered[iter], Ordered[iter_2]), sep = "_"),]
        }
      }
    }

  ## Ordered outputs:
  if(exists("SE") & class(SE)[1] != "try-error") SE <- estimated_parameters[,2]
  estimated_parameters <- as.numeric(estimated_parameters[,1])
  names(estimated_parameters) <- parameters

  #Output List:
  NFCP_list <- list(MLE = NFCP_output$value, estimated_parameters = estimated_parameters)
  if(exists("SE")) NFCP_list <- c(NFCP_list, list(standard_errors = SE))

  return(c(NFCP_list,
           NFCP::NFCP_Kalman_filter(parameter_values = estimated_parameters, parameter_names = parameters, log_futures = log_futures, futures_TTM = futures_TTM, ME_TTM = ME_TTM,
                                  dt = dt, verbose = TRUE)[-1], proc_time = list(proc.time() - time.0), list(genoud_value = NFCP_output)))
}
