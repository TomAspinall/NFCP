#' N-factor model parameter estimation through the Kalman filter and maximum likelihood estimation
#'
#'@description
#'\loadmathjax
#'The \code{NFCP_MLE} function performs parameter estimation of commodity pricing models under the N-factor framework of Cortazar and Naranjo (2006). It uses term structure futures data and estimates unknown parameters through maximum likelihood estimation.
#' \code{NFCP_MLE} allows for missing observations, a variable number of state variables, deterministic seasonality and a variable number of measurement error terms.
#'
#'@param log_futures \code{matrix}. The natural logarithm of observed futures prices. Each row must correspond to quoted futures prices at a particular date and every column must correspond to a unique futures contract.
#'NA values are allowed.
#'
#'@param dt \code{numeric}. Constant, discrete time step of observations, in years.
#'
#'@param futures_TTM \code{vector} or \code{matrix}.  The time-to-maturity of observed futures contracts, in years, at a given observation date. This time-to-maturity can either be constant (ie. class 'vector') or variable (ie. class 'matrix') across observations.
#'The number of columns of 'futures_TTM' must be identical to the number of columns of object 'log_futures'. The number of rows of object 'futures_TTM' must be either 1 or equal to the number of rows of object 'log_futures'.
#'
#'@param N_factors \code{numeric}. Number of state variables in the spot price process.
#'
#'@param N_season \code{numeric}. The number of deterministic, cyclical seasonal factors to include in the spot price process.
#'
#'@param N_ME \code{numeric}. The number of independent measuring errors of observable futures contracts to consider in the Kalman filter.
#'
#'@param ME_TTM \code{vector}. the time-to-maturity groupings to consider for observed futures prices. The length of \code{ME_TTM} must be equal to the number of 'ME' parameters specified in object 'parameter_names'. The maximum of 'ME_TTM' must be greater than the maximum value of 'futures_TTM'.
#'When the number of 'ME' parameter values is equal to one or the number of columns of object 'log_futures', this argument is ignored.
#'
#'@param GBM \code{logical}. When \code{TRUE}, factor 1 of the model is assumed to follow a Brownian Motion, inducing a unit-root in the spot price process.
#'
#'@param estimate_initial_state  \code{logical}. Should the initial state vector be specified as unknown parameters of the commodity pricing model? These are generally estimated with low precision (see \bold{details}).
#'
#'@param Domains \code{matrix}. An option matrix of two columns specifying the lower and upper bounds for parameter estimation. The 'NFCP_domains' function is recommended. When not specified, the default parameter bounds returned by the 'NFCP_domains' function are used.
#'
#'@param cluster \code{cluster}.	An optional object returned by one of the makeCluster commands in the \code{parallel} package to allow for parameter estimation
#'to be performed across multiple cluster nodes.
#'
#'@param ... additional arguments to be passed into the \code{genoud} genetic algorithm numeric optimization. These can highly influence the maximum likelihood estimation procedure. See \code{help(genoud)}
#'
#'@details
#'
#'The \code{NFCP_MLE} function facilitates parameter estimation of commodity pricing models under the N-factor framework through the Kalman filter and maximum likelihood estimation. \code{NFCP_MLE}
#'uses genetic algorithms through the \code{genoud} function of the \code{rgenoud} package to numerically optimize the log-likelihood score returned from the \code{NFCP_Kalman_filter} function.
#'
#'Parameter estimation of commodity pricing models can involve a large number of observations, state variables and unknown parameters. It also features an objective log-likelihood function that is nonlinear and
#'discontinuous with respect to model parameters. \code{NFCP_MLE} is designed to perform parameter estimation as efficiently as possible, maximizing the likelihood of attaining a global optimum.
#'
#'
#'Arguments passed to the \code{genoud} function can greatly influence estimated parameters as well as computation time and must be considered when performing parameter estimation. All arguments of the \code{genoud} function
#'may be passed through the \code{NFCP_MLE} function.
#'
#'
#'When \code{grad} is not specified, the \code{grad} function from the \code{numDeriv} package is called to
#'approximate the gradient within the \code{genoud} optimization algorithm through Richardsons extrapolation.
#'
#'Richardsons extrapolation is regarded for its ability to improve the approximation of estimation methods,
#'which may improve the likelihood of obtained a global maxmimum estimate of the log-likelihood.
#'
#'The population size can highly influence parameter estimates at the expense of increased computation time. For commodity pricing models with a large number of unknown parameters, large population sizes may be necessary to maximize the estimation process.
#'
#'\code{NFCP_MLE} by default performs boundary constrained optimization of log-likelihood scores and does not allow does not allow for out-of-bounds evaluations within
#'the \code{genoud} optimization process, preventing candidates from straying beyond the bounds provided by argument \code{Domains}.
#'
#'When \code{Domains} is not specified, the default bounds specified by the \code{NFCP_domains} function are used. The size of the search domains of unknown parameters can highly
#'influence the computation time of the \code{NFCP_MLE} function, however setting domains that are too restrictive may result in estimated parameters returned at the upper or lower bounds. Custom search domains can be used
#'through the \code{NFCP_domains} function and subsequently the \code{Domains} argument of this function.
#'
#'Finally, the maximum likelihood estimation process of parameters provides no in-built guarantee that the estimated parameters of commodity models are financially sensible results. When the commodity model has been over-parameterized
#'(i.e., the number of factors N specified is too high) or the optimization algorithm has failed to attain a global maximum likelihood estimate, estimated parameters may be irrational.
#'
#'Evidence of irrational parameter estimates include correlation coefficients that are extremely large (e.g., > 0.95 or < -0.95), risk-premiums or drift terms that are unrealistic, filtered state variables that are unrealistic and extremely large/small mean-reverting terms with associated large standard errors.
#'
#'Irrational parameter estimates may indicate that the number of stochastic factors (i.e., \code{N_factors}) of the model or number of seasonal factors (i.e., \code{N_season}) are too high.
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
#'\mjdeqn{log(S_{t}) = E + season(t) + \sum_{i=1}^N x_{i,t}}{log(S[t]) = E + season(t) + sum_{i=1}^n x[i,t]}
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
#'represented by the log-likelihood, to be directly compared with statistical testing possible through a chi-squared test. The AIC or BIC can also be used to compare models.
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
#'\code{NFCP_MLE} returns a \code{list} with 10 objects. 9 objects are returned when the user has specified not to calculate the hessian matrix at solution.
#'
#'\tabular{ll}{
#'
#'\code{MLE} \tab \code{numeric} The Maximum-Likelihood-Estimate of the solution. \cr
#'
#'\code{estimated_parameters} \tab \code{vector}. Estimated parameters. \cr
#'
#'\code{standard_errors} \tab \code{vector}. Standard error of the estimated parameters. Returned only when \code{hessian = T} is specified.  \cr
#'
#'\code{Information Criteria} \tab \code{vector}. The Akaikie and Bayesian Information Criterion. \cr
#'
#'\code{x_t} \tab \code{vector}. The final observation of the state vector.
#'
#'When deterministic seasonality is considered, it also returns the observation point
#'
#'along the deterministic curve. \cr
#'
#'\code{X} \tab \code{matrix}. Optimal one-step-ahead state vector.
#'
#'When deterministic seasonality is considered, it also returns the observation point
#'
#'along the deterministic curve. \cr
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
#'\code{Term Structure Volatility Fit} \tab \code{matrix}. Theoretical and empirical volatility of observed futures contract returns \cr
#'
#'\code{proc_time} \tab \code{list}. The real and CPU time (in seconds) the \code{NFCP_MLE} function has taken. \cr
#'
#'\code{genoud_value} \tab \code{list}. Outputs of \code{genoud}.
#'
#' }
#'
#'@references
#'
#'Hannan, E. J., et al. (1970). "The seasonal adjustment of economic time series." \emph{International economic review}, 11(1): 24-52.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Sørensen, C. (2002). "Modeling seasonality in agricultural commodity futures." \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products} 22(5): 393-426.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'Mebane, W. R., and J. S. Sekhon, (2011). Genetic Optimization Using Derivatives: The rgenoud Package for R.
#'\emph{Journal of Statistical Software}, 42(11), 1-26. URL http://www.jstatsoft.org/v42/i11/.
#'
#' @examples
#'# Estimate a 'one-factor' geometric Brownian motion model:
#'Oil_1F_estimated_model <- NFCP_MLE(
#'## Arguments
#'log_futures = log(SS_oil$contracts)[1:20,1:5],
#'dt = SS_oil$dt,
#'futures_TTM= SS_oil$contract_maturities[1:20,1:5],
#'N_factors = 1, N_ME = 1,
#'## Genoud arguments:
#'pop.size = 4, print.level = 0, gr = NULL,
#'max.generations = 0)
#'@export
NFCP_MLE <- function(log_futures, dt, futures_TTM, N_factors, N_season = 0, N_ME = 1, ME_TTM = NULL, GBM = TRUE, estimate_initial_state = FALSE, Domains = NULL, cluster = FALSE, ...){

  # ----------------------------------------------------------------------
  # Input checks:
  time_0 <- proc.time()

  ## Standardize format:
  log_futures <- as.matrix(log_futures)
  futures_TTM <- as.matrix(futures_TTM)

  ## "Contract" data or "Aggregate" Data?
  contract_data <- all(dim(futures_TTM)>1)
  if(contract_data && !all(dim(log_futures) == dim(futures_TTM))) stop("log_futures and futures_TTM have different dimensions")
  if(!contract_data && length(futures_TTM)!=ncol(log_futures)) stop("Aggregate futures data, however ncol(log_futures) and length(futures_TTM) have different dimensions")

  ME_TTM_used <- N_ME > 1 && N_ME < ncol(log_futures)
  ## Have enough ME maturity terms been specified?
  if(ME_TTM_used){
    if(is.null(ME_TTM)) stop("Multiple measurement error (ME) terms have been specified but the maturity terms of the measurement error (ME_TTM) have not.")
    if(length(ME_TTM) != N_ME) stop("Number of measurement error (ME) terms specified does not match the length of argument 'ME_TTM'")
    if(max(futures_TTM, na.rm = TRUE) > max(ME_TTM, na.rm = TRUE)) stop("Maximum observed contract maturity (futures_TTM) is greater than the max specified maturity grouping for the measurement error (ME_TTM)")
  }

  ## Unknown Parameters:
  parameters <- NFCP::NFCP_parameters(N_factors, GBM, estimate_initial_state, N_ME, N_season)

  # ----------------------------------------------------------------------
  # This section develops the arguments applied to the genoud function call.

  ## Custom arguments for NFCP specifically:
  solution.tolerance <- 0.1
  max <- TRUE
  boundary.enforcement <- 2
  hessian <- TRUE
  print.level <- 1
  gradient.check <- FALSE

  ##Gradient Function - Richardsons Extrapolation:
  Richardsons_extrapolation <- TRUE
  gr <- function(x,...) numDeriv::grad(func = NFCP::NFCP_Kalman_filter, x, parameter_names = parameters, log_futures = log_futures,
                                       dt = dt, futures_TTM = futures_TTM, ME_TTM = ME_TTM)

  ## Genoud standard parameter values:
  pop.size <- 1000
  starting.values <- project.path <- P9mix <- BFGSfn <- BFGShelp <- NULL
  output.append <- transform <- debug <- balance  <- lexical <- data.type.int <- FALSE
  hard.generation.limit <- MemoryMatrix <- BFGS <- TRUE
  wait.generations <- default.domains <- 10
  max.generations <- 100
  unif.seed <- round(stats::runif(1, 1, 2147483647L))
  int.seed <- round(stats::runif(1, 1, 2147483647L))
  share.type <- instance.number <- 0
  output.path <- "stdout"
  P1 <- P2 <- P3 <- P4 <- P5 <- P6 <- P7 <- P8 <- 50
  P9 <- 0
  BFGSburnin <- 0
  control <- list()
  optim.method <- ifelse(boundary.enforcement < 2, "BFGS","L-BFGS-B")

  fn_args <- match.call()

  # Assign arguments that have been specified in ...
  for(loop in 1:length(fn_args)){
    var <- names(fn_args)[loop]
    if(!var %in% c("log_futures", "dt", "futures_TTM", "N_factors", "N_season", "N_ME", 'ME_TTM', 'GBM','estimate_initial_state', 'cluster', "Domains")){
      # print(var)
      if(!is.null(fn_args[[var]])) assign(var, fn_args[[var]])
    }
  }

  # Base domains:
  if(is.null(Domains)) Domains <- NFCP::NFCP_domains(parameters)


  if(!Richardsons_extrapolation) gr <- NULL

  ##Parallel Processing?
  if(!any(class(cluster)=="cluster" | class(cluster)=="SOCKcluster")) cluster <- FALSE

  # ----------------------------------------------------------------------

  cat("----------------------------------------------------------------
Term Structure Estimation: \n")
  cat(paste(length(parameters), "unknown parameters \n"))
  cat(paste(nrow(log_futures), "observations \n"))
  cat(paste(ncol(log_futures)), "futures contracts \n")
  if(any(class(cluster)=="cluster" | class(cluster)=="SOCKcluster")) cat("parallel processing detected! \n")
  if(ME_TTM_used){
    cat(paste(length(ME_TTM), "measurement error maturity groupings \n"))
  } else {
    cat("\n")
  }

  ##Run the Genetic Algorithm Parameter Estimation:
  NFCP_output <- rgenoud::genoud(NFCP::NFCP_Kalman_filter, nvars = length(parameters), parameter_names = parameters,
                               log_futures = log_futures, dt = dt, futures_TTM = futures_TTM, ME_TTM = ME_TTM,

## Genoud Args:
max=max, pop.size=pop.size, max.generations=max.generations,
wait.generations=wait.generations, hard.generation.limit=hard.generation.limit, starting.values=starting.values,
MemoryMatrix=MemoryMatrix, Domains=Domains, default.domains=default.domains,
solution.tolerance=solution.tolerance, gr=gr, boundary.enforcement=boundary.enforcement, lexical=lexical,
gradient.check=gradient.check, BFGS=BFGS, data.type.int=data.type.int, hessian=hessian,
unif.seed=unif.seed,int.seed=int.seed,print.level=print.level, share.type=share.type,
instance.number=instance.number, output.path=output.path, output.append=output.append, project.path=project.path, P1=P1,
P2=P2, P3=P3, P4=P4, P5=P5, P6=P6, P7=P7, P8=P8, P9=P9, P9mix=P9mix, BFGSburnin=BFGSburnin, BFGSfn=BFGSfn, BFGShelp=BFGShelp,
control=control, optim.method=optim.method, transform=transform, debug=debug, cluster=cluster, balance=balance
)

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

  order_bool <- NFCP_output$generations != max.generations

  ## Sort the estimated parameters in terms of increasing Kappa's
  if(length(parameter_index) > 1 && order_bool){

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
  if(exists("SE")) if(!any(class(SE) == "try-error")) SE <- estimated_parameters[,2]
  estimated_parameters <- as.numeric(estimated_parameters[,1])
  names(estimated_parameters) <- parameters

  #Output List:
  NFCP_list <- list(MLE = NFCP_output$value, estimated_parameters = estimated_parameters)
  if(exists("SE")) NFCP_list <- c(NFCP_list, list(standard_errors = SE))

  return(c(NFCP_list,
           NFCP::NFCP_Kalman_filter(parameter_values = estimated_parameters, parameter_names = parameters, log_futures = log_futures, futures_TTM = futures_TTM, ME_TTM = ME_TTM,
                                  dt = dt, verbose = TRUE)[-1], proc_time = list(proc.time() - time_0), list(genoud_value = NFCP_output)))
}
