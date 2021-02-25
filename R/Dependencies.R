#' Specify parameters of N-factor model
#'
#' @description
#'\loadmathjax
#' the \code{NFCP.Parameters} function specifies the parameters of
#' a commodity pricing model under the N-factor framework first described by Cortazar and Naranjo (2006).
#' This function is a recommended starting position for the application of N-factor models within the \code{NFCP} package.
#'
#' @param N.factors \code{numeric}. Number of state variables in the spot price process.
#' @param GBM \code{logical}. If \code{GBM = T}, factor 1 of the model is assumed to follow a Brownian Motion, inducing a unit-root in the spot price process.
#' @param Initial.State \code{logical}. If \code{Initial.State = T}, the initial state vector is specified as unknown parameters of the commodity pricing model.
#' @param S.Constant \code{logical}. If \code{S.Constant = T}, the white noise of observable contracts are assumed and identical (and independent).
#' @param N.contracts \code{numeric}. The number of individual observation white noise terms when \code{S.Constant = F}. Optional when \code{S.Constant = T}. When \code{N.contracts = 0},
#' the value of \code{NFCP.Parameters} returns a vector without any "sigma.contract" or "sigma.contracts" elements.
#' @param verbose \code{logical}. If \code{verbose = T}, the specified N-factor model is printed when the function is called.
#'
#'
#'@details
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
#'When \code{N.contracts = 0}, "sigma.contract" parameters are not returned within the parameter vector.
#'
#'\bold{Diffuse Assumption}:
#'If \code{Initial.State = F}, a 'diffuse' assumption is made within Kalman filtering and parameter estimation functions (See \code{NFCP.MLE} or \code{NFCP.Kalman.filter} for more information)
#'
#'
#'@return A vector of parameter names for a specified N-factor spot price process. This vector is ideal for application within many other functions within the \code{NFCP} package
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'##Generate parameter of a Two-factor model Crude Oil model
#'##as first presented by Schwartz and Smith (2000):
#'SS.Oil.Two.Factor.parameters <- NFCP.Parameters(N.factors = 2,
#'                                                       GBM = TRUE,
#'                                                       Initial.State = FALSE,
#'                                                       S.Constant = FALSE,
#'                                                       N.contracts = 5)
#'print(SS.Oil.Two.Factor.parameters)
#'@export
NFCP.Parameters <- function(N.factors, GBM, Initial.State, S.Constant, N.contracts = NULL, verbose = TRUE){

  ##To start, specify the GBM growth factor:
  input_names <- "mu"

  ##Classic spaghetti code - Append the mean-reverting (kappa), risk-premiums (lambda) and volatilities of each factor:
  variables <- c("lambda", "kappa", "sigma")
  input_names <- c(input_names, sapply(variables, FUN = function(X) paste(X, 1:N.factors, sep ="_")))
  ##The correlations:
  for(i in 1:N.factors) for(j in i:N.factors) if(i != j)  input_names <- c(input_names, paste("rho", i, j, sep = "_"))
  ##Measurement Error (white noise of contracts):
  if(!S.Constant){
    sigma_e <- paste0("sigma.contract_", 1:N.contracts)
  } else {
    sigma_e <- "sigma.contracts"
  }
  input_names <- c(input_names, sigma_e)

  if(Initial.State)  input_names <- c(paste("X.0", 1:N.factors, sep = "_"), input_names)

  #####Model Initialization:
  ###If the first factor is a GBM Process, like in Schwartz and Smith (2000), etc.:
  if(GBM){
    input_names <- input_names[!input_names %in% "kappa_1"]
    ##Are we calculating risk premium directly (low estimation power) or the risk-adjusted growth rate (high estimation power)?
    input_names[input_names %in% "lambda_1"] = "mu_star"
  } else {
    #We don't estimate a growth rate, and in the KF function we'll set it to zero:
    input_names[input_names %in% "mu"] = "E"
  }

  if(verbose){
  #####Print the details:
  cat("\n----------------------------------------------------------------\n")
  cat(paste0(N.factors, ifelse(N.factors>1, " Factors", " Factor"), " model:", ifelse(GBM, " first factor is a GBM", " first factor is an MR"),  "\n\n"))
  cat("Risk Neutral SDE: \n\n")
  if(GBM) cat("log(S_t) = sum(x_t)\n\n") else cat("log(S_t) = E + sum(x_t)\n\n")
  cat("Where: \n")
  if(GBM){
    cat("d x1_t    = mu_star * dt  + sigma_1 * dW_1\n")
  }
  if(!GBM || (GBM && N.factors > 1)){
    for(i in ifelse(GBM, 2, 1): N.factors){
      cat(paste0("d x",i,"_t    = - (lambda_",i," + kappa_", i," * x",i,"_t) * dt  + sigma_",i," * dW_",i,"\n"))
    }
  }
  if(N.factors > 1){
    cat("\n And: \n\n")
    for(i in 1:N.factors) for(j in i:N.factors) if(i!=j) cat(paste0("E(dW_",i," * dW_",j,") = rho_",i,"_",j," * dt\n"))
  }
  }

  if(!is.null(N.contracts)) if(N.contracts == 0) input_names <- input_names[!grepl("contract", input_names)]

  return(c(input_names, use.names = FALSE))
}


#'Stitch Futures Contracts
#'@description Aggregate futures contract price data by stitching according to either approximate maturities and rollover frequency or contract number from closest maturity.
#'
#'@param Futures Contract futures price data. Each row of \code{Futures} should represent one observation of futures prices and each column should represent one quoted futures contract. NA's in \code{Futures} are allowed, representing missing observations.
#'@param TTM A \code{vector} of contract maturities to stitch
#'@param maturity.matrix The time-to-maturity (in years) for each contract at each given observation point. The dimensions of \code{maturity.matrix} should match those of \code{Futures}
#'@param rollover.frequency the frequency (in years) at which contracts should be rolled over
#'@param Contract.Numbers A \code{vector} of contract numbers offset from the closest-to-maturity contract at which to stitch contracts.
#'@param verbose \code{logical}. Should additional information be output? see \bold{details}
#'
#'@details
#'This function aggregates a set of futures contract data by stitching contract data over an observation period, resulting in a
#'set of futures observations that is 'complete' (ie. Does not feature missing observations). Aggregated futures
#'data benefit from several computational efficiencies compared to raw contract data, but results in the loss of futures price information.
#'
#'There are two methods of the \code{Stitch.Contracts} function that can be utilized the stitch contracts:
#'
#'\bold{Method 1}
#'
#'\code{Stitch.Contracts(Futures, Contract.Numbers, verbose = T)}
#'Futures data may be aggregated by stitching prices according to maturity matching. This method requires the inputs \code{TTM}, \code{maturity.matrix} and \code{rollover.frequency}.
#'This method stitched contracts by matching the observation prices according to which contract has the closest time-to-maturity of the desired maturity specified
#'in \code{TTM}. Contracts are rolled over at the frequency specified in \code{rollover.frequency}.
#'
#'\bold{Method 2}
#'
#'\code{Stitch.Contracts(Futures, TTM, maturity.matrix, rollover.frequency, verbose = T)}
#'Futures data may be stitched according to the contract numbers offset from the closest-to-maturity contract. This method requires only the
#'input \code{Contract.Numbers} specifying which contracts should be included. This method is most appropriate when the maturity of available
#'contracts are consistent (ie. contracts expire every month or three months).
#'
#'@return
#'\code{Stitch.Contracts} returns a matrix of stitched futures prices if \code{verbose = T} and a list with two or three objects otherwise (see below).
#'
#'\tabular{ll}{
#'
#' \code{Prices} \tab A data frame of Stitched futures prices. Each row represents an observation of the specified contracts. \cr
#'
#' \code{Maturities} \tab A data frame of the time-to-maturity of observed futures prices. Each row represents an observation of the
#'specified contracts. Returned only when \bold{Method 1} is used (see \bold{Details})  and \code{verbose = T}. \cr
#'
#' \code{Tickers} \tab  A data frame of the named columns of observed futures prices (e.g. contract tickers). Returned only when \code{Futures} or \code{maturity.matrix} have named columns and \code{verbose = T}. \cr
#' }
#'
#'@examples
#'##These examples approximately replicate the Crude Oil data utilized within the
#'##prominent work of Schwartz and Smith (2000):
#'
#'###Method 1 - Stitch crude oil contracts according to maturity matching:
#'SSOilStitched.M1 <- Stitch.Contracts(Futures = SS.Oil$Contracts,
#'TTM = c(1, 5, 9, 13, 17)/12, maturity.matrix = SS.Oil$Contract.Maturities,
#'rollover.frequency = 1/12, verbose = TRUE)
#'
#'###Method 2 - Stitch crude oil contracts according to nearest contract numbers:
#'SSOilStitched.M2 <- Stitch.Contracts(Futures = SS.Oil$Contracts,
#'Contract.Numbers = c(1, 5, 9, 13, 17), verbose = TRUE)
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@export
Stitch.Contracts <- function(Futures, TTM = NULL, maturity.matrix = NULL, rollover.frequency = NULL, Contract.Numbers = NULL, verbose = FALSE){

  Futures <- as.matrix(Futures)

  ###Rolling Over by trying to match maturities:
  if(is.null(Contract.Numbers)){
    TTM <- as.matrix(TTM)


    Stitched.Prices <- matrix(NA, nrow = nrow(Futures), ncol = length(TTM))
    colnames(Stitched.Prices) <- TTM
    rownames(Stitched.Prices) <- rownames(Futures)
    if(!is.null(colnames(Futures))){
      Stitched.Tickers <- matrix(NA, nrow = nrow(Futures), ncol = length(TTM))
      colnames(Stitched.Tickers) <- TTM
      rownames(Stitched.Tickers) <- rownames(Futures)
    }

    maturity.matrix <- as.matrix(maturity.matrix)
    ##Get rid of data we don't have Trading.Dates for:
    maturity.matrix[is.na(Futures)] <- NA


    Stitched.Maturities <- matrix(NA, nrow = nrow(Futures), ncol = length(TTM))
    colnames(Stitched.Maturities) <- TTM
    rownames(Stitched.Maturities) <- rownames(Futures)

    if(is.null(maturity.matrix) || is.null(rollover.frequency)) stop("Stitching by maturity matching requires both a maturity matrix and rollover frequency")

    ###We need a function that evaluates the time to maturities and picks the most optimal ones:
    i <- row <- start_row <- stitch_row <- 1
    stitch_fill_row <- 0

    while(stitch_fill_row < nrow(Futures)){

      ###Which columns/contracts?
      colnum <-  apply(data.frame(TTM), MARGIN = 1, FUN = function(x) which.min(abs(maturity.matrix[row,] - x)))

      ###Fill out observations until rollover:
      fill_row <- which.min(abs(maturity.matrix[,colnum[1]] - (maturity.matrix[row,colnum[1]] - rollover.frequency)))

      ###The row that you go to in the stitched data:
      stitch_fill_row <- min(nrow(Futures), fill_row - (start_row - 1))

      stitch_fill <- stitch_row:stitch_fill_row

      Stitched.Prices    [stitch_fill,] <- as.matrix(Futures[stitch_fill,colnum])
      Stitched.Maturities[stitch_fill,] <- as.matrix(maturity.matrix[stitch_fill,colnum])
      if(!is.null(colnames(Futures))) Stitched.Tickers[stitch_fill,] <- colnames(Futures)[colnum]

      ###We've stitched, now repeat
      ###Where we start from:
      row <- fill_row + 1

      ###Where you start in your stitch rows:
      stitch_row <- stitch_fill_row + 1
      ###We move to the next maturity contract:
    }
    if(anyNA(Stitched.Prices)) warning("NA's returned in Output Data")

    if(verbose){
      if(!is.null(colnames(Futures))) return(list(Prices = Stitched.Prices, Maturities = Stitched.Maturities, Tickers = Stitched.Tickers))
      return(list(Prices = Stitched.Prices, Maturities = Stitched.Maturities))
    } else {
      return(Stitched.Prices)
    }
  } else {


    Stitched.Prices <- matrix(NA, nrow = nrow(Futures), ncol = length(Contract.Numbers))
    colnames(Stitched.Prices) <- TTM
    rownames(Stitched.Prices) <- rownames(Futures)
    if(!is.null(colnames(Futures))){
      Stitched.Tickers <- matrix(NA, nrow = nrow(Futures), ncol = length(Contract.Numbers))
      colnames(Stitched.Tickers) <- TTM
      rownames(Stitched.Tickers) <- rownames(Futures)
    }

    ###We need a function that evaluates the time to maturities and picks the most optimal ones:
    start_row <- fill_row <- 1
    ###Which column do we start with?
    NA_col <- which(!is.na(Futures[1,]))[1]
    colnum <- Contract.Numbers + NA_col - 1

    for(t in 1:nrow(Futures)){
      Stitched.Prices[t,] <- Futures[t, which(!is.na(Futures[t,]))[colnum]]
      if(!is.null(colnames(Futures))) Stitched.Tickers[t,] <- colnames(Futures)[which(!is.na(Futures[t,]))[colnum]]
    }

    colnames(Stitched.Prices) <- paste0("F", Contract.Numbers)
    rownames(Stitched.Prices) <- rownames(Futures)
    Stitched.Prices <- as.data.frame(Stitched.Prices)

    if(!is.null(colnames(Futures))){
      colnames(Stitched.Tickers) <- paste0("F", Contract.Numbers)
      rownames(Stitched.Tickers) <- rownames(Futures)
      Stitched.Tickers <- as.data.frame(Stitched.Tickers)
    }

    if(anyNA(Stitched.Prices)) warning("NA's returned in Output Data")

    if(verbose){
      return(list(Prices = Stitched.Prices, Tickers = Stitched.Tickers))
    } else {
      return(Stitched.Prices)

    }
  }
}

#Applied within various functions throughout the rest of the package:

#Function for A(T):
#'Calculate \eqn{A(T)}
#'@description
#'Calculate the values of \eqn{A(T)} for a given N-factor model parameters and observations. Primarily purpose is for application within other functions of the \code{NFCP} package.
#'
#'@param parameters A named vector of parameters of an N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param Tt A vector or matrix of the time-to-maturity of observed futures prices
#'
#'@return A matrix of identical dimensions to \eqn{T} providing the values of function \eqn{A(T)} of a given N-factor model and observations.
#'
#'@details
#'
#'\loadmathjax
#'Under the assumption that Factor 1 follows a Brownian Motion, \eqn{A(T)} is given by:
#'\mjdeqn{A(T) = \mu^*T-\sum_{i=1}^N - \frac{1-e^{-\kappa_i T}\lambda_i}{\kappa_i}+\frac{1}{2}(\sigma_1^2T +
#'\sum_{i.j\neq 1} \sigma_i \sigma_j \rho_{i,j} \frac{1-e^{-(\kappa_i+\kappa_j)T}}{\kappa_i+\kappa_j})}{A(T) = mu^* * T - sum_{i=1}^N (1-e^(-kappa[i] T)lambda[i])/(kappa[i]) + 1/2 (sigma[1]^2 * T)
#' + sum_{i.j != 1} sigma[i] sigma[j] rho[i,j] (1 - e^(-(kappa[i] + kappa[j]) * T)) / (kappa[i] + kappa[j])}
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'##Calculate time homogeneous values of A(T) for the
#'##Schwartz and Smith (2000) two-factor model:
#'SS.Oil.A_T <- A_T(SS.Oil$Two.Factor, SS.Oil$Stitched.TTM)
#'
#'@export
A_T <- function(parameters, Tt){

  if(is.null(names(parameters))) stop("argument parameters must be a named vector")

  N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
  output <- diag(N.factors)

  if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0
  GBM <- parameters["kappa_1"] == 0

  if(GBM) output <- (parameters["mu_star"] + 0.5 * parameters["sigma_1"]^2) * Tt else output <- 0

  if(!(GBM && N.factors == 1)){
    for(i in ifelse(GBM, 2, 1):N.factors){
      output <- output - (1 - exp(-parameters[paste0("kappa_", i)] * Tt)) * (parameters[paste0("lambda_", i)]/parameters[paste0("kappa_", i)])
    }

    for(i in 1:N.factors){
      for(j in 1:N.factors){
        if(!(i == 1 && i == j)){
          kappa_sums <- sum(parameters[paste0("kappa_", j)], parameters[paste0("kappa_", i)])
          output <- output + 0.5 * parameters[paste0("sigma_",i)] * parameters[paste0("sigma_",j)] * ifelse(i==j, 1, parameters[paste("rho",min(i,j), max(i,j), sep = "_")]) * (1 - exp(-(kappa_sums * Tt)))/kappa_sums
        }
      }
    }}
  return(output)
}

#'N-factor Covariance:
#'
#'@description
#'\loadmathjax
#'Calculate the covariance matrix of state variables for a given N-factor model parameters and discrete time step.
#'
#'@param parameters a named vector of parameters of an N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param dt a discrete time step
#'
#'@details The primary purpose of the \code{model_covariance} function is to be called within other functions of the \code{NFCP} package. The covariance of an N-factor model is given by:
#'
#'
#'\mjdeqn{cov_{1,1}(x_{1,t},x_{1,t}) = \sigma_1^2t}{cov[1,1](x[1,t],x[1,t]) = sigma[1]^2 * t}
#'\mjdeqn{cov_{i,j}(x_{i,t},x_{j,t}) = \sigma_i\sigma_j\rho_{i,j}\frac{1-e^{-(\kappa_i+\kappa_j)t}}{\kappa_i+\kappa_j}}{cov[i,j](x[i,t],x[j,t]) = sigma[i] sigma[j] rho[i,j] (1-e^(-(kappa[i]+kappa[j])t ) / (kappa[i] + kappa[j])}
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@return A \code{matrix} object with dimensions \mjeqn{N \times N}{N X N}, where \eqn{N} is the number of factors of the specified N-factor model.
#'@examples
#'#Calculate the covariance matrix of a two-factor model over one discrete (weekly) time step:
#'SS.Oil.covariance <- cov_func(SS.Oil$Two.Factor, SS.Oil$dt)
#'@export
cov_func <- function(parameters, dt){

  if(is.null(names(parameters))) stop("argument parameters must be a named vector")

  N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
  output <- diag(N.factors)

  if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0
  GBM <- parameters["kappa_1"] == 0

  for(i in 1:N.factors){
    for(j in 1:N.factors){
      kappa_sums <- sum(parameters[paste0("kappa_", j)], parameters[paste0("kappa_", i)])
      output[i,j] <- parameters[paste0("sigma_",i)] * parameters[paste0("sigma_",j)] * ifelse(i==j, 1, parameters[paste("rho",min(i,j), max(i,j), sep = "_")]) * (1 - exp(-(kappa_sums * dt)))/kappa_sums
    }
  }
  if(GBM) output[1,1] <- parameters["sigma_1"]^2 * dt

  return(output)
}

#'N-Factor MLE Search Boundaries
#'@description Generate boundaries for the domain of parameters of the N-factor model for parameter estimation.
#'
#'@param parameters a vector of parameter names of an N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param kappa A vector of length two specifying the lower and upper bounds for the 'kappa' parameter
#'@param lambda A vector of length two specifying the lower and upper bounds for the 'lambda' parameter
#'@param sigma A vector of length two specifying the lower and upper bounds for the 'sigma' parameter
#'@param mu A vector of length two specifying the lower and upper bounds for the 'mu' parameter
#'@param mu_star A vector of length two specifying the lower and upper bounds for the 'mu_star' parameter
#'@param rho A vector of length two specifying the lower and upper bounds for the 'rho' parameter
#'@param contract A vector of length two specifying the lower and upper bounds for the 'contract' parameter
#'@param X.0 A vector of length two specifying the lower and upper bounds for the 'X.0' parameter
#'@param E A vector of length two specifying the lower and upper bounds for the 'E' parameter
#'
#'@details
#'The \code{NFCP.Domains} function generates lower and upper bounds for the parameter estimation procedure in the format required of the 'Domains' argument of the 'genoud' function. \code{NFCP.Domains}
#'allows easy setting of custom boundaries for parameter estimation, whilst also providing default domains of parameters.
#'
#'@return
#'A matrix of defaulted domains for the given unknown parameters. The first column corresponds to the lower bound of the
#'allowable search space for the parameter, whilst the second column corresponds to the upper bound. These values were set to allow for the
#''realistic' possible values of given parameters as well as restricting some parameters (such as variance and mean-reverting terms) from taking
#'negative values. The format of the returned matrix matches that required by the \code{Domains} argument of the \code{Genoud} function from the package \code{RGenoud}.
#'
#'@references
#'
#'Mebane, W. R., and J. S. Sekhon, (2011). Genetic Optimization Using Derivatives: The rgenoud Package for R.
#'\emph{Journal of Statistical Software}, 42(11), 1-26. URL http://www.jstatsoft.org/v42/i11/.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'
#'@examples
#'##Specify the Schwartz and Smith (2000) two-factor model
#'##with constant contract white noise:
#'SS.parameters <- NFCP.Parameters(N.factors = 2,
#'                                 GBM = TRUE,
#'                                 Initial.State = TRUE,
#'                                 S.Constant = TRUE)
#'
#'###Generate the default 'domains' argument of 'NFCP.MLE' function:
#'NFCP.MLE.Bounds <- NFCP.Domains(SS.parameters)
#'@export
NFCP.Domains <- function(parameters,
                        kappa = NULL,
                        lambda = NULL,
                        sigma = NULL,
                        mu = NULL,
                        mu_star = NULL,
                        rho = NULL,
                        contract = NULL,
                        X.0 = NULL,
                        E = NULL){

  if(is.null(kappa))    kappa    <- c(1e-5, 50)
  if(is.null(lambda))   lambda   <- c(-10, 10)
  if(is.null(sigma))    sigma    <- c(0, 10)
  if(is.null(mu))       mu       <- c(-10, 10)
  if(is.null(mu_star))  mu_star  <- c(-10, 10)
  if(is.null(rho))      rho      <- c(-1, 1)
  if(is.null(contract)) contract <- c(1e-10, 1)
  if(is.null(X.0))      X.0      <- c(-10, 10)
  if(is.null(E))        E        <- c(-10, 10)



  ##The Bounds:
  lower_bounds <- upper_bounds <- rep(0, length(parameters))

  ## Equilibrium Price (GBM = FALSE)
  lower_bounds[parameters == "E"] <- E[1]
  upper_bounds[parameters == "E"] <- E[2]

  ## Risk-premium (lambda)
  lower_bounds[grepl("lambda", parameters)] <- lambda[1]
  upper_bounds[grepl("lambda", parameters)] <- lambda[2]

  ## long-term growth rate
  lower_bounds[grepl("mu", parameters)] <- mu[1]
  upper_bounds[grepl("mu", parameters)] <- mu[2]

  ## kappa
  lower_bounds[grepl("kappa", parameters)] <- kappa[1]
  upper_bounds[grepl("kappa", parameters)] <- kappa[2]

  ## sigma
  lower_bounds[grepl("sigma", parameters)] <- sigma[1]
  upper_bounds[grepl("sigma", parameters)] <- sigma[2]


  ###Setting White noise bounds that are too low will result in a singular matrix (1e-5^2 == 1e-10):

  lower_bounds[grepl("contract", parameters)] <- contract[1]
  upper_bounds[grepl("contract", parameters)] <- contract[2]
  ##Correlation between -1 and 1:
  lower_bounds[grepl("rho", parameters)] <- rho[1]
  upper_bounds[grepl("rho", parameters)] <- rho[2]

  ###Initial Inputs:
  lower_bounds[grepl("X.0_", parameters)] <- X.0[1]
  upper_bounds[grepl("X.0_", parameters)] <- X.0[2]
  bounds <- cbind(lower_bounds, upper_bounds)
  colnames(bounds) <- c("Lower Bound", "Upper Bound")
  rownames(bounds) <- parameters

  return(bounds)
}
