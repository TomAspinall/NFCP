#' Specify the constant parameters of an N-factor model
#'
#' @description
#'\loadmathjax
#' the \code{NFCP_parameters} function specifies the parameters of
#' a commodity pricing model under the N-factor framework first described by Cortazar and Naranjo (2006).
#' This function is a recommended starting position for the application of N-factor models within the \code{NFCP} package.
#'
#' @param N_factors \code{numeric}. Number of state variables in the spot price process.
#' @param GBM \code{logical}. If \code{GBM = T}, factor 1 of the model is assumed to follow a Brownian Motion, inducing a unit-root in the spot price process.
#' @param initial_states \code{logical}. If \code{initial_states = T}, the initial state vector is specified as unknown parameters of the commodity pricing model.
#' @param N_ME \code{numeric}. The number of independent measuring errors of observable futures contracts to consider in the Kalman filter.
#' @param N_season \code{numeric}. The number of deterministic, cyclical seasonal factors to include in the spot price process.
#' @param verbose \code{logical}. If \code{verbose = T}, the stochastic differential equation of the spot price process is printed when the function is called.
#'
#'
#'@details
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
#'model's fit to observed prices, or as errors in the reporting of prices (Schwartz and Smith, 2000). These disturbances are typically assumed independent by the commodity pricing literature.
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
#'@return A vector of parameter names for a specified N-factor spot price process. This vector is ideal for application within many other functions within the \code{NFCP} package
#'
#'@references
#'
#'Hannan, E. J., et al. (1970). "The seasonal adjustment of economic time series." \emph{International economic review} 11(1): 24-52.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Sørensen, C. (2002). "Modeling seasonality in agricultural commodity futures." \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products} 22(5): 393-426.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'##Generate parameter of a Two-factor model Crude Oil model
#'##as first presented by Schwartz and Smith (2000):
#'two_factor_parameters <- NFCP_parameters(N_factors = 2,
#'                                         GBM = TRUE,
#'                                         initial_states = FALSE,
#'                                         N_ME = 5)
#'print(two_factor_parameters)
#'@export
NFCP_parameters <- function(N_factors, GBM, initial_states, N_ME, N_season = 0, verbose = TRUE){

  bool <- !c(is.numeric(N_factors), is.numeric(N_ME))
  if(any(bool)) stop(paste("parameter", paste(c("N_factors", "N_ME")[bool], collapse = ", "), "must be class 'numeric'!"))
  bool <- !c(is.logical(GBM), is.logical(initial_states), is.logical(verbose))
  if(any(bool)) stop(paste("parameter", paste(c("GBM", "initial_states", "verbose")[bool], collapse = ", "), "must be class 'logical'!"))

  if(N_factors < 1) stop("N_factors must be >= 1!")
  if(N_season < 0) stop("The number of seasonal factors must be >=0!")

  ## Hard set catch error:
  if(N_factors > 10) stop("N_factors probably shouldn't be > 10!")

  ##To start, specify the GBM growth factor:
  input_names <- c("mu", sapply(c("lambda", "kappa", "sigma"), FUN = function(X) paste(X, 1:N_factors, sep ="_")))

  ##The correlations:
  for(i in 1:N_factors) for(j in i:N_factors) if(i != j)  input_names <- c(input_names, paste("rho", i, j, sep = "_"))

  if(N_season > 0) for(i in 1:N_season) input_names <- c(input_names, paste("season", i, 1:2, sep = "_"))

  ##Measurement Error (white noise of contracts):
  if(N_ME > 0) input_names <- c(input_names, paste0("ME_", 1:N_ME))
  if(initial_states)  input_names <- c(paste("x_0", 1:N_factors, sep = "_"), input_names)

  #####Model Initialization:
  ###If the first factor is a GBM Process, like in Schwartz and Smith (2000), etc.:
  if(GBM){
    input_names <- input_names[!input_names %in% "kappa_1"]
    input_names[input_names %in% "lambda_1"] = "mu_rn"
  } else {
    #We don't estimate a growth rate, and in the KF function we'll set it to zero:
    input_names[input_names %in% "mu"] = "E"
  }

  if(verbose){
  #####Print the details:
  cat("\n----------------------------------------------------------------\n")
  cat(paste0(N_factors, ifelse(N_factors>1, " Factors", " Factor"), " model:", ifelse(GBM, " first factor is a GBM", " first factor is an MR"),  "\n\n"))
  cat("Risk Neutral SDE: \n\n")

  statement <- ifelse(GBM, "log(s_t) = ", "log(s_t) = E + ")
  if(N_season > 0) statement <- paste0(statement, "season_t")
  cat(paste0(statement, " + sum(x_t)\n\n"))

  cat("Where: \n")
  if(GBM){
    cat("d x1_t    = mu_rn * dt  + sigma_1 * dW_1\n")
  }
  if(!GBM || (GBM && N_factors > 1)){
    for(i in ifelse(GBM, 2, 1): N_factors){
      cat(paste0("d x",i,"_t    = - (lambda_",i," + kappa_", i," * x",i,"_t) * dt  + sigma_",i," * dW_",i,"\n"))
    }
  }
  if(N_factors > 1){
    cat("\n And: \n\n")
    for(i in 1:N_factors) for(j in i:N_factors) if(i!=j) cat(paste0("E(dW_",i," * dW_",j,") = rho_",i,"_",j," * dt\n"))
  }

  if(N_season > 0){
    cat("\nAdditionally, there is a deterministic seasonal component: \n\n")

    season_print <- "season_1_1 * cos(2 * PI * t) + season_1_2 * sin(2 * PI * t)"
    if(N_season > 1) for(i in 2:N_season) season_print <- paste0(season_print, " + season_", i, "_1 * cos(",2*i, " * PI * t) + season_", i, "_2 * sin(", 2*i, "* PI * t)" )
    cat(paste0("SEASON_t = ", season_print, "\n"))
  }
  }


  return(c(input_names, use.names = FALSE))
}


#'Stitch futures contracts
#'@description Aggregate futures contract price data by stitching according to either approximate maturities and rollover frequency or contract number from closest maturity.
#'
#'@param futures Contract futures price data. Each row of \code{Futures} should represent one observation of futures prices and each column should represent one quoted futures contract. NA's in \code{Futures} are allowed, representing missing observations.
#'@param futures_TTM A \code{vector} of contract maturities to stitch
#'@param maturity_matrix The time-to-maturity (in years) for each contract at each given observation point. The dimensions of \code{maturity_matrix} should match those of \code{Futures}
#'@param rollover_frequency the frequency (in years) at which contracts should be rolled over
#'@param contract_numbers A \code{vector} of contract numbers offset from the closest-to-maturity contract at which to stitch contracts.
#'@param verbose \code{logical}. Should additional information be output? see \bold{details}
#'
#'@details
#'This function aggregates a set of futures contract data by stitching contract data over an observation period, resulting in a
#'set of futures observations that is 'complete' (ie. Does not feature missing observations). Aggregated futures
#'data benefit from several computational efficiencies compared to raw contract data, but results in the loss of futures price information.
#'
#'There are two methods of the \code{stitch_contracts} function that can be utilized the stitch contracts:
#'
#'\bold{Method 1}
#'
#'\code{stitch_contracts(futures, contract_numbers, verbose = T)}
#'Futures data may be aggregated by stitching prices according to maturity matching. This method requires the inputs \code{futures_TTM}, \code{maturity_matrix} and \code{rollover_frequency}.
#'This method stitched contracts by matching the observation prices according to which contract has the closest time-to-maturity of the desired maturity specified
#'in \code{futures_TTM}. Contracts are rolled over at the frequency specified in \code{rollover_frequency}.
#'
#'\bold{Method 2}
#'
#'\code{stitch_contracts(futures, futures_TTM, maturity_matrix, rollover_frequency, verbose = T)}
#'Futures data may be stitched according to the contract numbers offset from the closest-to-maturity contract. This method requires only the
#'input \code{contract_numbers} specifying which contracts should be included. This method is most appropriate when the maturity of available
#'contracts are consistent (ie. contracts expire every month or three months).
#'
#'@return
#'\code{stitch_contracts} returns a matrix of stitched futures prices if \code{verbose = T} and a list with two or three objects otherwise (see below).
#'
#'\tabular{ll}{
#'
#' \code{prices} \tab A data frame of Stitched futures prices. Each row represents an observation of the specified contracts. \cr
#'
#' \code{maturities} \tab A data frame of the time-to-maturity of observed futures prices. Each row represents an observation of the
#'specified contracts. Returned only when \bold{Method 1} is used (see \bold{Details})  and \code{verbose = T}. \cr
#'
#' \code{tickers} \tab  A data frame of the named columns of observed futures prices (e.g. contract tickers). Returned only when \code{Futures} or \code{maturity_matrix} have named columns and \code{verbose = T}. \cr
#' }
#'
#'@examples
#'
#'
#'##These examples approximately replicate the Crude Oil data utilized within the
#'##prominent work of Schwartz and Smith (2000):
#'
#'###Method 1 - Stitch crude oil contracts according to maturity matching:
#'SS_stitched_M1 <- stitch_contracts(futures = SS_oil$contracts,
#'                                   futures_TTM = c(1, 5, 9, 13, 17)/12,
#'                                   maturity_matrix = SS_oil$contract_maturities,
#'                                   rollover_frequency = 1/12, verbose = TRUE)
#'
#'###Method 2 - Stitch crude oil contracts according to nearest contract numbers:
#'SS_stitched_M2 <- stitch_contracts(futures = SS_oil$contracts,
#'                                   contract_numbers = c(1, 5, 9, 13, 17), verbose = TRUE)
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@export
stitch_contracts <- function(futures, futures_TTM = NULL, maturity_matrix = NULL, rollover_frequency = NULL, contract_numbers = NULL, verbose = FALSE){

  futures <- as.matrix(futures)

  ###Rolling Over by trying to match maturities:
  if(is.null(contract_numbers)){
    futures_TTM <- as.matrix(futures_TTM)


    stitched_prices <- matrix(NA, nrow = nrow(futures), ncol = length(futures_TTM))
    colnames(stitched_prices) <- futures_TTM
    rownames(stitched_prices) <- rownames(futures)
    if(!is.null(colnames(futures))){
      stitched_tickers <- matrix(NA, nrow = nrow(futures), ncol = length(futures_TTM))
      colnames(stitched_tickers) <- futures_TTM
      rownames(stitched_tickers) <- rownames(futures)
    }

    maturity_matrix <- as.matrix(maturity_matrix)
    ##Get rid of data we don't have Trading.Dates for:
    maturity_matrix[is.na(futures)] <- NA


    stitched_maturities <- matrix(NA, nrow = nrow(futures), ncol = length(futures_TTM))
    colnames(stitched_maturities) <- futures_TTM
    rownames(stitched_maturities) <- rownames(futures)

    if(is.null(maturity_matrix) || is.null(rollover_frequency)) stop("Stitching by maturity matching requires both a maturity matrix and rollover frequency")

    ###We need a function that evaluates the time to maturities and picks the most optimal ones:
    i <- row <- start_row <- stitch_row <- 1
    stitch_fill_row <- 0

    while(stitch_fill_row < nrow(futures)){

      ###Which columns/contracts?
      colnum <-  apply(data.frame(futures_TTM), MARGIN = 1, FUN = function(x) which.min(abs(maturity_matrix[row,] - x)))

      ###Fill out observations until rollover:
      fill_row <- which.min(abs(maturity_matrix[,colnum[1]] - (maturity_matrix[row,colnum[1]] - rollover_frequency)))

      ###The row that you go to in the stitched data:
      stitch_fill_row <- min(nrow(futures), fill_row - (start_row - 1))

      stitch_fill <- stitch_row:stitch_fill_row

      stitched_prices    [stitch_fill,] <- as.matrix(futures[stitch_fill,colnum])
      stitched_maturities[stitch_fill,] <- as.matrix(maturity_matrix[stitch_fill,colnum])
      if(!is.null(colnames(futures))) stitched_tickers[stitch_fill,] <- colnames(futures)[colnum]

      ###We've stitched, now repeat
      ###Where we start from:
      row <- fill_row + 1

      ###Where you start in your stitch rows:
      stitch_row <- stitch_fill_row + 1
      ###We move to the next maturity contract:
    }
    if(anyNA(stitched_prices)) warning("NA's returned in Output Data")

    if(verbose){
      if(!is.null(colnames(futures))) return(list(prices = stitched_prices, maturities = stitched_maturities, tickers = stitched_tickers))
      return(list(prices = stitched_prices, maturities = stitched_maturities))
    } else {
      return(stitched_prices)
    }
  } else {


    stitched_prices <- matrix(NA, nrow = nrow(futures), ncol = length(contract_numbers))
    colnames(stitched_prices) <- futures_TTM
    rownames(stitched_prices) <- rownames(futures)
    if(!is.null(colnames(futures))){
      stitched_tickers <- matrix(NA, nrow = nrow(futures), ncol = length(contract_numbers))
      colnames(stitched_tickers) <- futures_TTM
      rownames(stitched_tickers) <- rownames(futures)
    }

    ###We need a function that evaluates the time to maturities and picks the most optimal ones:
    start_row <- fill_row <- 1
    ###Which column do we start with?
    NA_col <- which(!is.na(futures[1,]))[1]
    colnum <- contract_numbers + NA_col - 1

    for(t in 1:nrow(futures)){
      stitched_prices[t,] <- futures[t, which(!is.na(futures[t,]))[colnum]]
      if(!is.null(colnames(futures))) stitched_tickers[t,] <- colnames(futures)[which(!is.na(futures[t,]))[colnum]]
    }

    colnames(stitched_prices) <- paste0("F", contract_numbers)
    rownames(stitched_prices) <- rownames(futures)
    stitched_prices <- as.data.frame(stitched_prices)

    if(!is.null(colnames(futures))){
      colnames(stitched_tickers) <- paste0("F", contract_numbers)
      rownames(stitched_tickers) <- rownames(futures)
      stitched_tickers <- as.data.frame(stitched_tickers)
    }

    if(anyNA(stitched_prices)) warning("NA's returned in Output Data")

    if(verbose){
      return(list(prices = stitched_prices, tickers = stitched_tickers))
    } else {
      return(stitched_prices)

    }
  }
}

#Function for A(T):
#'Calculate \eqn{A(T)}
#'@description
#'Calculate the values of \eqn{A(T)} for a given N-factor model parameters and observations. Primarily purpose is for application within other functions of the \code{NFCP} package.
#'
#'@keywords internal
#'
#'@param parameters A named vector of parameters of an N-factor model. Function \code{NFCP_parameters} is recommended.
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
#'SS_oil_A_T <- A_T(SS_oil$two_factor, SS_oil$stitched_TTM)
#'
#'@keywords internal
#'@export
A_T <- function(parameters, Tt){

  if(is.null(names(parameters))) stop("argument parameters must be a named vector")

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
  output <- diag(N_factors)

  if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0
  GBM <- parameters["kappa_1"] == 0

  if(GBM) output <- (parameters["mu_rn"] + 0.5 * parameters["sigma_1"]^2) * Tt else output <- 0

  if(!(GBM && N_factors == 1)){
    for(i in ifelse(GBM, 2, 1):N_factors){
      output <- output - (1 - exp(-parameters[paste0("kappa_", i)] * Tt)) * (parameters[paste0("lambda_", i)]/parameters[paste0("kappa_", i)])
    }

    for(i in 1:N_factors){
      for(j in 1:N_factors){
        if(!(i == 1 && i == j)){
          kappa_sums <- sum(parameters[paste0("kappa_", j)], parameters[paste0("kappa_", i)])
          output <- output + 0.5 * parameters[paste0("sigma_",i)] * parameters[paste0("sigma_",j)] * ifelse(i==j, 1, parameters[paste("rho",min(i,j), max(i,j), sep = "_")]) * (1 - exp(-(kappa_sums * Tt)))/kappa_sums
        }
      }
    }}
  return(output)
}

#'N-factor model covariance:
#'
#'@description
#'\loadmathjax
#'Calculate the covariance matrix of state variables for a given N-factor model parameters and discrete time step.
#'
#' @keywords internal
#'
#'@param parameters a named vector of parameters of an N-factor model. Function \code{NFCP_parameters} is recommended.
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
#'SS_oil.covariance <- cov_func(SS_oil$two_factor, SS_oil$dt)
#'
#'@keywords internal
#'@export
cov_func <- function(parameters, dt){

  if(is.null(names(parameters))) stop("argument parameters must be a named vector")

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
  output <- diag(N_factors)

  if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0
  GBM <- parameters["kappa_1"] == 0

  for(i in 1:N_factors){
    for(j in 1:N_factors){
      kappa_sums <- sum(parameters[paste0("kappa_", j)], parameters[paste0("kappa_", i)])
      output[i,j] <- parameters[paste0("sigma_",i)] * parameters[paste0("sigma_",j)] * ifelse(i==j, 1, parameters[paste("rho",min(i,j), max(i,j), sep = "_")]) * (1 - exp(-(kappa_sums * dt)))/kappa_sums
    }
  }
  if(GBM) output[1,1] <- parameters["sigma_1"]^2 * dt

  return(output)
}

#'N-Factor MLE search boundaries
#'@description Generate boundaries for the domain of parameters of the N-factor model for parameter estimation.
#'
#'@param parameters a vector of parameter names of an N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param kappa A vector of length two specifying the lower and upper bounds for the 'kappa' parameter
#'@param lambda A vector of length two specifying the lower and upper bounds for the 'lambda' parameter
#'@param sigma A vector of length two specifying the lower and upper bounds for the 'sigma' parameter
#'@param mu A vector of length two specifying the lower and upper bounds for the 'mu' parameter
#'@param mu_rn A vector of length two specifying the lower and upper bounds for the 'mu_rn' parameter
#'@param rho A vector of length two specifying the lower and upper bounds for the 'rho' parameter
#'@param season A vector of length two specifying the lower and upper bounds for the 'season' parameter
#'@param ME A vector of length two specifying the lower and upper bounds for the 'ME' (i.e., measurement error) parameter
#'@param x_0 A vector of length two specifying the lower and upper bounds for the 'x_0' parameter
#'@param E A vector of length two specifying the lower and upper bounds for the 'E' parameter
#'
#'@details
#'The \code{NFCP_domains} function generates lower and upper bounds for the parameter estimation procedure in the format required of the 'Domains' argument of the 'genoud' function. \code{NFCP_domains}
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
#'##with fixed measurement error:
#'parameters_2F <- NFCP_parameters(N_factors = 2,
#'                                         GBM = TRUE,
#'                                         initial_states = TRUE,
#'                                         N_ME = 1)
#'
#'###Generate the default 'domains' argument of 'NFCP_MLE' function:
#'NFCP_MLE_bounds <- NFCP_domains(parameters_2F)
#'
#'@export
NFCP_domains <- function(parameters,
                        kappa = NULL,
                        lambda = NULL,
                        sigma = NULL,
                        mu = NULL,
                        mu_rn = NULL,
                        rho = NULL,
                        season = NULL,
                        ME = NULL,
                        x_0 = NULL,
                        E = NULL){

  if(is.null(kappa))    kappa    <- c(1e-5, 50)
  if(is.null(lambda))   lambda   <- c(-10, 10)
  if(is.null(sigma))    sigma    <- c(0, 10)
  if(is.null(mu))       mu       <- c(-10, 10)
  if(is.null(mu_rn))    mu_rn    <- c(-10, 10)
  if(is.null(rho))      rho      <- c(-1, 1)
  if(is.null(ME))       ME       <- c(1e-10, 1)
  if(is.null(x_0))      x_0      <- c(-10, 10)
  if(is.null(E))        E        <- c(-10, 10)
  if(is.null(season))       season       <- c(-1, 1)



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

  ## long-term risk-free growth rate
  lower_bounds[grepl("mu_rn", parameters)] <- mu_rn[1]
  upper_bounds[grepl("mu_rn", parameters)] <- mu_rn[2]

  ## kappa
  lower_bounds[grepl("kappa", parameters)] <- kappa[1]
  upper_bounds[grepl("kappa", parameters)] <- kappa[2]

  ## sigma
  lower_bounds[grepl("sigma", parameters)] <- sigma[1]
  upper_bounds[grepl("sigma", parameters)] <- sigma[2]

  ## Seasonality Factors:
  lower_bounds[grepl("season", parameters)] <- season[1]
  upper_bounds[grepl("season", parameters)] <- season[2]

  ###Setting White noise bounds that are too low will result in a singular matrix (1e-5^2 == 1e-10):

  lower_bounds[grepl("ME", parameters)] <- ME[1]
  upper_bounds[grepl("ME", parameters)] <- ME[2]
  ##Correlation between -1 and 1:
  lower_bounds[grepl("rho", parameters)] <- rho[1]
  upper_bounds[grepl("rho", parameters)] <- rho[2]

  ###Initial Inputs:
  lower_bounds[grepl("x_0_", parameters)] <- x_0[1]
  upper_bounds[grepl("x_0_", parameters)] <- x_0[2]
  bounds <- cbind(lower_bounds, upper_bounds)
  colnames(bounds) <- c("Lower Bound", "Upper Bound")
  rownames(bounds) <- parameters

  return(bounds)
}
