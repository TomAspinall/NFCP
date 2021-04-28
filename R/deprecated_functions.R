
#'Volatility Term Structure of futures returns
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'V_TSFit <- TSFit.Volatility(
#'  parameters = SS.Oil$Two.Factor,
#'  Futures = SS.Oil$Stitched.Futures,
#'  TTM = SS.Oil$Stitched.TTM,
#'  dt = SS.Oil$dt)
#'
#'## ->
#'
#'V_TSFit <- TSfit_volatility(
#'  parameters = SS_oil$two_factor,
#'  futures = SS_oil$stitched_futures,
#'  futures_TTM = SS_oil$stitched_TTM,
#'  dt = SS_oil$dt)
#'
#'@keywords internal
#'@export
TSFit.Volatility <- function(parameters, Futures, TTM, dt){

  #warning deprecation:
  .Deprecated(msg = "'TSFit.Volatility()' was deprecated in NFCP 0.3.0. \n Please use 'TSfit_volatility()' instead. \n In addition, the following arguments have been renamed: \n
'Futures' -> 'futures'
'TTM' -> 'futures_TTM'")

  # TSfit_volatility(parameters, Futures, TTM, dt)

  TTM <- c(TTM)
  if(is.null(parameters)) stop("'parameters' must be a named vector!")

  N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.numeric) & !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))

  if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0

  ###Theoretical Volatility Term Structure:
  VTS.theoretical <- rep(0, length(TTM))
  for(F_T in 1:length(TTM)){
    if(!is.na(TTM[F_T])){
      for(i in 1:N.factors){
        for(j in 1:N.factors){
          VTS.theoretical[F_T]  <- VTS.theoretical[F_T] +
            parameters[paste0("sigma_", i)] * parameters[paste0("sigma_", j)] * ifelse(i == j, 1, parameters[paste("rho", min(i,j), max(i,j), sep = "_")]) * exp(- (parameters[paste0("kappa_", i)] + parameters[paste0("kappa_", j)]) * TTM[F_T])
        }}}}

  N.obs <- nrow(Futures)

  ###Empirical Volatility Term Structure:
  VTS.empirical <- rep(0, length(TTM))
  for(F_T in 1:length(TTM)){
    if(!is.na(TTM[F_T])){

      dates <- which(!is.na(Futures[,F_T]))
      # dates = dates[min(N.obs, length(dates)):length(dates)]

      returns <- log(Futures[dates[-length(dates)],F_T]/Futures[dates[-1],F_T])
      mean.returns <- mean(returns, na.rm = T)
      VTS.empirical[F_T] <- sum((returns - mean.returns)^2) / (length(dates) * dt)

    }}
  Volatility_term_structure <- rbind(VTS.theoretical, VTS.empirical)
  Volatility_term_structure <- sqrt(Volatility_term_structure)
  Volatility_term_structure <- rbind(TTM, Volatility_term_structure)
  rownames(Volatility_term_structure) <- c("Maturity", "Theoretical Volatility", "Empirical Volatility")
  colnames(Volatility_term_structure) <- colnames(Futures)
  Volatility_term_structure <- Volatility_term_structure[,!is.na(Volatility_term_structure["Maturity",])]

  return(Volatility_term_structure)

}

#'N-Factor Model American Put Option Pricing
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#' @examples
#'
#'American.Option.Value(X.0 = c(3,0),
#'                      parameters = SS_oil$two_factor,
#'                      n = 1e2,
#'                      t = 1,
#'                      dt = 1/12,
#'                      K = 20,
#'                      r = 0.05,
#'                      verbose = FALSE,
#'                      orthogonal = "Power",
#'                      degree = 2)
#'
#'## ->
#'
#'American_option_value(x_0 = c(3,0),
#'                      parameters = SS_oil$two_factor,
#'                      N_simulations = 1e2,
#'                      t = 1,
#'                      dt = 1/12,
#'                      K = 20,
#'                      r = 0.05,
#'                      verbose = FALSE,
#'                      orthogonal = "Power",
#'                      degree = 2)
#'@keywords internal
#' @export
American.Option.Value <- function(X.0, parameters, n, t, dt, K, r, orthogonal = "Power", degree = 2, verbose = FALSE, debugging = FALSE){

  #warning deprecation:
  .Deprecated(msg = "'American.Option.Value()' was deprecated in NFCP 0.3.0. \n Please use 'American_option_value()' instead. \n In addition, the following arguments have been renamed: \n
'X.0' -> 'x_0'")

  American_option_value(x_0 = X.0, parameters = parameters, N_simulations = n, t = t, dt = dt, K = K, r = r, orthogonal = orthogonal, degree = degree, verbose = verbose, debugging = debugging)


}


#'N-Factor Model European Option Pricing
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#' @examples
#'
#'##Step 2 - Calculate 'call' option price:
#' European.Option.Value(X.0 = c(3,0),
#'                       parameters = SS.Oil$Two.Factor,
#'                       t = 1,
#'                       TTM = 1,
#'                       K = 20,
#'                       r = 0.05,
#'                       call = TRUE,
#'                       verbose = FALSE)
#'
#'## ->
#'
#'European_option_value(x_0 = c(3,0),
#'                      parameters = SS_oil$two_factor,
#'                      t = 1,
#'                      futures_maturity = 1,
#'                      option_maturity = 1,
#'                      K = 20,
#'                      r = 0.05,
#'                      call = TRUE,
#'                      verbose = FALSE)
#'
#'@keywords internal
#'@export
European.Option.Value <- function(X.0, parameters, t, TTM, K, r, call, verbose = FALSE){

  #warning deprecation:
  .Deprecated(msg = "'European.Option.Value()' was deprecated in NFCP 0.3.0. \n Please use 'European_option_value()' instead. \n In addition, the following arguments have been renamed: \n
'X.0' -> 'x_0'")

  # European_option_value(X.0, parameters, t, TTM, TTM, K, r, call, verbose)
  X.0 <- c(X.0, use.names = F)
  names(X.0) <- paste0("X.", 1:length(X.0))

  ###Inputs:
  parameters <- c(t = t, TTM = TTM, K = K, r = r, X.0, parameters)

  ##Remove unnecessary parameters:
  parameters <- parameters[!names(parameters) %in% c("mu")]
  parameters <- parameters[!grepl("contract", names(parameters))]

  ###Instantiate the actual function:
  European.Option.Calc <- function(X = 0, parameters, gradit = 0, call){

    if(gradit>0 && gradit<=length(parameters)) parameters[gradit] <- X

    if(parameters["TTM"] < parameters["t"]) stop("Contract Maturity must be greater than Option Expiration")

    ###Calculate the Variance:
    N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
    GBM <- !("kappa_1" %in% names(parameters))
    if(GBM) parameters["kappa_1"] <- 0

    ###Drift term:
    drift <- parameters["mu_star"] * parameters["t"] + NFCP::A_T(parameters, parameters["TTM"] - parameters["t"])
    for(i in 1:N.factors) drift <- drift + exp(-parameters[paste0("kappa_", i)] * parameters["TTM"]) * parameters[paste0("X.",i)]

    ###Volatility term:
    sigma_tT <- ifelse(GBM, parameters["sigma_1"]^2 * parameters["t"], 0)
    for(i in 1:N.factors){
      for(j in 1:N.factors){
        if(!(GBM && i == 1 && j ==1)){
          kappa_sums <- sum(parameters[paste0("kappa_", i)], parameters[paste0("kappa_", j)])

          sigma_tT <- sigma_tT + exp(-kappa_sums * (parameters["TTM"] - parameters["t"])) * parameters[paste0("sigma_",i)] * parameters[paste0("sigma_",j)] * ifelse(i==j, 1,
                                                                                                                                                                     parameters[paste("rho",min(i,j), max(i,j), sep = "_")]) * (1 - exp(-(kappa_sums * parameters["t"])))/kappa_sums
        }
      }}
    ###SD:
    sigma_tT <- c(sqrt(sigma_tT), use.names = FALSE)

    ###Expected Futures price:
    F_tT <- exp(drift + 0.5 * sigma_tT^2)
    # F_tT <- NFCP::Futures.Price.Forecast(X.0 = parameters[grepl("X.", names(parameters))],
    #                                      parameters = parameters, t = parameters["t"], TTM = parameters["TTM"])

    parameters["F_tT"] <- F_tT
    parameters["sigma_tT"] <- sigma_tT

    if(length(parameters) - gradit < 3) parameters[gradit] <- X

    d.first <-  log(parameters["F_tT"] / parameters["K"])/ (parameters["sigma_tT"])
    d.second <- 0.5 * parameters["sigma_tT"]

    d1 <-   d.first + d.second
    d2 <-   d.first - d.second
    d3 <- - d.first + d.second
    d4 <- - d.first - d.second

    if(call){
      value <- as.numeric(exp(-parameters["r"]*parameters["t"]) * (parameters["F_tT"] * stats::pnorm(d1) - parameters["K"] * stats::pnorm(d2)))
    } else {
      value <- as.numeric(exp(-parameters["r"]*parameters["t"]) * (parameters["K"] * stats::pnorm(d3) - parameters["F_tT"] * stats::pnorm(d4)))
    }

    return(value)
  }
  Value <- European.Option.Calc(parameters = parameters,call = call)
  if(!verbose){
    return(Value)
  } else {

    ##Calculate Initial Values:
    ###Calculate expected Futures price:
    F_tT <- NFCP::Futures.Price.Forecast(X.0, parameters, t, TTM)

    ###Calculate the Variance:
    N.factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters)))
    if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0
    GBM <- parameters["kappa_1"] == 0
    sigma_tT <- ifelse(GBM, parameters["sigma_1"]^2 * t, 0)
    for(i in 1:N.factors){
      for(j in 1:N.factors){
        if(!(GBM && i == 1 && j ==1)){
          kappa_sums <- sum(parameters[paste0("kappa_", i)], parameters[paste0("kappa_", j)])

          sigma_tT <- sigma_tT + exp(-kappa_sums * (TTM - t)) * parameters[paste0("sigma_",i)] * parameters[paste0("sigma_",j)] * ifelse(i==j, 1, parameters[paste("rho",min(i,j), max(i,j), sep = "_")]) * ((1 - exp(-(kappa_sums * t)))/kappa_sums)
        }}}
    ###SD:
    sigma_tT <- c(sqrt(sigma_tT), use.names = FALSE)
    if(GBM) parameters <- parameters[!names(parameters) %in% "kappa_1"]

    parameters <- c(parameters, F_tT = F_tT, sigma_tT = sigma_tT)

    ###Calculate Gradients:
    TheGreeks <- rep(0, length(parameters))
    names(TheGreeks) <- names(parameters)

    for(gradit in 1:length(parameters)){
      if(TTM - t < 0.02){
        if(names(parameters)[gradit] %in% c("t", "TTM")){
          if(names(parameters)[gradit] %in% "t")   TheGreeks[gradit] <- numDeriv::grad(func = European.Option.Calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call, side = -1)
          if(names(parameters)[gradit] %in% "TTM") TheGreeks[gradit] <- numDeriv::grad(func = European.Option.Calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call, side = 1)
        }
        else {
          TheGreeks[gradit] <- numDeriv::grad(func = European.Option.Calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call)
        }
      }
      else {
        TheGreeks[gradit] <- numDeriv::grad(func = European.Option.Calc, parameters = parameters, x = parameters[gradit], gradit = gradit, call = call)
      }
    }
    names(TheGreeks)[(length(TheGreeks)-1):length(TheGreeks)] <- c("Underlying Volatility", "Underlying Futures Price")


    final_output <- list(
      Value = Value,
      Annualized.Volatility = sigma_tT/sqrt(t),
      Parameter.Sensitivity = TheGreeks[!names(TheGreeks) %in% c("Underlying Volatility", "Underlying Futures Price", "t", "TTM", "K", "r")],
      Greeks = TheGreeks[c("Underlying Volatility", "Underlying Futures Price", "t", "TTM", "K", "r")]
    )
    return(final_output)
  }
}



#' N-Factor Model Parameter Estimation through the Kalman Filter and Maximum Likelihood Estimation
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'SS_2F <- NFCP.MLE(
#'  ####Arguments
#'  log.futures = log(SS_oil$stitched_futures)[1:5,1],
#'  dt = SS_oil$dt,
#'  TTM= SS_oil$stitched_TTM[1],
#'  S.Constant = FALSE, N.factors = 1, GBM = TRUE,
#'  ####Genoud arguments:
#'  hessian = TRUE,
#'  Richardsons.Extrapolation = FALSE,
#'  pop.size = 4, optim.method = "L-BFGS-B", print.level = 0,
#'  max.generations = 0, solution.tolerance = 10)
#'
#'## ->
#'NFCP_MLE(
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
#'
#'@keywords internal
#'@export
NFCP.MLE <- function(log.futures, dt, TTM, N.factors, GBM = TRUE, S.Constant = TRUE, Estimate.Initial.State = FALSE,
                     Richardsons.Extrapolation = TRUE, cluster = FALSE, Domains = NULL, ...){

  #warning deprecation:
  .Deprecated(msg = "'NFCP.MLE()' was deprecated in NFCP 0.3.0. \n Please use 'NFCP_MLE()' instead.\n In addition, the following arguments have been renamed:\n
'log.futures' -> 'log_futures'
'TTM' -> 'futures_TTM'
'N.factors' -> 'N_factors'
'S.Constant' -> 'fixed_ME'
'Estimate.Initial.State' -> 'estimate_initial_states'
'Richardsons.Extrapolation' -> 'Richardsons_extrapolation' \n
The following outputs have also been renamed: \n
'Estimated.Parameters' -> 'estimated_parameters'
'Standard.Errors' -> 'standard_errors'
'genoud.value' -> 'genoud_value'
")

# output <-  NFCP_MLE(log.futures, dt, TTM, N.factors, GBM, S.Constant, Estimate.Initial.State,
#            Richardsons.Extrapolation, cluster, Domains, ...)
# names(output) <- c("MLE", "Estimated.Parameters", "Standard.Errors", "X.t", "X", "Y", "V", "TSFit.Error", "TSFit.Volatility", "proc_time", "genoud.value")
# output

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


#'Forecast N-Factor Model Spot Prices
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'Spot.Price.Forecast(X.0 = c(3,0),
#'                   parameters = SS_oil$two_factor,
#'                   t = seq(0,9,1/12),
#'                   Percentiles = c(0.1, 0.9))
#'
#' ## ->
#'
#'spot_price_forecast(x_0 = c(3,0),
#'                    parameters = SS_oil$two_factor,
#'                    t = seq(0,9,1/12),
#'                    percentiles = c(0.1, 0.9))
#'
#'@keywords internal
#'@export
Spot.Price.Forecast <- function(X.0, parameters, t, Percentiles = NULL){

  #warning deprecation:
  .Deprecated(msg = "'Spot.Price.Forecast()' was deprecated in NFCP 0.3.0. \n Please use 'spot_price_forecast()' instead.\n In addition, the following arguments have been renamed: \n
'X.0' -> 'x_0'")

  spot_price_forecast(X.0, parameters, t, Percentiles)

}

#'Forecast N-Factor Model Futures Prices
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#' @examples
#'
#'Futures.Price.Forecast(X.0 = c(3,0),
#'                      parameters = SS.Oil$Two.Factor,
#'                      t = 0,
#'                      TTM = seq(0,9,1/12),
#'                      Percentiles = c(0.1, 0.9))
#'
#'## ->
#'
#'futures_price_forecast(x_0 = c(3,0),
#'                       parameters = SS_oil$two_factor,
#'                       t = 0,
#'                       futures_TTM = seq(0,9,1/12),
#'                       percentiles = c(0.1, 0.9))
#'
#'@keywords internal
#'@export
Futures.Price.Forecast <- function(X.0, parameters, t = 0, TTM = 1:10, Percentiles = NULL){

  #warning deprecation:
  .Deprecated(msg = "'Futures.Price.Forecast()' was deprecated in NFCP 0.3.0. \n Please use 'futures_price_forecast()' instead.\n In addition, the following arguments have been renamed: \n
'X.0' -> 'x_0'
'TTM' -> 'futures_TTM'")

  futures_price_forecast(X.0, parameters, t, TTM, Percentiles)

}

#'Simulate N-Factor Model Spot Prices
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'Spot.Price.Simulate(
#'X.0 = c(3,0),
#'parameters = SS.Oil$Two.Factor,
#'t = 1,
#'dt = 1/12,
#'n = 1e3,
#'antithetic = TRUE,
#'verbose = TRUE)
#'
#'## ->
#'
#'spot_price_simulate(
#'  x_0 = c(3,0),
#'  parameters = SS_oil$two_factor,
#'  t = 1,
#'  dt = 1/12,
#'  N_simulations = 1e3,
#'  antithetic = TRUE,
#'  verbose = TRUE)
#'
#'
#'@keywords internal
#'@export
Spot.Price.Simulate <- function(X.0, parameters, t = 1, dt = 1, n = 2, antithetic = TRUE, verbose = FALSE){


#warning deprecation:
.Deprecated(msg = "'Spot.Price.Simulate()' was deprecated in NFCP 0.3.0. \n Please use 'spot_price_simulate()' instead.\n In addition, the following arguments have been renamed: \n
'X.0' -> 'x_0'
'n' -> 'N_simulations' \n
The following outputs have also been renamed: \n
'State_variables' -> 'state_variables'
'Prices' -> 'prices'")

theoutput <-  spot_price_simulate(X.0, parameters, t, dt, n, antithetic, verbose)

if(!verbose) return(theoutput)

names(theoutput) <- c("State_variables", "Prices")
return(theoutput)

}

#'Simulate N-Factor Model Futures Prices
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'Futures.Price.Simulate(X.0 = c(log(SS.Oil$Spot[1,1]), 0),
#'                      parameters = SS.Oil.Two.Factor,
#'                      dt = SS.Oil$dt,
#'                      N.obs = nrow(SS.Oil$Contracts),
#'                      TTM = SS.Oil$Contract.Maturities)
#'
#'## ->
#'
#'futures_price_simulate(x_0 = c(log(SS_oil$spot[1,1]), 0),
#'                      parameters = SS_oil$two_factor,
#'                      dt = SS_oil$dt,
#'                      n_obs = nrow(SS_oil$stitched_futures),
#'                      futures_TTM = SS_oil$stitched_TTM)
#'@keywords internal
#'@export
Futures.Price.Simulate <- function(X.0, parameters, dt, N.obs, TTM, verbose = TRUE){

#warning deprecation:
.Deprecated(msg = "'Futures.Price.Simulate()' was deprecated in NFCP 0.3.0. \n Please use 'futures_price_simulate()' instead.\n In addition, the following arguments have been renamed: \n
'X.0' -> 'x_0'
'N.obs' -> 'n_obs'
'TTM' -> 'futures_TTM' \n
The following outputs have also been renamed: \n
'State_Vector' -> 'state_vector'
'Futures' -> 'futures'
'Spot' -> 'spot'")

  # output <- futures_price_simulate(X.0, parameters, dt, N.obs, TTM, verbose)
  # if(!verbose) return(output)
  # names(output) <- c("State_Vector", "Futures", "Spot")
  # return(output)

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


#' Specify parameters of N-factor model
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'NFCP.Parameters(N.factors = 2,
#'                GBM = TRUE,
#'                Initial.State = FALSE,
#'                S.Constant = FALSE,
#'                N.contracts = 5)
#'
#'## ->
#'
#'NFCP_parameters(N_factors = 2,
#'                GBM = TRUE,
#'                initial_states = FALSE,
#'                fixed_ME = FALSE,
#'                N_contracts = 5)
#'@keywords internal
#'@export
NFCP.Parameters <- function(N.factors, GBM, Initial.State, S.Constant, N.contracts = NULL, verbose = TRUE){

#warning deprecation:
.Deprecated(msg = "'NFCP.Parameters()' was deprecated in NFCP 0.3.0. \n Please use 'NFCP_parameters()' instead.\n In addition, the following arguments have been renamed: \n
'N.factors' -> 'N_factors'
'Initial.State' -> 'initial_states'  \n
'S.Constant' -> 'fixed_ME'
'N.contracts' -> 'N_contracts'
The following outputs have also been renamed: \n
'sigma.contract' -> 'ME'")


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
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'Stitch.Contracts(Futures = SS.Oil$Contracts,
#'TTM = c(1, 5, 9, 13, 17)/12, maturity.matrix = SS.Oil$Contract.Maturities,
#'rollover.frequency = 1/12, verbose = TRUE)
#'
#'## ->
#'
#'stitch_contracts(futures = SS_oil$contracts,
#'                futures_TTM = c(1, 5, 9, 13, 17)/12,
#'                maturity_matrix = SS_oil$contract_maturities,
#'                rollover_frequency = 1/12, verbose = TRUE)
#'
#'@keywords internal
#'@export
Stitch.Contracts <- function(Futures, TTM = NULL, maturity.matrix = NULL, rollover.frequency = NULL, Contract.Numbers = NULL, verbose = FALSE){

  #warning deprecation:
  .Deprecated(msg = "'Stitch.Contracts()' was deprecated in NFCP 0.3.0. \n Please use 'stitch_contracts()' instead.\n In addition, the following arguments have been renamed: \n
              'Futures' -> 'futures'
              'TTM' -> 'futures_TTM'")

#warning deprecation:
.Deprecated(msg = "'Stitch.Contracts()' was deprecated in NFCP 0.3.0. \n Please use 'stitch_contracts()' instead.\n In addition, the following arguments have been renamed: \n
'Futures' -> 'futures'
'maturity.matrix' -> 'maturity_matrix'  \n
'rollover.frequency' -> 'rollover_frequency'
'Contract.Numbers' -> 'contract_numbers'")

stitch_contracts(Futures, TTM, maturity.matrix, rollover.frequency, Contract.Numbers, verbose)

}



#'N-Factor MLE Search Boundaries
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'NFCP.Domains(names(SS_oil$two_factor))
#'
#'## ->
#'
#'NFCP_domains(names(SS_oil$two_factor))
#'
#'@keywords internal
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

  #warning deprecation:
  .Deprecated(msg = "'NFCP.Domains()' was deprecated in NFCP 0.3.0. \n Please use 'NFCP_domains()' instead.\n In addition, the following arguments have been renamed: \n
              'contract' -> 'ME'")
# NFCP_domains(parameters, kappa, lambda, sigma, mu, mu_star,rho, contract, X.0, E)
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


#'Crude Oil Term Structure Futures Data (1990 - 1995)
#'
#' @docType data
#'
#' @usage data(SS.Oil)
#'
#'
#' @keywords datasets internal
#'
#'@description \code{SS.Oil} has been deprecated as of version 0.3.0. Please use the \code{SS_oil} data object instead
#'
#'@format A \code{list} Containing eight objects:
#'
#'\describe{
#'  \item{Contracts}{A data frame with 268 rows and 82 columns. Each column represents a Crude Oil futures contract, and each row represents a closing
#'weekly price for that futures contract. Observation dates of the contract object are weekly in frequency from
#'\code{1990-02-06} to \code{1995-02-14}. Contracts without observations on a particular date are represented as \code{NA}.}
#'  \item{Stitched.Futures}{Schwartz and Smith (2000) applied stitched contract observation data to estimate commodity pricing models, which
#'are approximated within this object. The \code{Stitched.Futures} object was developed using the \code{Stitch.Contracts} function (see \code{Stitch.Contracts} examples for more details). Contracts were
#'stitched according to the contract numbers
#'specified within the object \code{Stitched.TTM}. \code{Stitched.Futures} is
#'identical to the futures data made available within the MATLAB program "SchwartzSmithModel" developed by Goodwin (2013).}
#'  \item{Spot}{A \code{data.frame} of spot prices of Crude Oil. weekly in frequency from
#'\code{1990-02-06} to \code{1995-02-14}.}
#'  \item{Final.Trading.Days}{Named vector listing the final trading days of each observed futures contract within the \code{Contracts} object. Each element of
#'\code{Final.Trading.Days} corresponds to a column of the \code{Contracts} object. The final
#'trading day of a futures contract is used to calculate the number of business
#'days from a given observation to the maturity of the contract (ie. a contract time to maturity).}
#'  \item{Contract.Maturities}{A data frame with identical dimensions to the \code{Contracts} data frame. This data frame
#'lists the time to maturity of a given futures contract in years at each observation point.
#'This is identical to the number of business days (in years) between the observed date and the final trading day of a particular futures contract.
#'The maturity matrix assumes 262 trading days a year. If the contract is not yet available or has expired, the \code{Contract.Maturities} element is \code{NA}.}
#'  \item{Stitched.TTM}{A vector corresponding to the constant time to maturities that was
#'assumed within the original study of Schwartz and Smith (2000).}
#'  \item{dt}{The discrete time step used to estimate parameters with this data. The time step is 5/262, which represents a
#'weekly frequency of observations where each weekday is a business day (ie. there are no business days on weekends).}
#'  \item{Two.Factor}{The crude oil two-factor commodity pricing model parameters presented within the work of Schwartz and Smith (2000).
#'These parameter estimates are prolific, benchmarked within several subsequent publications.}
#' }
#'
#'@references
#'
#'Dominice Goodwin (2013). Schwartz-Smith 2-factor model - Parameter estimation (https://www.mathworks.com/matlabcentral/fileexchange/43352-schwartz-smith-2-factor-model-parameter-estimation),
#'MATLAB Central File Exchange. Retrieved November 21, 2020.
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
"SS.Oil"



#'N-Factor Commodity Pricing Kalman Filter
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function was deprecated due to a change in the name of the function to adhere to the tidyverse style guide.
#'
#'@examples
#'
#'NFCP.Kalman.filter(
#'  parameter.values = SS_oil$two_factor,
#'  parameters = names(SS_oil$two_factor),
#'  log.futures = log(SS_oil$stitched_futures),
#'  TTM = SS_oil$stitched_TTM,
#'  dt = SS_oil$dt,
#'  verbose = FALSE)
#'
#'## ->
#'
#'NFCP_Kalman_filter(
#'parameter_values = SS_oil$two_factor,
#'parameter_names = names(SS_oil$two_factor),
#'log_futures = log(SS_oil$stitched_futures),
#'futures_TTM = SS_oil$stitched_TTM,
#'dt = SS_oil$dt,
#'verbose = FALSE)
#'
#'@keywords internal
#'@export
NFCP.Kalman.filter = function(parameter.values, parameters, log.futures, dt, TTM, verbose = FALSE, debugging = FALSE){

  #warning deprecation:
  .Deprecated(msg = "'NFCP.Kalman.filter()' was deprecated in NFCP 0.3.0. \n Please use 'NFCP_Kalman_filter()' instead.\n In addition, the following arguments have been renamed:\n
'parameter.values' -> 'parameter_values'
'parameters' -> 'parameter_names'
'TTM' -> 'futures_TTM'
'log.futures' -> 'log_futures' \n
The following outputs have also been renamed: \n
'X.t' -> 'x_t'
'Filtered.Error' -> 'filtered_error'
'TSFit.Error' -> 'TSfit_error'
'TSFit.Volatility' -> 'TSfit_volatility'
'd' -> 'd_t'
'Z' -> 'Z_t'
'H' -> 'H_t'
")

  # output <- NFCP_Kalman_filter(parameter.values, parameters, log.futures, dt, TTM, verbose, debugging)
  #
  # if(verbose || debugging){
  #   if(debugging){
  #     names(output) <- c("LL", "X.t", "X", "Y", "V", "Filtered.Error", "TSFit.Error", "TSFit.Volatility", "LL_t", "P_t", "F_t", "K_t", "d", "Z", "G_t", "c_t", "Q_t", "H")
  #   } else {
  #     names(output) <- c("LL", "X.t", "X", "Y", "V", "Filtered.Error", "TSFit.Error", "TSFit.Volatility")
  #   }
  #
  # }
  #
  # return(output)

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

}






