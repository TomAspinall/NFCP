
#'Calculate the volatility term structure of futures returns
#'
#'@description Estimate the theoretical and empirical volatility term structure of futures returns
#'
#'@param parameters \code{vector}. A named vector of parameters of an N-factor model. Function \code{NFCP_parameters} is recommended.
#'@param futures \code{matrix}. Historical observes futures price data. Each column must correspond to a listed futures contract and each row must correspond to a discrete observation of futures contracts. NA's are permitted.
#'@param futures_TTM \code{vector}. Each element of 'futures_TTM' must correspond to the time-to-maturity from the current observation point of futures contracts listed in object 'futures'.
#'@param dt \code{numeric}. Constant, discrete time step of observations, in years.
#'
#'@details
#'\loadmathjax
#'
#'The fit of an N-factor models theoretical volatility term structure of futures returns to those obtained directly from observed futures prices can be used as a measure of robustness for
#'the models ability to explain the behaviour of a commodities term structure.
#'
#'The theoretical model volatility term structure of futures returns is given by the following equation:
#'
#'\mjdeqn{\sigma_F(\tau) = \sum_{i=1}^N \sum_{j=1}^N \sigma_i \sigma_j \rho_{i,j} e^{-(\kappa_i + \kappa_j)\tau}}{sigma_F(tau) = sum_{i = 1, j = 1}^N sigma[i] sigma[j] rho[i,j] e^(-(kappa[i] + kappa[j]) tau)}
#'
#'Under the case that \mjeqn{\kappa_1 = 0}{kappa[1] = 0}, the model volatility term structure converges to \mjeqn{\sigma_1^2}{sigma[1]^2} as \mjeqn{\tau}{tau} grows large.
#'
#'The empirical volatility term structure of futures returns is given by:
#'
#'\mjdeqn{\hat\sigma_F^2(\tau) = \frac{1}{\Delta t}\sum_{i=1}^N(log(F(t_i,\tau)/F(t_i-\Delta t,\tau)) - \bar\mu)^2}{hat(sigma)[F^2](tau) = 1/(Delta * t) sum_{i=1}^N (log(F(t[i],tau) / F(t[i] - Delta t, tau)) - bar(mu))^2}
#'
#'According to Cortazar and Naranjo (2006): "A larger number of factors gives more flexibility to adjust first and second moments simultaneously, hence explaining why (a) four-factor (may) outperform (a) three-factor one in fitting the volatility term structure."
#'
#'@return
#'\code{TSfit_volatility} returns a matrix with the theoretical and empirical volatility term structure of futures returns, with the number of columns of this matrix coinciding with the number of input futures contracts.
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'# Test the volatility term structure fit of the Schwartz-Smith two-factor model on crude oil:
#'V_TSFit <- TSfit_volatility(
#'parameters = SS_oil$two_factor,
#'futures = SS_oil$stitched_futures,
#'futures_TTM = SS_oil$stitched_TTM,
#'dt = SS_oil$dt)
#'
#'@export
TSfit_volatility <- function(parameters, futures, futures_TTM, dt){

  futures_TTM <- c(futures_TTM)
  if(is.null(parameters)) stop("'parameters' must be a named vector!")
  if(length(futures_TTM) != ncol(futures)) stop("The length of 'futures_TTM' does not equal the number of columns of object 'futures'!")
  if(any(futures_TTM < 0, na.rm = TRUE)) stop("'futures_TTM' cannot have elements less than zero!")
  if(dt <= 0) stop("'dt' must be greater than zero!")

  N_factors <- max(which(paste0("sigma_", 1:length(parameters)) %in% names(parameters) & sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.numeric) & !sapply(parameters[paste0("sigma_",1:length(parameters))], FUN = is.na)))

  if(!"kappa_1" %in% names(parameters)) parameters["kappa_1"] <- 0

  ###Theoretical Volatility Term Structure:
  VTS_theoretical <- rep(0, length(futures_TTM))
  for(F_T in 1:length(futures_TTM)){
    if(!is.na(futures_TTM[F_T])){
      for(i in 1:N_factors){
        for(j in 1:N_factors){
          VTS_theoretical[F_T]  <- VTS_theoretical[F_T] +
            parameters[paste0("sigma_", i)] * parameters[paste0("sigma_", j)] * ifelse(i == j, 1, parameters[paste("rho", min(i,j), max(i,j), sep = "_")]) * exp(- (parameters[paste0("kappa_", i)] + parameters[paste0("kappa_", j)]) * futures_TTM[F_T])
        }}}}

  N.obs <- nrow(futures)

  ###Empirical Volatility Term Structure:
  VTS_empirical <- rep(0, length(futures_TTM))
  for(F_T in 1:length(futures_TTM)){
    if(!is.na(futures_TTM[F_T])){

    dates <- which(!is.na(futures[,F_T]))
    # dates = dates[min(N.obs, length(dates)):length(dates)]

     returns <- log(futures[dates[-length(dates)],F_T]/futures[dates[-1],F_T])
     mean.returns <- mean(returns, na.rm = T)
     VTS_empirical[F_T] <- sum((returns - mean.returns)^2) / (length(dates) * dt)

    }}
  volatility_term_structure <- rbind(VTS_theoretical, VTS_empirical)
  volatility_term_structure <- sqrt(volatility_term_structure)
  volatility_term_structure <- rbind(futures_TTM, volatility_term_structure)
  rownames(volatility_term_structure) <- c("Maturity", "Theoretical Volatility", "Empirical Volatility")
  colnames(volatility_term_structure) <- colnames(futures)
  volatility_term_structure <- volatility_term_structure[,!is.na(volatility_term_structure["Maturity",])]

  return(volatility_term_structure)
}



