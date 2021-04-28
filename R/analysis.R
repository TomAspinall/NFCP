
#'Calculate the volatility term structure of futures returns
#'
#'@description Estimate the Theoretical and Empirical Volatility Term Structure of futures returns
#'
#'@param parameters A named vector of parameters of an N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param futures A Matrix of futures price data. Each column corresponds to a given futures contract, and each row is an observation of the futures contracts.
#'@param futures_TTM A vector listing the Time to Maturities of each listed futures contract from the current observation point.
#'@param dt Numeric. The length of the discrete time step (years).
#'
#'@details
#'\loadmathjax
#'
#'The fit of the models theoretical volatility term structure of futures returns to those obtained directly from observed futures prices can be used as an additional measure of robustness for
#'the models ability to explain the behavior of a commodities term structure. A commodity pricing model should capture all dynamics of a commodities term structure,
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



