
#'Volatility Term Structure of futures returns
#'
#'@description Estimate the Theoretical and Empirical Volatility Term Structure of futures returns
#'
#'@param parameters A named vector of parameters of an N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param Futures A Matrix of futures price data. Each column corresponds to a given futures contract, and each row is an observation of the futures contracts.
#'@param TTM A vector listing the Time to Maturities of each listed Futures contract from the current observation point.
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
#'\code{TSFit.Volatility} returns a matrix with the theoretical and empirical volatility term structure of futures returns, with the number of columns of this matrix coinciding with the number of input futures contracts.
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'### Test the volatility term structure fit of the Schwartz-Smith two-factor model on crude oil:
#'V_TSFit <- TSFit.Volatility(
#'  parameters = SS.Oil$Two.Factor,
#'  Futures = SS.Oil$Stitched.Futures,
#'  TTM = SS.Oil$Stitched.TTM,
#'  dt = SS.Oil$dt)
#'
#'@export
TSFit.Volatility <- function(parameters, Futures, TTM, dt){

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


#'N-Factor Model European Option Pricing
#'@description Value European Option Put and Calls under the parameters of an N-factor model.
#'
#'@param X.0 Initial values of the state vector.
#'@param parameters Named vector of parameter values of a specified N-factor model. Function \code{NFCP.Parameters} is recommended.
#'@param t Time to expiration of the option
#'@param TTM Time to maturity of the Futures contract.
#'@param K Strike price of the European Option
#'@param r Risk-free interest rate.
#'@param call \code{logical} is the European option a call or put option?
#'@param verbose \code{logical}. Should additional information be output? see \bold{details}
#'
#'@details
#'
#'\loadmathjax
#'
#'The \code{European.Option.Value} function calculates analytic expressions of the value of European call and put options on futures contracts within the N-factor model. Under the assumption that future futures prices
#'are log-normally distributed under the risk-neutral process, there exist analytic expressions of the value of European call and put options on futures contracts. The following analytic expression follows from that presented by Schwartz and Smith (2000) extended to the N-factor framework. The value of a European option on a futures contract
#'is given by calculating its expected future value using the risk-neutral process and subsequently discounting at the risk-free rate.
#'
#'One can verify that under the risk-neutral process, the expected futures price at time \mjeqn{t}{t} is:
#'
#'\mjdeqn{E^*[F_{T,t}] = exp(\sum_{i=1}^Ne^{-\kappa_iT}x_{i}(0) + \mu^*t +  A(T-t) +
#'\frac{1}{2}(\sigma_1^2t+\sum_{i.j\neq1}{e^{-\left(\kappa_i+\kappa_j\right)\left(T-t\right)}\sigma_i\sigma_j\rho_{i,j}.\frac{1-e^{-\left(\kappa_i+\kappa_j\right)t}}{\kappa_i+\kappa_j}})) \equiv F_{T,0}}{
#'E^*[F[T,t]] = exp(sum_{i=1}^N (e^(-kappa[i]*T) * x[i](0) + mu^* * t +  A(T-t) +  1/2 (sigma[1]^2 + sum_[i.j != 1] (e^(- (kappa[i] + kappa[j])(T-t)) * sigma[i] * sigma[j] * rho[i,j] * (1 - e^(-(kappa[i] * kappa[j])t))
#'/(kappa[i] + kappa[j]))) equiv F[T,0]}
#'
#'This follows from the derivation provided within the vignette of the NFCP package as well as the details of the \code{Futures.Price.Forecast} package.
#'The equality of expected futures price at time \mjeqn{t}{t} being  equal to the time-\mjeqn{t}{t} current futures price \mjeqn{F_{T,0}}{F[T,0]} is proven by Futures prices
#'being given by expected spot prices under the risk-neutral process
#'\mjeqn{(F_{T,t}=E_t^\ast\left[S_T\right])}{F[T,t] = E[t]^*(S[T])} and the law of iterated expectations \mjeqn{\left(E^\ast\left[E_t^\ast\left[S_T\right]\right]=E^\ast\left[S_T\right]\right)}{E^*(E[t]^*(S[T])) = E^*(S[T])}
#'
#'Because future futures prices are log-normally distributed under the risk-neutral process, we can write a closed-form expression for valuing European put and call options on these futures. When \mjeqn{T=0}{T=0}
#'these are European options on the spot price of the commodity itself. The value of a European call option on a futures contract maturing at time
#'\mjeqn{T}{T}, with strike price \mjeqn{K}{K}, and with time \mjeqn{t}{t} until the option expires, is:
#'
#'\mjdeqn{e^{-rt}E^\ast\left[\max{\left(F_{T,t}-K,0\right)}\right]}{e^(-rt) * E^*(max(F[T,t] - K, 0))}
#'\mjdeqn{= e^{-rt}( F_{T,0}N(d) - KN(d-\sigma_\phi(t,T)))}{e^(-rt) * (F[T,0] * N(d) - K * N(d - sigma[phi](t,T)))}
#'
#'Where:
#'\mjdeqn{d=\frac{\ln(F/K)}{\sigma_\phi(t,T)}+\frac{1}{2}\sigma_\phi\left(t,T\right)}{d = ln(F/K) / sigma[phi](t,T) + 1/2 sigma[phi](t,T)}
#'
#'and:
#'
#'\mjdeqn{\sigma_\phi\left(t,T\right) = \sqrt(\sigma_1^2t+\sum_{i.j\neq1}{e^{-\left(\kappa_i+\kappa_j\right)\left(T-t\right)}\sigma_i\sigma_j\rho_{i,j}.
#'\frac{1-e^{-\left(\kappa_i+\kappa_j\right)t}}{\kappa_i+\kappa_j}})}{
#'sigma[phi](t,T) = sqrt( sigma[1]^2 + sum_[i.j != 1]( e^(-(kappa[i] + kappa[j])(T-t)) sigma[i] sigma[j] rho[i,j] * (1 - e^(-(kappa[i] + kappa[j])t))/(kappa[i] + kappa[j])))}
#'
#'Parameter \mjeqn{ N(d) }{N(d)} indicates cumulative probabilities for the standard normal distribution (i.e. \mjeqn{P(Z<d)}{P(Z<d)}).
#'
#'Similarly, the value of a European put with the same parameters is given by:
#'
#'\mjdeqn{e^{-rt} E^*[max(K-F_{T,t},0)]}{e^(-rt) E^*(max(K - F[T,t],0))}
#'
#'\mjdeqn{=e^{-rt}\left(-F_{T,0}N\left(-d\right)+KN\left(\sigma_\phi\left(t,T\right)-d\right)\right)}{e^(-rt) * (- F[T,0] * N(-d) + K * N(sigma[phi](t,T) - d))}
#'
#'
#'The presented option valuation formulas are analogous to the Black-Scholes formulas for valuing European options on stocks that do not pay dividends
#'
#'Under this terminology, the stock price corresponds to the present value of the futures commitment \mjeqn{(e^{-rt}F_{T,0})}{e^(-rt) F[T,0]}  and the equivalent annualized volatility would be \mjeqn{\sigma_\phi(t,T)/\sqrt t}{sigma[phi](t,T) / sqrt(t)}
#'
#'When \code{verbose = T}, the \code{European.Option.Value} function numerically calculates the sensitivity of option prices to the underlying parameters specified within the N-factor model, as well as some of the most common
#'"Greeks" related to European put and call option pricing. All gradients are calculated numerically by calling the \code{grad} function from the \code{numDeriv} package.
#'
#'@return
#'The \code{European.Option.Value} function returns a numeric value corresponding to the present value of an option when \code{verbose = F}.
#'When \code{verbose = T}, \code{European.Option.Value} returns a list with three objects:
#'
#'\tabular{ll}{
#'
#' \code{Value} \tab Present value of the option. \cr
#'
#'\code{Annualized.Volatility} \tab Annualized volatility of the option. \cr
#'
#'\code{Parameter.Sensitivity} \tab Sensitivity of the option value to each parameter of the N-factor model. \cr
#'
#'\code{Greeks} \tab Sensitivity of the option value to different option parameters. \cr
#'
#' }
#'
#'@references
#'
#'Schwartz, E. S., and J. E. Smith, (2000). Short-Term Variations and Long-Term Dynamics in Commodity Prices. \emph{Manage. Sci.}, 46, 893-911.
#'
#'Cortazar, G., and L. Naranjo, (2006). An N-factor Gaussian model of oil futures prices. \emph{Journal of Futures Markets: Futures, Options, and Other Derivative Products}, 26(3), 243-268.
#'
#'@examples
#'##Example 1 - A European 'put' option on a stock following 'GBM'
#'##growing at the risk-free rate:
#'
#'### Risk-free rate:
#'rf <- 0.05
#'### Stock price:
#'S0 <- 20
#'### Stock volatility:
#'S.sigma <- 0.2
#'### Option maturity:
#'Tt <- 1
#'### Exercise price:
#'K <- 20
#'### Calculate 'put' option price:
#'European.Option.Value(X.0 = log(S0), parameters = c(mu_star = rf, sigma_1 = S.sigma),
#'t = Tt, TTM = Tt, K = K, r = rf, call = FALSE)
#'
#'##Example 2 - A European call option under a two-factor crude oil model:
#'
#'##Step 1 - Obtain current (i.e. most recent) state vector by filtering the
#'##two-factor oil model:
#'Schwartz.Smith.Oil <- NFCP.Kalman.filter(parameter.values = SS.Oil$Two.Factor,
#'                                       parameters = names(SS.Oil$Two.Factor),
#'                                       log.futures = log(SS.Oil$Stitched.Futures),
#'                                       dt = SS.Oil$dt,
#'                                       TTM = SS.Oil$Stitched.TTM,
#'                                       verbose = TRUE)
#'
#'##Step 2 - Calculate 'call' option price:
#' European.Option.Value(X.0 = Schwartz.Smith.Oil$X.t,
#'                       parameters = SS.Oil$Two.Factor,
#'                       t = 1,
#'                       TTM = 1,
#'                       K = 20,
#'                       r = 0.05,
#'                       call = TRUE,
#'                       verbose = FALSE)
#'@export
European.Option.Value <- function(X.0, parameters, t, TTM, K, r, call, verbose = FALSE){

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

