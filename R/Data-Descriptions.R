
#'Crude Oil Term Structure Futures Data (1990 - 1995)
#'
#' @docType data
#'
#' @usage data(SS.Oil)
#'
#'
#' @keywords datasets
#'
#'@description The \code{SS.Oil} \code{list} object features the approximate weekly observations of Crude Oil (WTI) futures contracts used to develop a two-factor
#'commodity pricing model within the prominent work of Schwartz and Smith (2000) titled: "Short-Term Variations and long-Term Dynamics in Commodity Prices".
#'The two-factor commodity pricing model presented within this study is also included. The \code{SS.Oil} list object is used extensively within the
#'\code{NFCP} package to provide working examples and showcase the features of the package.
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
