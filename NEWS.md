# NFCP 1.1.0

- The 'NFCP' package now allows for the consideration of deterministic seasonality within commodity pricing models throughout the package.
- The 'American_option_value' function now calculates American option on futures contracts of commodities, which is more reflective of actual options in commodity markets.
- The vignette of 'NFCP' has been updated, with more detail throughout the document.
- Deprecated function names removed
- Fixed a bug in 'NFCP_MLE' when the ... argument 'hessian' is FALSE.
- Removed the 'richardsons_extrapolation' argument of the 'NFCP_MLE' function, however it is still supported within the function.
- Minor documentation changes

# NFCP 1.0.1

- Minor bug fix in the 'NFCP_Kalman_filter' function

# NFCP 1.0.0

- Original function names of the NFCP package are now deprecated and will be removed in a future version

- Extensive function and parameter name changes to better fit the tidyverse style guide

- Functions throughout the NFCP package now support the application of measurement errors that group contracts by set maturities. This is particularly useful for model estimation that uses full contract data, as it allows for the relaxing of the assumption that measurement errors are identical across observations.

- Bug fixes for European option pricing

- Other bug fixes and changes

# NFCP 0.2.1

- Minor bug fix to Ornstein-Uhlenbeck simulation of the 'Spot.Price.Simulate' function

# NFCP 0.2.0

- NFCP now depends on the 'LSMRealOptions' package
- American put option pricing is now supported through the 'AmericanOptionValue' function and the 'LSMRealOptions' package
- Minor adjustments to parameters of 'TSFit.Volatility' function
- 'NFCP.bounds' has been renamed to 'NFCP.Domains' and now allows for further customization of the search applied in 'NFCP.MLE'
- Added American option pricing example to vignette
- The 'NFCP.MLE' function now orders outputs in terms of increasing mean-reversion rates
- The 'NFCP.Kalman.filter' function now returns the bias and RMSE of input parameters when 'verbose = TRUE' through the 'Filtered.Error' list object
- The 'NFCP.MLE' function now returns the user, system and elapsed time through the 'proc_time' list object.
- Minor bug fixes, particularly to the 'Spot.Price.Simulate' function


# NFCP 0.1.0

* Release of NFCP
