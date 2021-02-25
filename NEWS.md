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
