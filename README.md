# quant

VEGA quant risk package

Relies on gonum.org

Current set-up:
- misc package for various basic numerical calculations that are not problem-specific
- riskmeasures package that calculates risk measures for various distributions as well as empirical data
- bsformula all things related to the Black-Scholes formula (call / put prices, greeks)
- riskmodelsbs the risk model for Forwards and European calls / puts based on the Black-Scholes model i.e. log-normal distributions of future prices
