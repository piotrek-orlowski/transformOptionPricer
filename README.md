# transformOptionPricer

Characteristic-function based option pricer. 

The `cosTransform` pricer uses the method described in http://ta.twi.tudelft.nl/mf/users/oosterle/oosterlee/COS.pdf. 

The `gaussLaguerre` pricer evaluates the CF integrals from http://www.math.nyu.edu/research/carrp/papers/pdf/jcfpub.pdf with the use of Gauss-Laguerre quadrature of the difference between the calculated CF and known auxiliary model (Black-Scholes).

# Installation

To install the package first load the `devtools` library and type
```
install_github(repo= "piotrek-orlowski/transformOptionPricer")
```

# Usage

The package allows for calculating vanilla option prices via numerical integration of (a transform of) the characteristic function of the log stock return. Thus, it allows for calculation prices in the Black-Scholes model, any model based on a Levy process, and in affine and quadratic-affine jump-diffusion models.

By default, the package loads `affineModelR` (https://github.com/piotrek-orlowski/affineModelR), which allows for calculating CF values in affine jump diffusion models.

# Examples

The file https://github.com/piotrek-orlowski/transformOptionPricer/blob/master/inst/tests/test-pricers.R contains example calculations and some accuracy checks.

# Authors

The package is being developed by Piotr Orłowski from a codebase started together with Andras Sali (https://github.com/andrewsali/).
