---
bibliography: references.bib
---

# ImbRegSamp

Sampling methods for imbalanced regression. Implemented resampling methods are:

1- [SMOTER](https://link.springer.com/chapter/10.1007/978-3-642-40669-0_33) [@10.1007/978-3-642-40669-0_33]

2- [SMOGN](https://proceedings.mlr.press/v74/branco17a) [@pmlr-v74-branco17a]

Implemented metrics are:

1- WSSE (weighted sum of squared error)

2- WMSE (weighted mean squared error)

3- WMAD (weighted mean absolute deviation)

4- WMAPE (weighted mean absolute percentage error)

Implemented relevance functions are:

1- PCHIP (Piecewise Cubic Hermite Interpolating Polynomials) [@torgo2007utility]

2- [Inverse Kernel Density Estimation](https://link.springer.com/article/10.1007/s10994-021-06023-5) [@steininger2021]

More methods will be added.

# R installation

devtools::install_github("<https://github.com/fatihsaglam/ImbRegSamp>")

# References
