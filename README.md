# ImbRegSamp

Sampling methods for imbalanced regression. Implemented resampling
methods are:

1-
[SMOTER](https://link.springer.com/chapter/10.1007/978-3-642-40669-0_33)
(Luís Torgo et al. 2013)

2- [SMOGN](https://proceedings.mlr.press/v74/branco17a) (Branco, Torgo,
and Ribeiro 2017)

3- Random over-under sampling

4- WERSC (Branco, Torgo, and Ribeiro 2019)

5- GNO (Gaussian Noise Oversampling) (Branco, Torgo, and Ribeiro 2019)

6- GSMOTER (Geometric SMOTE for Regression) (Camacho, Douzas, and Bacao
2022)

Implemented metrics are:

1- WSSE (weighted sum of squared error)

2- WMSE (weighted mean squared error)

3- WMAD (weighted mean absolute deviation)

4- WMAPE (weighted mean absolute percentage error)

Implemented relevance functions are:

1- PCHIP (Piecewise Cubic Hermite Interpolating Polynomials) (Luis Torgo
and Ribeiro 2007)

2- [Inverse Kernel Density
Estimation](https://link.springer.com/article/10.1007/s10994-021-06023-5)
(Steininger et al. 2021)

More methods will be added.

# R installation

devtools::install_github(“<https://github.com/fatihsaglam/ImbRegSamp>”)

# References

Branco, Paula, Luis Torgo, and Rita P. Ribeiro. 2019. “Pre-Processing
Approaches for Imbalanced Distributions in Regression.” *Neurocomputing*
343: 76–99.
https://doi.org/<https://doi.org/10.1016/j.neucom.2018.11.100>.

Branco, Paula, Luís Torgo, and Rita P. Ribeiro. 2017. “SMOGN: A
Pre-Processing Approach for Imbalanced Regression.” In *Proceedings of
the First International Workshop on Learning with Imbalanced Domains:
Theory and Applications*, edited by Paula Branco Luís Torgo and Nuno
Moniz, 74:36–50. Proceedings of Machine Learning Research. PMLR.
<https://proceedings.mlr.press/v74/branco17a.html>.

Camacho, Luís, Georgios Douzas, and Fernando Bacao. 2022. “Geometric
SMOTE for Regression.” *Expert Systems with Applications* 193: 116387.
https://doi.org/<https://doi.org/10.1016/j.eswa.2021.116387>.

Steininger, Michael, Konstantin Kobs, Padraig Davidson, Anna Krause, and
Andreas Hotho. 2021. “Density-Based Weighting for Imbalanced
Regression.” *Machine Learning* 110 (8): 2187–2211.
<https://doi.org/10.1007/s10994-021-06023-5>.

Torgo, Luis, and Rita Ribeiro. 2007. “Utility-Based Regression.” In
*Knowledge Discovery in Databases: PKDD 2007: 11th European Conference
on Principles and Practice of Knowledge Discovery in Databases, Warsaw,
Poland, September 17-21, 2007. Proceedings 11*, 597–604. Springer.

Torgo, Luís, Rita P. Ribeiro, Bernhard Pfahringer, and Paula Branco.
2013. “SMOTE for Regression.” In *Progress in Artificial Intelligence*,
edited by Luís Correia, Luís Paulo Reis, and José Cascalho, 378–89.
Berlin, Heidelberg: Springer Berlin Heidelberg.
