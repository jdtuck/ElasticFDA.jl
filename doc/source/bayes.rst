Bayesian Alignment
==================

The following functions align functional data in the srsf framework using a
Bayesian approach. These functions are *experimental* and results are not
fully tested

Alignment
---------
.. function:: pair_warping_baye(f1, f2; iter=15000, times=5, powera=1, showplot=true)

    Compute pair warping between two functions using Bayesian method

    + ``f1, f2`` vectors describing functions
    + ``iter`` number of iterations
    + ``times`` MCMC parameter
    + ``powera`` MCMC parameter
    + ``showplot`` show plots

    Returns Dict containing:

    + ``f1`` function f1,
    + ``f2_q`` srsf registration,
    + ``gam_q`` warping funtion,
    + ``f2a`` registered f2,
    + ``gam`` warping function,
    + ``dist_collect`` distance,
    + ``best_match`` best match,

.. function:: group_warping bayes(f; iter=20000, times=5, powera=1, showplot=true)

    Group alignment of functions using Bayesian method

    + ``f`` array (M,N) of N functions,
    + ``iter`` number of MCMC iterations,
    + ``times`` time slicing,
    + ``powera`` MCMC parameter,
    + ``showplot`` show plots,

