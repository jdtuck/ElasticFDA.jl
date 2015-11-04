Curve Alignment
===============

These functions are for processing of N-D curves using the *square-root
velocity framework (srvf)*

SRVF Functions
--------------
.. function:: curve_to_q(beta)

    Convert curve to square-root velocity function (srvf)

    ``beta`` is an array of shape (n,T) describing the curve, where n is the dimension and T is the number of sample points

.. function:: q_to_curve(q)

    Convert srvf to curve

    ``q`` is an array of shape (n,T) describing the srvf, where n is the dimension and T is the number of sample points

.. function:: optimum_reparam(beta1, beta2, lam, method="DP", w=0.01, rotated=true, isclosed=false))

    Calculates the optimum reparamertization (warping) between two curves beta1 and beta2, using the srvf framework

    + ``beta1`` array (n,T) describing curve 1
    + ``beta2`` array (n,T) describing curve 2
    + ``lam`` control amount of warping (default=0.0)
    + ``method`` optimization method to find warping, default is Dynamic Programming ("DP"). Other options are Coordinate Descent ("DP2"), Riemanain BFGS ("LRBFGS").
    + ``w`` Controls LRBFGS (default = 0.01)
    + ``rotated`` calculate rotation (default = true)
    + ``isclosed`` closed curve (default = false)

    Returns:

    + ``gam`` warping function
    + ``R`` rotation matrix
    + ``tau`` seed value

.. function:: calc_shape_dist(beta1, beta2)

    Calculate elastic shape distance between two curves beta1 and beta2

    ``beta1`` and ``beta2`` are arrays of shape (n,T) describing the curve, where n is the dimension and T is the number of sample points

.. function:: resamplecurve(x, N=100)

    Resmaples Curve

    + ``x`` array describing curve (n,T)
    + ``N`` Number of samples to re-sample curve, N usually is > T

Alignment and Statistics
------------------------
.. function:: curve_srvf_align(beta; mode='O', maxit=20)

    Aligns a collection of curves using the elastic square-root velocity (srvf) framework.

    + ``beta`` array (n,T,N) for N number of curves
    + ``mode`` Open ('O') or Closed ('C') curves
    + ``maxit`` maximum number of iterations

    Returns:

    + ``betan`` aligned curves
    + ``qn`` aligned srvfs
    + ``betamean`` mean curve
    + ``q_mu`` mean srvf

.. function:: curve_karcher_mean(beta; mode='O', maxit=20)

    Calculates Karcher mean of a collection of curves using the elastic square-root velocity (srvf) framework.

    + ``beta`` array (n,T,N) for N number of curves
    + ``mode`` Open ('O') or Closed ('C') curves
    + ``maxit`` maximum number of iterations

    Returns:

    + ``mu`` mean srvf
    + ``betamean`` mean curve
    + ``v`` shooting vectors
    + ``q`` array of srvfs

.. function:: curve_karcher_cov(betamean, beta; mode='O')

    Calculate Karcher Covariance of a set of curves

    + ``betamean`` array (n,T) of mean curve
    + ``beta`` array (n,T,N) for N number of curves
    + ``mode`` Open ('O') or Closed ('C') curves

    Returns:

    + ``K`` covariance matrix

.. function:: curve_principal_directions(betamean, mu, K; mode='O', no=3, N=5)

    Calculate principal directions of a set of curves

    + ``betamean`` array (n,T) of mean curve
    + ``mu`` array (n,T) of mean srvf
    + ``K`` array (T,T) covariance matrix
    + ``mode`` Open ('O') or Closed ('C') curve
    + ``no`` number of components
    + ``N`` number of samples on each side of mean

    Returns:

    + ``pd`` array describing principal directions

.. function:: sample_shapes(mu, K; mode='O', no=3, numSamp=10)

    Sample shapes from model

    + ``mu`` array (n,T) mean srvf
    + ``K`` array (T,T) covariance matrix
    + ``mode`` Open ('O') or Closed ('C') curves
    + ``no`` number of principal components
    + ``numSamp`` number of samples

    Return:

    + ``samples`` array (n,T,numSamp) of sample curves
