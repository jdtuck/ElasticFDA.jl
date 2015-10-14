Functional Principal Component Analysis
=======================================

These functions are for computing functional principal component anlaysis
(fPCA) on aligned data and generating random samples

fPCA Functions
--------------
.. function:: vert_fPCA(fn, timet, qn; no=1)

    Calculates vertical functional principal component analysis on aligned data

    + ``fn`` array of shape (M,N) of N aligned functions with M samples
    + ``timet`` vector of size M describing the sample points
    + ``qn`` array of shape (M,N) of N aligned SRSF with M samples
    + ``no`` number of components to extract (default = 1)

    Returns Dict containing:

    + ``q_pca`` srsf principal directions
    + ``f_pca`` functional principal directions
    + ``latent`` latent values
    + ``coef`` coefficients
    + ``U`` eigenvectors

.. function:: horiz_fPCA(gam, timet; no=1)

    Calculates horizontal functional principal component analysis on aligned data

    + ``gam`` array of shape (M,N) of N warping functions with M samples
    + ``timet`` vector of size M describing the sample points
    + ``no`` number of components to extract (default = 1)

    Returns Dict containing:

    + ``gam_pca`` warping principal directions
    + ``psi_pca`` srsf functional principal directions
    + ``latent`` latent values
    + ``U`` eigenvectors
    + ``gam_mu`` mean warping function

.. function:: gauss_model(fn, timet, qn, gam; n=1, sort_samples=false)

    Computes random samples of functions from aligned data using Gaussian model

    + ``fn`` aligned functions (M,N)
    + ``timet`` vector (M) describing time
    + ``qn`` aligned srvfs (M,N)
    + ``gam`` warping functions (M,N)
    + ``n`` number of samples
    + ``sort_samples`` sort samples

    Returns Dict containing:

    + ``fs`` random aligned functions
    + ``gams`` random warping functions
    + ``ft`` random functions

