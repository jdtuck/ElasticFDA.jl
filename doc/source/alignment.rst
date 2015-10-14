Functional Alignment
====================

The main functions deal with the alignment of functional data using the
*square-root slope (srsf)* framework. Where an input into a function is
expecting an array of functions the shape is assumed to be ``(M,N)`` with ``M``
being the number of sample points and ``N`` being the number of functions.

SRSF Functions
-------------
.. function:: f_to_srsf(f::Array, timet=0, smooth=false)

    Convert function to square-root slope (srsf) functions

    ``f`` is an array of shape ``(M,N)`` as described above. By default the
    function will generate timing information, otherwise ``timet`` should be
    vector of length ``M`` describing the timing information. If
    ``smooth=true`` the input data will be smoothed first using smoothing
    splines.

.. function:: srsf_to_f(q::Array, timet, f0=0.0)

    Convert srsf to function space

    ``q`` is an array with the standard shape. ``timet`` is a vector of timing
    infomration. ``f0`` is the initial value of the function in f-space, this
    is required to make the transformation a bijection.

.. function:: smooth_data(f::Array, sparam=10)

    Smooth functional data using a box filter

    ``f`` is an array with the standard shape. ``sparam`` is the number of
    times to run the filter.

.. function:: smooth_data!(f::Array, sparam=10)

    same as smooth_data, except the smoothing is done in-place

.. function:: trapz(x::Vector, y::Array, dim=1)

    Trapezodial Integration

    ``x`` is a vector of time samples. ``y`` is the reponse and ``dim`` is the
    dimension to integrate along.

.. function:: optimum_reparam(q1, timet, q2, lam=0.0, method="SIMUL", w=0.01, f1o=0.0, f2o=0.0)

    Calculates the optimum reparamertization (warping) between two srsfs q1 and
    q2.

    ``q1`` and ``q2`` can be vectors or arrays of the standard shape. ``timet``
    is a vector describing the time samples. ``lam`` controls the amount of
    warlping. ``method`` is the optimization method to find the warping. The
    default is Simultaneous Alignment ("SIMUL"). Other options are Dynamic
    Programming ("DP" or "DP2") and Riemannian BFGS ("LRBFGS").

.. function:: warp_f_gamma(time::Vector, f::Vector, gam::Vector)

    Warp function f by warping function gamma

.. function:: warp_q_gamma(time::Vector, q::Vector, gam::Vector)

    Warp srsf q by warping function gamma

.. function:: elastic_distance(f1::Vector, f2::Vector, timet::Vector)

    Caclulates the elastic distance between two functions and returns the
    amplitude distance ``da`` and phase distance ``dp``.

.. function:: rgam(N, sigma, num)

    Generate random warping functions of length ``N``. ``sigma`` controls the
    standard deviation across the random samples and ``num`` is the number of
    random samples.

Alignment
---------
.. function:: srsf_align(f, timet; method="mean", smooth=false, sparam=10, lam=0.0, optim="SIMUL")


    Aligns a collection of functions using the elastic square-root slope (srsf)
    framework.

    ``f`` is and array of shape (M,N) of N functions with M samples
    ``timet`` is a vector of size M describing the sample points
    ``method`` (string) calculate Karcher Mean or Median (options = "mean" or "median") (default="mean")
    ``smooth`` Smooth the data using a box filter (default = false)
    ``sparam`` Number of times to run smoothing filter (default 10)
    ``lam`` controls the elasticity (default = 0)
    ``optim`` optimization method to find warping, default is Simultaneous Alignment ("SIMUL"). Other options are Dynamic Programming ("DP2"), Riemanain BFGS ("LRBFGS")

    Returns Dict containing
    ``fn`` aligned functions - array of shape (M,N) of N functions with M samples
    ``qn`` aligned srsfs - similar structure to fn
    ``q0`` original srsfs - similar structure to fn
    ``fmean`` function mean or median - vector of length N
    ``mqn`` srvf mean or median - vector of length N
    ``gam`` warping functions - similar structure to fn
    ``orig_var`` Original Variance of Functions
    ``amp_var`` Amplitude Variance
    ``phase_var`` Phase Variance

.. function:: align_fPCA(f, timet; num_comp=3, smooth=false, sparam=10)

    Aligns a collection of functions while extracting principal components.
    The functions are aligned to the principal components

    ``f`` array of shape (M,N) of N functions with M samples
    ``timet`` vector of size M describing the sample points
    ``num_comp`` Number of componets (default = 3)
    ``smooth`` Smooth the data using a box filter (default = false)
    ``sparam`` Number of times to run smoothing filter (default 10)

    Returns Dict containing
    ``fn`` aligned functions - array of shape (M,N) of N functions with M samples
    ``qn`` aligned srvfs - similar structure to fn
    ``q0`` original srvf - similar structure to fn
    ``mqn`` srvf mean or median - vector of length M
    ``gam`` warping functions - similar structure to fn
    ``q_pca`` srsf principal directions
    ``f_pca`` functional principal directions
    ``latent`` latent values
    ``coef`` coefficients
    ``U`` eigenvectors
    ``orig_var`` Original Variance of Functions
    ``amp_var`` Amplitude Variance
    ``phase_var`` Phase Variance

