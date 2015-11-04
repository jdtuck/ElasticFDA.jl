Image Alignment
===============

These functions are for processing of images using the *q-map framework*

Alignment
---------
.. function:: pair_align_image(I1, I2; M=5, ortho=true, basis_type="t",
                               resizei=true, N=64, stepsize=1e-5, itermax=1000)

    Pairwise align two images

    + ``I1`` reference image
    + ``I2`` image to warp
    + ``M`` number of basis elements
    + ``ortho`` orthonormalize basis
    + ``basis_type`` type of basis ("t", "s", "i", "o")
    + ``resizei`` resize image
    + ``N`` size of resized image
    + ``stepsize`` gradient stepsize
    + ``itermax`` maximum number of iterations

    Returns:

    + ``I2_new`` aligned I2
    + ``gam`` warping function

