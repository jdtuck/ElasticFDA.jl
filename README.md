# ElasticFDA

[![ElasticFDA](http://pkg.julialang.org/badges/ElasticFDA_release.svg)](http://pkg.julialang.org/?pkg=ElasticFDA&ver=release)
[![Build Status](https://travis-ci.org/jdtuck/ElasticFDA.jl.svg?branch=master)](https://travis-ci.org/jdtuck/ElasticFDA.jl)
[![Coverage Status](https://coveralls.io/repos/jdtuck/ElasticFDA.jl/badge.png?branch=master)](https://coveralls.io/r/jdtuck/ElasticFDA.jl?branch=master)

A julia package for functional data analysis using the square root slope framework
and curves using the square root velocity framework which performs pair-wise and
group-wise alignment as well as modeling using functional component analysis and
regression.

### Installation
This package can be installed using and is only currently supported on linux

    Pkg.add("ElasticFDA")

This package relies on two c/cpp optimization routines which will either compile
with icc or g++. One of the libraries relies LAPACK and BLAS. The makefile will
detect if icc is installed and use it, otherwise it will default to g++. If icc
is detected it will use MKL as the BLAS and LAPACK implementation. Otherwise
OpenBLAS is used/required.

We are working on an compiling in Windows and OSX, though time is limited and
any help would be appreciated.

Documentation is in the works and will hopefully be done soon.

### References
Tucker, J. D. 2014, Functional Component Analysis and Regression using Elastic
Methods. Ph.D. Thesis, Florida State University.

Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds with
Applications. Ph.D. Thesis, Florida State University.

Srivastava, A., Wu, W., Kurtek, S., Klassen, E. and Marron, J. S. (2011). Registration of Functional Data Using Fisher-Rao Metric. arXiv:1103.3817v2 [math.ST].

Tucker, J. D., Wu, W. and Srivastava, A. (2013). Generative models for functional data using phase and amplitude separation. Computational Statistics and Data Analysis 61, 50-66.

Joshi, S.H., Srivastava, A., Klassen, E. and Jermyn, I. (2007). A Novel Representation for Computing Geodesics Between n-Dimensional Elastic Curves. IEEE Conference on computer Vision and Pattern Recognition (CVPR), Minneapolis, MN.

Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.

Wen Huang, Kyle A. Gallivan, Anuj Srivastava, Pierre-Antoine Absil. "Riemannian
Optimization for Elastic Shape Analysis", Short version, The 21st International
Symposium on Mathematical Theory of Networks and Systems (MTNS 2014).

