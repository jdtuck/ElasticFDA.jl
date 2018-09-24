# ElasticFDA
*Julia library for elastic functional data analysis*

[![Build Status](https://travis-ci.org/jdtuck/ElasticFDA.jl.svg?branch=master)](https://travis-ci.org/jdtuck/ElasticFDA.jl)
[![Build status](https://img.shields.io/appveyor/ci/jdtuck/elasticfda-jl.svg?style=flat-square&label=windows)](https://ci.appveyor.com/project/jdtuck/elasticfda-jl/branch/master)
[![Coverage Status](http://img.shields.io/coveralls/jdtuck/ElasticFDA.jl.svg?style=flat-square)](https://coveralls.io/r/jdtuck/ElasticFDA.jl?branch=master)

[![ElasticFDA](http://pkg.julialang.org/badges/ElasticFDA_0.7.svg)](http://pkg.julialang.org/?pkg=ElasticFDA)
[![ElasticFDA](http://pkg.julialang.org/badges/ElasticFDA_1.0.svg)](http://pkg.julialang.org/?pkg=ElasticFDA)

A Julia package for functional data analysis using the square root slope framework
and curves using the square root velocity framework which performs pair-wise and
group-wise alignment as well as modeling using functional component analysis and
regression.

### Installation
This package can be installed using and is only currently supported on linux

    (v1.0) add ElasticFDA

This package relies on two c/cpp optimization routines which will either compile
with icc or g++. One of the libraries relies LAPACK and BLAS. The makefile will
detect if icc is installed and use it, otherwise it will default to g++. If icc
is detected it will use MKL as the BLAS and LAPACK implementation. Otherwise
OpenBLAS is used/required.

### Doumentation
<http://elasticfdajl.readthedocs.org/en/latest/>

### References
Tucker, J. D. 2014, Functional Component Analysis and Regression using Elastic
Methods. Ph.D. Thesis, Florida State University.

Robinson, D. T. 2012, Function Data Analysis and Partial Shape Matching in the
Square Root Velocity Framework. Ph.D. Thesis, Florida State University.

Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds with
Applications. Ph.D. Thesis, Florida State University.

Srivastava, A., Wu, W., Kurtek, S., Klassen, E. and Marron, J. S. (2011).
Registration of Functional Data Using Fisher-Rao Metric. arXiv:1103.3817v2
[math.ST].

Tucker, J. D., Wu, W. and Srivastava, A. (2013). Generative models for
functional data using phase and amplitude separation. Computational Statistics
and Data Analysis 61, 50-66.

J. D. Tucker, W. Wu, and A. Srivastava, ``Phase-Amplitude Separation of
Proteomics Data Using Extended Fisher-Rao Metric," Electronic Journal of
Statistics, Vol 8, no. 2. pp 1724-1733, 2014.

J. D. Tucker, W. Wu, and A. Srivastava, "Analysis of signals under compositional
noise With applications to SONAR data," IEEE Journal of Oceanic Engineering, Vol
29, no. 2. pp 318-330, Apr 2014.

S. Kurtek, A. Srivastava, and W. Wu. Signal estimation under random
time-warpings and nonlinear signal alignment. In Proceedings of Neural
Information Processing Systems (NIPS), 2011.

Joshi, S.H., Srivastava, A., Klassen, E. and Jermyn, I. (2007). A Novel
Representation for Computing Geodesics Between n-Dimensional Elastic Curves.
IEEE Conference on computer Vision and Pattern Recognition (CVPR), Minneapolis, MN.

Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of
elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence,
IEEE Transactions on 33 (7), 1415-1428.

Wen Huang, Kyle A. Gallivan, Anuj Srivastava, Pierre-Antoine Absil. "Riemannian
Optimization for Elastic Shape Analysis", Short version, The 21st International
Symposium on Mathematical Theory of Networks and Systems (MTNS 2014).

Q. Xie, S. Kurtek, E. Klassen, G. E. Christensen and A. Srivastava. Metric-based
pairwise and multiple image registration. IEEE European Conference on Computer
Vision (ECCV), September, 2014

Cheng, W., Dryden, I. L., & Huang, X. (2016). Bayesian registration of functions
and curves. Bayesian Analysis, 11(2), 447–475.

Y. Lu, R. Herbei and S. Kurtek (2017). "Bayesian Registration of Functions with a Gaussian Process Prior." Journal of Computational and Graphical Statistics: in press: 1-34
