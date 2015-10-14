Getting Started
===============

Installation
------------
The ElasticFDA package is available through the Julia package system by running
``Pkg.add("ElasticFDA")``. Throughout, we have assume that you have installed
the package.

This package relies on two c/cpp optimization routines which will either compile
with icc or g++. One of the libraries relies LAPACK and BLAS. The makefile will
detect if icc is installed and use it, otherwise it will default to g++. If icc
is detected it will use MKL as the BLAS and LAPACK implementation. Otherwise
OpenBLAS is used/required.

References
----------
This package is based on code from the following publications:

* Tucker, J. D. 2014, Functional Component Analysis and Regression using Elastic Methods. Ph.D. Thesis, Florida State University.

* Robinson, D. T. 2012, Function Data Analysis and Partial Shape Matching in the Square Root Velocity Framework. Ph.D. Thesis, Florida State University.

* Huang, W. 2014, Optimization Algorithms on Riemannian Manifolds with Applications. Ph.D. Thesis, Florida State University.

* Srivastava, A., Wu, W., Kurtek, S., Klassen, E. and Marron, J. S. (2011). Registration of Functional Data Using Fisher-Rao Metric. arXiv:1103.3817v2 [math.ST].

* Tucker, J. D., Wu, W. and Srivastava, A. (2013). Generative models for functional data using phase and amplitude separation. Computational Statistics and Data Analysis 61, 50-66.

* Tucker, J.D., Wu, W. and Srivastava, A. (2014). Phase-Amplitude Separation of Proteomics Data Using Extended Fisher-Rao Metric. Electronic Journal of Statistics 8 (2), 1724-1733.

* Tucker, J.D., Wu, W. and Srivastava, A. (2014). Analysis of signals under compositional noise With applications to SONAR data. IEEE Journal of Oceanic Engineering 29 (2),  318-330.

* Kurtek, S., Srivastava, A. and Wu, W. (2011). Signal estimation under random time-warpings and nonlinear signal alignment. In Proceedings of Neural Information Processing Systems (NIPS).

* Joshi, S.H., Srivastava, A., Klassen, E. and Jermyn, I. (2007). A Novel Representation for Computing Geodesics Between n-Dimensional Elastic Curves. IEEE Conference on computer Vision and Pattern Recognition (CVPR), Minneapolis, MN.

* Srivastava, A., Klassen, E., Joshi, S., Jermyn, I., (2011). Shape analysis of elastic curves in euclidean spaces. Pattern Analysis and Machine Intelligence, IEEE Transactions on 33 (7), 1415-1428.

* Wen Huang, Kyle A. Gallivan, Anuj Srivastava, Pierre-Antoine Absil. "Riemannian Optimization for Elastic Shape Analysis", Short version, The 21st International Symposium on Mathematical Theory of Networks and Systems (MTNS 2014).

