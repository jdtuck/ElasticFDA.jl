Elastic Functional Regression
=============================

These functions compute elastic standard, logistic, and m-logistic regression
models. This code is experimental and results are not guaranteed

Regression Models and Prediction
--------------------------------
.. function:: elastic_regression(f, y, timet; B=None, lambda=0, df=20, max_itr=20, smooth=false)

    Calculate elastic regression from function data f, for response y

    ``f`` array (M,N) of N functions
    ``y`` vector (N) of responses
    ``timet`` vector (N) describing time samples
    ``B`` matrix describing basis functions (M,N) (default=None generates a B-spline basis
    ``lambda`` regularization parameter
    ``df`` degree of freedom of basis
    ``max_itr`` maximum number of iterations
    ``smooth`` smooth data

    Returns Dict describing regression
    ``alpha`` intercept
    ``beta`` regression function
    ``fn`` aligned functions
    ``qn`` aligned srsfs
    ``gamma`` warping functions
    ``q`` original srsfs
    ``B`` basis functions
    ``type`` type of regression
    ``b`` coefficients
    ``SSE`` sum of squared error

.. function:: elastic_logistic(f, y, timet; B=None, df=20, max_itr=20, smooth=false)

    Calculate elastic logistic regression from function data f, for response y

    ``f`` array (M,N) of N functions
    ``y`` vector (N) of responses
    ``timet`` vector (N) describing time samples
    ``B`` matrix describing basis functions (M,N) (default=None generates a B-spline basis
    ``df`` degree of freedom of basis
    ``max_itr`` maximum number of iterations
    ``smooth`` smooth data

    Returns Dict describing regression
    ``alpha`` intercept
    ``beta`` regression function
    ``fn`` aligned functions
    ``qn`` aligned srsfs
    ``gamma`` warping functions
    ``q`` original srsfs
    ``B`` basis functions
    ``type`` type of regression
    ``b`` coefficients
    ``LL`` logistic loss

.. function:: elastic_mlogistic(f, y, timet; B=None, df=20, max_itr=20, smooth=false)

    Calculate elastic m-logistic regression from function data f, for response y

    ``f: array (M,N) of N functions
    ``y: vector (N) of responses
    ``timet: vector (N) describing time samples
    ``B: matrix describing basis functions (M,N) (default=None generates a B-spline basis
    ``df: degree of freedom of basis
    ``max_itr: maximum number of iterations
    ``smooth: smooth data

    Returns Dict describing regression
    ``alpha`` intercept
    ``beta`` regression function
    ``fn`` aligned functions
    ``qn`` aligned srsfs
    ``gamma`` warping functions
    ``q`` original srsfs
    ``B`` basis functions
    ``type`` type of regression
    ``b`` coefficients
    ``n_classes`` number of classes
    ``LL`` logistic loss

.. function:: elastic_prediction(f, timet, model; y=None, smooth=false)

    Prediction from elastic regression model

    ``f`` functions to predict
    ``timet`` vector describing time samples
    ``model`` calculated model (regression, logistic, mlogistic)
    ``y`` true responses (default = None)
    ``smooth`` smooth data (default = false)

    Returns
    ``y_pred`` predicted value
    ``y_labels`` labels of predicted value
    ``Perf`` Performance metric if truth is supplied

