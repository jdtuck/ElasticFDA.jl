"""
Optimization function to calculate warping for elastic regression

    regression_warp(beta, timet, q, y, alpha)
    :param beta: regression function
    :param timet: vector describing time samples
    :param q: vector describing srsf
    :param y: response value
    :param alpha: intercept

    Returns
    :return gamma_new: new gamma
"""
function regression_warp(beta::Vector, timet::Vector, q::Vector, y::Float64,
                         alpha::Float64)
    gam_M = optimum_reparam(beta, timet, q, method="DP2");
    qM = warp_q_gamma(timet, q, gam_M);
    y_M = trapz(timet, qM.*beta);

    gam_m = optimum_reparam(-1.0*beta, timet, q, method="DP2");
    qm = warp_q_gamma(timet, q, gam_m);
    y_m = trapz(timet, qm.*beta);

    if y > alpha+y_M
        gamma_new = gam_M;
    elseif y < alpha + y_m
        gamma_new = gam_m;
    else
        gamma_new = zero_crossing(y-alpha, q, beta, timet, y_M, y_m,
                                  gam_M, gam_m);
    end

    return gamma_new

end


"""
Calculate warping for logistic regression

    logistic_warp(beta, timet, q, y)
    :param beta: regression function
    :param timet: time samples
    :param q: srsf
    :param y: response

    Returns
    :return gamma: new gamma
"""
function logistic_warp(beta::Vector, timet::Vector, q::Array, y)
    if y == 1
        gamma = optimum_reparam(beta, timet, q, method="DP");
    elseif y == -1
        gamma = optimum_reparam(-1.0*beta, timet, q, method="DP");
    end

    return gamma

end


"""
Calculate logistic optimization function

    logit_optm(x::Vector, grad::Vector, Phi, y)
    :param x: samples
    :param grad: gradient
    :param Phi: coefficient matrix
    :param y: response
"""
function logit_optm(x::Vector, grad::Vector, Phi, y)
    if length(grad) > 0
        logit_gradient!(x, grad, Phi, y);
    end
    out = logit_loss(x, Phi, y);
    return out
end


"""
Logistic function

    phi(t)

    Returns
    :return out: phi(t)
"""
function phi(t)
    idx = t .> 0;
    out = Array{Float64}(undef, size(t));
    if sum(idx) > 0
        out[idx] = 1.0/(1 .+ exp.(-1*t[idx]));
    end
    exp_t = exp.(t[.~idx]);
    out[.~idx] = exp_t ./ (1 .+ exp_t);
    return out
end


"""
Calculate logistic loss function

    logit_loss(b, X, y)
    :param b: coefficients
    :param X: matrix
    :param y: response

    Returns
    :return out: loss function
"""
function logit_loss(b, X, y)
    z = X * b;
    yz = y .* z;
    idx = yz .> 0;
    out = zeros(size(yz));
    out[idx] = log.(1 .+ exp.(-1.0*yz[idx]));
    out[.~idx] = (-yz[.~idx] + log.(1 .+ exp.(yz[.~idx])));
    out = sum(out);
    return out
end


"""
Calculate gradient of logistic optimization in place

    logit_gradient!(b, grad, X, y)
    :param b: coefficients
    :param grad: gradient
    :param X: matrix
    :param y: response
"""
function logit_gradient!(b, grad, X, y)
    z = X * b;
    z = phi(y.*z);
    z0 = (z-1) .* y;
    grad[:] = transpose(X) * z0;
end


"""
Calculate Hessian of logistic optimization

    logit_hessian(s, b, X, y)
    :param s:
    :param b: coefficients
    :param X: matrix
    :param y: response
"""
function logit_hessian(s, b, X, y)
    z = X*b;
    z = phi(y.*z);
    d = z*(1-z);
    wa = d.*(X * s);
    Hs = transpose(X) * wa;
    return Hs
end


"""
Calculate m-logistic warping using gradient method

    mlogit_warp_grad(alpha, beta, timet, q, y; max_itr=8000, tol=1e-10,
                     delt=0.008, display=0)
    :param alpha: intercept
    :param beta: regression function
    :param timet: vector describing time samples
    :param q: srsf
    :param y: response
    :param max_itr: maximum number of iterations
    :param tol: stopping tolerance
    :param delt: gradient step size
    :param display: display optimization iterations
"""
function mlogit_warp_grad(alpha, beta, timet, q, y; max_itr=8000,
                          tol=1e-10, delt=0.008, display=0)
    m1 = length(timet);
    m = size(beta,2);
    q /= norm(q);
    alpha /= norm(alpha);
    for i in 1:m
        beta[:, i] = beta[:, i] / norm(beta[:, i]);
    end
    gam1 = collect(LinRange(0, 1, m1));
    gamout = zeros(m1);
    beta1 = reshape(beta, m1*m, 1);

    ccall((:mlogit_warp_grad, libfdasrsf), Cvoid,
          (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ref{Int32}, Ref{Float64},
          Ref{Float64},  Ref{Int32}, Ptr{Float64}),
          m1, m, alpha, beta1, timet, gam1, q, y, max_itr, tol, delt,
          display, gamout)

    return gamout
end


"""
Calculate warping for m-logistic elastic regression

    mlogit_optim(x, grad, Phi, Y)
    :param x: sample
    :param grad: gradient
    :param Phi: matrix
    :param Y: response matrix
"""
function mlogit_optm(x::Vector, grad::Vector, Phi, Y)
    if length(grad) > 0
        mlogit_gradient!(x, grad, Phi, Y);
    end
    out = mlogit_loss(x, Phi, Y);
    return out
end


"""
Calculate loss for m-logistic regression

    mlogit_loss(b, X, Y)
    :param b: coefficients
    :param X: matrix
    :param Y: response matrix
"""
function mlogit_loss(b, X, Y)
    N, m = size(Y);  # n_samples, n_classes
    M = size(X,2);   # n_features
    B = reshape(b, M, m);
    Yhat = X * B;
    Yhat -= repeat(minimum(Yhat, dims=2),1,size(Yhat,2));
    Yhat = exp.(-1.0*Yhat);
    # l1-normalize
    Yhat = Yhat./repeat(sum(Yhat, dims=2),1,size(Yhat,2));

    Yhat = Yhat .* Y;
    nll = sum(log.(sum(Yhat,dims=2)));
    nll /= -N;

    return nll
end


"""
Calculate m-logistic elastic regression loss function gradient in place

    mlogit_gradient!(b, grad, X, Y)
    :param b: coefficients
    :param grad: gradient
    :param X: matrix
    :param Y: response matrix
"""
function mlogit_gradient!(b, grad, X, Y)
    N, m = size(Y);  # n_samples, n_classes
    M = size(X,2); # n_features
    B = reshape(b, M, m);
    Yhat = X * B;
    Yhat -= repeat(minimum(Yhat, 2),1,size(Yhat,2));
    Yhat = exp(-1.0*Yhat);
    # l1-normalize
    Yhat = Yhat./repeat(sum(Yhat, 2),1,size(Yhat,2));

    _Yhat = Yhat .* Y;
    _Yhat = _Yhat./repeat(sum(_Yhat,2),1,size(_Yhat,2));
    Yhat -= _Yhat;
    grad1 = transpose(X) * Yhat;
    grad1 /= -N;
    grad[:] = reshape(grad1, M*m, 1);
end
