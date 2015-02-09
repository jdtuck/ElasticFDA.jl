function regression_warp(beta::Vector, timet::Vector, q::Vector, y::Float64,
                         alpha::Float64)
    gam_M = optimum_reparam(beta, timet, q);
    qM = warp_q_gamma(timet, q, gam_M);
    y_M = trapz(timet, qM.*beta);

    gam_m = optimum_reparam(-1.*beta, timet, q);
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


function logistic_warp(beta::Vector, timet::Vector, q::Array, y)
    if y == 1
        gamma = optimum_reparam(beta, timet, q);
    elseif y == -1
        gamma = optimum_reparam(-1.*beta, timet, q);
    end

    return gamma

end


function logit_optm(x::Vector, grad::Vector, Phi, y)
    if length(grad) > 0
        logit_gradient!(x, grad, Phi, y);
    end
    out = logit_loss(x, Phi, y);
    return out
end


function phi(t)
    idx = t .> 0;
    out = Array(Float64, length(t));
    if sum(idx) > 0
        out[idx] = 1./(1+exp(-1*t[idx]));
    end
    exp_t = exp(t[~idx]);
    out[~idx] = exp_t ./ (1 + exp_t);
    return out
end


function logit_loss(b, X, y)
    z = X * b;
    yz = y .* z;
    idx = yz .> 0;
    out = zeros(yz);
    out[idx] = log(1+exp(-1.*yz[idx]));
    out[~idx] = (-yz[~idx] + log(1+exp(yz[~idx])));
    out = sum(out);
    return out
end


function logit_gradient!(b, grad, X, y)
    z = X * b;
    z = phi(y.*z);
    z0 = (z-1) .* y;
    grad[:] = X.' * z0;
end


function logit_hessian(s, b, X, y)
    z = X*b;
    z = phi(y.*z);
    d = z*(1-z);
    wa = d.*(X * s);
    Hs = X.' * wa;
    return Hs
end


function mlogit_warp_grad(alpha, beta, timet, q, y; max_itr=8000,
                          tol=1e-10, delt=0.008, display=0)
    m1 = length(timet);
    m = size(beta,2);
    q /= norm(q);
    alpha /= norm(alpha);
    for i in 1:m
        beta[:, i] = beta[:, i] / norm(beta[:, i]);
    end
    gam1 = linspace(0, 1, m1);
    gamout = zeros(m1);
    beta1 = reshape(beta, m1*m, 1);

    ccall((:mlogit_warp_grad, libfdasrsf), Void,
          (Ptr{Int32}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64},
          Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Int32}, Ptr{Float64},
          Ptr{Float64},  Ptr{Int32}, Ptr{Float64}),
          &m1, &m, alpha, beta1, timet, gam1, q, y, &max_itr, &tol, &delt,
          &display, gamout)

    return gamout
end


function mlogit_optm(x::Vector, grad::Vector, Phi, Y)
    if length(grad) > 0
        mlogit_gradient!(x, grad, Phi, Y);
    end
    out = mlogit_loss(x, Phi, Y);
    return out
end


function mlogit_loss(b, X, Y)
    N, m = size(Y);  # n_samples, n_classes
    M = size(X,2); # n_features
    B = reshape(b, M, m);
    Yhat = X * B;
    Yhat -= repmat(minimum(Yhat, 2),1,size(Yhat,2));
    Yhat = exp(-1.*Yhat);
    # l1-normalize
    Yhat = Yhat./repmat(sum(Yhat, 2),1,size(Yhat,2));

    Yhat = Yhat .* Y;
    nll = sum(log(sum(Yhat,2)));
    nll /= -N;

    return nll
end


function mlogit_gradient!(b, grad, X, Y)
    N, m = size(Y);  # n_samples, n_classes
    M = size(X,2); # n_features
    B = reshape(b, M, m);
    Yhat = X * B;
    Yhat -= repmat(minimum(Yhat, 2),1,size(Yhat,2));
    Yhat = exp(-1.*Yhat);
    # l1-normalize
    Yhat = Yhat./repmat(sum(Yhat, 2),1,size(Yhat,2));

    _Yhat = Yhat .* Y;
    _Yhat = _Yhat./repmat(sum(_Yhat,2),1,size(_Yhat,2));
    Yhat -= _Yhat;
    grad1 = X.' * Yhat;
    grad1 /= -N;
    grad[:] = reshape(grad1, M*m, 1);
end

