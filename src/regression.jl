"""
Calculate elastic regression from function data f, for response y

    elastic_regression(f, y, timet; B=None, lambda=0, df=20, max_itr=20,
                       smooth=false)
    :param f: array (M,N) of N functions
    :param y: vector (N) of responses
    :param timet: vector (N) describing time samples
    :param B: matrix describing basis functions (M,N) (default=None generates a
              B-spline basis
    :param lambda: regularization parameter
    :param df: degree of freedom of basis
    :param max_itr: maximum number of iterations
    :param smooth: smooth data

    Returns Dict describing regression
    :return alpha: intercept
    :return beta: regression function
    :return fn: aligned functions
    :return qn: aligned srsfs
    :return gamma: warping functions
    :return q: original srsfs
    :retrun B: basis functions
    :return type: type of regressions
    :return b: coefficients
    :return SSE: sum of squared error
"""
function elastic_regression(f::Array, y::Vector, timet::Vector; B=Union{}, lambda=0, df=20,
                            max_itr=20, smooth=false)
    M, N = size(f);

    if M > 500
        parallel = true;
    elseif N > 100
        parallel = true;
    else
        parallel = false;
    end

    binsize = mean(diff(timet));

    # Create B-spline basis if none provided
    if B == Union{}
        B = bs(timet, df, 4);
    end
    Nb = size(B,2);

    # second derivative for regularization
    Bdiff = zeros(M, Nb);
    for ii in 1:Nb
        Bdiff[:, ii] = gradient(gradient(B[:,ii], binsize), binsize);
    end

    q = f_to_srsf(f, timet, smooth);

    gamma = repeat(linspace(0,1,M), 1, N);
    fn = zeros(M, N);
    qn = zeros(M, N);

    itr = 1;
    SSE = zeros(max_itr);
    alpha = 0;
    b = zeros(Nb+1);
    while itr <= max_itr
        @printf("Iteration: %d\n", itr)
        # align data
        fn = zeros(M, N);
        qn = zeros(M, N);
        for ii in 1:N
            xout = (timet[end] - timet[1]) * gamma[:, ii] + timet[1];
            fn[:,ii] = approx(timet, f[:, ii], xout);
            qn[:,ii] = warp_q_gamma(timet, q[:,ii], gamma[:,ii]);
        end

        # OLS using basis
        Phi = ones(N, Nb+1);
        for ii in 1:N
            for jj in 2:Nb+1
                Phi[ii,jj] = trapz(timet, qn[:,ii].*B[:,jj-1]);
            end
        end

        R = zeros(Nb+1,Nb+1);
        for ii in 2:Nb+1
            for jj in 2:Nb+1
                R[ii,jj] = trapz(timet, Bdiff[:,ii-1].*Bdiff[:,jj-1]);
            end
        end

        xx = transpose(Phi) * Phi;
        inv_xx = inv(xx + lambda * R);
        xy = transpose(Phi) * y;
        b = inv_xx * xy;

        alpha = b[1];
        beta = B * b[2:Nb+1];

        # compute the SSE
        int_X = zeros(N);
        for ii in 1:N
            int_X[ii] = trapz(timet, qn[:,ii].*beta);
        end

        SSE[itr] = sum((y-alpha-int_X).^2);

        # find gamma
        gamma_new = zeros(M, N);
        if parallel
            gamma_new = @distributed (hcat) for i=1:N
                regression_warp(beta, timet, q[:, i], y[i], alpha);
            end
        else

            for i = 1:N
                gamma_new[:,i] = regression_warp(beta, timet, q[:, i], y[i],
                                                 alpha);
            end
        end

        if norm(gamma-gamma_new) < 1e-5
            break
        else
            gamma = copy(gamma_new);
        end

        itr += 1;
    end

    # Last step with centering of gam
    gamI = sqrt_mean_inverse(gamma);
    gamI_dev = gradient(gamI, 1/(M-1));
    timet0 = (timet[end] -timet[1]) .* gamI  + timet[1];
    beta = approx(timet, beta, timet0) .* sqrt(gamI_dev);

    for k in 1:N
        qn[:, k] = approx(timet, qn[:, k], timet0) .* sqrt(gamI_dev);
        fn[:, k] = approx(timet, fn[:, k], timet0);
        gamma[:, k] = approx(timet, gamma[:, k], timet0);
    end

    out = Dict("alpha" => alpha, "beta" => beta, "fn" => fn, "qn" => qn,
               "gamma" => gamma, "q" => q, "B" => B, "type" => "linear",
               "b" => b[2:end], "SSE" => SSE[1:itr-1]);
    return out
end


"""
Calculate elastic logistic regression from function data f, for response y

    elastic_logistic(f, y, timet; B=None, df=20, max_itr=20, smooth=false)
    :param f: array (M,N) of N functions
    :param y: vector (N) of responses
    :param timet: vector (N) describing time samples
    :param B: matrix describing basis functions (M,N) (default=None generates a
              B-spline basis
    :param df: degree of freedom of basis
    :param max_itr: maximum number of iterations
    :param smooth: smooth data

    Returns Dict describing regression
    :return alpha: intercept
    :return beta: regression function
    :return fn: aligned functions
    :return qn: aligned srsfs
    :return gamma: warping functions
    :return q: original srsfs
    :retrun B: basis functions
    :return type: type of regressions
    :return b: coefficients
    :return LL: logistic loss
"""
function elastic_logistic(f, y, timet; B=Union{}, df=20, max_itr=20,
                          smooth=false)

    M, N = size(f);

    if M > 500
        parallel = true;
    elseif N > 100
        parallel = true;
    else
        parallel = false;
    end

    binsize = mean(diff(timet));

    # Create B-spline basis if none provided
    if B == Union{}
        B = bs(timet, df, 4);
    end
    Nb = size(B,2);

    q = f_to_srsf(f, timet, smooth);

    gamma = repeat(linspace(0,1,M), 1, N);
    fn = zeros(M, N);
    qn = zeros(M, N);

    itr = 1;
    LL = zeros(max_itr);
    b = zeros(Nb+1);
    alpha = 0.0;
    beta = zeros(length(timet));
    while itr <= max_itr
        @printf("Iteration: %d\n", itr)
        # align data
        fn = zeros(M, N);
        qn = zeros(M, N);
        for ii in 1:N
            xout = (timet[end] - timet[1]) * gamma[:, ii] + timet[1];
            fn[:,ii] = approx(timet, f[:, ii], xout);
            qn[:,ii] = warp_q_gamma(timet, q[:,ii], gamma[:,ii]);
        end

        Phi = ones(N, Nb+1);
        for ii in 1:N
            for jj in 2:Nb+1
                Phi[ii,jj] = trapz(timet, qn[:,ii].*B[:,jj-1]);
            end
        end

        # find alpha and beta using l_bfgs
        b0 = zeros(Nb+1);
        opt = Opt(:LD_LBFGS, Nb+1);
        maxeval!(opt, 250);
        xtol_abs!(opt, 1e-12);
        ftol_abs!(opt, 1e-30);
        min_objective!(opt, (x,grad) -> logit_optm(x,grad,Phi,y));
        optf, b, ret = optimize(opt, b0);

        alpha = b[1];
        beta = B * b[2:Nb+1];

        # compute the logistic loss
        LL[itr] = logit_loss(b, Phi, y);

        # find gamma
        gamma_new = zeros(M, N);
        if parallel
            gamma_new = @distributed (hcat) for i=1:N
                logistic_warp(beta, timet, q[:, i], y[i]);
            end
        else

            for i = 1:N
                gamma_new[:,i] = logistic_warp(beta, timet, q[:, i], y[i]);
            end
        end

        if norm(gamma-gamma_new) < 1e-5
            break
        else
            gamma = copy(gamma_new);
        end

        itr += 1;
    end

    out = Dict("alpha" => alpha, "beta" => beta, "fn" => fn, "qn" => qn,
               "gamma" => gamma, "q" => q, "B" => B, "type" => "logistic",
               "b" => b[2:end], "Loss" => LL[1:itr-1]);
    return out
end


"""
Calculate elastic m-logistic regression from function data f, for response y

    elastic_mlogistic(f, y, timet; B=None, df=20, max_itr=20, smooth=false)
    :param f: array (M,N) of N functions
    :param y: vector (N) of responses
    :param timet: vector (N) describing time samples
    :param B: matrix describing basis functions (M,N) (default=None generates a
              B-spline basis
    :param df: degree of freedom of basis
    :param max_itr: maximum number of iterations
    :param smooth: smooth data

    Returns Dict describing regression
    :return alpha: intercept
    :return beta: regression function
    :return fn: aligned functions
    :return qn: aligned srsfs
    :return gamma: warping functions
    :return q: original srsfs
    :retrun B: basis functions
    :return type: type of regressions
    :return b: coefficients
    :return n_classes: number of classes
    :return LL: logistic loss
"""
function elastic_mlogistic(f, y, timet; B=Union{}, df=20, max_itr=20,
                           delta=.01, smooth=false)

    M, N = size(f);
    # Code labels
    m = Int(maximum(y));
    Y = zeros(Int32, N, m);
    for ii in 1:N
        Y[ii, y[ii]] = 1;
    end

    if M > 500
        parallel = true;
    elseif N > 100
        parallel = true;
    else
        parallel = false;
    end

    binsize = mean(diff(timet));

    # Create B-spline basis if none provided
    if B == Union{}
        B = bs(timet, df, 4);
    end
    Nb = size(B,2);

    q = f_to_srsf(f, timet, smooth);

    gamma = repeat(linspace(0,1,M), 1, N);
    fn = zeros(M, N);
    qn = zeros(M, N);

    itr = 1;
    LL = zeros(max_itr);
    b = zeros(m*(Nb+1));
    alpha = zeros(m);
    beta = zeros(M, m);
    while itr <= max_itr
        @printf("Iteration: %d\n", itr)
        # align data
        fn = zeros(M, N);
        qn = zeros(M, N);
        for ii in 1:N
            xout = (timet[end] - timet[1]) * gamma[:, ii] + timet[1];
            fn[:,ii] = approx(timet, f[:, ii], xout);
            qn[:,ii] = warp_q_gamma(timet, q[:,ii], gamma[:,ii]);
        end

        Phi = ones(N, Nb+1);
        for ii in 1:N
            for jj in 2:Nb+1
                Phi[ii,jj] = trapz(timet, qn[:,ii].*B[:,jj-1]);
            end
        end

        # find alpha and beta using l_bfgs
        b0 = zeros(m * (Nb+1));
        opt = Opt(:LD_LBFGS, m * (Nb+1));
        maxeval!(opt, 250);
        xtol_abs!(opt, 1e-12);
        ftol_abs!(opt, 1e-30);
        min_objective!(opt, (x,grad) -> mlogit_optm(x,grad,Phi,Y));
        optf, b, ret = optimize(opt, b0);

        B0 = reshape(b, Nb+1, m);
        alpha = B0[1,:];
        beta = zeros(M, m);
        for i in 1:m
            beta[:, i] = B * B0[2:Nb+1, i];
        end

        # compute the logistic loss
        LL[itr] = mlogit_loss(b, Phi, Y);

        # find gamma
        gamma_new = zeros(M, N);
        if parallel
            gamma_new = @distributed (hcat) for i=1:N
                mlogit_warp_grad(alpha, beta, timet, q[:, i], Y[i, :],
                                 delta=delta);
            end
        else

            for i = 1:N
                gamma_new[:,i] = mlogit_warp_grad(alpha, beta, timet, q[:, i],
                                                  Y[i, :], delt=delta);
            end
        end

        if norm(gamma-gamma_new) < 1e-5
            break
        else
            gamma = copy(gamma_new);
        end

        itr += 1;
    end

    out = Dict("alpha" => alpha, "beta" => beta, "fn" => fn, "qn" => qn,
               "gamma" => gamma, "q" => q, "B" => B, "type" => "mlogistic",
               "b" => b[2:end], "n_classes" => m, "Loss" => LL[1:itr-1]);
    return out
end


"""
Prediction from elastic regression model

    elastic_prediction(f, timet, model; y=None, smooth=false)
    :param f: functions to predict
    :param timet: vector describing time samples
    :param model: calculated model (regression, musicologist, mlogistic)
    :param y: true response (default = None)
    :param smooth: smooth data (default = false)

    Returns
    :return y_pred: predicted value
    :return y_labels: labels of predicted value
    :return Perf: Performance metric if truth is supplied
"""
function elastic_prediction(f, timet, model::Dict; y=Union{}, smooth=false)
    q = f_to_srsf(f, timet, smooth);
    n = size(q, 2);

    if model["type"] == "linear" || model["type"] == "logistic"
        y_pred = zeros(n);
    elseif model["type"] == "mlogistic"
        m = model["n_classes"];
        y_pred = zeros(n, m);
    end

    for ii in 1:n
        diff = model["q"] - repeat(q[:, ii], 1, size(model["q"], 2));
        dist = sum(abs(diff).^2, 1).^(1/2);
        q_tmp = warp_q_gamma(timet, q[:, ii],
                             model["gamma"][:, argmin(dist)]);
        if model["type"] == "linear"
            y_pred[ii] = model["alpha"] + trapz(timet, q_tmp.*model["beta"]);
        elseif model["type"] == "logistic"
            y_pred[ii] = model["alpha"] + trapz(timet, q_tmp.*model["beta"]);
        elseif model["type"] == "mlogistic"
            for jj in 1:m
                y_pred[ii,jj] = model["alpha"][jj] + trapz(timet, q_tmp.*model["beta"][:, jj]);
            end
        end
    end

    if y==Union{}
        if model["type"] == "linear"
            Perf = Union{};
            y_labels =  Array(Integer, 0);
        elseif model["type"] == "logistic"
            y_pred = phi(y_pred);
            y_labels = ones(n);
            y_labels[y_pred .< 0.5] = -1;
            Perf = Union{};
        elseif model["type"] == "mlogisistic"
            y_pred = phi(reshape(y_pred,n*m,1));
            y_pred = reshpae(y_pred,n,m);
            max_val, y_labels = findmax(y_pred, 2);
            Perf = Union{};
        end
    else
        if model["type"] == "linear"
            Perf = sum((y-y_pred).^2);
            y_labels = Array(Integer, 0);
        elseif model["type"] == "logistic"
            y_pred = phi(y_pred);
            y_labels = ones(n);
            y_labels[y_pred .< 0.5] = -1;
            TP = sum(y[y_labels .== 1] .== 1);
            FP = sum(y[y_labels .== -1] .== 1);
            TN = sum(y[y_labels .== -1] .== -1);
            FN = sum(y[y_labels .== 1] .== -1);
            Perf = (TP+TN)/(TP+FP+FN+TN);
        elseif model["type"] == "mlogistic"
            y_pred = phi(reshape(y_pred,n*m,1));
            y_pred = reshape(y_pred,n,m);
            max_val, y_labels = findmax(y_pred, 2);
            temp = ind2sub((n,m),vec(y_labels));
            y_labels = temp[2];
            Perf = zeros(m);
            cls_set = collect(1:m+1);
            for ii in 1:m
                cls_sub = setdiff(cls_set, ii);
                TP = sum(y[y_labels .== ii] .== ii);
                FP = sum(y[indexin(y_labels, cls_sub).>0] .== ii)
                TN = sum(y[indexin(y_labels,cls_sub).>0] .==
                         y_labels[indexin(y_labels,cls_sub).>0]);
                FN = sum(indexin(y[y_labels.==ii], cls_sub).>0);
                Perf[ii] = (TP+TN)/(TP+FP+FN+TN);
            end
            Perf = sum(y .==y_labels)/length(y_labels);
        end
    end

    return y_pred, y_labels, Perf
end
