function srsf_align(f, timett; method="mean", smooth=false, sparam=10, lam=0.0,
                    optim="DP")
    ####################################################################
    # This function aligns a collection of functions using the elastic
    # square-root slope (srsf) framework.

    # :param f:darray of shape (M,N) of M functions with N samples
    # :param time: vector of size N describing the sample points
    # :param method: (string) warp calculate Karcher Mean or Median
    # (options = "mean" or "median") (default="mean")
    # :param showplot: Shows plots of results using matplotlib (default = T)
    # :param smoothdata: Smooth the data using a box filter (default = F)
    # :param lam: controls the elasticity (default = 0)
    # :param optim: optimization method to find warping, default is
    #               Dynamic Programming ("DP"). Other options are
    #               Coordiante Descent ("DP2") and Riemanain BFGS
    #               ("LRBFGS")

    # :return fn: aligned functions - array of shape (M,N) of M
    # functions with N samples
    # :return qn: aligned srvfs - similar structure to fn
    # :return q0: original srvf - similar structure to fn
    # :return fmean: function mean or median - vector of length N
    # :return mqn: srvf mean or median - vector of length N
    # :return gam: warping functions - similar structure to fn
    # :return orig_var: Original Variance of Functions
    # :return amp_var: Amplitude Variance
    # :return phase_var: Phase Variance
    ####################################################################

    M, N = size(f);
    if smooth
        smooth_data!(f, sparam);
    end

    if M > 500
        parallel = true;
    elseif N > 100
        parallel = true;
    else
        parallel = false;
    end

    epsilon = eps(Float64);
    f0 = copy(f);

    methods_ = ["mean", "median"];

    method_ = "";
    for ii in methods_
        if method == ii
            method_ = ii;
        end
    end
    if method_ == ""
        warn("improper method, defaulting to mean")
        method = "mean";
    end

    # Compute SRSF
    binsize = mean(diff(timet));
    q = Array(Float64, M, N);
    for ii in 1:N
        q[:, ii] = f_to_srsf(f[:, ii], timet);
    end
    fo = vec(f[1,:]);
    println("Initializing...")
    mnq = mean(q,2);
    d1 = repmat(mnq,1,N);
    d = (q - d1).^2;
    dqq = sqrt(sum(d,1));
    min_ind = indmax(dqq);
    mq = q[:, min_ind];
    mf = f[:, min_ind];

    if parallel
        gam = @parallel (hcat) for i=1:N
            optimum_reparam(mq, timet, q[:, i], lam, method=optim,
                            f1o=mf[1], f2o=fo[i]);
        end
    else
        gam = optimum_reparam(mq, timet, q, lam, method=optim,
                              f1o=mf[1], f2o=fo);
    end

    gamI = sqrt_mean_inverse(gam);
    xout = (timet[end] -timet[1]) .* gamI + timet[1];
    mf = approx(timet, mf, xout);
    mq = f_to_srsf(mf ,timet);

    # Compute Karcher Mean
    if method == "mean"
        @printf("Compute Karcher Mean of %d functions in SRSF space..\n",N)
    elseif method == "median"
        @printf("Compute Karcher Median of %d functions in SRSF space..\n",N)
    end

    MaxItr = 20;
    ds = zeros(MaxItr+2);
    ds[1] = Inf;
    qun = zeros(MaxItr+1);
    tmp = zeros(M, MaxItr+2);
    tmp[:, 1] = mq;
    mq = copy(tmp);
    tmp = zeros(M, MaxItr+2);
    tmp[:, 1] = mf;
    mf = copy(tmp);
    tmp = zeros(M, N, MaxItr+2);
    tmp[:, :, 1] = f;
    f = copy(tmp);
    tmp = zeros(M, N, MaxItr+2);
    tmp[:, :, 1] = q;
    q = copy(tmp);
    r1 = 0;

    for r in 1:MaxItr
        @printf("updating step: r=%d\n", r)
        if r == MaxItr
            println("maximal number of iterations reached")
        end

        # Matching Step
        if parallel
            gam = @parallel (hcat) for i=1:N
                optimum_reparam(mq[:,r], timet, q[:, i, 1], lam, method=optim,
                                f1o=mf[1,r], f2o=fo[i]);
            end
        else
            gam = optimum_reparam(mq[:,r], timet, q[:,:,1], lam, method=optim,
                                  f1o=mf[1,r], f2o=fo);
        end

        gam_dev = zeros(M,N);
        for k in 1:N
            xout = (timet[end]-timet[1]) .* gam[:, k] + timet[1];
            f[:, k, r+1] = approx(timet, f[:, k, 1], xout);
            q[:, k, r+1] = f_to_srsf(f[:, k, r+1], timet);
            gam_dev = gradient(gam[:, k], 1/(M-1));
        end

        mqt = mq[:, r];
        a = repmat(mqt,1,N);
        d = (q[:, :, r+1] - a).^2;
        if method == "mean"
            d1 = sum(trapz(timet, d));
            d2 = sum(trapz(timet, (1-sqrt(gam_dev)).^2));
            ds_tmp = d1 + lam * d2;
            ds[r+1] = ds_tmp;

            # Minimization Step
            # compute the mean of the matched function
            qtemp = q[:, :, r+1];
            mq[:, r+1] = mean(qtemp, 2);
            mf[:, r+1] = mean(f[:, :, r+1], 2);

            qun[r] = norm(mq[:,r+1] - mq[:,r])/norm(mq[:,r]);
        end

        if method == "median"
            d1 = sqrt(sum(trapz(timet, d)));
            d2 = sum(trapz(timet, (1-sqrt(gam_dev)).^2));
            ds_tmp = d1 * lam * d2;
            ds[r+1] = ds_tmp;

            # Minimization Step
            # compute the median of the matched functions
            dist_iinv = ds[r+1].^(-1);
            qtemp = q[:, :, r+1]/ds[r+1];
            mq[:, r+1] = sum(qtemp,2) .* dist_iinv;
            mf[:, r+1] = sum(f[:,:,r+1]/ds[r+1],2) .* dist_iinv;

            qun[r] = norm(mq[:, r+1] - mq[:,r])/ norm(mq[:,r]);
        end

        r1 = r;
        if (qun[r] < 1e-1) | r >= MaxItr
            break
        end
    end

    # Last Step with Centering of gam
    r = r1 + 1;
    if parallel
        gam = @parallel (hcat) for i=1:N
            optimum_reparam(mq[:,r], timet, q[:, i, 1], lam, method=optim,
                            f1o=mf[1,r], f2o=fo[i]);
        end
    else
        gam = optimum_reparam(mq[:,r], timet, q[:,:,1], lam, method=optim,
                              f1o=mf[1,r], f2o=fo);
    end

    gam_dev = zeros(M,N);
    for k in 1:N
        gam_dev = gradient(gam[:, k], 1/(M-1));
    end

    gamI = sqrt_mean_inverse(gam);
    gamI_dev = gradient(gamI, 1/(M-1));
    timet0 = (timet[end]-timet[1]) .* gamI  + timet[1];
    mq[:, r+1] = approx(timet, mq[:, r], timet0) .* sqrt(gamI_dev);

    for k in 1:N
        q[:, k, r+1] = approx(timet, q[:, k, r+1], timet0) .* sqrt(gamI_dev);
        f[:, k, r+1] = approx(timet, f[:, k, r+1], timet0);
        gam[:, k] = approx(timet, gam[:, k], timet0);
    end

    # Aligned data & stats
    fn = f[:, :, r];
    qn = q[:, :, r];
    q0 = q[:, :, 1];
    mean_f0 = mean(f0, 2);
    std_f0 = std(f0, 2);
    mean_fn = mean(fn, 2);
    std_fn = std(fn, 2);
    mqn = mq[:, r];
    fmean = mean(f0[1,:]) + cumtrapz(timet, mqn .* abs(mqn));

    fgam = zeros(M, N);
    for ii in 1:N
        xout = (timet[end] -timet[1]) .* gam[:, ii] + timet[1];
        fgam[:, ii] = approx(timet, fmean, xout);
    end
    var_fgam = var(fgam,2);

    orig_var = trapz(timet, std_f0.^2);
    amp_var = trapz(timet, std_fn.^2);
    phase_var = trapz(timet, var_fgam);

    out = ["fn" => fn, "qn" => qn, "q0" => q0, "fmean" => fmean,
           "mqn" => mqn, "gam" => gam, "orig_var" => orig_var,
           "amp_var" => amp_var, "phase_var" => phase_var];
    return out
end


function align_fPCA(f, timett; num_comp=3, smooth=false, sparam=10)
    ####################################################################
    # aligns a collection of functions while extracting principal components.
    # The functions are aligned to the principal components

    # :param f: array of shape (M,N) of M functions with N samples
    # :param time: vector of size N describing the sample points
    # :param num_comp: number of fPCA components
    # :param showplot: Shows plots of results using matplotlib (default = T)
    # :param smooth_data: Smooth the data using a box filter (default = F)
    # :param sparam: Number of times to run box filter (default = 25)

    # :return fn: aligned functions - numpy ndarray of shape (M,N) of M
    #             functions with N samples
    # :return qn: aligned srvfs - similar structure to fn
    # :return q0: original srvf - similar structure to fn
    # :return mqn: srvf mean or median - vector of length N
    # :return gam: warping functions - similar structure to fn
    # :return q_pca: srsf principal directions
    # :return f_pca: functional principal directions
    # :return latent: latent values
    # :return coef: coefficients
    # :return U: eigenvectors
    # :return orig_var: Original Variance of Functions
    # :return amp_var: Amplitude Variance
    # :return phase_var: Phase Variance
    ####################################################################
    lam = 0.0;
    MaxItr = 50;
    coef = [-2:2];
    Nstd = length(coef);
    M, N = size(f);

    if smooth
        smooth_data!(f, sparam);
    end

    if M > 500
        parallel = true;
    elseif N > 100
        parallel = true;
    else
        parallel = false;
    end

    epsilon = eps(Float64);
    f0 = copy(f);

    # Compute SRSF
    binsize = mean(diff(timet));
    q = Array(Float64, M, N);
    for ii in 1:N
        q[:, ii] = f_to_srsf(f[:, ii], timet);
    end
    println("Initializing...")
    mnq = mean(q,2);
    d1 = repmat(mnq,1,N);
    d = (q - d1).^2;
    dqq = sqrt(sum(d,1));
    min_ind = indmax(dqq);
    @printf("Compute %d functions in SRSF space to %d fPCA components..\n",N,num_comp)

    gamI = sqrt_mean_inverse(gam);
    xout = (timet[end] -timet[1]) .* gamI + timet[1];
    mf = approx(timet, mf, xout);
    mq = f_to_srsf(mf ,timet);

    # Compute Karcher Mean
    itr = 1;
    mq = zeros(M, MaxItr+1);
    mq[:, 1] = q[:, min_ind];
    fi = zeros(M, N, MaxItr+1);
    fi[:, :, 1] = f;
    qi = zeros(M, N, MaxItr+1);
    qi[:, :, 1] = q;
    gam = zeros(M, N, MaxItr+1);
    cost = zeros(MaxItr+1);
    itrf = 0;

    while itr<=MaxItr
        @printf("updating step: r=%d\n", itr)
        if r == MaxItr
            println("maximal number of iterations reached")
        end

        # PCA Step
        a = repmat(mq[:, itr], 1, N);
        qhat_cent = qi[:, :, itr] - a;
        K = cov(qi[:, :, itr], vardim=2);
        U, s, V = svd(K);

        alpha_i = zeros(num_comp, N);
        for ii in 1:num_comp
            for jj in 1:N
                alpha_i[ii,jj] = trapz(timet, qhat_cent[:, jj] .* U[:, ii]);
            end
        end

        U1 = U[:,1:num_comp];
        tmp = U1 * alpha_i;
        qhat = a + tmp;

        # Matching Step
        if parallel
            gam_t = @parallel (hcat) for i=1:N
                optimum_reparam(qhat[:,i], timet, qi[:,i,itr], lam);
            end
            gam[:, :, itr] = gam_t;
        else
            gam[:, :, itr] = optimum_reparam(qhat, timet, qi[:,:,itr], lam);
        end

        for k in 1:N
            xout = (timet[end] -timet[1]) .* gam[:, k] + timet[1];
            fi[:, k, itr+1] = approx(timet, fi[:, k, itr], xout);
            qi[:, k, itr+1] = f_to_srsf(fi[:, k, itr+1], timet);
        end

        qtemp = qi[:, :, itr+1];
        mq[:, itr+1] = mean(qtemp, 2);

        cost_temp = zeros(N);
        for ii in 1:N
            cost_temp[ii] = norm(qtemp[:,ii] - qhat[:,ii]).^2;
        end

        cost[itr+1] = mean(cost_temp);

        if abs(cost[itr+1] - cost[itr] < 1e-6)
            break
        end

        itr += 1;
        itrf = itr;
    end

    if itrf >= MaxItr
        itrf = MaxItr;
    else
        itrf += 1;
    end
    cost = cost[1:itrf];

    # Aligned data & stats
    fn = f[:, :, itrf];
    qn = q[:, :, itrf];
    q0 = q[:, :, 1];
    mean_f0 = mean(f0, 2);
    std_f0 = std(f0, 2);
    mqn = mq[:, r];
    gamf = gam[:, :, 0];
    for k in 1:itrf-1
        gam_k = gam[:, :, k];
        for l in 1:N
            time0 = (timet[end] -timet[1]) .* gam_k[:, l] + timet[1];
            gamf[:,l] = approx(time, gamf[:,l], time0);
        end
    end

    # Center Mean
    gamI = sqrt_mean_inverse(gamf);
    gamI_dev = gradient(gamI, 1/(M-1));
    timet0 = (timet[end] -timet[1]) .* gamI  + timet[1];
    mqn = approx(timet, mqn, timet0) .* sqrt(gamI_dev);
    for k in 1:N
        qn[:, k] = approx(timet, qn[:, k], timet0) .* sqrt(gamI_dev);
        fn[:, k] = approx(timet, fn[:, k], timet0);
        gamf[:, k] = approx(timet, gamf[:, k], timet0);
    end
    mean_fn = mean(fn, 2);
    std_fn = std(fn, 2);

    # Get Final PCA
    out = vert_fPCA(fn, timet, qn, num_comp);

    fmean = mean(f0[1,:]) + cumtrapz(timet, mqn .* abs(mqn));

    fgam = zeros(M, N);
    for ii in 1:N
        xout = (timet[end] -timet[1]) .* gamf[:, ii] + timet[1];
        fgam[:, ii] = approx(timet, fmean, xout);
    end
    var_fgam = var(fgam,2);

    orig_var = trapz(timet, std_f0.^2);
    amp_var = trapz(timet, std_fn.^2);
    phase_var = trapz(timet, var_fgam);

    out2 = ["fn" => fn, "qn" => qn, "q0" => q0, "fmean" => fmean,
           "mqn" => mqn, "gam" => gamf, "q_pca" => out["q_pca"],
           "f_pca" => out["f_pca"], "latent" => out["latent"],
           "coef" => out["coef"], "U" => out["U"], "orig_var" => orig_var,
           "amp_var" => amp_var, "phase_var" => phase_var];
    return out2
end
