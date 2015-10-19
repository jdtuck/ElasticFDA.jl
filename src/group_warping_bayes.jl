"""
Group alignment of functions using Bayesian method

    group_warping_bayes(f; iter=20000, times=5, powera=1)
    :param f: array (M,N) of N functions
    :param iter: number of MCMC iterations
    :param times: time slicing
    :param powera: MCMC parameter

    Returns Dict containing
    :return f_q: registered srvfs
    :return gam_q: warping function
    :return f_a: registered functions
    :return gam: warping functions
"""
function group_warping_bayes(f; iter=20000, times=5, powera=1)

    tau = ceil(times*.4);
    gp = collect(1:size(f,2));
    # Default setting shall work for many situations.
    # If convergence issues arise then adjust proposal variance tau.
    if (times == 2)
        warn("Small times may lead to convergence issues")
    end

    burnin = NaN;
    kappa = 1000;
    alpha = 1;
    beta = 0.001;
    var_const = 10000;
    thin = 1;
    scale_data = true;
    cut = 5*times;

    m = size(f, 1)-1;
    if mod(m, times) != 0
        error(@sprintf("Number of points on q function = %d is not a mulitple of times = %d", m, times))
    end
    timet = collect(linspace(0, 1, m+1));
    n = size(f, 2)
    qt_matrix = zeros(m, n);
    qt_fitted_matrix = zeros(m, n);

    for j = 1:n
        qt_matrix[:, j] = sign(diff(f[:,j])).*sqrt(abs(diff(f[:,j]))./diff(timet));
        if scale_data
            rescale = sqrt(m/sum(qt_matrix[:, j].^2));
            qt_matrix[:, j] = rescale.*qt_matrix[:, j];
        end
    end

    row = collect(1:times:m);
    L = length(row);

    res_dp = DP_mean(f, times);
    mu_5 = copy(res_dp["estimator2"]);
    mu_5 = vec(mu_5);
    match_matrix = copy(res_dp["match_matrix"]);

    MAP = copy(mu_5);
    best_match_matrix = copy(match_matrix);
    dist_vec = 100.*ones(n);
    best_vec = copy(dist_vec);
    sumdist = zeros(iter);
    kappa_collect = zeros(iter);
    log_collect = zeros(iter);
    logmax = 0;
    mu_prior  = ones(m);
    cov_prior = diagm(var_const.*ones(m));
    mu_q = zeros(round(Integer, iter/thin), m);
    mu_q_standard = zeros(mu_q);

    burnin = round(Integer, 0.5*iter/thin);
    AVG = length(collect(burnin:(round(Integer, iter/thin))));
    res = itermatch(iter, n, m, mu_5, match_matrix, qt_matrix,
                    qt_fitted_matrix, L, tau, times, kappa, alpha, beta,
                    powera, best_vec, dist_vec, best_match_matrix, mu_prior,
                    var_const, sumdist, thin, mu_q, mu_q_standard, logmax,
                    burnin, AVG);

    mu_q_standard = res["mu_q_standard"];
    mu_q = res["mu_q"];
    MAP = res["MAP"];
    best_match_matrix = res["best_match_matrix"];
    kappafamily = res["kappafamily"];
    sumdist = res["sumdist"];
    log_posterior = res["log_posterior"];
    dist = res["dist"]

    mu_est = mean(mu_q_standard[burnin:round(Integer, iter/thin),:], 1);
    rescale = sqrt(m/sum(mu_est.^2));
    mu_est *= rescale;
    mu_est2 = mean(mu_q[burnin:round(Integer, iter/thin),:], 1);
    rescale = sqrt(m/sum(mu_est2.^2));
    mu_est2 *= rescale;

    bayes_warps = res["bayes_warps"];
    gam_q = zeros(m+1, n);
    gam_a = zeros(m+1, n);
    f_a = zeros(m+1, n);
    f_q = zeros(m+1, n);

    for t = 1:n
         f_i = InterpIrregular([row; m+1], best_match_matrix[:, t], BCnearest, InterpLinear);
         gam_q[:, t] = f_i[1:(m+1)];
         f_o = InterpGrid(f[:, t], BCnil, InterpCubic);
         tmp = f_o[linspace(1,m,times*(m+1)-1)];
         f_q[:, t] = tmp[round(Integer, (gam_q[:, t]-1)*times+1)];
         f_i = InterpIrregular([row; m+1], bayes_warps[:, t], BCnearest, InterpLinear);
         gam_a[:, t] = f_i[1:(m+1)];
         f_a[:, t] = tmp[round(Integer, (gam_a[:, t]-1)*times+1)];
    end

    out = Dict("f_q" => f_q, "gam_q" => gam_q, "gam_a" => gam_a, "f_a" => f_a);

end
