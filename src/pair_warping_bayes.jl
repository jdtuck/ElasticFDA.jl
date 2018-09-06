"""
Compute pair warping between two functions using Bayesian method

    pair_warping_bayes(f1, f2; iter=15000, times=5, powera=1)
    :param f1, f2: vectors describing functions
    :param iter: number of iterations
    :param times: MCMC parameter
    :param powera: MCMC parameter

    Returns Dict containing
    :return f1:
    :return f2_q: srsf registration
    :return gam_q: warping function
    :return f2a: registered f2
    :return gam: warping function
    :return dist_collect: distance
    :return best_match: best match
"""
function pair_warping_bayes(f1, f2; iter=15000, times=5, powera=1)

    tau = ceil(times*.4);

    # Default setting shall work for many situations.
    # If convergence issues arise then adjust proposal variance tau.
    if (times == 2)
        warn("Small times may lead to convergence issues")
    end

    burnin = NaN;
    kappa = 1000;
    thin = 1;
    cut = 5*times;
    alpha = 1;
    beta = 0.001;
    scale_data = true;

    # Convert to SRSF
    q1 = f_to_srsf(f1);
    q2 = f_to_srsf(f2);

    times = 5;
    L = round(Integer,length(q1)/times);
    row = times*collect(0:L-1)+1;
    p = length(q1);

    if (mod(p,times) != 0)
        error(@sprintf("Number of points on q function = %d is not a multiple of times = %d", p, times))
    end

    if (scale_data)
        rescale = sqrt(p/sum(q1.^2));
        q1 = rescale.*q1;
        rescale = sqrt(p/sum(q2.^2));
        q2 = rescale.*q2;
    end

    MatchIn2, NDist, q2LL = dp_bayes(q1[row], q1, q2, times, cut);
    match = [MatchIn2; p+1];
    match_collect = Array(Float64, round(Integer, iter/thin), round(Integer,L+1));
    best_match = copy(match);
    dist = NaN;
    dist_collect = zeros(iter+1);
    f_i = interpolate(([row;p+1],), float(match), Gridded(Linear()))
    idy = round(Integer, f_i[1:p]);
    idy[idy .> p] = p;
    scale1 = sqrt(diff(match)*(1/times));
    scalevec = kron(scale1, ones(times));
    dist = Enorm(q1-scalevec.*q2[idy])^2/p;
    dist_collect[1] = dist;
    dist_min = dist;

    res = simuiter(iter, p, q1, q2, L, tau, times, kappa, alpha, beta,
                   powera, dist, dist_min, best_match, match, thin, cut);

    best_match = res["best_match"];
    match_collect = res["match_collect"];
    dist_min = res["dist_min"];
    log_posterior = res["log_posterior"];
    dist_collect = res["dist_collect"];
    kappafamily = res["kappafamily"];
    f_i = interpolate(([row; p+1],), float(best_match), Gridded(Linear()))
    bestidy = f_i[1:p];
    bestidy[bestidy .> p] = p;
    push!(bestidy,p+1);
    burnin = round(Integer, 0.5*iter/thin);
    LowerP = Array(Float64, L+1);
    UpperP = Array(Float64, L+1);
    MeanP = Array(Float64, L+1);
    idx_sub = round(Integer, iter/thin);
    for i = 1:L+1
        LowerP[i] = quantile(match_collect[burnin:idx_sub,i], 0.025);
        UpperP[i] = quantile(match_collect[burnin:idx_sub,i], 0.975);
        MeanP[i] = mean(match_collect[burnin:idx_sub,i]);
    end

    f_i = interpolate(([row; p+1],), float(MeanP), Gridded(Linear()))
    Meanidy = f_i[1:p];
    Meanidy[Meanidy .> p] = p;
    push!(Meanidy,p+1);

    f_i = interpolate((collect(1:p),), f2, Gridded(Linear()))
    reg_q = f_i[LinRange(1,p,times*(p+1)-1)];
    reg_q = reg_q[round(Integer, (bestidy-1)*times+1)];
    reg_a = f_i[LinRange(1,p,times*(p+1)-1)];
    reg_a = reg_a[floor(Integer, (Meanidy-1)*times+1)];

    out = Dict("f1" => f1, "f2_q" => reg_q, "gam_q" => (bestidy-1)/p,
               "f2a" => reg_a, "gam_a" => (Meanidy-1)/p,
               "dist_collect" => dist_collect, "best_match" => best_match);
    return out
end
