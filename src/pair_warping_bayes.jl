function pair_warping_bayes(f1, f2; iter=15000, times=5, powera=1,
                            showplot=true)

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
    L = int(length(q1)/times);
    row = times*[0:L-1]+1;
    p = length(q1);

    if (mod(p,times) != 0)
        error(@sprintf("Number of points on q function = %d is not a mulitple of times = %d", p, times))
    end

    if (scale_data)
        rescale = sqrt(p/sum(q1.^2));
        q1 = rescale.*q1;
        rescale = sqrt(p/sum(q2.^2));
        q2 = rescale.*q2;
    end

    MatchIn2, NDist, q2LL = dp_bayes(q1[row], q1, q2, times, cut);
    match = [MatchIn2, p+1];
    match_collect = Array(Float64,int(iter/thin), int(L+1));
    best_match = copy(match);
    dist = NaN;
    dist_collect = zeros(iter+1);
    f_i = InterpIrregular([row, p+1], float(match), BCnearest, InterpLinear);
    idy = round(f_i[1:p]);
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
    f_i = InterpIrregular([row, p+1], float(best_match), BCnearest, InterpLinear);
    bestidy = f_i[1:p];
    bestidy[bestidy .> p] = p;
    push!(bestidy,p+1);
    burnin = round(0.5*iter/thin);
    LowerP = Array(Float64, L+1);
    UpperP = Array(Float64, L+1);
    MeanP = Array(Float64, L+1);
    for i = 1:L+1
        LowerP[i] = quantile(match_collect[burnin:(iter/thin),i], 0.025);
        UpperP[i] = quantile(match_collect[burnin:(iter/thin),i], 0.975);
        MeanP[i] = mean(match_collect[burnin:(iter/thin),i]);
    end

    f_i = InterpIrregular([row, p+1], float(MeanP), BCnearest, InterpLinear);
    Meanidy = f_i[1:p];
    Meanidy[Meanidy .> p] = p;
    push!(Meanidy,p+1);

    f_i = InterpIrregular([0:p], f2, BCnil, InterpLinear);
    reg_q = f_i[linspace(1,p,times*(p+1)-1)];
    reg_q = reg_q[int((bestidy-1)*times+1)];
    reg_a = f_i[linspace(0,p,times*(p+1)-1)];
    reg_a = reg_a[floor((Meanidy-1)*times+1)];

    if (showplot)
        timet = linspace(0, 1, length(f1));
        f11 = qtocurve(q1, timet);
        f21 = qtocurve(q2, timet);
        ranges = maximum(f11) - minimum(f11);
        curve1 = f11 - mean(f11);
        curve2 = f21 - mean(f21)+1*ranges;
        figure()
        plot(timet,curve1)
        oplot(timet,curve2,"b")
        for n in 1:length(best_match)
            oplot([timet[times*(n-1)+1] timet[best_match[n]]], [curve1[times*(n-1)+1] curve2[best_match[n]]],"r")
        end
        xlabel("t")
        title("Correspondence between 2 functions")

        figure()
        time0 = ([1:(L+1)]-1)./L;
        temp = (best_match-1)/p;
        plot(time0, temp,"b")
        temp = (LowerP-1)./p;
        oplot(time0, temp,"r")
        temp = (UpperP-1)./p;
        oplot(time0, temp,"r")
        temp = (MeanP-1)./p;
        oplot(time0, temp,"k")

        figure()
        plot(timet, f1, "k")
        oplot(timet, f2, "b")
        title("Original Functions")

        figure()
        plot(timet, f1, "k")
        oplot(timet, reg_q, "b")
        title("Registration by Quotient estimate")

        figure()
        plot(timet, f1, "k")
        oplot(timet, reg_a, "b")
        title("Registration by Bayesian estimate")

        figure()
        plot(kappafamily[burnin:iter])
        title("Traceplot of kappa after burn-in")

        figure()
        plot(log_posterior[burnin:iter])
        title("Traceplot of log posterior after burn-in")
    end

    out = Dict("f1" => f1, "f2_q" => reg_q, "gam_q" => (bestidy-1)/p,
               "f2a" => reg_a, "gam_a" => (Meanidy-1)/p,
               "dist_collect" => dist_collect, "best_match" => best_match);
    return out
end
