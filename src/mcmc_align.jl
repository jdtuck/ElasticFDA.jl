"""
MCMC iteration for pairwise alignment
"""
function simuiter(iter, p, q1, q2, L, tau, times, kappa, alpha, beta,
                  powera, dist, dist_min, best_match, match, thin, cut)

    match -= 1;
    best_match -= 1;

    Ixout = Array(Int64, 2*times);
    Ioriginal = Array(Int64, L+1);
    increment = 0.0; n_dist = 0.0; o_dist = 0.0; adjustcon = 0.0;
    ration = 0.0; prob = 0.0; u = 0.0; logpost = 0.0;
    increment_int = 0; newj = 0; tempnew = 0; tempold = 0; tempx = 0;
    newmatchvec = Array(Float64, 3);
    oldmatchvec = Array(Float64, 3);
    idenmatchvec = Array(Float64, 3);
    pnewvec = Array(Float64, 2);
    poldvec = Array(Float64, 2);
    internew = Array(Float64, 2*times);
    interold = Array(Float64, 2*times);
    qt1_5_interx = Array(Float64, 2*times);
    qt2_5_internew = Array(Float64, 2*times);
    qt2_5_interold = Array(Float64, 2*times);
    diff_ynew = Array(Float64, 2*times);
    diff_yold = Array(Float64, 2*times);
    original = Array(Float64, L+1);
    idy = Array(Float64, p);
    scalevec = Array(Float64, p);
    qt_5_fitted = Array(Float64, p);
    kappa_collect = Array(Float64, iter);
    log_collect = Array(Float64, iter);
    dist_collect = Array(Float64, iter);
    Irow = [1:p]-1;
    row = float(Irow);
    scale1 = Array(Float64, L);
    match_collect = Array(Float64, int(iter/thin), L+1);

    q2LLlen = (p-1)*times+1;
    q2LL = Array(Float64, q2LLlen);
    tempspan1 = [1:q2LLlen]-1;
    temp_span2 = float(tempspan1);
    timesf = float(times);
    q2LL_time = (temp_span2*(1/timesf));
    q2L_time1 = [1:p]-1;
    q2L_time2 = float(q2L_time1);
    q2LL = approx(q2L_time2, q2, q2LL_time);

    bar = Progress(iter, 1)   # minumum update interval: 1 second
    for j = 1:iter
        for i = 2:L
            if ((match[i+1]-match[i-1])>2)
                increment = rand(Normal(0, tau));
                increment_int = int(round(increment));
                if (increment_int == 0)
                    increment_int = (increment>0)?(1):(-1);
                end
                newj = match[i] + increment_int;
                if (newj < match[i+1] && newj > match[i-1])
                    newmatchvec[1] = match[i-1];
                    newmatchvec[2] = newj;
                    newmatchvec[3] = match[i+1];
                    oldmatchvec[1] = match[i-1];
                    oldmatchvec[2] = match[i];
                    oldmatchvec[3] = match[i+1];
                    idenmatchvec[1] = times*(i-2);
                    idenmatchvec[2] = times*(i-1);
                    idenmatchvec[3] = times*i;
                    Ixout = [1:2*times]+times*(i-2)-1;
                    xout = float(Ixout);
                    internew = approx(idenmatchvec, newmatchvec, xout);
                    interold = approx(idenmatchvec, oldmatchvec, xout);
                    interx = copy(xout);
                    interynew = copy(internew);
                    push!(interynew, match[i+1]);
                    interyold = copy(interold);
                    push!(interyold, match[i+1]);
                    diff_ynew = interynew[2:2*times+1]-interynew[1:2*times];
                    diff_yold = interyold[2:2*times+1]-interyold[1:2*times];
                    for ll = 1:2*times
                        tempx = round(interx[ll]);
                        qt1_5_interx[ll] = q1[tempx+1];
                        internew[ll] = (internew[ll]>(p-1))?(p-1):(internew[ll]);
                        tempnew = round(times*internew[ll]+1);
                        qt2_5_internew[ll] = q2LL[tempnew]*sqrt(diff_ynew[ll]);
                        interold[ll] = (interold[ll]>(p-1))?(p-1):(interold[ll]);
                        tempold = round(times*interold[ll]+1);
                        qt2_5_interold[ll] = q2LL[tempold]*sqrt(diff_yold[ll]);
                    end
                    n_dist = norm(qt1_5_interx - qt2_5_internew)^2/p;
                    o_dist = norm(qt1_5_interx - qt2_5_interold)^2/p;
                    pnewvec = (newmatchvec[2:3] - newmatchvec[1:2])/p;
                    poldvec = (oldmatchvec[2:3] - oldmatchvec[1:2])/p;
                    adjustcon = exp((powera-1)*(sum(log(pnewvec))-
                                    sum(log(poldvec))));
                    ratio = adjustcon*exp(kappa*o_dist-kappa*n_dist);
                    prob = (ratio < 1)?(ratio):(1);
                    u = rand(Uniform(0,1));
                    match[i] = (u<prob)?(newj):(match[i]);
                end
            end
        end

        Ioriginal = ([1:L+1]-1)*times;
        original = float(Ioriginal);
        idy = round(approx(original, match, row));
        for ii = 1:L
            scale1[ii] = sqrt((match[ii+1]-match[ii])/times);
        end
        for kk = 1:p
            idy[kk] = (idy[kk]<p)?(idy[kk]):(p-1);
            scalevec[kk] = scale1[ceil(kk/times)];
            qt_5_fitted[kk] = scalevec[kk]*q2[idy[kk]+1];
        end
        dist = norm(q1-qt_5_fitted)^2/p;
        dist_collect[j] = dist;
        if (dist < dist_min)
            best_match = copy(match);
            dist_min = dist;
        end
        if (mod(j,thin)==0)
            match_collect[ceil(j/thin),:] = match;
        end
        kappa = rand(Gamma(p/2+alpha, 1/(dist+beta)));
        kappa_collect[j] = kappa;
        logpost = (p/2+alpha)*log(kappa)-kappa*(beta+dist);
        log_collect[j] = logpost;
        next!(bar);
    end
    match_collect += 1;
    best_match += 1;

    out = Dict("best_match" => best_match, "match_collect" => match_collect,
           "dist_collect" => dist_collect, "kappafamily" => kappa_collect,
           "log_posterior" => log_collect, "dist_min" => dist_min);
    return out
end


"""
MCMC iteration for group alignment
"""
function itermatch(iter, n, m, mu_5, match_matrix, qt_matrix,
                   qt_fitted_matrix, L, tau, times, kappa, alpha, beta,
                   powera, best_vec, dist_vec, best_match_matrix, mu_prior,
                   var_const, sumdist, thin, mu_q, mu_q_standard, logmax,
                   burnin, AVG)

    burnin -= 1;
    match_matrix -= 1;
    best_match_matrix -= 1;

    scale1 = Array(Float64, L);
    increment = 0.0; adjustcon = 0.0; ratio = 0.0; prob = 0.0; u = 0.0;
    n_dist = 0.0; o_dist = 0.0; dist = 0.0; logpost = 0.0; SigmaVar = 0.0;
    rescale = 0.0;
    increment_int = 0; newj = 0; tempnew = 0; tempold = 0; tempx = 0;
    newmatchvec = Array(Float64, 3);
    oldmatchvec = Array(Float64, 3);
    idenmatchvec = Array(Float64, 3);
    qt2_5 = Array(Float64, m);
    match = Array(Float64, L+1);
    interynew = Array(Float64, 2*times);
    interyold = Array(Float64, 2*times);
    interx = Array(Float64, 2*times);
    scalevec = Array(Float64, m);
    qt_5_fitted = Array(Float64, m);
    Ixout = Array(Int64, 2*times);
    Ioriginal = Array(Int64, L+1);
    xout = Array(Float64, 2*times);
    internew = Array(Float64, 2*times);
    interold = Array(Float64, 2*times);
    mu_5_interx = Array(Float64, 2*times);
    qt2_5_internew = Array(Float64, 2*times);
    qt2_5_interold = Array(Float64, 2*times);
    diff_ynew = Array(Float64, 2*times);
    diff_yold = Array(Float64, 2*times);
    pnewvec = Array(Float64, 2);
    poldvec = Array(Float64, 2);
    original = Array(Float64, L+1);
    idy = Array(Float64, m);
    Irow = [1:m]-1;
    row = float(Irow);
    karcher_res = Array(Float64, m+1);
    tempkarcher_res = Array(Float64, m+1);
    revscalevec = Array(Float64, L*times);
    kappa_collect = Array(Float64, iter);
    log_collect = Array(Float64, iter);
    MAP = Array(Float64, m);
    Mean = Array(Float64, m);
    qt2_5_warp = Array(Float64, m);
    mu_5_warp = Array(Float64, m);
    tempmatch_matrix = Array(Float64, L+1, n);
    bayes_warps = zeros(L+1, n);

    q2LLlen = (m-1)*times+1;
    q2LL_mat = Array(Float64, q2LLlen, n);
    tempspan1 = [1:q2LLlen]-1;
    temp_span2 = float(tempspan1);
    timesf = float(times);
    q2LL_time = (temp_span2*(1/timesf));
    q2L_time1 = [1:m]-1;
    q2L_time2 = float(q2L_time1);

    for LLL in 1:n
        q2LL_mat[:, LLL] = approx(q2L_time2, qt_matrix[:, LLL], q2LL_time);
    end

    bar = Progress(iter, 1);
    for j = 1:iter
        for t = 1:n
            qt2_5 = qt_matrix[:, t];
            match = match_matrix[:, t];
            for i = 2:L
                if ((match[i+1]-match[i-1])>2)
                    increment = rand(Normal(0, tau));
                    increment_int = int(round(increment));
                    if (increment_int == 0)
                        break
                        # increment_int = (increment>0)?(1):(-1);
                    end
                    newj = match[i] + increment_int;
                    if (newj < match[i+1] && newj > match[i-1])
                        newmatchvec[1] = match[i-1];
                        newmatchvec[2] = newj;
                        newmatchvec[3] = match[i+1];
                        oldmatchvec[1] = match[i-1];
                        oldmatchvec[2] = match[i];
                        oldmatchvec[3] = match[i+1];
                        idenmatchvec[1] = times*(i-2);
                        idenmatchvec[2] = times*(i-1);
                        idenmatchvec[3] = times*i;
                        Ixout = [1:2*times]+times*(i-2)-1;
                        xout = float(Ixout);
                        internew = approx(idenmatchvec, newmatchvec, xout);
                        interold = approx(idenmatchvec, oldmatchvec, xout);
                        interx = copy(xout);
                        interynew = copy(internew);
                        push!(interynew, match[i+1]);
                        interyold = copy(interold);
                        push!(interyold, match[i+1]);
                        diff_ynew = interynew[2:(2*times+1)]-interynew[1:2*times];
                        diff_yold = interyold[2:(2*times+1)]-interyold[1:2*times];
                        for ll = 1:2*times
                            tempx = round(interx[ll]);
                            mu_5_interx[ll] = mu_5[tempx+1];
                            internew[ll] = (internew[ll]>(m-1))?(m-1):(internew[ll]);
                            tempnew = round(times*internew[ll]+1);
                            tempnew = (tempnew > q2LLlen)?(q2LLlen):(tempnew);
                            qt2_5_internew[ll] = q2LL_mat[tempnew, t] *
                                sqrt(diff_ynew[ll]);
                            interold[ll] = (interold[ll]>(m-1))?(m-1):(interold[ll]);
                            tempold = round(times*interold[ll]+1);
                            tempold = (tempold > q2LLlen)?(q2LLlen):(tempold);
                            qt2_5_interold[ll] = q2LL_mat[tempold, t] *
                                sqrt(diff_yold[ll]);
                        end
                        n_dist = norm(mu_5_interx - qt2_5_internew)^2/m;
                        o_dist = norm(mu_5_interx - qt2_5_interold)^2/m;
                        pnewvec = (newmatchvec[2:3] - newmatchvec[1:2])/m;
                        poldvec = (oldmatchvec[2:3] - oldmatchvec[1:2])/m;
                        adjustcon = exp((powera-1)*(sum(log(pnewvec))-
                                        sum(log(poldvec))));
                        ratio = adjustcon*exp(kappa*o_dist-kappa*n_dist);
                        prob = (ratio < 1)?(ratio):(1);
                        u = rand(Uniform(0,1));
                        match[i] = (u<prob)?(newj):(match[i]);
                    end
                end
            end
            Ioriginal = ([1:L+1]-1)*times;
            original = float(Ioriginal);
            idy = approx(original, match, row);
            qt2_5_warp = approx(q2L_time2, qt2_5, idy);

            for ii = 1:L
                scale1[ii] = sqrt((match[ii+1]-match[ii])/times);
            end

            for kk = 1:m
                scalevec[kk] = scale1[ceil(kk/times)];
                qt_5_fitted[kk] = scalevec[kk]*qt2_5_warp[kk];
            end

            dist = norm(mu_5-qt_5_fitted)^2/m;
            match_matrix[:, t] = match;
            dist_vec[t] = dist;
            qt_fitted_matrix[:, t] = qt_5_fitted;
            if ((j/thin)>=burnin)
                bayes_warps[:, t] = bayes_warps[:, t] + match./AVG;
            end
        end
        sumdist[j] = sum(dist_vec);
        if (sumdist[j] < sum(best_vec))
            best_vec = copy(dist_vec);
            best_match_matrix = copy(match_matrix);
        end
        if (mod(j,thin)==0)
            mu_q[ceil(j/thin),:] = mu_5;
            karcher_res, revscalevec = findkarcherinv(match_matrix+1, times);
            # karcher_res = karcher_res[1:m];
            # karcher_res[karcher_res .>= m] = m;
            mu_5_warp = approx(q2L_time2, vec(mu_5), karcher_res-1);
            for qq in 1:m
                mu_q_standard[ceil(j/thin), qq] = revscalevec[qq]*mu_5_warp[qq];
            end
        end
        kappa = rand(Gamma(n*m/2+alpha, 1/(sum(dist_vec)+beta)));
        kappa_collect[j] = kappa;
        logpost = (n*m/2+alpha)*log(kappa)-kappa*(beta+sum(dist_vec));
        log_collect[j] = logpost;
        if (logpost > logmax)
            logmax = logpost;
            MAP = copy(mu_5);
        end
        Mean = (var_const*m)/(m+2*kappa*n*var_const)*(2*kappa*n*mean(qt_fitted_matrix,2)/m + mu_prior/var_const);
        SigmaVar = (var_const*m)/(m+2*kappa*n*var_const);
        mu_5 = Mean + rand(Normal(0, sqrt(SigmaVar)), m);
        rescale = sqrt(m/sum(mu_5.^2));
        mu_5 *= rescale;
        next!(bar);
    end

    best_match_matrix += 1;
    bayes_warps += 1;

    out = Dict("mu_q_standard" => mu_q_standard, "mu_q" => mu_q,
               "best_match_matrix" => best_match_matrix,
               "sumdist" => sumdist, "kappafamily" => kappa_collect,
               "log_posterior" => log_collect, "dist" => sum(best_vec),
               "MAP" => MAP, "bayes_warps" => bayes_warps);
    return out
end
