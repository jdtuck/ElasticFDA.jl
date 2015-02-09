function dp_bayes(q1, q1L, q2L, times, cut)

    colnum = length(q1L);
    rownum = div(colnum,times);
    q2LLlen = (colnum-1)*times+1;
    q2LL = zeros(q2LLlen);
    tempspan1 = [1:q2LLlen]-1;
    temp_span2 = float(tempspan1);
    timesf = float(times);
    q2LL_time = (temp_span2*(1/timesf))+1;
    q2L_time1 = [1:colnum];
    q2L_time2 = float(q2L_time1);
    q2LL_i = InterpGrid(q2L, BCnil, InterpLinear);
    q2LL = q2LL_i[q2LL_time];
    ID = zeros(rownum+1, colnum+1);
    S = zeros(ID);
    jend = 0; start = 0; end1 = 0; k = 0; index =0;
    interx = zeros(Int64, times-1);
    intery2 = zeros(times-1);
    intery = zeros(Int64, times-1);
    span1 = [1:(times-1)];
    span2 = float(span1);
    q1x = Array(Float64, times-1);
    q2y = Array(Float64, times-1);

    for i in 1:(rownum-1)
        if (i*cut+1) < colnum
            jend = i*cut+1;
        else
            jend = colnum;
        end

        for j in i:jend-1
            start = max(i-1,j-cut);
            end1 = min(j-1, cut*(i-1));
            n = start-1+[1:(end1-start+1)];
            k = length(n);
            interx = times*(i-1)+span1;
            Energy = Array(Float64, k);

            for m in 1:k
                jf = j;
                divisor = (jf-n[m])/times;
                intery2 = divisor.*span2+n[m];
                for l in 1:(times-1)
                    intery[l] = round(times*intery2[l]);
                    q1x[l] = q1L[interx[l]+1];
                    q2y[l] = q2LL[intery[l]+1];
                end
                Energy[m] = S[i,n[m]+1]+(q1[i]-sqrt(divisor)*q2L[n[m]+1])^2 +
                    (norm(q1x-sqrt(divisor)*q2y))^2;
            end
            min_val, index = findmin(Energy);
            loc = n[index];
            S[i+1,j+1] = min_val;
            ID[i+1,j+1] = loc;
        end
    end

    i = rownum;
    j = colnum;
    start = max(i-1,j-cut);
    end1 = min(j-1,cut*(i-1));
    n = start-1+[1:(end1-start+1)];
    k = length(n);
    interx = times*(i-1)+span1;
    span2 = float(span1);
    Energy = Array(Float64,k);

    for m in 1:k
        jf = j;
        divisor = (jf-n[m])/times;
        intery2 = divisor.*span2+n[m];
        for l in 1:times-1
            intery[l] = round(times*intery2[l]);
            q1x[l] = q1L[interx[l]+1];
            if intery[l]+1 < q2LLlen
                q2y[l] = q2LL[intery[l]+1];
            else
                q2y[l] = 0.0;
            end
        end
        Energy[m] = S[i,n[m]+1]+(q1[i]-sqrt(divisor)*q2L[n[m]+1])^2 +
                    (norm(q1x-sqrt(divisor)*q2y))^2;
    end

    min_val, index = findmin(Energy);
    loc = n[index];
    S[i+1,j+1] = min_val;
    ID[i+1,j+1] = loc;

    path = Array(Int64,rownum);
    count = ID[i+1,j+1];
    oldcount = 0;
    path[i] = count;

    while count>1
        i -= 1;
        oldcount = count;
        count = ID[i+1,oldcount+1];
        path[i] = count;
    end

    MatchIn2 = path+1;
    NDist = S[rownum+1,colnum+1]/colnum;

    return MatchIn2, NDist, q2LL
end


function DP_mean(f, times=5, fig=false)
    cut = 5*times;
    iter = 20;
    timet = linspace(0,1,size(f,1));
    group = [1:size(f, 2)];
    m = length(timet) - 1;
    n = size(f, 2);
    qt_matrix = zeros(m, n);

    for j = 1:n
        qt_matrix[:, j] = curvetoq(f[:, j], timet);
        rescale = sqrt(m/sum(qt_matrix[:, j].^2));
        qt_matrix[:, j] = rescale.*qt_matrix[:, j];
    end
    row = [1:times:m];
    qt_fitted_matrix = zeros(m, n);

    search = Array(Float64, n);
    meanq = mean(qt_matrix, 2);
    for j = 1:n
        search[j] = Enorm(vec(qt_matrix[:, j] - meanq));
    end
    position = indmin(search);

    mu_5 = qt_matrix[:, position];
    mu_curve = zeros(iter, m+1);
    dist_matrix = zeros(iter+1, n);
    for j = 1:n
        dist_matrix[1, j] = (Enorm(vec(mu_5 - qt_matrix[:, j])))^2/m;
    end
    rtmatrix = zeros(m+1, n);
    match_matrix = zeros(length(row)+1, n);

    i = 1;
    diffdist = 1000.;
    while ((i<iter) & (diffdist > 0.001))
        for j = 1:n
            MatchIn2, NDist, q2LL = dp_bayes(mu_5[row], mu_5, qt_matrix[:, j],
                                             times, cut);
            match = [MatchIn2, m+1];
            f_i = InterpIrregular([row, m+1], float(match), BCnearest,
                                  InterpLinear);
            idy = int(f_i[1:m]);
            idy[idy .> m] = m;
            scale1 = sqrt(diff(match)*(1/times));
            scalevec = kron(scale1, ones(times));
            tmpidy = copy(idy) - 1;
            extended_idy = tmpidy[1:m]*times+1;
            q_fitted = scalevec.*q2LL[extended_idy];
            qt_fitted_matrix[:, j] = q_fitted;
            dist_matrix[i+1, j] = NDist;
            rtmatrix[:, j] = [idy, m+1];
            match_matrix[:, j] = match;
        end
        diffdist = abs(sum(dist_matrix[i+1, :]) - sum(dist_matrix[i, :]));
        mu_5 = mean(qt_fitted_matrix, 2);
        rescale = sqrt(m/sum(mu_5.^2));
        mu_5 *= rescale;
        i += 1;
    end

    estimator2 = copy(mu_5);
    invidy, revscalevec = findkarcherinv(match_matrix, times);
    invidy = invidy[1:m];
    invidy[invidy .>= m] = m;
    f_i = InterpIrregular([1:m], vec(mu_5), BCnearest, InterpLinear);
    mu_5 = revscalevec.*f_i[invidy];
    rescale = sqrt(m/sum(mu_5.^2));
    estimator = rescale.*copy(mu_5);
    reg_curve = zeros(m+1, n);
    for j = 1:n
        f_i = InterpIrregular([0:m], f[:, j], BCnearest, InterpLinear);
        tmp = f_i[linspace(0, m, times*(m+1)-1)];
        reg_curve[:, j] = tmp[(rtmatrix[:,j] - 1)*times+1];
    end
    crossmean = mean(reg_curve, 2);
    sumdist = sum(reg_curve, 2);

    if (fig)
        plotl = min(f);
        plotu = max(f);

        figure()
        plot(timet, reg_curve[:, 1])
        for j in 2:n
            oplot(timet, reg_curve[:, n])
        end
        title("registered curves")

        figure()
        plot(timet, crossmean, "r")
        title("Cross sectional mean")
        for j in 1:n
            oplot(timet, reg_curve[:, j])
        end
    end

    out = ["distfamily" => dist_matrix, "match_matrix" => match_matrix,
           "position" => position, "mu_5" => mu_5,
           "rtmatrix" => rtmatrix, "sumdist" => sumdist,
           "qt_fiited" => qt_fitted_matrix, "estimator" => estimator,
           "estimator2" => estimator2, "regcurve" => reg_curve];
    return out
end
