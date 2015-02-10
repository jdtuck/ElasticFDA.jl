function Enorm(x::Array{Float64,1})
    n = sqrt(sum(diagm(x'*x)));
    n = n[1];
    return n
end


function Enorm(x::Array{Complex64,1})
    n = sqrt(real(x'*x));
    return n[1]
end


function smooth_data!(f::Array{Float64,2}, sparam=10)
    M, N = size(f);
    for r in 1:sparam
        for i in 1:N
            f[2:(M-1), i] = (f[1:(M-2),i]+2*f[2:(M-1),i] + f[3:M,i])/4;
        end
    end
end


function smooth_data!(f::Array{Float64,1}, sparam=10)
    M = length(f);
    for r in 1:sparam
        f[2:(M-1)] = (f[1:(M-2)]+2*f[2:(M-1)] + f[3:M])/4;
    end
end


function smooth_data(f::Array{Float64,1}, sparam=10)
    M = length(f);
    g = zeros(Float64, M);
    g = copy(f);
    for r in 1:sparam
        g[2:(M-1)] = (g[1:(M-2)]+2*g[2:(M-1)] + g[3:M])/4;
    end

    return g
end


function gradient_spline(timet::Vector, f, smooth=false)
    M = size(f,1);
    spar = 0.0;

    if ndims(f) > 1
        N = size(f,2);
        f0 = zeros(M, N);
        g = zeros(f0);
        g2 = zeros(f0);
        for k in 1:N
            if smooth
                spar = length(timet) * (0.25 * abs(maximum(f[:,k])))^2;
            end
            spl = Spline1D(timet, f[:,k]; s=spar);
            f0[:, k] = evaluate(spl, timet);
            g[:, k] = derivative(spl, timet; nu=1);
            g2[:, k] = derivative(spl, timet; nu=2);
        end
    else
        if smooth
            spar = length(timet) * (0.25 * abs(maximum(f)))^2;
        end
        spl = Spline1D(timet, f; s=spar);
        f0 = evaluate(spl, timet);
        g = derivative(spl, timet; nu=1);
        g2 = derivative(spl, timet; nu=2);
    end

    return f0, g, g2
end


function f_to_srsf(f::Array, timet=0, smooth=false)
    if (timet == 0)
        timet = linspace(0,1,length(f));
    end
    f0, g, g2 = gradient_spline(timet, f, smooth);
    q = g ./ sqrt(abs(g)+eps(Float64));
    return q
end


function optimum_reparam(q1::Array{Float64,1}, timet::Array{Float64,1},
                       q2::Array{Float64,1}, lam::Float64=0.0)
    q1 = q1./norm(q1);
    M = length(q2);
    n1 = 1;
    G = zeros(M);
    T = zeros(M);
    sizei = Cdouble[0];
    ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
           Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Float64), q1, timet, q2, timet,
           n1, M, M, timet, timet, M, M, G, T, sizei, lam)
    G = G[1:int(sizei[1])];
    T = T[1:int(sizei[1])];
    yi = InterpIrregular(T, G, BCnil, InterpLinear);
    gam = yi[timet];
    gam = (gam-gam[1]) ./ (gam[end] - gam[1]);
    return gam
end


function optimum_reparam(q1::Array{Float64,1}, time1::Array{Float64,1},
                         q2::Array{Float64,1}, time2::Array{Float64,1},
                         lam::Float64=0.0)
    q1 = q1./norm(q1);
    M1 = length(q1);
    M2 = length(q2);
    n1 = 1;
    G = zeros(M1);
    T = zeros(M1);
    sizei = Cdouble[0];
    ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
          (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
           Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32, Ptr{Float64},
           Ptr{Float64}, Ptr{Float64}, Float64), q1, time1, q2, time2,
           n1, M1, M2, time1, time2, M1, M2, G, T, sizei, lam)
    G = G[1:int(sizei[1])];
    T = T[1:int(sizei[1])];
    yi = InterpIrregular(T, G, BCnil, InterpLinear);
    gam = yi[time1];
    gam = (gam-gam[1]) ./ (gam[end] - gam[1]);
    return gam
end


function optimum_reparam(q1::Array{Float64,1}, timet::Array{Float64,1},
                       q2::Array{Float64,2}, lam::Float64=0.0)
    q1 = q1./norm(q1);
    M, N = size(q2);
    n1 = 1;
    sizei = Cdouble[0];
    gam = zeros(M, N);
    for ii in 1:N
        G = zeros(M);
        T = zeros(M);
        qi = q2[:, ii];
        ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
              Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32,
              Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), q1, timet,
              qi, timet, n1, M, M, timet, timet, M, M, G, T, sizei, lam)
        G = G[1:int(sizei[1])];
        T = T[1:int(sizei[1])];
        yi = InterpIrregular(T, G, BCnil, InterpLinear);
        gam0 = yi[timet];
        gam[:, ii] = (gam0-gam0[1]) ./ (gam0[end] - gam0[1]);
    end

    return gam
end


function optimum_reparam(q1::Array{Float64,2}, timet::Array{Float64,1},
                       q2::Array{Float64,2}, lam::Float64=0.0)
    M, N = size(q1);
    n1 = 1;
    sizei = Cdouble[0];
    gam = zeros(M, N);
    for ii in 1:N
        q1i = q1[:, ii] ./ norm(q1[:, ii]);
        q2i = q2[:, ii] ./ norm(q2[:, ii]);
        G = zeros(M);
        T = zeros(M);
        ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
              Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32,
              Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), q1i, timet,
              q2i, timet, n1, M, M, timet, timet, M, M, G, T, sizei, lam)
        G = G[1:int(sizei[1])];
        T = T[1:int(sizei[1])];
        yi = InterpIrregular(T, G, BCnil, InterpLinear);
        gam0 = yi[timet];
        gam[:, ii] = (gam0-gam0[1]) ./ (gam0[end] - gam0[1]);
    end

    return gam
end


function warp_q_gamma(time::Vector, q::Vector, gam::Vector)
    M = length(gam);
    gam_dev = gradient(gam, 1/(M-1));
    tmp = InterpIrregular(time, q, BCnil, InterpLinear);
    xout = (time[end] - time[1]) .* gam + time[1];
    q_temp = tmp[xout] .* sqrt(gam_dev);

    return q_temp
end


function rgam(N, sigma, num)
    gam = zeros(N, num);

    TT = N-1;
    timet = linspace(0,1,TT);
    mu = sqrt(ones(N-1)*TT/(N-1));
    omega = (2*pi)/TT;
    for k in 1:num
        alpha_i = randn(1)*sigma;
        v = alpha_i .* ones(TT);
        cnt = 1;
        for l in 2:10
            alpha_i = randn(1)*sigma;
            if (mod(l,2) != 0)
                v = v+alpha_i.*sqrt(2).*cos(cnt*omega*timet);
                cnt += 1;
            elseif (mod(l,2) == 0)
                v = v+alpha_i.*sqrt(2).*sin(cnt*omega*timet);
            end
        end

        v = v - (mu * transpose(v)) * mu/TT;
        vn = norm(v)/sqrt(TT);
        psi = cos(vn)*mu + sin(vn) * v/vn;
        gam[2:end, k] = cumsum(psi.*psi)./TT;
        gam[:, k] = (gam[:,k] - gam[1,k]) ./ (gam[end,k] - gam[1,k]);
    end

    return gam
end


function random_gamma(gam, num)
    mu, gam_mu, psi, vec1 = sqrt_mean(gam);
    K = cov(vec1, vardim=2);

    U, s, V = svd(K);
    n = 5;
    TT = size(vec1, 1) + 1;
    vm = mean(vec1, 2);
    rgam = zeros(TT, num);
    for k in 1:num
        a = randn(n)
        v = zeros(size(vm,1));
        for i in 1:n
            v = v + a[i] * sqrt(s[i]) * U[:, i];
        end

        vn = norm(v) / sqrt(TT);
        psi = cos(vn) .* mu + sin(vn) .* v./vn;
        tmp = zeros(TT);
        tmp[2:TT] = cumsum(psi.*psi)/TT;
        rgam[:, k] = (tmp - tmp[1]) ./ (tmp[end] - tmp[1]);
    end
    return rgam
end


function invert_gamma(gam)
    N = length(gam);
    x = [1:N]/N;
    gamI = approx(gam,x,x);
    gamI = (gamI - gamI[1]) ./ (gamI[end] - gamI[1]);
    return gamI
end


function qtocurve(q, timet=0)
    m = length(q);
    if (timet == 0)
        timet = linspace(0, 1, m+1);
    end
    curve = zeros(m+1);
    for i = 2:(m+1)
        curve[i] = q[i-1]*abs(q[i-1])*(timet[i]-timet[i-1])+curve[i-1]
    end
    return curve
end


function sqrt_mean_inverse(gam)
    eps1 = eps(Float64);
    T1, n = size(gam);
    dt = 1.0/(T1-1);
    psi = Array(Float64,T1-1,n);
    for k = 1:n
        psi[:,k] = sqrt(diff(gam[:,k]) / dt + eps1);
    end

    # Find Direction
    mnpsi = mean(psi, 2);
    d1 = repmat(mnpsi, 1, n);
    d = (psi - d1).^2;
    dqq = sqrt(sum(d,1));
    min_ind = indmin(dqq);
    mu = psi[:, min_ind];
    maxiter = 20;
    tt = 1;
    lvm = zeros(maxiter);
    vec = zeros(T1-1, n);
    for itr in 1:maxiter
        for k in 1:n
            dot = trapz(linspace(0,1,T1-1), mu.*psi[:,k]);
            if dot > 1
                dot = 1;
            elseif dot < (-1)
                dot = -1;
            end
            leng = acos(dot);
            if leng > 0.0001
                vec[:, k] = (leng/sin(leng)) * (psi[:,k] - cos(leng) * mu);
            else
                vec[:, k] = zeros(T1-1);
            end
        end
        vm = mean(vec,2);
        vm1 = vm.*vm;
        lvm[itr] = sqrt(sum(vm1) * dt);
        if lvm[itr] == 0
            break
        end
        mu = cos(tt * lvm[itr]) * mu + (sin(tt * lvm[itr])/lvm[itr]) * vm;
        if (lvm[itr] < 1e-6) | (itr >= maxiter)
            break
        end
    end
    tmp = mu .* mu;
    gam_mu = zeros(T1)
    gam_mu[2:end] = cumsum(tmp)./T1;
    gam_mu = (gam_mu - minimum(gam_mu)) ./ (maximum(gam_mu) - minimum(gam_mu));
    gamI = invert_gamma(gam_mu);
    return gamI
end


function zero_crossing(Y, q, bt, timet, y_max, y_min, gmax, gmin)
    max_itr = 100;
    a = zeros(max_itr);
    a[1] = 1;
    f = zeros(max_itr);
    f[1] = y_max - Y;
    f[2] = y_min - Y;
    mrp = f[1];
    mrn = f[2];
    mrp_ind = 1;  # most recent positive index
    mrn_ind = 2;  # most recent negative index
    out_ii = 1;

    for ii in 3:max_itr
        x1 = a[mrp_ind];
        x2 = a[mrn_ind];
        y1 = mrp;
        y2 = mrn;
        a[ii] = (x1*y2 - x2*y1) / (y2-y1);

        gam_m = a[ii] *gmax + (1-a[ii]) * gmin;
        qtmp = warp_q_gamma(timet, q, gam_m);
        f[ii] = trapz(timet, qtmp.*bt) - Y;

        if abs(f[ii]) < 1e-5
            break
            out_ii = ii;
        elseif f[ii] > 0
            mrp = f[ii];
            mrp_ind = ii;
        else
            mrn = f[ii];
            mrn_ind = ii;
        end
        out_ii = ii;
    end

    gamma = a[out_ii] * gmax + (1-a[out_ii]) * gmin;

    return gamma

end


function sqrt_mean(gam)
    eps1 = eps(Float64);
    TT, n = size(gam);
    psi = Array(Float64,TT-1,n);
    for k = 1:n
        psi[:,k] = sqrt(diff(gam[:,k]) * TT);
    end

    # Find Direction
    mnpsi = mean(psi, 2);
    d1 = repmat(mnpsi, 1, n);
    d = (psi - d1).^2;
    dqq = sqrt(sum(d,1));
    min_ind = indmin(dqq);
    mu = psi[:, min_ind];
    maxiter = 20;
    tt = 1;
    lvm = zeros(maxiter);
    vec = zeros(TT-1, n);
    for itr in 1:maxiter
        for k in 1:n
            dot = trapz(linspace(0,1,TT-1), mu.*psi[:,k]);
            if dot > 1
                dot = 1;
            elseif dot < (-1)
                dot = -1;
            end
            leng = acos(dot);
            if leng > 0.0001
                vec[:, k] = (leng/sin(leng)) * (psi[:,k] - cos(leng) * mu);
            else
                vec[:, k] = zeros(TT-1);
            end
        end
        vm = mean(vec,2);
        vm1 = vm.*vm;
        lvm[itr] = sqrt(sum(vm1) / TT);
        if lvm[itr] == 0
            break
        end
        mu = cos(tt * lvm[itr]) * mu + (sin(tt * lvm[itr])/lvm[itr]) * vm;
        if (lvm[itr] < 1e-6) | (itr >= maxiter)
            break
        end
    end
    tmp = mu .* mu;
    gam_mu = zeros(TT)
    gam_mu[2:end] = cumsum(tmp)./TT;
    gam_mu = (gam_mu - minimum(gam_mu)) ./ (maximum(gam_mu) - minimum(gam_mu));
    return mu, gam_mu, psi, vec
end


function findkarcherinv(warps, times, roundi=false)
    m, n = size(warps);
    psi_m = zeros(m-1, n);
    for j = 1:n
        psi_m[:, j] = sqrt(diff(warps[:, j])./times);
    end
    w = mean(psi_m, 2);
    mupsi = w./sqrt(sum(w.^2/(m-1)));
    v_m = zeros(m-1, n);
    mupsi_update = zeros(m);
    check = 1.;
    while (check > 0.01)
        for i = 1:n
            theta = acos(sum(mupsi.*psi_m[:, i]./(m-1)));
            v_m[:, i] = theta/sin(theta) * (psi_m[:, i]-cos(theta)*mupsi);
        end
        vbar = vec(mean(v_m, 2));
        check = Enorm(vbar)/sqrt(m-1);
        if (check > 0)
            mupsi_update = cos(0.01*check)*mupsi+sin(0.01*check)*vbar/check;
        else
            mupsi_update = cos(0.01*check)*mupsi;
        end
    end
    karcher_s = 1 + [0, cumsum(mupsi_update.^2)*times];
    f_i = InterpIrregular(vec(karcher_s), float([1:times:(m-1)*times+1]),
                          BCnil, InterpLinear);
    if (roundi)
        invidy = [round(f_i[1:((m-1)*times)]), (m-1)*times+1];
    else
        invidy = [f_i[1:((m-1)*times)], (m-1)*times+1];
    end
    revscalevec = sqrt(diff(invidy));

    return invidy, revscalevec
end
