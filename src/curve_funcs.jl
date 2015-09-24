function resamplecurve(x, N=100)
    n, T = size(x);
    xn = zeros(n, T);

    delta = zeros(T);

    for r in 2:T
        delta[r] = norm(x[:,r]-x[:,r-1]);
    end

    cumdel = cumsum(delta)/sum(delta);
    newdel = linspace(0,1,N);

    for r in 1:n
        s = Spline1D(cumdel,x[r,:]);
        xn[r,:] = evaluate(s,newdel);
    end

    return xn
end


function calculatecentroid(beta)
    n, T = size(beta);
    betadot = zeros(n,T);
    for i = 1:n
        betadot[n,:] = gradient(beta[n,:], 1.0/(T-1));
    end
    normbetadot = zeros(T);
    integrand = zeros(n,T);
    for i = 1:T
        normbetadot[i] = norm(betadot[:,i]);
        integrand[:,i] = beta[:,i]*normbetadot[i];
    end
    scale = trapz(linspace(0,1,T), normbetadot);
    centroid = trapz(linspace(0,1,T), integrand, 2)./scale;

    return centroid
end


function curve_to_q(beta)
    n, T = size(beta);
    v = zeros(n,T);
    for i = 1:n
        v[n,:] = gradient(beta[n,:], 1.0/(T-1));
    end

    len = sum(sqrt(sum(v.*v)))/T;
    v = v./len;
    q = zeros(n,T);
    for i = 1:T
        L = sqrt(norm(v[:,i]));
        if (L > 0.0001)
            q[:,i] = v[:,i]./L;
        else
            q[:,i] = v[:,i]*0.0001;
        end
    end

    q = q./sqrt(innerprod_q2(q, q));

    return q
end


function q_to_curve(q)
    T = size(q,2);
    qnorm = zeros(T);
    for i = 1:T
        qnorm[i] = norm(q[:,i]);
    end

    integrand = zeros(2,T);
    integrand[1,:] = q[1,:].*qnorm;
    integrand[2,:] = q[2,:].*qnorm;
    beta = cumtrapz([1:T], integrand, 2)./T;

    return beta
end


function optimum_reparam(beta1::Array{Float64,2}, beta2::Array{Float64,2},
                         lam::Float64=0.0; method::ASCIIString="DP", w=0.01,
                         rotated::Bool=true, isclosed::Bool=false)
    n1, M = size(beta2);
    timet = linspace(0,1,M);
    skipm = 4;
    auto = 2;
    tau = 0;
    if (method == "DP")
        # Optimze over SO(n) x Gamma
        q1 = curve_to_q(beta1);

        # Optimzie over SO(n)
        beta2, R, tau = find_rotation_seed_coord(beta1, beta2);
        q2 = curve_to_q(beta2);

        # Optimzie over Gamma
        q1i = reshape(q1, M*n1, 1);
        q2i = reshape(q2, M*n1, 1);
        G = zeros(M);
        T = zeros(M);
        sizei = Cdouble[0];
        ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
            (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
            Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Float64), q1i, timet, q2i, timet,
            n1, M, M, timet, timet, M, M, G, T, sizei, lam)
        G = G[1:int(sizei[1])];
        T = T[1:int(sizei[1])];
        yi = InterpIrregular(T, G, BCnil, InterpLinear);
        gam = yi[timet];

    elseif (method == "DP2")
        c1 = reshape(beta1', M*n1, 1);
        c2 = reshape(beta2', M*n1, 1);
        opt = zeros(M+n1*n1+1);
        swap = false;
        fopts = zeros(5);
        comtime = zeros(5);
        @cpp ccall((:optimum_reparam, libgropt), Void,
                   (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                    Ptr{Float64}, Bool, Bool, Bool, Int32, Int32, Ptr{Float64},
                    Ptr{Float64}), c1, c2, M, n1, 0.0, true, rotated, isclosed,
                    skipm, auto, opt, swap, fopts, comtime)

        gam = opt[1:end-5];
        R = reshape(opt[end-4:end-1],2,2);

        if swap
            gam = invertGamma(gam);
            R = R';
        end

    else
        c1 = reshape(beta1', M*n1, 1);
        c2 = reshape(beta2', M*n1, 1);
        opt = zeros(M+n1*n1+1);
        swap = false;
        fopts = zeros(5);
        comtime = zeros(5);
        @cpp ccall((:optimum_reparam, libgropt), Void,
                   (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                    Ptr{Float64}, Bool, Bool, Bool, Int32, Int32, Ptr{Float64},
                    Ptr{Float64}), c1, c2, M, n1, w, false, rotated, isclosed,
                    skipm, auto, opt, swap, fopts, comtime)

        if fopts[1] == 1000
            @cpp ccall((:optimum_reparam, libgropt), Void,
                       (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                        Ptr{Float64}, Bool, Bool, Bool, Int32, Int32,
                        Ptr{Float64}, Ptr{Float64}), c1, c2, M, n1, 0.0, true,
                        rotated, isclosed, skipm, auto, opt, swap, fopts,
                        comtime)
        end

        gam = opt[1:end-5];
        R = reshape(opt[end-4:end-1],2,2);

        if swap
            gam = invertGamma(gam);
            R = R';
        end
    end

    gam = (gam-gam[1]) ./ (gam[end] - gam[1]);

    return gam, R, tau
end


function find_best_rotation(q1, q2)
    epsilon = eps(Float64);
    n, T = size(q1);
    A = q1*q2';
    U, S, V = svd(A);
    if (abs(det(U) * det(V) - 1) < 10 * eps)
        S = eye(n)
    else
        S = eye(n)
        S[:,end] = -S[:,end];
    end
    R = U*S*V';
    q2new = R*q2;

    return q2new, R
end


function calculate_variance(beta)
    n, T = size(beta);
    betadot = zeros(n, T);
    for i = 1:n
        betadot[i,:] = gradient(beta[i,:], 1.0/(T-1));
    end
    normbetadot = zeros(T);
    centroid = calculatecentroid(beta);
    integrand = zeros(n, n, T);
    t = linspace(0,1,T);
    for i = 1:T
        normbetadot[i] = norm(betadot[:,i]);
        a1 = beta[:,i] - centroid;
        integrand[:, :, i] = a1 * a1' * normbetadot[i];
    end
    l = trapz(t, normbetadot);
    variance = trapz(t, integrand, 3);
    variance /= l;

    return variance
end


function psi(x, a, q)
    T = size(q,2);
    covmat = calculate_variance(x+repmat(a,1,T));
    psi1 = covmat[1,1] - covmat[2,2];
    psi2 = covmat[1,2];
    psi3 = x[1,end];
    psi4 = x[2,end];

    return psi1, psi2, psi3, psi4
end


function find_basis_normal(q)
    n, T = size(q);

    f1 = zeros(n,T);
    f2 = zeros(n,T);
    for i = 1:T
        f1[:,i] = q[1,i]*q[:,i]/norm(q[:,i])+[norm(q[:,i]);0];
        f2[:,i] = q[2,i]*q[:,i]/norm(q[:,i])+[0;norm(q[:,i])];
    end
    h3 = copy(f1);
    h4 = copy(f2);
    integrandb3 = zeros(T);
    integrandb4 = zeros(T);
    for i = 1:T
        integrandb3[i] = q[:,i]'*h3[:,i];
        integrandb4[i] = q[:,i]'*h4[:,i];
    end
    b3 = h3 - q*trapz(linspace(0,1,T),integrandb3);
    b4 = h4 - q*trapz(linspace(0,1,T),integrandb4);

    basis = Array(any,2);
    basis[1] = b3;
    basis[2] = b4;

    return basis
end


function calc_j(basis)
    b1 = basis[1];
    b2 = basis[2];
    T = size(b1,2);

    integrand11 = zeros(T);
    integrand12 = zeros(T);
    integrand22 = zeros(T);

    for i = 1:T
        integrand11[i] = b1[:,i]'*b1[:,i];
        integrand12[i] = b1[:,i]'*b2[:,i];
        integrand22[i] = b2[:,i]'*b2[:,i];
    end

    j = zeros(2,2);
    j[1,1] = trapz(linspace(0,1,T), integrand11);
    j[1,2] = trapz(linspace(0,1,T), integrand12);
    j[2,2] = trapz(linspace(0,1,T), integrand22);
    j[2,1] = j[1,2];

    return j
end


function shift_f(f, tau)
    n, T = size(f);
    fn = circshift(f[:,1:T-1],[0 tau]);
    fn[:,T] = fn[:,1];

    return fn
end


function find_rotation_seed_coord(beta1,beta2)
    n,T = size(beta1);
    q1 = curve_to_q(beta1);
    Ltwo = zeros(T);
    Rlist = zeros(n, n, T);
    for ctr = 1:T
        beta2n = shift_f(beta2, ctr);
        beta2new, R = find_best_rotation(beta1, beta2n);
        q2new = curve_to_q(beta2new);
        Ltwo[ctr] = innerprod_q2(q1-q2new, q1-q2new);
        Rlist[:,:,ctr] = R;
    end

    tau = indmin(Ltwo);
    O_hat = Rlist[:,:,tau];
    beta2new = shift_f(beta2, tau);
    beta2new = O_hat * beta2new;

    return beta2new, O_hat, tau
end


function find_rotation_and_seed_q(q1,q2)
    n, T = size(q1);
    Ltwo = zeros(T);
    Rlist = zeros(n, n, T);
    for ctr = 1:T
        q2n = shift_f(q2, ctr);
        q2new, R = find_best_rotation(q1, q2n);
        Ltwo[ctr] = innerprod_q2(q1-q2new, q1-q2new);
        Rlist[:, :, ctr] = R;
    end

    tau = indmin(Ltwo);
    O_hat = Rlist[:,:,tau];
    q2new = shift_f(q2,tau);
    q2new = O_hat * q2new;

    return q2new, O_hat, tau
end


function group_action_by_gamma(q, gamma)
    n, T = size(q);
    gammadot = gradient(gamma, 1.0/T);
    qn = zeros(n, T);

    for j = 1:n
        s = Spline1D(linspace(0, 1, T), q[j, :]);
        qn[j, :] = evaluate(s, gamma) * sqrt(gammadot);
    end

    qn /= sqrt(innerprod_q2(qn,qn));

    return qn
end


function group_action_by_gamma_coord(f, gamma)
    n, T = size(f);
    fn = zeros(n,T);

    for j = 1:n
        s = Spline1D(linspace(0,1,T), f[j, :]);
        fn[j,:] = evaluate(s, gamma);
    end

    return fn
end


function project_curve(q)
    T = size(q,2);
    tol = 1e-5;
    maxit = 200;
    itr = 0;
    delta = 0.5;
    x = q_to_curve(q);
    a = -1 * calculatecentroid(x);

    psi1, psi2, psi3, psi4 = psi(x,a,q);

    r = [psi3,psi4];
    rnorm = zeros(maxit+1);
    rnomr[1] = norm(r);

    while itr <= maxit
        basis = find_basis_normal(q);

        # calculate Jacobian
        J = calc_j(basis);

        # Newton-Raphson step to update q
        y = J\(-r);
        dq = delta * (y[1]*basis[1] + y[2]*basis[2]);
        normdq = sqrt(innerprod_q2(dq,dq));
        q = cos(normdq) * q + sin(normdq) * dq/normdq;
        q /= sqrt(innerprod_q2(q,q));

        # update x and a from the new q
        beta_new = q_to_curve(q);
        x = copy(beta_new);
        a = -1 * calculatecentroid(x);
        beta_new = x + repmat(a,1,T);

        # calculate the new value of psi
        psi1, psi2, psi3, psi4 = psi(x,a,q);
        r = [psi3,psi4];
        rnomr[iter+1] = norm(r);

        if norm(r) < tol
            break
        end

        itr += 1;

    end

    rnorm = rnorm[1:itr-1];

    return q
end


function pre_proc_curve(beta, T=100)
    beta = resamplecurve(beta, T);
    q = curve_to_q(beta);
    qnew = project_curve(q);
    x = q_to_curve(qnew);
    a = -1 * calculatecentroid(x);
    betanew = x + repmat(a,1,T);
    A = eye(2)

    return betanew, qnew, A
end


function inverse_exp_coord(beta1, beta2)
    T = size(beta1,2);
    centroid1 = calculatecentroid(beta1);
    beta1 -= repamat(centroid1,1,T);
    centroid2 = calculatecentroid(beta2);
    beta2 -= repmat(centroid2,1,T);

    q1 = curve_to_q(beta1);

    # Iteratively optimize over SO(n) x Gamma using old DP
    gam, R, tau = optimum_reparam(beta1, beta2);
    beta2 = R * shift_f(beta2, tau);
    gamI = invertGamma(gam);
    beta2 = group_action_by_gamma_coord(beta2, gamI);
    beta2, R, tau = find_rotation_seed_coord(beta1, beta2);
    q2n = curve_to_q(beta2);

    # Compute geodesic distance
    q1dotq2 = innerprod_q2(q1, q2n);
    dist = acos(q1dotq2);

    # Compute shooting vector
    if q1dotq2 > 1
        q1dotq2 = 1.;
    end

    u = q2n - q1dotq2 * q1;
    normu = sqrt(innerprod_q2(u, u));

    if normu > 1e-4
        v = u * acos(q1dotq2)/normu;
    else
        v = zeros(2, T);
    end

    return v, dist
end


function inverse_exp(q1, q2, beta2)
    T = size(q1,2);
    centroid1 = calculatecentroid(beta2);
    beta2 -= repmat(centroid1,1,T);

    # Optimize over SO(n) x Gamma
    beta1 = q_to_curve(q1);
    gam, R, tau = optimum_reparam(beta1, beta2);
    gamI = invertGamma(gam);
    beta2 = R * shift_f(beta2, tau);

    # Applying optimal re-parameterization to the second curve
    beta2 = group_action_by_gamma_coord(beta2, gamI);
    q2 = curve_to_q(beta2);

    # Optimize over SO(n)
    q2, R2, tau = find_rotation_and_seed_q(q1,q2);

    # Compute geodesic distance
    q1dotq2 = innerprod_q2(q1, q2);
    dist = acos(q1dotq2);

    # Compute shooting vector
    if q1dotq2 > 1
        q1dotq2 = 1.;
    end

    u = q2 - q1dotq2 * q1;
    normu = sqrt(innerprod_q2(u, u));

    if normu > 1e-4
        v = u * acos(q1dotq2) / normu;
    else
        v = zeros(2, T);
    end

    return v
end


function gram_schmidt(basis)
    b1 = basis[1];
    b2 = basis[2];

    basis1 = b1 / sqrt(innerprod_q2(b1, b1));
    b2 = b2 - innerprod_q2(basis1,b2)*basis1;
    basis2 = b2 / sqrt(innerprod_q2(b2, b2));

    basis_o = Array(any,2);
    basis_o[1] = basis1;
    basis_o[2] = basis2;

    return basis_o
end


function project_tangent(w, q, basis)
    w = w - innerprod_q2(w, q) * q;
    bo = gram_schmidt(basis);

    wproj = w - innerprod_q2(w, bo[1])*bo[1] - innerprod_q2(w,bo[1])*bo[1];

    return wproj
end


function scale_curve(beta)
    n, T = size(beta);
    normbetadot = zeros(T);
    betadot = zeros(n,T);
    for i = 1:n
        betadot[n,:] = gradient(beta[n,:], 1.0/(T-1));
    end
    for i = 1:T
        normbetadot[i] = norm(betadot[:,i]);
    end
    scale = trapz(linspace(0,1,T), normbetadot);
    beta_scaled = beta / scale;

    return beta_scaled, scale
end


function parallel_translate(w, q1, q2, basis, mode='O')
    wtilde = w - 2*innerprod_q2(w,q2) / innerprod_q2(q1+q2, q1+q2) * (q1+q2);
    l = sqrt(innerprod_q2(wtilde, wtilde));

    if mode == 'C'
        wbar = project_tangent(wtilde, q2, basis);
        normwbar = sqrt(innerprod_q2(wbar, wbar));
        if normwbar > 1e-4
            wbar = wbar * l / normwbar;
        end
    else
        wbar = wtilde;
    end

    return wbar
end


function curve_zero_crossing(Y, beta, bt, y_max, y_min, gmax, gmin)
    T = size(beta,2);
    betanu = q_to_curve(bt);
    max_itr = 100;
    a = zeros(max_itr);
    a[1] = 1.;
    f = zeros(max_itr);
    f[1] = y_max - Y;
    f[2] = y_min - Y;
    mrp = f[1];
    mrn = f[2];
    mrp_ind = 1;  # most recent positive index
    mrn_ind = 2;  # most recent negative index

    for ii = 3:max_itr
        x1 = a[mrp_ind];
        x2 = a[mrn_ind];
        y1 = mrp;
        y2 = mrn;
        a[ii] = (x1 * y1 - x2 * y2) / float(y2-y1);

        beta1, O_hat, tau = find_rotation_seed_coord(betanu, beta);
        gamma = a[ii] * gmax + (1-a[ii]) * gmin;

        beta1 = group_action_by_gamma_coord(beta1, gamma);
        beta1, O_hat1, tau = find_rotation_seed_coord(betanu, beta1);
        q1 = curve_to_q(beta1);
        f[ii] = innerprod_q2(q1, bt) - Y;

        if abs(f[ii]) < 1e-5
            break
        elseif f[ii] > 0
            mrp = f[ii];
            mrp_ind = ii;
        else
            mrn = f[ii];
            mrn_ind = ii;
        end
    end

   gamma = a[ii] * gmax + (1-a[ii]) * gmin;

   beta1, O_hat, tau = find_rotation_seed_coord(betanu, beta1);
   beta1 = group_action_by_gamma_coord(beta1, gamma);
   beta1, O_hat1, tau = find_rotation_seed_coord(betanu, beta1);
   O_hat = O_hat*O_hat1;

  return gamma, O_hat, tau
end


function rot_mat(theta)
    O = [cos(theta) -1*sin(theta); sin(theta) cos(theta)];

    return O
end


function innerprod_q2(q1, q2)
    T = size(q1,2);
    val = sum(sum(q1.*q2))/T;

    return val
end

