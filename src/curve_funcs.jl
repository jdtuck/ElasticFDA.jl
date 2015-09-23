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
    centroid = zeros(n);
    for i = 1:n
        centroid[n] = trapz(linspace(0,1,T), integrand[n,:]);
    end

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
    beta = zeros(2,T);
    for i = 1:2
        beta[i,:] = cumtrapz([1:T], integrand[i,:]) / T;
    end

    return beta
end


function optimum_reparam(q1::Array{Float64,2}, q2::Array{Float64,2},
                         lam::Float64=0.0; method::ASCIIString="DP", w=0.01,
                         beta1::Array{Float64,2}, beta2::Array{Float64,2},
                         rotated::Bool=true, isclosed::Bool=false)
    n1, M = size(q2);
    q1i = reshape(q1, M*n1, 1);
    q2i = reshape(q2, M*n1, 1);
    timet = linspace(0,1,M);
    skipm = 4;
    auto = 2;
    if (method == "DP")

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

    return gam, R
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
    variance = zeors(n, n);
    for i = 1:n
        for j = 1:n
            variance(i,j) = trapz(t, integrand(i,j,:));
        end
    end
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


function innerprod_q2(q1, q2)
    T = size(q1,2);
    val = sum(sum(q1.*q2))/T

    return val
end


