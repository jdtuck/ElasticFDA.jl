"""
Calculate normal Energy on an Array

    Enorm(x::Array{Float64,1})
"""
function Enorm(x::Array{Float64,1})
    n = sqrt(sum(diagm(x'*x)));
    n = n[1];
    return n
end

"""
    Enorm(x::Array{Complex64,1})
"""
function Enorm(x::Array{Complex64,1})
    n = sqrt(real(x'*x));
    return n[1]
end


"""
Smooth functional data using box filter in place

    smooth_data!(f::Array{Float64,2}, sparam=10)
"""
function smooth_data!(f::Array{Float64,2}, sparam=10)
    M, N = size(f);
    for r in 1:sparam
        for i in 1:N
            f[2:(M-1), i] = (f[1:(M-2),i]+2*f[2:(M-1),i] + f[3:M,i])/4;
        end
    end
end


"""
    smooth_data!(f::Array{Float64,1}, sparam=10)
    :param sparam: Number of times to run filter (default = 10)
"""
function smooth_data!(f::Array{Float64,1}, sparam=10)
    M = length(f);
    for r in 1:sparam
        f[2:(M-1)] = (f[1:(M-2)]+2*f[2:(M-1)] + f[3:M])/4;
    end
end


"""
Smooth functional data using box filter

    smooth_data(f::Array{Float64,1}, sparam=10)
    :param sparam: Number of times to run filter (default = 10)
"""
function smooth_data(f::Array{Float64,1}, sparam=10)
    M = length(f);
    g = zeros(Float64, M);
    g = copy(f);
    for r in 1:sparam
        g[2:(M-1)] = (g[1:(M-2)]+2*g[2:(M-1)] + g[3:M])/4;
    end

    return g
end


"""
Calculate gradient of function using B-splines

    gradient_spline(timet::Vector, f, smooth=false)
    :param: timet: Vector describing time samples
    :param: f: Vector or Array (M,N) describing functions of M samples
    :param: smooth: smooth data (default = false)
"""
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


"""
Convert function to square-root slope function (srsf)

    f_to_srsf(f::Array, timet=0, smooth=false)
    :param f: array of shape (M,N) describing N functions of M samples
    :param timet: vector describing time samples (default = 0) will generate
                  linearly spaced time vector of length M
    :param smooth: smooth data (default = false)
"""
function f_to_srsf(f::Array, timet=0, smooth=false)
    if (timet == 0)
        timet = collect(linspace(0,1,length(f)));
    end
    f0, g, g2 = gradient_spline(timet, f, smooth);
    q = g ./ sqrt(abs(g)+eps(Float64));
    return q
end


"""
Convert square-root slope function (srsf) to f

    srsf_to_f(q::Array, time, f0=0.0)
    :param q: array of shape (M,N) describing N srsf of M samples
    :param time: vector describing time samples of length M
    :param f0: initial value of f
"""
function srsf_to_f(q::Array, time, f0=0.0)
    M = size(q,1);
    if ndims(q) > 1
        N = size(q,2);
        f = zeros(M,N);
        for i = 1:N
            qnorm = abs(q[:,i]);
            integrand = q[:,i].*qnorm;
            f[:,i] = f0[i] + cumtrapz(time, integrand);
        end
    else
        integrand = q.*abs(q);
        f = f0+cumtrapz(time, integrand);
    end

    return f
end


"""
Calculate elastic distance between two functions

    elastic_distance(f1::Vector, f2::Vector, timet::Vector,
                     method::AbstractString="SIMUL")
    :param f1: vector of function 1 samples
    :param f2: vector of function 2 samples
    :param timet: vector of time samples
    :param method: optimization method to find warping, default is
                   Dynamic Programming ("DP"). Other options are
                   Coordinate Descent ("DP2"), Riemannian BFGS
                   ("RBFGS"), Simultaneous Alignment ("SIMUL")

    :return da: amplitude distance
    :return dp: phase distance
"""
function elastic_distance(f1::Vector, f2::Vector, timet::Vector,
                          method::AbstractString="SIMUL")
    q1 = f_to_srsf(f1, timet);
    q2 = f_to_srsf(f2, timet);

    gam = optimum_reparam(q1, timet, q2, 0.0, method=method);

    q2a = warp_q_gamma(timet, q2, gam);

    da = sqrt(sum(trapz(timet, (q1-q2a).^2)));

    M = length(gam);
    psi = sqrt(diff(gam)*(M-1));
    mu = ones(M-1);
    dp = real(acos(sum(mu.*psi)/(M-1)));

    return da, dp
end


function optimum_reparam(q1::Array{Float64,1}, timet::Array{Float64,1},
                         q2::Array{Float64,1}, lam::Float64=0.0;
                         method::AbstractString="DP", w=0.01, f1o::Float64=0.0,
                         f2o::Float64=0.0)
    q1 = q1./norm(q1);
    q2 = q2./norm(q2);
    c1 = srsf_to_f(q1,timet,f1o);
    c2 = srsf_to_f(q2,timet,f2o);
    M = length(q2);
    rotated = false;
    isclosed = false;
    skipm = 0;
    auto = 0;
    n1 = 1;
    if (method == "DP")
        G = zeros(M);
        T = zeros(M);
        sizei = Cdouble[0];
        ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
            (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
            Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Float64), q1, timet, q2, timet,
            n1, M, M, timet, timet, M, M, G, T, sizei, lam)
        G = G[1:round(Integer,sizei[1])];
        T = T[1:round(Integer,sizei[1])];
        yi = interpolate((T,), G, Gridded(Linear()))
        gam = yi[timet];
    elseif (method == "SIMUL")
        s1,s2,g1,g2,ext1,ext2,mpath = simul_align(c1,c2);
        u = linspace(0,1,length(g1));
        tmin = minimum(timet);
        tmax = maximum(timet);
        timet2 = copy(timet);
        timet2 = (timet2-tmin)/(tmax-tmin);
        gam = simul_gam(collect(u),g1,g2,timet2,s1,s2,timet2);
    elseif (method == "DP2")
        opt = zeros(M+n1*n1+1);
        swap = false;
        fopts = zeros(5);
        comtime = zeros(5);
        @cpp ccall((:optimum_reparam, libgropt), Void,
                   (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                    Bool, Bool, Int32, Int32, Ptr{Float64}, Bool, Ptr{Float64},
                    Ptr{Float64}), c1, c2, M, n1, 0.0, true, rotated, isclosed,
                    skipm, auto, opt, swap, fopts, comtime)

        gam = opt[1:end-2];

        if swap
            gam = invertGamma(gam);
        end

    else
        opt = zeros(M+n1*n1+1);
        swap = false;
        fopts = zeros(5);
        comtime = zeros(5);
        @cpp ccall((:optimum_reparam, libgropt), Void,
                   (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                    Bool, Bool, Int32, Int32, Ptr{Float64}, Bool, Ptr{Float64},
                    Ptr{Float64}), c1, c2, M, n1, w, false, rotated, isclosed,
                    skipm, auto, opt, swap, fopts, comtime)

        if fopts[1] == 1000
            @cpp ccall((:optimum_reparam, libgropt), Void,
                       (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                        Bool, Bool, Int32, Int32, Ptr{Float64}, Bool,
                        Ptr{Float64}, Ptr{Float64}), c1, c2, M, n1, 0.0, true,
                        rotated, isclosed, skipm, auto, opt, swap, fopts,
                        comtime)
        end

        gam = opt[1:end-2];

        if swap
            gam = invertGamma(gam);
        end
    end

    gam = (gam-gam[1]) ./ (gam[end] - gam[1]);

    return gam
end


"""
Calculate optimum parameterization (warping of q2 to q1)

    optimum_reparam(q1, timet, q2, lam=0.0, method="SIMUL", w=0.01, f1o=0.0,
                    f2o=0.0)
    :param q1: array (M,N) or vector (M) describing srsf set 1
    :param timet: vector describing time samples of length M
    :param q2: array (M,N) or vector (M) describing srsf of set 2
    :param lam: control amount of warping (default=0.0)
    :param method: optimization method to find warping, default is
                   Dynamic Programming ("DP"). Other options are
                   Coordinate Descent ("DP2"), Riemannian BFGS
                   ("RBFGS"), Simultaneous Alignment ("SIMUL")
    :param w: Controls RBFGS (default = 0.01)
    :param f1o: initial value of f1, vector or scalar depending on q1, defaults
                to zero
    :param f2o: initial value of f2, vector or scalar depending on q1, defaults
                to zero

    optimum_reparam(q1, time1, q2, time2, lam=0.0, method="DP", w=0.01, f1o=0.0,
                    f2o=0.0)
    same as above, but different timing for q1 and q2

    optimum_reparam(beta1, beta2, lam, method="DP", w=0.01, rotated=true,
                    isclosed=false)
    :param beta1: array (n,T) describing curve 1
    :param beta2: array (n,T) describing curve 2
    :param lam: control amount of warping (default=0.0)
    :param method: optimization method to find warping, default is
                   Dynamic Programming ("DP"). Other options are
                   Coordinate Descent ("DP2"), Riemanain BFGS
                   ("RBFGS")
    :param w: Controls RBFGS (default = 0.01)
    :param rotated: calculate rotation (default = true)
    :param isclosed: closed curve (default = false)

    :return gam: warping function
    :return R: rotation matrix
    :return tau: seed value
"""
function optimum_reparam(q1::Array{Float64,1}, time1::Array{Float64,1},
                         q2::Array{Float64,1}, time2::Array{Float64,1},
                         lam::Float64=0.0; method::AbstractString="DP", w = 0.01,
                         f1o::Float64=0.0, f2o::Float64=0.0)
    q1 = q1./norm(q1);
    q2 = q2./norm(q2);
    c1 = srsf_to_f(q1,time1,f1o);
    c2 = srsf_to_f(q2,time2,f2o);
    M1 = length(q1);
    M2 = length(q2);
    n1 = 1;
    rotated = false;
    isclosed = false;
    skipm = 0;
    auto = 0;
    if (M1 != M2)
        method = "DP";
    end
    if (method == "DP")
        G = zeros(M1);
        T = zeros(M1);
        sizei = Cdouble[0];
        ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
            (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
            Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32, Ptr{Float64},
            Ptr{Float64}, Ptr{Float64}, Float64), q1, time1, q2, time2,
            n1, M1, M2, time1, time2, M1, M2, G, T, sizei, lam)
        G = G[1:round(Integer,sizei[1])];
        T = T[1:round(Integer,sizei[1])];
        yi = interpolate((T,), G, Gridded(Linear()))
        gam = yi[time1];
    elseif (method == "SIMUL")
        s1,s2,g1,g2,ext1,ext2,mpath = simul_align(c1,c2);
        u = linspace(0,1,length(g1));
        tmin = minimum(time1);
        tmax = maximum(time1);
        timet1 = copy(time1);
        timet1 = (timet1-tmin)/(tmax-tmin);
        gam = simul_gam(collect(u),g1,g2,timet1,s1,s2,timet1);
    elseif (method == "DP2")
        opt = zeros(M1+n1*n1+1);
        swap = false;
        fopts = zeros(5);
        comtime = zeros(5);
        @cpp ccall((:optimum_reparam, libgropt), Void,
                   (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                    Bool, Bool, Int32, Int32, Ptr{Float64}, Bool,  Ptr{Float64},
                    Ptr{Float64}), c1, c2, M1, n1, 0.0, true, rotated, isclosed,
                    skipm, auto, opt, swap, fopts, comtime)

        gam = opt[1:end-2];

        if swap
            gam = invertGamma(gam);
        end

    else
        opt = zeros(M1+n1*n1+1);
        swap = false;
        fopts = zeros(5);
        comtime = zeros(5);
        @cpp ccall((:optimum_reparam, libgropt), Void,
                   (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                    Bool, Bool, Int32, Int32, Ptr{Float64}, Bool, Ptr{Float64},
                    Ptr{Float64}), c1, c2, M1, n1, w, false, rotated, isclosed,
                    skipm, auto, opt, swap, fopts, comtime)

        if fopts[1] == 1000
            @cpp ccall((:optimum_reparam, libgropt), Void,
                       (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                        Bool, Bool, Int32, Int32, Ptr{Float64}, Bool,
                        Ptr{Float64}, Ptr{Float64}), c1, c2, M1, n1, 0.0, true,
                        rotated, isclosed, skipm, auto, opt, swap, fopts,
                        comtime)
        end

        gam = opt[1:end-2];

        if swap
            gam = invertGamma(gam);
        end
    end

    gam = (gam-gam[1]) ./ (gam[end] - gam[1]);

    return gam
end


function optimum_reparam(q1::Array{Float64,1}, timet::Array{Float64,1},
                         q2::Array{Float64,2}, lam::Float64=0.0;
                         method::AbstractString="DP", w=0.01, f1o::Float64=0.0,
                         f2o::Array{Float64,1}=zeros(length(q2)))
    q1 = q1./norm(q1);
    c1 = srsf_to_f(q1,timet,f1o);
    M, N = size(q2);
    n1 = 1;
    gam = zeros(M, N);
    rotated = false;
    isclosed = false;
    skipm = 0;
    auto = 0;
    for ii in 1:N
        qi = q2[:, ii];
        qi = qi./norm(qi);
        ci = srsf_to_f(qi,timet,f2o[ii]);
        if (method == "DP")
            G = zeros(M);
            T = zeros(M);
            sizei = Cdouble[0];
            ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
              Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32,
              Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), q1, timet,
              qi, timet, n1, M, M, timet, timet, M, M, G, T, sizei, lam)
            G = G[1:round(Integer,sizei[1])];
            T = T[1:round(Integer,sizei[1])];
            yi = interpolate((T,), G, Gridded(Linear()))
            gam0 = yi[timet];
        elseif (method == "SIMUL")
            s1,s2,g1,g2,ext1,ext2,mpath = simul_align(c1,ci);
            u = linspace(0,1,length(g1));
            tmin = minimum(timet);
            tmax = maximum(timet);
            timet2 = copy(timet);
            timet2 = (timet2-tmin)/(tmax-tmin);
            gam0 = simul_gam(collect(u),g1,g2,timet2,s1,s2,timet2);
        elseif (method == "DP2")
            opt = zeros(M+n1*n1+1);
            swap = false;
            fopts = zeros(5);
            comtime = zeros(5);
            @cpp ccall((:optimum_reparam, libgropt), Void,
                    (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                     Bool, Bool, Int32, Int32, Ptr{Float64}, Bool, Ptr{Float64},
                     Ptr{Float64}), c1, ci, M, n1, 0.0, true, rotated, isclosed,
                     skipm, auto, opt, swap, fopts, comtime)

            gam0 = opt[1:end-2];

            if swap
                gam0 = invertGamma(gam0);
            end

        else
            opt = zeros(M+n1*n1+1);
            swap = false;
            fopts = zeros(5);
            comtime = zeros(5);
            @cpp ccall((:optimum_reparam, libgropt), Void,
                    (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                     Bool, Bool, Int32, Int32, Ptr{Float64}, Bool,
                     Ptr{Float64}, Ptr{Float64}), c1, ci, M, n1, w, false,
                     rotated, isclosed, skipm, auto, opt, swap, fopts, comtime)

            if fopts[1] == 1000
                @cpp ccall((:optimum_reparam, libgropt), Void,
                        (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                         Bool, Bool, Int32, Int32, Ptr{Float64}, Bool,
                         Ptr{Float64}, Ptr{Float64}), c1, ci, M, n1, 0.0, true,
                         rotated, isclosed, skipm, auto, opt, swap, fopts,
                         comtime)
            end

            gam0 = opt[1:end-2];

            if swap
                gam0 = invertGamma(gam0);
            end
        end

        gam[:, ii] = (gam0-gam0[1]) ./ (gam0[end] - gam0[1]);
    end

    return gam
end


function optimum_reparam(q1::Array{Float64,2}, timet::Array{Float64,1},
                         q2::Array{Float64,2}, lam::Float64=0.0;
                         method::AbstractString="DP", w=0.01,
                         f1o::Array{Float64,1}=zeros(length(q1)),
                         f2o::Array{Float64,1}=zeros(length(q2)))
    M, N = size(q1);
    n1 = 1;
    gam = zeros(M, N);
    rotated = false;
    isclosed = false;
    skipm = 0;
    auto = 0;
    for ii in 1:N
        q1i = q1[:, ii] ./ norm(q1[:, ii]);
        q2i = q2[:, ii] ./ norm(q2[:, ii]);
        c1i = srsf_to_f(q1i, timet, f1o[ii]);
        c2i = srsf_to_f(q2i, timet, f2o[ii]);
        if (method == "DP")
            G = zeros(M);
            T = zeros(M);
            sizei = Cdouble[0];
            ccall((:DynamicProgrammingQ2, libfdasrsf), Void,
              (Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Int32,
              Int32, Int32, Ptr{Float64},Ptr{Float64}, Int32, Int32,
              Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Float64), q1i, timet,
              q2i, timet, n1, M, M, timet, timet, M, M, G, T, sizei, lam)
            G = G[1:round(Integer,sizei[1])];
            T = T[1:round(Integer,sizei[1])];
            yi = interpolate((T,), G, Gridded(Linear()))
            gam0 = yi[timet];
            sizei = Cdouble[0];
        elseif (method == "SIMUL")
            s1,s2,g1,g2,ext1,ext2,mpath = simul_align(c1i,c2i);
            u = linspace(0,1,length(g1));
            tmin = minimum(timet);
            tmax = maximum(timet);
            timet2 = copy(timet);
            timet2 = (timet2-tmin)/(tmax-tmin);
            gam0 = simul_gam(collect(u),g1,g2,timet2,s1,s2,timet2);
        elseif (method == "DP2")
            opt = zeros(M+n1*n1+1);
            swap = false;
            fopts = zeros(5);
            comtime = zeros(5);
            @cpp ccall((:optimum_reparam, libgropt), Void,
                    (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                     Bool, Bool, Int32, Int32, Ptr{Float64}, Bool, Ptr{Float64},
                     Ptr{Float64}), c1i, c2i, M, n1, 0.0, true, rotated,
                     isclosed, skipm, auto, opt, swap, fopts, comtime)

            gam0 = opt[1:end-2];

            if swap
                gam0 = invertGamma(gam0);
            end

        else
            opt = zeros(M+n1*n1+1);
            swap = false;
            fopts = zeros(5);
            comtime = zeros(5);
            @cpp ccall((:optimum_reparam, libgropt), Void,
                    (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                     Bool, Bool, Int32, Int32, Ptr{Float64}, Bool, Ptr{Float64},
                     Ptr{Float64}), c1i, c2i, M, n1, w, false, rotated,
                     isclosed, skipm, auto, opt, swap, fopts, comtime)

            if fopts[1] == 1000
                @cpp ccall((:optimum_reparam, libgropt), Void,
                        (Ptr{Float64}, Ptr{Float64}, Int32, Int32, Float64, Bool,
                         Bool, Bool, Int32, Int32, Ptr{Float64}, Bool,
                         Ptr{Float64}, Ptr{Float64}), c11, c2i, M, n1, 0.0,
                         true, rotated, isclosed, skipm, auto, opt, swap, fopts,
                         comtime)
            end

            gam0 = opt[1:end-2];

            if swap
                gam0 = invertGamma(gam0);
            end
        end

        gam[:, ii] = (gam0-gam0[1]) ./ (gam0[end] - gam0[1]);
    end

    return gam
end


"""
Warp srsf by gamma

    warp_q_gamma(time::Vector, q::Vector, gam::Vector)
    :param time: describes time samples
    :param q: describes srsf
    :param gam: describes warping function
"""
function warp_q_gamma(time::Vector, q::Vector, gam::Vector)
    M = length(gam);
    gam_dev = gradient(gam, 1/(M-1));
    tmp = interpolate((time,), q, Gridded(Linear()))
    xout = (time[end] - time[1]) .* gam + time[1];
    q_temp = tmp[xout] .* sqrt(gam_dev);

    return q_temp
end


"""
Warp function by gamma

    warp_f_gamma(time::Vector, f::Vector, gam::Vector)
    :param time: describes time samples
    :param f: describes function
    :param gam: describes warping function
"""
function warp_f_gamma(time::Vector, f::Vector, gam::Vector)
    M = length(gam);
    tmp = interpolate((time,), f, Gridded(Linear()))
    xout = (time[end] - time[1]) .* gam + time[1];
    f_temp = tmp[xout];

    return f_temp
end


"""
Generate random warping functions

    rgam(N, sigma, num)
    :param N: number of time points
    :param sigma: standard deviation across samples
    :param num: number of random warpings
"""
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


"""
Generates random warping functions based on gam

    random_gamma(gam, num)
    :param gam: array (M,N) describing warping functions
    :param num: number of functions to generate
"""
function random_gamma(gam, num)
    mu, gam_mu, psi, vec1 = sqrt_mean(gam);
    K = Base.covm(vec1, mean(vec1,2), 2);

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


"""
Invert warping function

    invert_gamma(gam)
    :param gam: vector describing warping function
"""
function invert_gamma(gam::Vector)
    N = length(gam);
    x = collect(1:N)/N;
    gamI = approx(gam,x,x);
    gamI = (gamI - gamI[1]) ./ (gamI[end] - gamI[1]);
    return gamI
end


"""
Bayesian qtocurve function

    qtocurve(q, timet=0)
"""
function qtocurve(q, timet=0)
    m = length(q);
    if (timet == 0)
        timet = linspace(0, 1, m+1);
    end
    curve = zeros(m+1);
    for i = 2:(m+1)
        curve[i] = q[i-1]*abs(q[i-1])*(timet[i]-timet[i-1])+curve[i-1];
    end

    return curve
end


"""
Calculate sqrt mean inverse of warping function

    sqrt_mean_inverse(gam)
    :param gam: array (M,N) describing warping functions
"""
function sqrt_mean_inverse(gam::Array)
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
            dot1 = trapz(collect(linspace(0,1,T1-1)), mu.*psi[:,k]);
            if dot1 > 1
                dot1 = 1;
            elseif dot1 < (-1)
                dot1 = -1;
            end
            leng = acos(dot1);
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


"""
Calculate zero crossing of optimal warping function

    zero_crossing(Y,q,bt,timet,y_max,y_min,gmax,gmin)
"""
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


"""
Calculate sqrt mean of warping functions

    sqrt_mean(gam)
    :param gam: array (M,N) describing warping functions
"""
function sqrt_mean(gam::Array)
    eps1 = eps(Float64);
    TT, n = size(gam);
    psi = Array(Float64,TT-1,n);
    for k = 1:n
        psi[:,k] = sqrt(diff(gam[:,k]) * TT + eps1);
    end

    # Find Direction
    mnpsi = mean(psi, 2);
    w = mean(psi, 2);
    mu = w./sqrt(sum(w.^2/(TT-1)));
    maxiter = 500;
    tt = 1;
    lvm = zeros(maxiter);
    vec1 = zeros(TT-1, n);
    for itr in 1:maxiter
        for k in 1:n
            dot1 = sum(mu.*psi[:,k]./(TT-1));
            if dot1 > 1
                dot1 = 1;
            elseif dot1 < (-1)
                dot1 = -1;
            end
            leng = dot1
            if leng > 0.0001
                vec1[:, k] = (leng/sin(leng)) * (psi[:,k] - cos(leng) * mu);
            else
                vec1[:, k] = zeros(TT-1);
            end
        end
        vm = mean(vec1,2);
        vm1 = vm.*vm;
        lvm[itr] = Enorm(vm[:])/sqrt(TT);
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
    return mu, gam_mu, psi, vec1
end


"""
Find Karcher inverse of warping functions

    findkarcherinv(warps, times, roundi=false)
"""
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
    karcher_s = 1 + [0; cumsum(mupsi_update.^2)*times]
    f_i = interpolate((vec(karcher_s),), float(collect(1:times:(m-1)*times+1)),
                      Gridded(Linear()))
    if (roundi)
        invidy = [round(f_i[1:((m-1)*times)]); (m-1)*times+1];
    else
        invidy = [f_i[1:((m-1)*times)]; (m-1)*times+1];
    end
    revscalevec = sqrt(diff(invidy));

    return invidy, revscalevec
end


"""
Simultaneous alignment between two functions

    simul_align(f1::Vector, f2::Vector)
"""
function simul_align(f1::Vector,f2::Vector)
    # parameterize by arc-length
    s1 = arclength(f1);
    s2 = arclength(f2);

    len1 = maximum(s1);
    len2 = maximum(s2);

    f1 /= len1;
    s1 /= len1;
    f2 /= len2;
    s2 /= len2;

    # get srvf (should be +/-1)
    q1 = diff(f1)./diff(s1);
    q1[diff(s1).==0] = 0;
    q2 = diff(f2)./diff(s2);
    q2[diff(s2).==0] = 0;

    # get extreme points
    ext1, d1 = extrema_1s(s1,q1);
    ext2, d2 = extrema_1s(s2,q2);

    D, P, mpath = match_ext(s1,ext1,d1,s2,ext2,d2);

    te1 = s1[ext1];
    te2 = s2[ext2];

    g1, g2 = simul_reparam(te1,te2,mpath);

    return s1,s2,g1,g2,ext1,ext2,mpath
end


"""
Calculate arc-length parametrization of function

    arclength(f::Vector)
"""
function arclength(f::Vector)
    t1 = zeros(length(f));
    t1[2:end] = abs(diff(f));
    t1 = cumsum(t1);

    return t1
end


"""
Find location of change of sign of srsf that is arclength parameterized

    extrema_1s(t::Vector, q::Vector)
"""
function extrema_1s(t::Vector, q::Vector)
    q = round(q);

    if q[1] !=0
        d = -q[1];
    else
        d = q[q.!=0];
        d = d[1];
    end

    ext = find(diff(q))+1;

    ext2 = zeros(Integer, length(ext)+2);
    ext2[1] = 1;
    ext2[2:end-1] = round(ext);
    ext2[end] = length(t);

    return ext2, d
end


"""
Find matching between two extremas

    match_ext(t1, ext1, d1, t2, ext2, d2)
"""
function match_ext(t1,ext1::Array{Integer,1},d1,t2,ext2::Array{Integer,1},d2)
    te1 = t1[ext1];
    te2 = t2[ext2];

    # We'll pad each sequence to start on a 'peak' and end on a 'valley'
    pad1 = zeros(Integer, 2);
    pad2 = zeros(Integer, 2);
    if d1==-1
        te1a = zeros(length(te1)+1);
        te1a[2:end] = te1
        te1a[1] = te1a[2];
        te1 = copy(te1a);
        pad1[1] = 1;
    end

    if mod(length(te1),2)==1
        te1a = zeros(length(te1)+1);
        te1a[1:end-1] = te1;
        te1a[end] = te1[end];
        te1 = copy(te1a);
        pad1[2] = 1;
    end

    if d2 == -1
        te2a = zeros(length(te2)+1);
        te2a[2:end] = te2;
        te2a[1] = te2a[2];
        te2 = copy(te2a);
        pad2[1] = 1;
    end

    if mod(length(te2),2) == 1
        te2a = zeros(length(te2)+1);
        te2a[1:end-1] = te2;
        te2a[end] = te2[end];
        te2 = copy(te2a);
        pad2[2] = 1;
    end

    n1 = length(te1);
    n2 = length(te2);

    # initialize weight and path matrices
    D = zeros(n1,n2);
    P = zeros(n1,n2,2);

    for i = 1:n1
        for j = 1:n2
            if mod(i+j,2)==0
                for ib = (i-1):-2:(1+mod(i,2))
                    for jb = (j-1):-2:(1+mod(j,2))
                        icurr = ib+1:2:i;
                        jcurr = jb+1:2:j;
                        W = sqrt(sum(te1[icurr]-te1[icurr-1])) *
                            sqrt(sum(te2[jcurr]-te2[jcurr-1]));
                        Dcurr = D[ib,jb] + W;
                        if Dcurr > D[i,j]
                            D[i,j] = Dcurr;
                            P[i,j,:] = [ib,jb];
                        end
                    end
                end
            end
        end
    end

    D = D[1+pad1[1]:end-pad1[2],1+pad2[1]:end-pad2[2]];
    P = P[1+pad1[1]:end-pad1[2],1+pad2[1]:end-pad2[2],:];
    P[:,:,1] -= pad1[1];
    P[:,:,2] -= pad2[1];

    # Retrieve Best Path
    if pad1[2] == pad2[2]
        mpath = collect(size(D));
    elseif D[end-1,end] > D[end,end-1]
        mpath = collect(size(D)) - [1,0];
    else
        mpath = collect(size(D)) - [0,1];
    end
    mpath = round(Integer, mpath);
    P = round(Integer, P);
    prev_vert = reshape(P[mpath[1],mpath[2],:],1,2);

    mpath = mpath';

    while any(prev_vert.>0)
        mpath = [prev_vert;mpath];
        prev_vert = reshape(P[mpath[1,1],mpath[1,2],:],1,2);
    end

    return D, P, mpath
end


"""
Find simultaneous re-parametrization

    simul_reparam(te1, te2, mpath)
"""
function simul_reparam(te1, te2, mpath)
    g1 = [0.];
    g2 = [0.];

    if mpath[1,2] == 2
        g1 = [g1;0.];
        g2 = [g2; te2[2]];
    elseif mpath[1,1] == 2
        g1 = [g1;te1[2]];
        g2 = [g2; 0.];
    end

    m = size(mpath,1);
    for i = 1:m-1
        gg1, gg2 = simul_reparam_segment(vec(mpath[i,:]),vec(mpath[i+1,:]),
                                         te1,te2);

        g1 = [g1;gg1];
        g2 = [g2;gg2];
    end

    n1 = length(te1);
    n2 = length(te2);
    if (mpath[end,1] == n1-1) || (mpath[end,2] == n2-1)
        g1 = [g1; 1.];
        g2 = [g2; 1.];
    end

    return g1, g2
end


function simul_reparam_segment(src, tgt, te1, te2)
    i1 = src[1]+1:2:tgt[1];
    i2 = src[2]+1:2:tgt[2];

    v1 = sum(te1[i1]-te1[i1-1]);
    v2 = sum(te2[i2]-te2[i2-1]);
    R = v2/v1;

    a1 = src[1];
    a2 = src[2];
    t1 = te1[a1];
    t2 = te2[a2];
    u1 = 0.;
    u2 = 0.;

    gg1 = Array(Float64,0);
    gg2 = Array(Float64,0);

    while (a1<tgt[1]) && (a2<tgt[2])
        if a1==tgt[1]-1 && a2==tgt[2]-1
            a1 = tgt[1];
            a2 = tgt[2];
            gg1 = [gg1; te1[a1]];
            gg2 = [gg2; te2[a2]];
        else
            p1 = (u1+te1[a1+1]-t1)/v1;
            p2 = (u2+te2[a2+1]-t2)/v2;
            if p1<p2
                lam = t2 + R*(te1[a1+1]-t1);
                gg1 = [gg1; te1[a1+1]; te1[a1+2]];
                gg2 = [gg2; lam; lam];
                u1 = u1 + te1[a1+1] - t1;
                u2 = u2 + lam - t2;
                t1 = te1[a1+2];
                t2 = copy(lam);
                a1 += 2;
            else
                lam = t1 + (1./R)*(te2[a2+1]-t2);
                gg1 = [gg1; lam; lam];
                gg2 = [gg2; te2[a2+1]; te2[a2+2]];
                u1 = u1 + lam - t1;
                u2 = u2 + te2[a2+1] - t2;
                t1 = copy(lam);
                t2 = te2[a2+2];
                a2 += 2;
            end
        end
    end

    return gg1, gg2
end


"""
Calculate warping from q2 to q2 from simultaneous warping

    simul_gam(u, g1,g2,t,s1,s2,tt)
"""
function simul_gam(u::Array{Float64,1},g1,g2,t::Array{Float64,1},s1,s2,
                   tt::Array{Float64,1})
    gs1 = interp1_flat(u,g1,tt);
    gs2 = interp1_flat(u,g2,tt);

    gt1 = interp1_flat(s1,t,gs1);
    gt2 = interp1_flat(s2,t,gs2);

    gam = interp1_flat(gt1,gt2,tt);

    return gam
end


function old_dp(q1::Vector, q2::Vector,lam)
    N = length(q1);
    M = 5*N;
    spl = Spline1D(collect(1.:N)/N, q2);
    q2L = spl(collect(1.:M)/M);

    Nbrs = [1 1; 1 2; 2 1; 2 3; 3 2; 1 3; 3 1; 1 4; 3 4; 4 3; 4 1; 1 5; 2 5; 3 5; 4 5; 5 4; 5 3; 5 2; 5 1];

    E = zeros(N,N);
    E[1,:] = Inf;
    E[:,1] = Inf;
    E[1,1] = 0.0;
    Path = zeros(N,N,2);
    for i = 2:N
        for j = 2:N
            CandE = 100000*ones(size(Nbrs,1));
            for Num = 1:size(Nbrs,1)
                k = i - Nbrs[Num,1];
                l = j - Nbrs[Num,2];
                if (k>0 && l>0)
                    CandE[Num] = E[k,l] + cost_fn(q1,q2,q2L,k,l,i,j,N,lam);
                end
                E[i,j], idx = findmin(CandE);
                Path[i,j,1] = i - Nbrs[idx,1];
                Path[i,j,2] = j - Nbrs[idx,2];
            end
        end
    end

    x = zeros(N);
    y = zeros(N);
    x[1] = N;
    y[1] = N;
    cnt = 1;
    while (x[cnt]>1)
        yi = round(Integer,y[cnt]);
        xi = round(Integer,x[cnt]);
        y[cnt+1] = Path[yi,xi,1];
        x[cnt+1] = Path[yi,xi,2];
        cnt += 1;
    end
    x = x[1:cnt];
    y = y[1:cnt];
    idx = sortperm(x);
    x = x[idx];
    y = y[idx];
    yy = zeros(N);
    for i = 1:N
        F = abs(i-x);
        idx = indmin(F);
        if x[idx] == i
            yy[i] = y[idx];
        else
            if x[idx] > i
                a = x[idx] - i;
                b = i - x[idx-1];
                yy[i] = (a*y[idx-1] + b*y[idx])/(a+b);
            else
                a = i - x[idx];
                b = x[idx+1] - i;
                yy[i] = (a*y[idx+1] + b*y[idx])/(a+b);
            end
        end
    end
    gam = yy/N;

    return gam
end


function cost_fn(q1,q2,q2L,k,l,i,j,N,lam)
    M = length(q2L);
    x = collect(k:1:i);
    m = (j-l)/(i-k);
    y = (x-k)*m + l;
    idx = round(Integer, y*M/N);
    vec = sqrt(m)*q2L[idx];
    E = norm(q1[x]-vec)^2/N;

    return E
end
