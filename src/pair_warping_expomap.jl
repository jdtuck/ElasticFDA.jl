"""
Compute pair warping between two functions using Bayesian method

    pair_warping_expomap(f1, f2, timet; iter=20000, burnin=min(5000,iter/2),
                         alpha0=0.1, beta0=0.1, pbetas=[0.5,0.05,0.005,0.0001],
                         probs=[0.1,0.1,0.7,0.1],propvar=1.0,
                         init_coef=zeros(20),npoints=200,extrainfo=false)

    :param f1, f2: vectors describing functions
    :param timet: vector describing timing
    :param iter: number of MCMC iterations
    :param burnin: number of MCMC burnin iterations
    :param alpha0, beta0: IG parameters for prior of sigma1
    :param pbetas: mixture ratios for pCN
    :param probs: mixcture probabilities for zpCN
    :param propvar: variance of proposal distribution
    :param init_coef: initial g coefficients
    :param npoits: number of sample points to use during alignment
    :param extrainfo: T/F whether additional information is returned (accept, logl, gamma_mat, gamma_stats, xdist, ydist)

    Returns mcmc struct containing
    :return f2_warped: warped f2
    :return gamma: warping function
    :return g_coef: g_coeficients
    :return sigma1: variance samples
    :return accept: accept samples
    :return logl: log-likelihood
    :return gamma_mat: posterior warping function samples
    :return gamma_stats: posterior warping function samples 95% credible intervals
    :return xdist: phase distance of posterior warping functions
    :return ydist: amplitude distance of posterior warping functions
"""
function pair_warping_expomap(f1, f2, timet; iter=20000, burnin=min(5000,iter/2),
                              alpha0=0.1, beta0=0.1, pbetas=[0.8,0.1,0.005,0.0001],
                              probs=[0.1,0.1,0.7,0.1],propvar=1.0,init_coef=zeros(20),
                              npoints=200,extrainfo=false)

    # check inputs
    if (length(f1)!=length(f2))
        error("Length of f1 and f2 must be equal")
    end
    if (length(f1)!=length(timet))
        error("Length of f1 and f2 must be equal")
    end
    if (length(pbetas)!=length(probs))
        error("In zpcn, betas must equal length of probs")
    end
    if ((length(init_coef) % 2) !=0)
        error("Length of init_cef must be even")
    end

    SIG_GAM = 13
    Mg = length(init_coef)
    sigma1_ini = 1.0
    g_coef_ini = init_coef
    nbasis = Int(Mg/2)

    # normalize time
    timet = collect(LinRange(0,1,length(timet)))
    f1 = func(timet,f1)
    f2 = func(timet,f2)

    g_basis = basis_fourier(collect(LinRange(0,1,npoints)), nbasis, 1)

    function propose_g_coef(g_coef_curr)
        probm = [0; cumsum(probs)]
        z = rand()
        d = MvNormal(propvar./sort(repeat(1:nbasis,2)))
        ind = 0
        prop = g_coef_curr
        for i in 1:length(pbetas)
            if (z <= probm[i+1] && z > probm[i])
                g_coef_new = rand(d)
                prop = sqrt(1 - pbetas[i] ^ 2) * g_coef_curr + pbetas[i] * g_coef_new
                ind = i
            end
        end
        return(prop,ind)
    end

    # srsf transformation
    q1 = f_to_srsf(f1)
    q2 = f_to_srsf(f2)

    obs_domain = q1.x
    sim_domain = collect(LinRange(0,1,npoints))
    valid_index = Int(burnin):iter

    g_temp = f_basistofunction(g_basis.x,g_coef_ini,g_basis)
    tmp = f_exp1(g_temp)
    if (minimum(tmp.y)<0)
        error("Invalid initial value of g")
    end

    # result variables
    g_coef = Array{Float64}(undef,length(g_coef_ini),iter)
    sigma1 = Vector{Float64}(undef,iter)
    logl = Vector{Float64}(undef,iter)
    SSE = Vector{Float64}(undef,iter)
    accept = Vector{Bool}(undef,iter).*false
    accept_betas = Vector{Int64}(undef,iter)

    g_coef_curr = g_coef_ini
    sigma1_curr = sigma1_ini

    g = f_basistofunction(g_basis.x, g_coef_curr, g_basis)
    SSE_curr = f_SSEg_pw(g, q1, q2)
    logl_curr = f_logl_pw(g, sigma1_curr^2, q1, q2, SSEg=SSE_curr)

    g_coef[:,1] = g_coef_curr
    sigma1[1] = sigma1_curr
    SSE[1] = SSE_curr
    logl[1] = logl_curr
    accept[1] = true

    @showprogress 1 "Computing MCMC..." for m in 2:iter
        # update g
        g_coef_curr, SSE_curr, logl_curr, accepti = f_updateg_pw(g_coef_curr, g_basis,
                                                                sigma1_curr^2, q1, q2,
                                                                SSE_curr, propose_g_coef)

        # update sigma1
        newshape = length(q1.x) / 2 + alpha0
        newscale = 1 / 2 * SSE_curr + beta0
        d = InverseGamma(newshape, newscale)
        sigma1_curr = sqrt(rand(d))
        g1 = f_basistofunction(g_basis.x, g_coef_curr, g_basis)
        logl_curr = f_logl_pw(g1, sigma1_curr^2, q1, q2, SSEg=SSE_curr)

        g_coef[:,m] = g_coef_curr
        sigma1[m] = sigma1_curr
        SSE[m] = SSE_curr
        if (extrainfo)
            logl[m] = logl_curr
            accept[m] = accepti
        end
    end

    # calculate posterior mean of psi
    pw_sim_est_psi_matrix = Array{Float64}(undef,npoints,length(valid_index))
    cnt = 1
    for k in valid_index
        g_temp = f_basistofunction(sim_domain, g_coef[:,k], g_basis)
        psi_temp = f_exp1(g_temp)
        pw_sim_est_psi_matrix[:,cnt] = psi_temp.y
        cnt += 1
    end

    result_posterior_psi_simDomain = f_psimean(sim_domain,pw_sim_est_psi_matrix)

    # resample to same number of points
    result_posterior_psi = approx(sim_domain, result_posterior_psi_simDomain.y, q1.x)

    # transform to γ
    psi_tmp = func(obs_domain, result_posterior_psi)
    result_posterior_gamma = f_phiinv(psi_tmp)
    gam0 = result_posterior_gamma.y
    result_posterior_gamma = func(obs_domain,norm_gam(gam0))

    # warped f_2
    f2_warped = warp_f_gamma(result_posterior_gamma.x, f2.y, result_posterior_gamma.y)

    Dx = []
    Dy = []
    gamma_mat = []
    gamma_stats = []
    if (extrainfo)
        gamma_mat = zeros(length(obs_domain),size(pw_sim_est_psi_matrix,2))
        Dx = Vector{Float64}(undef, size(pw_sim_est_psi_matrix,2))
        Dy = Vector{Float64}(undef, size(pw_sim_est_psi_matrix,2))
        one_v = ones(size(pw_sim_est_psi_matrix,1))
        for i = 1:size(pw_sim_est_psi_matrix,2)
            tmp = approx(sim_domain, pw_sim_est_psi_matrix[:,i], obs_domain)
            tmp = func(obs_domain,tmp)
            gam0 = f_phiinv(tmp)
            gamma_mat[:,i] = norm_gam(gam0.y)

            v = inv_exp_map(one_v, pw_sim_est_psi_matrix[:,i])
            Dx[i] = sqrt(trapz(sim_domain, v.^2))
            q2warp = warp_q_gamma(q2.x, q2.y, gamma_mat[:,i])
            Dy[i] = sqrt(trapz(q2.x, (q1.y-q2warp).^2))
        end
        gamma_stats = statsFun(gamma_mat)
    end

    # return object
    out = mcmc_results(f2_warped, result_posterior_gamma.y, g_coef[:,valid_index], result_posterior_psi,
                       sigma1[valid_index], accept[valid_index], logl[valid_index], gamma_mat, gamma_stats, Dx, Dy)

    return out
end

function f_SSEg_pw(g::func, q1::func, q2::func)
    par_domain = g.x
    obs_domain = q1.x
    f = f_exp1(g)
    exp1g_temp = f_predictfunction(f, obs_domain)
    pt = zeros(length(obs_domain))
    for i in 2:length(pt)
        tmp = func(obs_domain[1:i],exp1g_temp.y[1:i])
        pt[i] = f_L2norm(tmp)^2
    end
    temp = f_predictfunction(q2,pt)
    vec1 = (q1.y - temp.y .* exp1g_temp.y).^2

    return sum(vec1)
end

function f_logl_pw(g::func, var1, q1::func, q2::func; SSEg=0.0)
    if (SSEg == 0)
        SSEg = f_SSEg_pw(g, q1, q2)
    end
    n = length(q1.y)
    out = n*log(1/sqrt(2*pi)) - n*log(sqrt(var1)) - SSEg/(2*var1)
    return out
end

function f_updateg_pw(g_coef_curr, g_basis, var1_curr, q1::func, q2::func,
                      SSE_curr, propose_g_coef::Function)

    g_coef_prop, ind = propose_g_coef(g_coef_curr)
    g_prop = f_basistofunction(g_basis.x, g_coef_prop, g_basis)
    tmp = f_exp1(g_prop)
    while (minimum(tmp.y)<0)
        g_coef_prop, ind = propose_g_coef(g_coef_curr)
        g_prop = f_basistofunction(g_basis.x, g_coef_prop, g_basis)
        tmp = f_exp1(g_prop)
    end

    g_curr = f_basistofunction(g_basis.x, g_coef_curr, g_basis)
    if (SSE_curr == 0)
        SSE_curr = f_SSEg_pw(g_curr, q1, q2)
    end

    SSE_prop = f_SSEg_pw(g_prop, q1, q2)

    logl_curr = f_logl_pw(g_curr, var1_curr, q1, q2, SSEg=SSE_curr)
    logl_prop = f_logl_pw(g_prop, var1_curr, q1, q2, SSEg=SSE_prop)

    ratio = min(1, exp(logl_prop-logl_curr))
    u = rand()
    if (u <= ratio)
        accept = true
        return g_coef_prop, SSE_prop, logl_prop, accept
    else
        accept = false
        return g_coef_curr, SSE_curr, logl_curr, accept
    end

end

function statsFun(vec)
    a = zeros(size(vec,1))
    b = zeros(size(vec,1))
    for i in 1:size(vec,1)
        a[i] = quantile(vec[i,:],0.025)
        b[i] = quantile(vec[i,:],0.975)
    end
    out = [a b]
    return out
end
