"""
Compute pair warping between two functions using Bayesian method

    pair_warping_expomap(f1, f2, timet; iter=15000, times=5, powera=1)
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
function pair_warping_expomap(f1, f2, timet; iter=20000, burnin=min(5000,iter/2),
                              alpha0=0.1, beta0=0.1, pbetas=[0.5,0.05,0.005,0.0001],
                              probs=[0.1,0.1,0.7,0.1],propvar=1,init_coef=zeros(20),
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
    Mg = length(init.coef)
    sigma1_ini = 1

    # normalize time
    timet = collect(LinRange(0,1,length(timet)))
    f1 = func(x=timet,y=f1)
    f2 = func(x=timet,y=f2)

    g_basis = basis_fourier(LinRange(0,1,npoints), Mg/2, 1)

    function propose_g_coef(g_coef_curr)
        probm = [0, cumsum(probs)]
        z = rand()
        d = MvNormal(propvar./sort(repeat(1:10,2)))
        for i in 1:length(pbetas)
            if (z <= probm[i+1] && z > probm[i])
                g_coef_new = randn(d)
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
    valid_index = burnin:iter

    tmp = f_exp1(f_basistofunction(g_basis.x,g_coef_ini,g_basis))
    if (minimum(tmp.y)<0)
        error("Invalid initial value of g")
    end

    # result variables
    g_coef = Array{Float64}(undef,length(g_coef_ini),iter)
    sigma1 = Vector{Float64}(undef,iter)
    logl = Vector{Float64}(undef,iter)
    SSE = Vector{Float64}(undef,iter)
    accept = Vector{Bool}(undef,iter)
    accept_betas = Vector{Int64}(undef,iter)

    g_coef_curr = g_coef_ini
    sigma1_curr = sigma1_ini

    g = f_basistofunction(g_basis.x, g_coef_curr, g_basis)
    SSE_curr = f_SSEg_pw(g, q1, q2)
    logl_curr = f_logl_pw(g, sigma1_curr^2, q1, q2, SSE_curr)

    g_coef[:,1] = g_coef_curr
    sigma1[1] = sigma1_curr
    SSE[1] = SSE_curr
    logl[1] = logl_curr
    accept[1] = true

    for m in 2:iter
        # update g
        g_coef_curr, SSE_curr, logl_curr, accept = f_updateg_pw(g_coef_curr, g_basis,
                                                                sigma1_curr^2, q1, q2,
                                                                SSE_curr, propose_g_coef)

        # update sigma1
        newshape = length(q1.x) / 2 + alpha0
        newscale = 1 / 2 * SSE_curr + beta0
        d = InverseGamma(newshape, newscale)
        sigma1_curr = rand(d)
        g = f_basistofunction(g_basis.x, g_coef_curr, g_basis)
        logl_curr = f_logl_pw(g, sigma1_curr^2, q1, q2, SSE_curr)

        g_coef[:,m] = g_coef_curr
        sigma1[m] = sigma1_curr
        SSE[m] = SSE_curr
        if (extrainfo)
            logl[m] = logl_curr
            accept[m] = accept
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
    result_posterior_gamma = f_phiinv(result_posterior_psi)
    gam0 = result_posterior_gamma.y
    result_posterior_gamma.y = norm_gam(gam0)

    # warped f_2
    f2_warped = warp_f_gamma(result_posterior_gamma.x, f2, result_posterior_gamma.y)

    Dx = []
    Dy = [];
    gamma_mat = [];
    gamma_stats = [];
    if (extrainfo)
        gamma_mat = zeros(size(pw_sim_est_psi_matrix))
        Dx = Vector{Float64}(undef, size(pw_sim_est_psi_matrix,2))
        Dy = Vector{Float64}(undef, size(pw_sim_est_psi_matrix,2))
        one_v = ones(size(pw_sim_est_psi_matrix,1))
        for i = 1:size(pw_sim_est_psi_matrix,2)
            tmp = approx(sim_domain, pw_sim_est_psi_matrix, q1.x)
            tmp = func(q1.x,tmp)
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
    out = mcmc_results(f2_warped, result_posterior_gamma.y, g_coef, result_posterior_psi.y,
                       sigma1, accept, logl, gamma_mat, gamma_stats, Dx, Dy)

    return out
end

function f_SSEg_pw(g::func, q1::func, q2::func)
    par_domain = g.x
    obs_domain = q1.x
    f = f_exp1(g)
    exp1g_temp = f_predictfunction(f, obs_domain)
    pt = zeros(length(obs_domain))
    for i in 2:length(pt)
        tmp = func(obs.domain[1:i],exp1g_temp.y[1:i])
        pt[i] = f_L2norm(tmp)^2
    end
    temp = f_predictfunction(q2,pt)
    vec = (q1.y - temp.y * exp1g_temp.y).^2

    return sum(vec)
end

function f_logl_pw(g::func, q1::func, q2::func, var1; SSEg=0)
    if (SSEg == 0)
        SSEg = f_SSEg_pw(g, q1, q2)
    end
    n = length(q1.y)
    out = n*log(1/sqrt(2*pi)) - n*log(sqrt(var1)) - SSEg/(2*var1)
    return out
end

function f_updateg_pw(g_coef_curr, g_basis, var1_curr, q1::func, q2::func,
                      SSEg_curr, propose_g_coef::Function)

    g_coef_prop = propose_g_coef(g_coef_curr)

    g_curr = f_basistofunction(g_basis.x, g_coef_curr, g_basis)
    if (SSE_curr == 0)
        SSE_curr = f_SSEg_pw(g_curr, q1, q2)
    end

    g_prop = f_basistofunction(g_basis.x, g_coef_prop, g_basis)
    SSE_prop = f_SSEg_pw(g_prop, q1, q2)

    logl_curr = f_logl_pw(g_curr, var1_curr, q1, q2, SSE_curr)
    logl_prop = f_logl_pw(g_prop, var1_curr, q1, q2, SSE_prop)

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
