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
function pair_warping_bayes(f1, f2, timet; iter=20000, burnin=min(5000,iter/2), alpha0=0.1, beta0=0.1, pbetas=[0.5,0.05,0.005,0.0001],probs=[0.1,0.1,0.7,0.1],propvar=1,init_coef=zeros(20),npoints=200,extrainfo=false)

    # check inputs
    if (length(f1)!=length(f2))
        error('Length of f1 and f2 must be equal')
    end
    if (length(f1)!=length(timet))
        error('Length of f1 and f2 must be equal')
    end
    if (length(pbetas)!=length(probs))
        error('In zpcn, betas must equal length of probs')
    end
    if ((length(init_coef) % 2) !=0)
        error('Length of init_cef must be even')
    end

    SIG_GAM = 13
    Mg = length(init.coef)

    # normalize time
    timet = collect(LinRange(0,1,length(timet)))
    f1 = func(x=timet,y=f1)
    f2 = func(x=timet,y=f2)

    g_basis = basis_fourier(LinRange(0,1,npoints), Mg/2, 1)

    function popose_g_coef(g_coef_curr)
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

    

    return
end
