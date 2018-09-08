"""
Computes random samples of functions from aligned data using Gaussian model

    gauss_model(warp_data; n=1, sort_samples=false)
    :param warp_data: fdawarp type from srsf_align of aligned data
    :param n: number of samples
    :param sort_samples: sort samples

    Returns warp_data containing
    :return fs: random aligned functions
    :return gams: random warping functions
    :return ft: random functions
"""
function gauss_model(warp_data::fdawarp; n=1, sort_samples=false)
    fn = warp_data.fn
    timet = warp_data.time
    qn = warp_data.qn
    gam = warp_data.gam
    # Parameters
    eps1 = eps(Float64);
    binsize = mean(diff(timet));
    M = length(timet);

    # compute mean and covariance in q-domain
    mq_new = mean(qn,dims=2);
    mididx = round(Integer, M/2);
    m_new = sign.(fn[mididx,:]) .* sqrt.(abs.(fn[mididx,:]));
    mqn = vec([mq_new; mean(m_new)]);
    qn2 = vcat(qn, m_new');
    C = Statistics.covm(qn2, mean(qn2,dims=2), 2);

    tmp = cholesky(C,check=false)
    if (!issuccess(tmp))
        for i = 1:M+1
            C[i,i] = C[i,i] + 1e-7
        end
    end
    di = MvNormal(mqn, C);
    q_s = rand(di, n);

    # compute the correspondence to the orignal function domain
    fs = zeros(M, n);
    for k in 1:n
        fs[:, k] = cumtrapzmid(timet, q_s[1:M,k].*abs.(q_s[1:M,k]),
                               sign(q_s[M+1,k])*(q_s[M+1,k]^2), mididx);
    end
    fbar = mean(fn,dims=2)
    fsbar = mean(fs,dims=2)
    err = repeat(fbar-fsbar,1,n)
    fs = fs + err

    # random warping generation
    rgam = random_gamma(gam, n);
    gams = zeros(M, n);
    for k in 1:n
        gams[:,k] = invert_gamma(rgam[:, k]);
    end

    if sort_samples
        mx = maximum(fs,1);
        seq1 = sortperm(vec(mx));

        # compute the psi-function
        psi = zeros(size(rgam));
        for ii in 1:size(rgam, 2)
            fy = gradient(rgam[:, ii], binsize);
            psi[:,ii] = fy ./ sqrt(abs(fy) + eps1);
        end
        ip = zeros(n);
        len = zeros(n);

        for i in 1:n
            tmp = ones(1,M);
            psi_tmp = (psi[:,i]/M);
            ip[i] = (tmp * psi_tmp)[1];
            len[i] = acos(ip[i]);
        end

        seq2 = sortperm(len);

        # combine x-variability and y-variability
        ft = zeros(M, n);
        time0 = [0:(M-1)]./(M-1);
        for k in 1:n
            ft[:,k] = approx(time0, fs[:,seq1[k]], gams[:,seq2[k]]);
            tmp = isnan(ft[:, k]);
            while any(tmp)
                rgam2 = random_gamma(gam, 1);
                gam_tmp = invert_gamma(rgam2);
                ft[:, k] = approx(time0, fs[:,seq1[k]], gam_tmp);
            end
        end
    else
        # combine x-variability and y-variability
        ft = zeros(M, n);
        time0 = collect(0:(M-1))./(M-1);
        for k in 1:n
            ft[:,k] = approx(time0, fs[:,k], gams[:,k]);
            tmp = isnan.(ft[:, k]);
            while any(tmp)
                rgam2 = random_gamma(gam, 1);
                gam_tmp = invert_gamma(rgam2);
                ft[:, k] = approx(time0, fs[:,k], gam_tmp);
            end
        end
    end

    out = fdawarp(warp_data.f, timet, fn, qn, warp_data.q0, warp_data.fmean,
                  warp_data.mqn, gam, warp_data.orig_var, warp_data.amp_var,
                  warp_data.phase_var, warp_data.cost, warp_data.lambda,
                  warp_data.method, warp_data.omethod, warp_data.gamI, true,
                  fs, rgam, ft, q_s)

    return out

end
