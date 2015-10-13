"""
Computes random samples of functions from aligned data using Gaussian model

    gauss_model(fn, timet, qn, gam; n=1, sort_samples=false)
    :param fn: aligned functions (M,N)
    :param timet: vector (M) describing time
    :param qn: aligned srvfs (M,N)
    :param gam: warping functions (M,N)
    :param n: number of samples
    :param sort_samples: sort samples

    Returns Dict containing
    :return fs: random aligned functions
    :return gams: random warping functions
    :return ft: random functions
"""
function gauss_model(fn, timet, qn, gam; n=1, sort_samples=false)
    # Parameters
    eps1 = eps(Float64);
    binsize = mean(diff(timet));
    M = length(timet);

    # compute mean and covariance in q-domain
    mq_new = mean(qn,2);
    mididx = round(Integer, M/2);
    m_new = sign(fn[mididx,:]) .* sqrt(abs(fn[mididx,:]));
    mqn = vec([mq_new; mean(m_new)]);
    qn2 = vcat(qn, m_new);
    C = cov(qn2, vardim=2);

    q_s = mvnrand(mqn, C, n);

    # compute the correspondence to the orignal function domain
    fs = zeros(M, n);
    for k in 1:n
        fs[:, k] = cumtrapzmid(timet, q_s[1:M,k].*abs(q_s[1:M,k]),
                               sign(q_s[M+1,k])*(q_s[M+1,k]^2));
    end

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
            tmp = isnan(ft[:, k]);
            while any(tmp)
                rgam2 = random_gamma(gam, 1);
                gam_tmp = invert_gamma(rgam2);
                ft[:, k] = approx(time0, fs[:,k], gam_tmp);
            end
        end
    end

    out = Dict("fs" => fs, "gams" => rgam, "ft" => ft);
    return out

end
