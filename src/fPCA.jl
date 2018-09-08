"""
Calculates vertical functional principal component analysis on aligned data

    vert_fPCA(warp_data, qn; no=1)
    :param warp_data: fdawarp type from srsf_align of aligned data
    :param no: number of components to extract (default = 1)

    Returns vfpca type containing
    :return q_pca: srsf principal directions
    :return f_pca: functional principal directions
    :return latent: latent values
    :return coef: scores
    :return U: eigenvectors
    :return id: point used for f(0)
"""
function vert_fPCA(warp_data::fdawarp; no=1)
    qn = warp_data.qn
    fn = warp_data.fn
    timet = warp_data.time
    coef = collect(-2:3);
    Nstd = length(coef);

    # fPCA
    mq_new = mean(qn, dims=2);
    N = length(mq_new);
    mididx = round(Integer, length(timet)/2);
    m_new = sign.(fn[mididx, :]) .* sqrt.(abs.(fn[mididx, :]));
    mqn = [mq_new; mean(m_new)];
    qn2 = vcat(qn, m_new');
    K = Statistics.covm(qn2, mean(qn2,dims=2), 2);

    U, s, V = svd(K);
    stdS = sqrt.(s);

    # compute the PCA in the q domain
    q_pca = Array{Float64}(undef, (N+1, Nstd, no));
    for k in 1:no
        for l in 1:Nstd
            q_pca[:, l, k] = mqn + coef[l] * stdS[k] * U[:, k];
        end
    end

    # compute the correspondence in the f domain
    f_pca = Array{Float64}(undef, (N, Nstd, no));
    for k in 1:no
        for l in 1:Nstd
            f_pca[:, l, k] = cumtrapzmid(timet, q_pca[1:N,l,k].*
                                         abs.(q_pca[1:N,l,k]),
                                         sign.(q_pca[N+1,l,k]*q_pca[N,l,k]^2),
                                         mididx);
        end
        fbar = mean(fn,dims=2)
        fsbar = mean(f_pca[:, :, k],dims=2)
        err = repeat(fbar-fsbar,1,Nstd)
        f_pca[:, :, k] = f_pca[:, :, k] + err
    end

    N2 = size(qn,2);
    c = zeros(N2, no);
    for k in 1:no
        for l in 1:N2
            c[l,k] = sum(vec((qn2[:,l] - mqn).* U[:,k]));
        end
    end

    out = vfpca(q_pca, f_pca, s[1:no], c, U[:,1:no], mididx, mqn[:], timet);

    return out
end


"""
Calculates horizontal functional principal component analysis on aligned data

    horiz_fPCA(warp_data; no=1)
    :param warp_data: fdawarp type from srsf_align of aligned data
    :param no: number of components to extract (default = 1)

    Returns hfpca containing
    :return gam_pca: warping principal directions
    :return psi_pca: srsf functional principal directions
    :return latent: latent values
    :return U: eigenvectors
    :return coef: scores
    :return vec: shooting vectors
    :return gam_mu: mean warping function
"""
function horiz_fPCA(warp_data::fdawarp; no=1)
    gam = warp_data.gam;
    mu, gam_mu, psi, vec1 = sqrt_mean(gam);
    tau = collect(1:5);
    m = length(mu);
    n = length(tau);
    TT = size(vec1,1)

    # TFPCA
    K = Statistics.covm(vec1, mean(vec1,dims=2), 2);

    U, s, V = svd(K);
    vm = mean(vec1, dims=2);

    gam_pca = Array{Float64}(undef, (m, n, no));
    psi_pca = Array{Float64}(undef, (m, n, no));
    for j in 1:no
        for k in tau
            v = (k-3)*sqrt(s[j]).*U[:,j];
            vn = norm(v) / sqrt(TT);
            if vn < 0.0001
                psi_pca[:, k, j] = mu;
            else
                psi_pca[:, k, j] = cos(vn).*mu + sin(vn).*v/vn;
            end

            gam0 = cumtrapz(collect(LinRange(0,1,TT)), psi_pca[:,k,j].*psi_pca[:,k,j]);
            gam_pca[:,k,j] = norm_gam(gam0)
        end
    end

    N2 = size(gam,2);
    c = zeros(N2, no);
    for k in 1:no
        for l in 1:N2
            c[l,k] = sum((vec1[:,l]-vm).* U[:,k]);
        end
    end

    out = hfpca(gam_pca, psi_pca, s[1:no], U[:,1:no], c, vec1, gam_mu)
    return out
end
