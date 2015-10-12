"""
Calculates vertical functional principal component analysis on aligned data

    vert_fPCA(fn, timet, qn; no=1)
    :param fn: array of shape (M,N) of N aligned functions with M samples
    :param timet: vector of size M describing the sample points
    :param qn: array of shape (M,N) of N aligned SRSF with M samples
    :param no: number of components to extract (default = 1)

    Returns Dict containing
    :return q_pca: srsf principal directions
    :return f_pca: functional principal directions
    :return latent: latent values
    :return coef: coefficients
    :return U: eigenvectors
"""
function vert_fPCA(fn, timet, qn; no=1)
    coef = [-2:3];
    Nstd = length(coef);

    # fPCA
    mq_new = mean(qn, 2);
    N = length(mq_new);
    mididx = round(length(timet)/2);
    m_new = sign(fn[mididx, :]) .* sqrt(abs(fn[mididx, :]));
    mqn = [mq_new, mean(m_new)];
    qn2 = vcat(qn, m_new);
    K = cov(qn2, vardim=2);

    U, s, V = svd(K);
    stdS = sqrt(s);

    # compute the PCA in the q domain
    q_pca = Array(Float64, N+1, Nstd, no);
    for k in 1:no
        for l in 1:Nstd
            q_pca[:, l, k] = mqn + coef[l] * stdS[k] * U[:, k];
        end
    end

    # compute the correspondence in the f domain
    f_pca = Array(Float64, N, Nstd, no);
    for k in 1:no
        for l in 1:Nstd
            f_pca[:, l, k] = cumtrapzmid(timet, q_pca[1:N,l,k].*
                                         abs(q_pca[1:N,l,k]),
                                         sign(q_pca[N+1,l,k]*q_pca[N,l,k]^2));
        end
    end

    N2 = size(qn,2);
    c = zeros(N2, no);
    for k in 1:no
        for l in 1:N2
            c[l,k] = sum([qn[:,l], m_new[l]] .* U[:,k]);
        end
    end

    out = Dict("q_pca" => q_pca, "f_pca" => f_pca, "latent" => s,
               "coef" => c, "U" => U);
    return out
end


"""
Calculates horizontal functional principal component analysis on aligned data

    horiz_fPCA(gam, timet; no=1)
    :param gam: array of shape (M,N) of N warping functions with M samples
    :param timet: vector of size M describing the sample points
    :param no: number of components to extract (default = 1)

    Returns Dict containing
    :return gam_pca: warping principal directions
    :return psi_pca: srsf functional principal directions
    :return latent: latent values
    :return U: eigenvectors
    :return gam_mu: mean warping function
"""
function horiz_fPCA(gam, timet; no=1)
    mu, gam_mu, psi, vec1 = sqrt_mean(gam);
    tau = [1:6];
    TT = length(timet);
    n = length(tau);
    m = length(mu);

    # TFPCA
    K = cov(vec1, vardim=2);

    U, s, V = svd(K);
    vm = mean(vec1, 2);

    gam_pca = Array(Float64, n, m+1, no);
    psi_pca = Array(Float64, n, m, no);
    for j in 1:no
        for k in tau
            v = (k-3)*sqrt(s[j])*U[:,j];
            vn = norm(v) / sqrt(TT);
            if vn < 0.0001
                psi_pca[k, :, j] = mu;
            else
                psi_pca[k, :, j] = cos(vn).*mu + sin(vn).*v/vn;
            end

            tmp = zeros(TT);
            tmp[2:TT] = cumsum(psi_pca[k,:,j] .* psi_pca[k,:,j],2);
            gam_pca[k,:,j] = (tmp - tmp[1]) ./ (tmp[end] - tmp[1]);
        end
    end

    out = Dict("gam_pca" => gam_pca, "psi_pca" => psi_pca, "latent" => s,
               "U" => U, "gam_mu" => gam_mu);
    return out
end
