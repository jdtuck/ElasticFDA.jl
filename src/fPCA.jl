function vert_fPCA(fn, timet, qn; no=1)
    ####################################################################
    # This function calculates vertical functional principal component analysis
    # on aligned data

    # :param fn: numpy ndarray of shape (M,N) of M aligned functions with N
    #            samples
    # :param time: vector of size N describing the sample points
    # :param qn: numpy ndarray of shape (M,N) of M aligned SRSF with N samples
    # :param no: number of components to extract (default = 1)
    # :param showplot: Shows plots of results using matplotlib (default = T)
    # :type showplot: bool
    # :type no: int

    # :return q_pca: srsf principal directions
    # :return f_pca: functional principal directions
    # :return latent: latent values
    # :return coef: coefficients
    # :return U: eigenvectors
    ####################################################################
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

    out = ["q_pca" => q_pca, "f_pca" => f_pca, "latent" => s,
           "coef" => c, "U" => U];
    return out
end


function horiz_fPCA(fn, timet, qn; no=1)
    ####################################################################
    # This function calculates horizontal functional principal component
    # analysis on aligned data

    # :param gam: numpy ndarray of shape (M,N) of M warping functions
    # :param time: vector of size N describing the sample points
    # :param no: number of components to extract (default = 1)
    # :param showplot: Shows plots of results using matplotlib (default = T)
    # :type showplot: bool
    # :type no: int

    # :return q_pca: srsf principal directions
    # :return f_pca: functional principal directions
    # :return latent: latent values
    # :return coef: coefficients
    # :return U: eigenvectors
    ####################################################################

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

    out = ["gam_pca" => gam_pca, "psi_pca" => psi_pca, "latent" => s,
           "U" => U, "gam_mu" => gam_mu];
    return out
end
