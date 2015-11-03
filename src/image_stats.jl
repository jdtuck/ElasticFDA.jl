"""
Find optimum reparam between two images

    reparam_image(It::Array,Im::Array,gam::Array,b::Array;
                  stepsize::Float64=1e-5, itermax::Int32=1000,
                  lmark::Bool=false)
    :param It: Template Image
    :param Im: Test Image
    :param gam: initial warping
    :param b: basis matrix
    :param stepsize: gradient stepsize
    :param itermax: maximum number of iterations
    :param lmark: use landmarks

    :return gamnew: final warping
    :return Inew: parameterized image
    :return H: energy
    :return stepsize: final stepsize
"""
function reparam_image(It::Array, Im::Array, gam::Array, b::Array;
                       stepsize::Float64=1e-5, itermax::Integer=1000,
                       lmark::Bool=false)
    m = size(It,1);
    n = size(It,2);
    gamid = makediffeoid(m,n);

    # Main Loop
    H = zeros(itermax+1);
    Iold = apply_gam_to_imag(Im,gam);
    Iold -= minimum(Iold);
    Iold /= maximum(Iold);
    gamold = copy(gam);
    qt = imag_to_q(It);
    qm = imag_to_q(Iold);

    gamnew = copy(gamold);
    Inew = copy(Iold);
    iter = 1;
    H = zeros(itermax+1);
    H[iter] = comp_energy(qt,qm);
    @printf("Iteration %d, energy %f\n",iter-1,H[iter])

    gamupdate = update_gam(qt,qm,b);
    cutoff = 1e-3;
    for iter = 2:(itermax+1)
        gaminc = gamid + stepsize.*gamupdate;
        G = check_crossing(gamnew);
        if !G
            @printf("Possible Crossing!")
            gamnew = copy(gamold);
            stepsize = 0.67*stepsize;
            H[iter] = H[iter=1];
            continue
        else
            gamnew = apply_gam_to_gam(gamnew, gaminc);
        end

        Inew = apply_gam_to_imag(Im, gamnew);
        Inew -= minimum(Inew);
        Inew /= maximum(Inew);
        qm = imag_to_q(Inew);
        H[iter] = comp_energy(qt,qm);
        @printf("Iteration %d, energy %f\n",iter-1,H[iter])

        if (iter > 4)
            hstop = 1.;
            for i = 1:4
                hstop = hstop*(H[iter] >= H[iter-1]);
            end
            if (hstop!=0)
                @printf("Warning: energy constantly increasing")
                break
            end
        end

        if (iter > 4) && (H[iter]>=H[iter-1]) && (H[iter-1]>=H[iter-2]) &&
            (H[iter-2] >= H[iter-3])
            @printf("Warning: energy is not changing")
        end

        if ((iter>1) && (H[iter]>H[iter-1])) || ((iter>3) &&
            ((H[iter-1]<=H[iter-2]) && (H[iter-2]>H[iter-3])))
            stepsize = 0.9*stepsize;
            gamnew = copy(gamold);
            H[iter] = H[iter-1];
            continue
        end

        gamold = copy(gamnew);
        gamupdate = update_gam(qt, qm, b);
    end

    H = H[1:iter];

    return gamnew, Inew, H, stepsize
end


"""
Pairwise align two images

    pair_align_image(I1, I2; M=5, ortho=true, basis_type="t", resizei=true,
                     N=64, stepsize=1e-5, itermax=1e3)
    :param I1: reference image
    :param I2: image to warp
    :param M: number of basis elements
    :param ortho: orthonormalize basis
    :param basis_type: type of basis ("t", "s", "i", "o")
    :param resizei: resize image
    :param N: size of resized image
    :param stepsize: gradient stepsize
    :param itermax: maximum number of iterations

    :return I2_new: aligned I2
    :return gam: warping function
"""
function pair_align_image(I1, I2; M=5, ortho=true, basis_type="t", resizei=true,
                          N=64, stepsize=1e-5, itermax=1000)
    m,n = size(I1);
    F1 = zeros(m,n,2);
    m1,n1 = size(I2);
    F2 = zeros(m,n,2);

    # Take Gradient
    fu, fv = gradient2(I1,1./(m-1),1./(n-1));
    F1[:,:,1] = fu;
    F1[:,:,2] = fv;
    fu, fv = gradient2(I2,1./(m1-1),1./(n1-1));
    F2[:,:,1] = fu;
    F2[:,:,2] = fv;

    # Resize Data and Center
    if resizei
        if N > m || N > n
            @printf("Not resizing, N is larger than image size")
        else
            m_n = linspace(1,m,N);
            n_n = linspace(1,n,N);
            F1a = zeros(N,N,2);
            spl = Spline2D(collect(1:m),collect(1:n),F1[:,:,1]);
            F1a[:,:,1] = evalgrid(spl,m_n,n_n);
            spl = Spline2D(collect(1:m),collect(1:n),F1[:,:,2]);
            F1a[:,:,2] = evalgrid(spl,m_n,n_n);
            F1 = copy(F1a);

            m_n = linspace(1,m1,N);
            n_n = linspace(1,n1,N);
            F2a = zeros(N,N,2);
            spl = Spline2D(collect(1:m1),collect(1:n1),F2[:,:,1]);
            F2a[:,:,1] = evalgrid(spl,m_n,n_n);
            spl = Spline2D(collect(1:m1),collect(1:n1),F2[:,:,2]);
            F2a[:,:,2] = evalgrid(spl,m_n,n_n);
            F2 = copy(F2a);
        end
    end
    F1 -= minimum(F1);
    F1 /= maximum(F1);
    F2 -= minimum(F2);
    F2 /= maximum(F2);

    # Generate basis
    b, gamid = run_basis(F1, M, basis_type, ortho);
    gamp = copy(gamid);
    gam, F2_new, H, stepsize = reparam_image(F1, F2, gamp, b,
                                             stepsize=stepsize,
                                             itermax=itermax);

    I2_new = apply_gam_to_imag(I2, gam);

    return I2_new, gam
end

