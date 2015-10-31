"""
Find optimum reparam between two images

    reparam_image(It::Array,Im::Array,gam::Array,b::Array;
    stepsize::Float64=1e-5, itermax::Int32=1000,lmark::Bool=false)
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
    stepsize::Float64=1e-5, itermax::Int32=1000,lmark::Bool=false)
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

