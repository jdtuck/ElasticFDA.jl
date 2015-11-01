"""
Compose to Image Reparameterizations

    apply_gam_to_gam(gamnew::Array, gam::Array)
    :param gamnew: 3-D Array describing new gamma
    :param gam:: 3-D Array describing current gamma
"""
function gamcum = apply_gam_to_gam(gamnew::Array, gam::Array)
    # gam \circ gam0
    m, n, D = size(gam);
    md = 8;
    mt = md*m;
    nt = md*n;
    U = linspace(0,1,n);
    V = linspace(0,1,m);
    Ut = linspace(0,1,nt);
    Vt = linspace(0,1,mt);
    gam_tmp = zeros(mt, nt, D);
    gam_new_tmp = zeros(mt, nt, D);
    for i = 1:D
        spl = Spline2D(U,V,gam[:,:,i],kx=1,ky=1);
        gam_tmp[:,:,i] = evalgrid(spl,Ut,Vt);
        spl = Spline2D(U,V,gamnew[:,:,i],kx=1,ky=1);
        gam_new_tmp[:,:,i] = evalgrid(spl,Ut,Vt);
    end

    gam_cum_tmp = zeros(mt,nt,D);
    for i = 1:D
        spl = Spline2D(Ut,Vt,gam_new_tmp[:,:,i],kx=1,ky=1);
        for j = 1:mt
            gam_cum_tmp[j,:,i] = spl(gam_tmp[j,:,1],gam_tmp[j,:,2]);
        end
    end

    gam_cum = zeros(m,n,D);
    for i = 1:D
        spl = Spline2D(Ut,Vt,gam_cum_tmp[:,:,i],kx=1,ky=1);
        gam_cum[:,:,i] = evalgrid(spl,U,V);
    end

    return gam_cum
end


"""
Warp img by gam

    apply_gam_to_imag(img::Array, gam::Array)
    :param img: 2-D or 3-D Array
    :param gam: 3-D Array
"""
function apply_gam_to_imag(img::Array, gam::Array)
    if ndims(img) == 3
        m,n,d = size(img);
        img_new = zeros(m,n,d);
    elseif ndims(img) == 2
        m,n = size(img);
        d = 1;
        img_new = zeros(m,n);
    end

    U = linspace(0,1,n);

    if d == 1
        spl = Spline2D(U,V,img,kx=1,ky=1);
        for j = 1:m
            img_new[j,:] = spl(gam[j,:,1],gam[j,:,2]);
        end
    else
        for i = 1:d
            spl = Spline2D(U,V, img[:,:,i],kx=1,ky=1);
            for j = 1:m
                img_new[j,:,i] = spl(gam[j,:,1],gam[j,:,2]);
            end
        end
    end

    return img_new
end


"""
Compose Image Reparameterization with Identity

    apply_gam_gamid(gamid::Array,gaminc::Array)
    :param gamid: 3-D Array
    :param gaminc: 3-D Array
"""
function apply_gam_gamid(gamid::Array,gaminc::Array)
    m,n,d = size(gamid);
    U = linspace(0,1,n);
    V = linspace(0,1,m);

    gam_cum = zeros(m,n,d);
    for j = 1:d
        spl = Spline2D(U,V,gamid[:,:,j]);
        for i = 1:m
            gam_cum[i,:,j] = spl(gaminc[i,:,1],gaminc[i,:,2]);
        end
    end

    return gam_cum
end


"""
Find 2-D gradient

    compgrad2D(f::Array)
    :param f: Array

    :return dfdu: Array
    :return dfdv: Array
"""
function compgrad2D(f)
    dims = ndims(f);
    if dims < 3
        n,t = size(f);
        d = 1;
        dfdu = zeros(n,t);
        dfdv = zeros(n,t);
    else
        n,t,d = size(f);
        dfdu = zeros(n,t,d);
        dfdv = zeros(n,t,d);
    end

    @cpp ccall((:findgrad2D, libfdaqmap), Void, (Ptr{Float64}, Ptr{Float64},
               Ptr{Float64}, Int32, Int32, Int32),
               dfdu, dfdv, f, n, t, d)

    return dfdu, dfdv
end


"""
Compute Gram Schmidt of basis using c library gradient

    gram_schmidt_c(b)
    :param b: basis
"""
function gram_schmidt_c(b)
    m,n,tmp,N = size(b);
    ds = 1./((m-1)*(n-1));

    cnt = 1;
    G = zeros(m,n,tmp,N);
    G[:,:,:,cnt] = b[:,:,:,cnt];

    dvx, dvy = compgrad2D(G[:,:,:,cnt]);
    l = vec(dvx)'*vec(dvx)*ds + vec(dvy)'*vec(dvy)*ds;

    for i = 2:N
        G[:,:,:,i] = b[:,:,:,i];
        dv1x, dv1y = compgrad2D(G[:,:,:,i]);
        for j = 1:i-1
            dv2x, dv2y = compgrad2D(G[:,:,:,j]);
            t = vec(dv1x)'*vec(dv2x)*ds + vec(dv1y)'*vec(dv2y)*ds;
            G[:,:,:,i] = G[:,:,:,i]-t.*G[:,:,:,j];
        end

        v = G[:,:,:,i];
        l = vec(v)'*vec(v)*ds;
        if l>0
            cnt += 1;
            G[:,:,:,cnt] = G[:,:,:,cnt]./sqrt(l);
        end

    return G
end


"""
Find jacobian of image F

    jacob_imag(F::Array)
    :param F: Array
"""
function jacob_imag(F::Array)
    m,n,d = size(F);

    dfdu, dfdv = compgrad2D(F);

    mult_factor = zeros(m,n);
    if (d==2)
        mult_factor = dfdu[:,:,1]*dfdv[:,:,2] - dfdu[:,:,2]*dfdv[:,:,1];
        mult_factor = abs(mult_factor);
    elseif (d==3)
        mult_factor = (dfdu[:,:,2]*dfdv[:,:,3] - dfdu[:,:,3]*dfdv[:,:,2]).^2
            + (dfdu[:,:,1]*dfdv[:,:,3] - dfdu[:,:,3]*dfdv[:,:,1]).^2
            + (dfdu[:,:,1]*dfdv[:,:,2] - dfdu[:,:,2]*dfdv[:,:,1]).^2;
        mult_factor = sqrt(mult_factor);
    end

    return mult_factor
end


"""
Form identity image warping

    makediffeoid(nrow::Integer,ncol::Integer)
    :param nrow: number of rows
    :param ncol: number of columns
"""
function makediffeoid(nrow::Integer,ncol::Integer)
    D = 2;
    gamid = zeros(nrow,ncol,D);
    U, V = meshgrid(linspace(0,1,ncol),linspace(0,1,nrow));

    gamid[:,:,1] = U;
    gamid[:,:,2] = V;

    return gamid
end


"""
Convert image to q-map

    imag_to_q(F::Array)
    :param F: Array of Image
"""
function imag_to_q(F::Array)
    dims = ndims(F);
    if dims < 3
        error("Data dimension is wrong!")
    end
    d = size(F,3);
    if d < 2
        error("Data dimension is wrong!")
    end
    q = copy(F);

    sqrtmultfact = sqrt(jacob_imag(F));
    for i = 1:d
        q[:,:,i] = sqrtmultfact.*F[:,:,i];
    end

    return q
end


"""
Compute distance between to q-maps

    comp_dist(q1::Array, q2::Array)
    :param q1: q-map1
    :param q2: q-map2
"""
function comp_dist(q1::Array, q2::Array)
    m = size(q1,1);
    n = size(q1,2);
    ds = 1./((m-1)*(n-1));

    tmp = q1-q2;

    H = sum(sqrt(vec(tmp).*vec(tmp))*ds);

    return H
end


"""
Compute energy of gradient

    comp_energy(q1::Array,q2::Array)
    :param q1: q-map1
    :param q2: q-map2
"""
function comp_energy(q1::Array,q2::Array)
    m,n,d = size(q1);
    ds = 1./(m-1)/(n-1);

    tmp = q1-q2;
    H = vec(tmp)'*vec(tmp)*ds;

    return H
end


"""
Update warping function

    update_gam(qt, qm, b)
"""
function update_gam(qt, qm, b)
    v = qt - qm;
    w = findphistar(qt,b);
    gamupdate, innp = findupdategam(v,w,b);

    return gamupdate
end


"""
Find phi*

    findphistar(q,b)
"""
function findphistar(q, b)
    m,n,tmp,K = size(b);
    d = size(q,3);
    dbxdu = zeros(m,n,K);
    dbydv = zeros(m,n,K);
    expr1 = zeros(m,n,d,K);
    expr2 = zeros(m,n,d,K);

    for k = 1:K
        d1, d2 = compgrad2D(b[:,:,1,k]);
        dbxdu[:,:,k] = d1;
        d1, d2 = compgrad2D(b[:,:,1,k]);
        dbydv[:,:,k] = d2;
    end
    divb = dbxdu + dbydv;

    dqdu, dqdv = compgrad2D(q);

    for k = 1:K
        for j = 1:d
            expr1[:,:,j,k] = divb[:,:,k].*q[:,:,j];
            expr2[:,:,j,k] = dqdu[:,:,j].*b[:,:,1,k] + dqdv[:,:,j].*b[:,:,2,k];
        end
    end

    w = 0.5.*expr1+expr2;

    return w
end


"""
Find updated warping function

    findupdategam(v, w, b)
"""
function findupdategam(v, w, b)
    m,n,D,K = size(b);
    ds = 1./((m-1)*(n-1));

    innp = zeros(K);

    gamupdate = zeros(m,n,D);

    for k = 1:K
        vt = w[:,:,:,k];
        innp[k] = vec(v)'*vec(vt)*ds;

        gamupdate = gamupdate + innp[k].*b[:,:,:,k];
    end

    return gamupdate, innp
end


"""
Check if warping function is a diffeo

    check_crossing(F)
"""
function check_crossing(f)
    n,t,D = size(f);

    if D!=2
        error("Third dimension of first argument expected to be 2")
    end

    diffeo = @cpp ccall((:check_crossing, libfdaqmap), Int32, (Ptr{Float64},
                        Int32, Int32, Int32), f, n, t, D);

    if diffeo == 0
        is_diffeo = false;
    elseif diffeo == 1
        is_diffeo = true;
    else
        error("check_crossing error")
    end

    return is_diffeo
end


"""
Compute Basis Tid

    formbasisTid(M,m,n;basis_type="tso")
    :param M: number of basis elements
    :param m:
    :param n:
    :param basis_type: type of basis (options are "t","s","o","i")
"""
function formbasisTid(M,m,n;basis_type="t")
    U, V = meshgrid(linspace(0,1,n), linspace(0,1,m));

    idx = 1;

    if basis_type=="t"
        b = zeros(m,n,2,2*M);
        for s = 1:M
            c = sqrt(2)*pi*s;
            sPI2 = 2*pi*s;

            b[:,:,1,idx] = zeros(m,n);
            b[:,:,2,idx] = (cos(sPI2*V)-1)/c;

            b[:,:,1,idx+1] = (cos(sPI2*U)-1)/c;
            b[:,:,2,idx+1] = zeros(m,n);

            idx += 2;
        end
    end

    if basis_type == "s"
        b = zeros(m,n,2,2*M);
        for s = 1:M
            c = sqrt(2)*pi*s;
            sPI2 = 2*pi*2;

            b[:,:,1,idx] = zeros(m,n);
            b[:,:,2,idx] = sin(sPI2*U)/c;

            b[:,:,1,idx+1] = sin(sPI2*U)/c;
            b[:,:,2,idx+1] = zero(m,n);

            idx += 2;
        end
    end

    if basis_type == "i"
        b = zeros(m,n,2,M*M*8);
        for s1 = 1:M
            s1PI2 = 2*pi*s1;
            for s2 = 1:M
                s2PI2 = 2*pi*s2;
                c = pi * sqrt(s1^2+3*s2^2);
                b[:,:,1,idx] = (cos(s1PI2*U)-1).*(cos(s2PI2*V))/c;
                b[:,:,2,idx] = zeros(m,n);
                b[:,:,1,idx+2] = ((cos(s1PI2*U)-1).*sin(s2PI2*V))/c;
                b[:,:,2,idx+2] = zeros(m,n);
                c = pi*sqrt(s1^2+s2^2);
                b[:,:,1,idx+4] = sin(s1PI2*U).*(cos(s2PI2*V))/c;
                b[:,:,2,idx+4] = zeros(m,n);
                b[:,:,1,idx+6] = (sin(s1PI2*U).*sin(s2PI2*V))/c;
                b[:,:,2,idx+6] = zeros(m,n);

                c = pi*sqrt(s1^2+3*s2^2);
                b[:,:,1,idx+1] = zeros(m,n);
                b[:,:,2,idx+1] = (cos(s1PI2*V)-1).*(cos(s2PI2*U))/c;
                b[:,:,1,idx+3] = zeros(m,n);
                b[:,:,2,idx+3] = ((cos(s1PI2*V)-1).*sin(s2PI2*U))/c;
                c = pi*sqrt(s1^2+s2^2);
                b[:,:,1,idx+5] = zeros(m,n);
                b[:,:,2,idx+5] = sin(s1PI2*V).*(cos(s2PI2*U))/c;
                b[:,:,1,idx+7] = zeros(m,n);
                b[:,:,2,idx+7] = (sin(s1PI2*V).*sin(s2PI2*U))/c;
                idx += 8;
            end
        end
    end

    if basis_type == "o"
        b = zeros(m,n,2,M*4);
        for s = 1:M
            c = sqrt((4*pi^2*s^2+9)/6);
            sPI2 = 2*pi*s;
            b[:,:,1,idx] = (cos(sPI2*U)-1)*V/c;
            b[:,:,2,idx] = zeros(m,n);
            b[:,:,2,idx+1] = (cos(sPI2*V)-1).*U/c;
            b[:,:,1,idx+1] = zeros(m,n);
            b[:,:,1,idx+2] = sin(sPI2*U).*V/c;
            b[:,:,2,idx+2] = zeros(m,n);
            b[:,:,2,idx+3] = sin(sPI2*V).*U/c;
            b[:,:,1,idx+3] = zeros(m,n);
            idx += 4;
        end
    end

    return b, U, V
end


"""
Generate basis elements on T_id(\gamma)

    run_basis(Ft, M=10, basis_type="t", is_orthon=true)
    :param Ft: Image Array
    :param M: number of basis elements
    :param basis_type: type of basis (option "t", "s", "i", "o")
    :param is_orthon: make basis orthonormal

    :return b: basis elements
    :retrun gamid: identity diffeo
"""
function run_basis(Ft, M=10, basis_type="t", is_orthon=true)
    m = size(Ft,1);
    n = size(Ft,2);

    gamid = makediffeoid(m,n);

    b = formbasisTid(M, m, n, basis_type);
    if is_orthon
        b = gram_schmidt_c(b);
    end

    return b, gamid
end

