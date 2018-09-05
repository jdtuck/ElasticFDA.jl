"""
1-D Gradient

    gradient(a::Array, dx::Float64=1)
    :param f: Vector
    :param dx: stepsize

    :return g: gradient
"""
function gradient(f::Array{Float64, 1}, dx::Float64=1.)
    n = length(f)
    g = zeros(n)
    h = dx*(1:n)
    # Take forward differences on left and right edges
    g[1] = (f[2] - f[1])/(h[2]-h[1])
    g[n] = (f[n] - f[n-1])/(h[end]-h[end-1])

    # Take centered differences on interior points
    h = h[3:n]-h[1:(n-2)]
    g[2:(n-1)] = (f[3:n]-f[1:(n-2)])/h[1]

    return g
end


"""
2-D Gradient

    gradient2(a::Array, dx::Float64=1, dy::Float=1)
    :param a: matrix
    :param dx: stepsize
    :param dy: stepsize

    :return dxdu: derivatives along x
    :return dydv: derivatives along y
"""
function gradient2(a::Array, dx::Float64=1., dy::Float64=1.)
    m,n = size(a);
    dxdu = zeros(m,n);
    dydv = zeros(m,n);

    for i=1:m
        dxdu[i,:] = gradient(vec(a[i,:]), dx)
    end

    for i=1:n
        dydv[:,i] = gradient(a[:,i], dy)
    end

    return dxdu, dydv
end


"""
Creates Rectangular Grid in 2-D space

    meshgrid(a::LinRange,b::LinRange)
"""
function meshgrid(a::LinRange,b::LinRange)
    grid_a = [i for i in a, j in b]';
    grid_b = [j for i in a, j in b]';

    return grid_a, grid_b
end


"""
Linear interpolation

    approx(xd, yd, xi)
    :param xd: x samples
    :param yd: response samples
    :param xi: new x samples
"""
function approx(xd, yd, xi)
    nd = length(xd);
    ni = length(xi);

    yi = zeros(ni);
    for i in 1:ni
        if (xi[i] <= xd[1])
            t = (xi[i]-xd[1]) / (xd[2] - xd[1]);
            yi[i] = (1.0 - t) * yd[1] + t * yd[2];
        elseif (xd[nd] <= xi[i])
            t = (xi[i] - xd[nd-1]) / (xd[nd] - xd[nd-1]);
            yi[i] = (1.0 - t) * yd[nd-1] + t * yd[nd];
        else
            for k in 2:nd
                if (xd[k-1] <= xi[i] && xi[i] <= xd[k])
                    t = (xi[i] - xd[k-1]) / (xd[k] - xd[k-1]);
                    yi[i] = (1.0 - t) * yd[k-1] + t * yd[k];
                    break
                end
            end
        end
    end

    return yi
end


"""
Trapezoidal Integration

    trapz(x, y, dim=1)
    :param x: vector of time samples
    :param y: array of response samples
    :param dim: dimension along which to integrate
"""
function trapz(x::Array{Float64, 1}, y::Array{Float64}, dim::Integer=1)
    perm = [dim:max(ndims(y),dim); 1:dim-1];
    y = permutedims(y, perm);
    if ndims(y) == 1
        m = 1;
    else
        m = size(y,1);
    end

    if m == 1
        M = length(y);
        out = sum(diff(x).*(y[1:M-1] + y[2:M])/2.0);
    else
        out = transpose(diff(x)) * (y[1:m-1,:] + y[2:m,:])/2.0;
        siz = size(y); siz = collect(siz); siz[1] = 1;
        out = reshape(out, tuple(siz...));
        out = ipermutedims(out, perm);
        ind = find(collect(size(out)).==1);
        out = squeeze(out,ind[1]);
        if length(out) == 1;
            out = out[1];
        end
    end

    return out
end


"""
Cumulative Trapezoidal Integration

    cumtrapz(x, y, dim=1)
    :param x: vector describing time samples
    :param y: array describing response
    :param dim: dimension to integrate over
"""
function cumtrapz(x::Array{Float64, 1}, y::Array{Float64}, dim::Integer=1)
    perm = [dim:max(length(size(y)),dim); 1:dim-1];
    y = permutedims(y, perm);
    if ndims(y) == 1
        n = 1;
        m = length(y);
    else
        m, n = size(y);
    end

    if n == 1
        dt = diff(x)/2.0;
        z = [0; cumsum(dt.*(y[1:(m-1)] + y[2:m]))];
    else
        dt = repmat(diff(x)/2.0,1,n);
        z = [zeros(1,n); cumsum(dt.*(y[1:(m-1), :] + y[2:m, :]),1)];
        z = ipermutedims(z, perm);
    end

    return z
end


"""
Cumulative Trapezoidal Integration using midpoint

    cumtrapzmid(x, y, c)
    :param x: time samples
    :param y: resposne samples
    :param c: midpoint
    :param mid: midpoint location
"""
function cumtrapzmid(x, y, c, mid)
    a = length(x);

    # case < mid
    fn = zeros(a);
    tmpx = x[(mid-1):-1:1];
    tmpy = y[(mid-1):-1:1];
    tmp = c + cumtrapz(tmpx, tmpy);
    fn[1:(mid-1)] = reverse(tmp);

    # case >= mid
    fn[mid:a] = c + cumtrapz(x[mid:a],y[mid:a]);

    return fn

end


"""
Multivariate Normal random number generation

    mvnrand(mu, C, n)
    :param mu: mean vector
    :param C: covariance matrix
    :param n: number of samples
"""
function mvnrand(mu, C, n)
    tmp = cholfact(C, :U, Val{true});
    R = tmp[:U];
    R = Array(R);
    R = R[:, tmp.piv];
    retval = randn(n, size(R,1)) * R;
    retval += transpose(repmat(mu, 1, n));
    return transpose(retval)
end


"""
Linear interpolation when response contains flat regions

    interp1_flat(x, y, xx)
    :param x: time samples
    :param y: response samples
    :param xx: new time samples
"""
function interp1_flat(x,y,xx)
    flat = find(diff(x).<=0);
    n = length(flat);

    if n==0
        tmp = interpolate((x,), y, Gridded(Linear()))
        yy = tmp[xx];
    else
        yy = zeros(size(xx));
        i1 = 1;
        if flat[1] == 1
            i2 = 1;
            j = xx.==x[i2];
            yy[j] = minimum(y[i2:i2+1]);
        else
            i2 = flat[1];
            j = (xx.>=x[i1]) & (xx.<=x[i2]);
            tmp = interpolate((x[i1:i2],), y[i1:i2], Gridded(Linear()))
            yy[j] = tmp[xx[j]];
            i1 = copy(i2);
        end
        for k = 2:n
            i2 = flat[k];
            if i2 > i1+1
                j = (xx.>=x[i1]) & (xx.<=x[i2]);
                yi = interpolate((x[i1+1:i2],), y[i1+1:i2], Gridded(Linear()))
                yy[j] = tmp[xx[j]];
            end
            j = xx.==x[i2];
            yy[j] = minimum(y[i2:i2+1]);
            i1 = copy(i2);
        end
        i2 = length(x);
        j = (xx.>=x[i1]) & (xx.<=x[i2]);
        if i1+1 == i2
            yy[j] = y[i2];
        else
            tmp = interpolate((x[i1+1:i2],), y[i1+1:i2], Gridded(Linear()))
            yy[j] = tmp[xx[j]];
        end
    end

    return yy
end
