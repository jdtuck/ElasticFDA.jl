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


function trapz(x, y)
    N = size(y,2);
    if N == 1
        M = length(y);
        out = sum(diff(x).*(y[1:M-1] + y[2:M])/2);
    else
        M = size(y,1);
        out = zeros(N);
        for i in 1:N
            out[i] = sum(diff(x).*(y[1:M-1,i]+y[2:M,i])/2);
        end
    end
    return out
end


function cumtrapz(x, y)
    m = length(y);
    dt = diff(x)/2.0;
    z = [0, cumsum(dt.*(y[1:(m-1)] + y[2:m]))];

    return z
end


function cumtrapzmid(x, y, c)
    a = length(x);
    mid = round(a/2);

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


function mvnrand(mu, C, n)
    tmp = cholfact(C, :U, pivot=true);
    R = tmp[:U];
    R = R[:, tmp.piv];
    retval = randn(n, size(R,1)) * R;
    retval += transpose(repmat(mu, 1, n));
    return transpose(retval)
end
