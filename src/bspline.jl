"""
Compute B-Spline basis
    bs(x, df, norder, nderiv=0)
    :param x: time samples
    :param df: degree of freedom
    :param norder: order of splines
    :param nderiv: derivative number
"""
function bs(x::Vector, df, norder, nderiv=0)
    n = length(x);

    nbreaks = df - norder + 2;
    breaks = linspace(minimum(x), maximum(x), nbreaks);
    # Check Inputs
    if df <= 0
        error("Degrees of freedom <= 0.")
    end

    if norder <= 0
        error("order of splines <= 0.")
    end

    if minimum(diff(x)) < 0
        isrt = sortperm(x);
        sort!(x);
        sortwrd = true;
    else
        sortwrd = false;
    end

    if x[1] - breaks[1] < -1e-10 || x[n] - breaks[nbreaks] > 1e-10
        println([x[1], x[n]])
        error("Argument values out of range")
    end

    # set some abbreviations
    k = norder;  # order of splines
    km1 = k-1;
    nb = length(breaks);  # number of break points
    nx = length(x);  # number of argument values
    nd = nderiv+1;  # ND is order of derivative plus one
    ns = nb-2+k;  # number of splines to compute
    if ns < 1
        println("There are no B-splines for the given input")
        bsplinemat = [];
        return bsplinemat
    end
    onenx = ones(nx);
    onenb = ones(k);
    onens = ones(ns);

    # augment break sequence to get knots by addinga  K-fold knot at
    # each end.

    knots = [breaks[1].*ones(km1); breaks; breaks[nb].*ones(km1)];
    nbasis = length(knots) - k;

    # For each i, determine left(i) so that K<= left(i) < nbasis+1, and
    # within that restriction, knots(left(i)) <= pts(i) < knot(left(i)+1)

    knotslower = knots[1:nbasis];
    index = sortperm([knotslower; x]);
    point = find(index .> nbasis) - collect(1:length(x));
    left = maximum(hcat(point, k*onenx),2);

    # compute bspline values and derivatives, if needed:

    # initialize the b array.
    temp = transpose([1; zeros(km1)]);
    b = repmat(temp, nd*nx, 1);
    nxs = nd.*collect(1:nx);

    # run the recurrence simultaneously for all x(i)

    # First, bring it upt to the intended level:

    for j in 1:(k-nd)
        saved = zeros(nx);
        for r in 1:j
            leftpr = round(Integer, left + r);
            tr = knots[leftpr] - x;
            tl = x - knots[leftpr-j];
            term = b[nxs,r]./(tr+tl);
            b[nxs,r] = saved + tr.*term;
            saved = tl.*term;
        end
        b[nxs,j+1] = saved;
    end

    # save the B-spline values in successive blocks in b.

    for jj in 1:(nd-1)
        j = k - nd + jj;
        saved = zeros(nx);
        nxn = nxs - 1;
        for r in 1:j
            leftpr = left+r;
            tr = knots[leftpr] - x;
            tl = x-knots[leftpr-j];
            term = b[nxs,r]./(tr+tl);
            b[nxn,r] = saved + tr.*term;
            saved = tl.*term;
        end
        b[nxn, j+1] = saved;
        nxs = nxn;
    end

    # now use the fact that derivative values can be obtained by differencing

    for jj in (nd-1):-1:1
        j = k - jj;
        temp = collect(jj:(nd-1)).*onenx + ones(nd-jj)*nxn;
        nxs = reshape(temp, (nd-1-jj+1)*nx, 1);
        for r in j:-1:1
            leftpr = left + r;
            temp = ones(nd-jj) * (knots[leftpr] - knots[leftpr-j])/j;
            b[nxs,r] = -1.*b[nxs,r]./temp;
            b[nxs,r+1] = b[nxs,r+1] - b[nxs,r];
        end
    end

    # Finally, zero out all rows of b corresponding to x outside the basic
    # interval, [breaks[1] ... breaks[nb]]

    index = find((x.<breaks[1]) | (x .> breaks[nb]));
    if !isempty(index)
        temp = [(1-nd):0].*ones(length(index))+nd*ones(nd).*index;
        b[temp, :] = zeros(nd*length(index),k);
    end

    # set up output matrix bsplinemat
    width = maximum([ns, nbasis]) + km1 + km1;
    cc = zeros(nx*width);
    index = collect((1-nx):0)*onenb.' + nx * (left*onenb.' + onenx*collect(-km1:0).');
    index = round(Integer, index);
    cc[index] = b[nd*collect(1:nx), :];
    # (This uses the fact that for a column bector v and a matrix A,
    #  v(A)(i,j) = v(A(i,j)), all i, j.)
    index2 = round(Integer, collect((1-nx):0)*onens.' + nx * onenx*collect(1:ns).');
    bsplinemat = reshape(cc[index2],nx,ns);

    if sortwrd
        temp = copy(bsplinemat);
        bsplinemat[isrt,:] = temp;
    end

    return bsplinemat
end
