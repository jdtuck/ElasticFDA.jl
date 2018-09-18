function exp_map(psi::Vector, v::Vector)
    v_norm = l2_norm(v)
    expgam = cos(v_norm) .* psi + sin(v_norm) .* v / v_norm

    return expgam
end

function inv_exp_map(Psi::Vector, psi::Vector)
    dot = inner_product(Psi,psi)
    if dot > 1
        dot = 1
    end
    if dot < -1
        dot = -1
    end

    theta = acos(dot)

    if (theta < 1e-10)
        exp_inv = zeros(length(psi))
    else
        exp_inv = theta / sin(theta) .* (psi .- cos(theta) .* Psi)
    end

    return exp_inv
end

function l2_norm(psi::Vector; timet=LinRange(0,1,length(psi)))
    timet = collect(timet)
    l2norm = sqrt(trapz(timet, psi.*psi))
    return l2norm
end

function inner_product(psi1, psi2; timet=LinRange(0,1,length(psi1)))
    timet = collect(timet)
    ip = trapz(timet, psi1.*psi2)
    return ip
end

function f_exp1(g::func)
    area = f_L2norm(g)
    y = cos(area) .+ sin(area)/area*g.y
    if (area==0)
        y = ones(length(g.x))
    end
    out = func(g.x,y)
    return out
end

function f_L2norm(f::func)
    y = f.y
    x = f.x
    ind = sortperm(x)
    area = trapz(x[ind], (y[ind]).^2)
    return sqrt(area)
end

function f_psimean(x, y; e1=0.001, e2=0.3)
    rmy = mean(y,dims=2)
    tmp = func(x, rmy)
    result = rmy / f_L2norm(tmp)
    out = func(x, result)
    return out
end

function f_phiinv(psi::func)
    f_domain = psi.x

    result = zeros(length(f_domain))
    for i in 2:length(result)
        tmp = func(f_domain[1:i], psi.y[1:i])
        result[i] = f_L2norm(tmp)^2
    end
    out = func(f_domain,result)
    return out
end
