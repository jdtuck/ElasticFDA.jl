function exp_map(psi::Vector, v::Vector)
    v_norm = l2_norm(v)
    expgam = cos(v_norm) .* psi + sin(v_norm) .* v / v_norm

    return expgam
end

function inv_exp_map(Psi::Vector, psi::Vector)
    theta = acos(inner_product(Psi,psi))

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
