struct warp_data
    f::Array{Float64}
    time::Vector{Float64}
    fn::Array{Float64}
    qn::Array{Float64}
    q0::Array{Float64}
    fmean::Vector{Float64}
    mqn::Vector{Float64}
    gam::Array{Float64}
    orig_var::Float64
    amp_var::Float64
    phase_var::Float64
    cost::Vector{Float64}
    lambda::Float64
    method::String
    omethod::String
    gamI::Vector{Float64}
    rsamps::Bool
end

struct vfpca
    q_pca::Array{Float64}
    f_pca::Array{Float64}
    latent::Vector{Float64}
    coef::Array{Float64}
    U::Array{Float64}
    id::Float64
    mqn::Vector{Float64}
    time::Vector{Float64}
end

struct hfpca
    gam_pca::Array{Float64}
    psi_pca::Array{Float64}
    latent::Vector{Float64}
    U::Array{Float64}
    coef::Array{Float64}
    vec::Array{Float64}
    gam_mu::Vector{Float64}
end
