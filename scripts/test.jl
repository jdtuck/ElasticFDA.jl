using JLD2, Distributions, ProgressMeter, Interpolations, Plots
@load "test/simu_data.jld2"
include("src/types.jl")
include("src/misc_funcs.jl")
include("src/srsf_funcs.jl")
include("src/geometry.jl")
include("src/pair_warping_expomap.jl")
f1 = f[:,1]
f2 = f[:,9]
iter=20000
burnin=min(5000,iter/2)
alpha0=0.1
beta0=0.1
pbetas=[0.5,0.05,0.005,0.0001]
probs=[0.1,0.1,0.7,0.1]
propvar=1.0
init_coef=zeros(20)
npoints=200
extrainfo=false

out = pair_warping_expomap(f1,f2,timet,extrainfo=true)
