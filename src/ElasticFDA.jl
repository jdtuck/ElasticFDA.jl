module ElasticFDA

## Requirements
using Grid
using Dierckx
using NLopt
using ProgressMeter
using Distributions
using Cpp
using Winston

export
    smooth_data!,
    smooth_data,
    f_to_srsf,
    srsf_to_f,
    optimum_reparam,
    warp_q_gamma,
    warp_f_gamma,
    vert_fPCA,
    horiz_fPCA,
    gauss_model,
    group_warping_bayes,
    pair_warping_bayes,
    elastic_distance,
    elastic_regression,
    elastic_logistic,
    elastic_mlogistic,
    elastic_prediction,
    srsf_align,
    align_fPCA,
    trapz,
    rgam,
    curve_karcher_mean,
    curve_srvf_align,
    curve_principal_directions,
    curve_karcher_cov,
    sample_shapes,
    curve_to_q,
    q_to_curve,
    calc_shape_dist

# load fdasrsf library
unixpath = "../deps/src/fdasrsf/fdasrsf"
winpath = "../deps/bin$WORD_SIZE/fdasrsf"
const libfdasrsf = joinpath(dirname(@__FILE__), @unix? unixpath : winpath)

# Ensure library is available.
if (Libdl.dlopen_e(libfdasrsf) == C_NULL)
    error("libfdasrsf not properly installed. Run Pkg.build(\"ElasticFDA\")")
end

# load gropt library
unixpath1 = "../deps/src/gropt/gropt"
winpath1 = "../deps/bin$WORD_SIZE/gropt"
const libgropt = joinpath(dirname(@__FILE__), @unix? unixpath1 : winpath1)

# Ensure library is available.
if (Libdl.dlopen_e(libgropt) == C_NULL)
    error("libgropt not properly installed. Run Pkg.build(\"ElasticFDA\")")
end

### source files
include("curve_funcs.jl")
include("curve_stats.jl")
include("srsf_funcs.jl")
include("misc_funcs.jl")
include("fPCA.jl")
include("bspline.jl")
include("regression_funcs.jl")
include("dp_bayes.jl")
include("mcmc_align.jl")
include("gauss_model.jl")
include("group_warping_bayes.jl")
include("pair_warping_bayes.jl")
include("regression.jl")
include("regression_funcs.jl")
include("time_warping.jl")

end # module
