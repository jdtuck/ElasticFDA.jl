module ElasticFDA

## Requirements
using Interpolations
using Dierckx
using NLopt
using ProgressMeter
using Distributions
using Printf
using Distributed
using Libdl

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
    calc_shape_dist,
    curve_pair_align,
    curve_geodesic,
    resamplecurve,
    bs,
    project_curve,
    inverse_exp_coord,
    inverse_exp,
    project_tangent

# load fdasrsf library
if (Sys.isunix())
    path = "../deps/src/fdasrsf/fdasrsf"
else
    path = "..\\deps\\fdasrsf"
end
const libfdasrsf = joinpath(dirname(@__FILE__), path)

# load gropt library
# if (Sys.isunix())
#     path1 = "../deps/src/gropt/gropt"
# else
#     path1 = "..\\deps\\gropt"
# end
# const libgropt = joinpath(dirname(@__FILE__), path1)

# load fdaqmap library
# if (Sys.isunix())
#     path2 = "../deps/src/fdaqmap/fdaqmap"
# else
#     path2 = "..\\deps\\fdaqmap"
# end
# const libfdaqmap = joinpath(dirname(@__FILE__), path2)


# Ensure library is available.
function __init__()
    if (dlopen_e(libfdasrsf) == C_NULL)
        error("libfdasrsf not properly installed. Run build(\"ElasticFDA\")")
    end

  #  if (dlopen_e(libgropt) == C_NULL)
  #       error("libgropt not properly installed. Run build(\"ElasticFDA\")")
  #   end

  #  if (dlopen_e(libfdaqmap) == C_NULL)
  #      error("libfdaqmap not properly installed. Run build(\"ElasticFDA\")")
  #  end

end

### source files
include("curve_funcs.jl")
include("curve_stats.jl")
include("srsf_funcs.jl")
include("misc_funcs.jl")
include("fPCA.jl")
include("bspline.jl")
include("dp_bayes.jl")
include("mcmc_align.jl")
include("gauss_model.jl")
include("group_warping_bayes.jl")
include("pair_warping_bayes.jl")
include("regression.jl")
include("regression_funcs.jl")
include("time_warping.jl")
#include("image_funcs.jl")
#include("image_stats.jl")

end # module
