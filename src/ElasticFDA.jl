module ElasticFDA

## Requirements
using Grid
using Dierckx
using NLopt
using ProgressMeter
using Distributions

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
    reparam_image,
    pair_align_image,
    bs,
    @cpp

# load fdasrsf library
if (is_unix())
    path = "../deps/src/fdasrsf/fdasrsf"
else
    path = "..\\deps\\fdasrsf"
end
const libfdasrsf = joinpath(dirname(@__FILE__), path)

# load gropt library
if (is_unix())
    path1 = "../deps/src/gropt/gropt"
else
    path1 = "..\\deps\\gropt"
end
const libgropt = joinpath(dirname(@__FILE__), path1)

# load fdaqmap library
if (is_unix())
    path2 = "../deps/src/fdaqmap/fdaqmap"
else
    path2 = "..\\deps\\fdaqmap"
end
const libfdaqmap = joinpath(dirname(@__FILE__), path2)


# Ensure library is available.
function __init__()
    if (Libdl.dlopen_e(libfdasrsf) == C_NULL)
        error("libfdasrsf not properly installed. Run Pkg.build(\"ElasticFDA\")")
    end

    if (Libdl.dlopen_e(libgropt) == C_NULL)
        error("libgropt not properly installed. Run Pkg.build(\"ElasticFDA\")")
    end

    if (Libdl.dlopen_e(libfdaqmap) == C_NULL)
        error("libfdaqmap not properly installed. Run Pkg.build(\"ElasticFDA\")")
    end

end

# CPP.jl workaround not compiling on windows
# the following macro was copied from https://github.com/timholy/Cpp.jl with the
# following modifications:
# 1) Update const for new julia v0.4 warnings
# 2) Formating
# the licensce for the code is included in the package directory
# Useful references:
# http://www.agner.org/optimize/calling_conventions.pdf
# http://mentorembedded.github.io/cxx-abi/abi.html#mangling

const SUBSTITUTION = [""; map(string, collect(0:9)); map(string, collect('A':'Z'))];

# Allow calling of functions in C++ shared libraries.
# Usage: a = @cpp ccall((:mysymbol,mylib),...)
macro cpp(ex)
    # If you get "undefined symbol" errors, use nm or readelf to view
    # the names of the symbols in the library. Then use the resulting
    # information to help improve this macro!
    # Note: for global objects without class or namespace qualifiers
    # (e.g., global consts), you do not need to use this macro (it
    # will work with a plain ccall, even though it is a C++ library)
    msg = "@cpp requires a ccall((:mysymbol, mylib),...)  expression"
    if !isa(ex,Expr) || ex.head != :ccall
        error(msg)
    end

    # Parse the library symbol's name
    exlib = ex.args[1]
    if !(isa(exlib,Expr) && exlib.head == :tuple)
        error(msg)
    end
    sym = exlib.args[1]
    fstr = string(eval(sym))
    #GNU3-4 ABI
    fstr = string("_Z",length(fstr),fstr)

    # Parse the arguments to ccall and construct the parameter type string
    exargtypes = ex.args[3]
    if exargtypes.head != :tuple
        error(msg)
    end
    exargs = exargtypes.args
    pstr = ""
    symtable = (:Void,:Bool,:Cchar,:Char,:String,:Int,:Int8,:Uint8,:Int16,:Uint16,:Int32,:Cint,:Uint32,:Int64,:Uint64,:Float32,:Float64)
    # GNU3-4 ABI v.3 and v.4
    ptable =   ('v', 'b', 'c', 'w', "Pc", 'i', 'a', 'h', 's', 't', 'i', 'i', 'j', 'l', 'm', 'f', 'd')
    msub = String[]
    for iarg = 1:length(exargs)
        thisarg = exargs[iarg]
        thisargm = ""
        while isa(thisarg,Expr) && thisarg.head == :curly && thisarg.args[1] == :Ptr
            thisargm = string(thisargm,'P')
            thisarg = thisarg.args[2]
        end
        matched = false
        for isym = 1:length(symtable)
            if thisarg == symtable[isym]
                matched = true
                thisargm = string(thisargm, ptable[isym])
                # Cchar is a special notation just for name mangling,
                # convert back to :Int8
                if thisarg == :Cchar
                    ex.args[3].args[iarg] = :Int8
                end
                break
            end
        end
        if matched
            if length(thisargm) > 1
                # Use substitution
                idx = indexin([thisargm], msub)[1]
                if idx != 0
                    pstr = string(pstr, "S"*SUBSTITUTION[idx]*"_")
                else
                    push!(msub, thisargm)
                    pstr = string(pstr, thisargm)
                end
            else
                pstr = string(pstr, thisargm)
            end
        else
            println(thisarg)
            error("@cpp: argument not recognized")
        end
    end
    ex.args[1].args[1] = Expr(:quote, Symbol(string(fstr,pstr)))
    ex.args[1].args[2] = esc(ex.args[1].args[2])
    for ival = 4:length(ex.args)
        ex.args[ival] = esc(ex.args[ival])
    end
    ex
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
#include("regression.jl")
#include("regression_funcs.jl")
include("time_warping.jl")
include("image_funcs.jl")
include("image_stats.jl")

end # module
