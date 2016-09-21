if is_unix()
    cd(joinpath(dirname(@__FILE__), "src", "fdasrsf"))
    suffix = @osx? "dylib" : "so"
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
    cd(joinpath(dirname(@__FILE__), "src", "gropt"))
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
    cd(joinpath(dirname(@__FILE__), "src", "fdaqmap"))
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
else # Windows
    using BinDeps
    # these binaries were cross-compiled from Cygwin for x86_64 only using
    # the Makefile_win in the corresponding src directories and the windows bin
    # release of openblas
    BinDeps.download_cmd("https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.4.0/fdaqmap.7z")
    BinDeps.download_cmd("https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.4.0/fdasrsf.7z")
    BinDeps.download_cmd("https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.4.0/gropt.7z")
    run(`7z x -y fdasrsf.7z`)
    run(`7z x -y gropt.7z`)
    run(`7z x -y fdaqmap.7z`)
end
