@unix_only begin
    cd(joinpath(dirname(@__FILE__), "src", "fdasrsf"))
    suffix = @osx? "dylib" : "so"
    run(`make SUFFIX=$suffix`)
    cd(joinpath(dirname(@__FILE__), "src", "gropt"))
    run(`make SUFFIX=$suffix`)
    cd(joinpath(dirname(@__FILE__), "src", "fdaqmap"))
    run(`make SUFFIC=$suffix`)
end

@windows_only begin
    # these binaries were cross-compiled from Cygwin for x86_64 only using
    # the Makefile_win in the corresponding src directories and the windows bin
    # release of openblas
    run(`curl -LO https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.3.0/gropt.7z`)
    run(`curl -LO https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.3.0/fdasrsf.7z`)
    run(`curl -LO https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.3.0/fdaqmap.7z`)
    run(`7z x -y fdasrsf.7z`)
    run(`7z x -y gropt.7z`)
    run(`7z x -y fdaqmap.7z`)
end
