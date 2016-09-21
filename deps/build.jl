if is_unix()
    cd(joinpath(dirname(@__FILE__), "src", "fdasrsf"))
    suffix = is_apple() ? "dylib" : "so"
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
    cd(joinpath(dirname(@__FILE__), "src", "gropt"))
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
    cd(joinpath(dirname(@__FILE__), "src", "fdaqmap"))
    run(`make clean`)
    run(`make SUFFIX=$suffix`)
else # Windows
    # these binaries were cross-compiled from Cygwin for x86_64 only using
    # the Makefile_win in the corresponding src directories and the windows bin
    # release of openblas
    url1 = "https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.4.0/fdaqmap.7z"
    url2 = "https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.4.0/fdasrsf.7z"
    url3 = "https://github.com/jdtuck/ElasticFDA.jl/releases/download/v0.4.0/gropt.7z"
    try
        run(`curl -LO $url1`)
        run(`curl -LO $url2`)
        run(`curl -LO $url3`)
    catch
        run(`powershell -Command "(new-object net.webclient).DownloadFile(\"$url1\", \"fdaqmap.7z\")"`)
        run(`powershell -Command "(new-object net.webclient).DownloadFile(\"$url2\", \"fdasrsf.7z\")"`)
        run(`powershell -Command "(new-object net.webclient).DownloadFile(\"$url3\", \"gropt.7z\")"`)
    end
    run(`7z x -y fdaqmap.7z`)
    run(`7z x -y fdasrsf.7z`)
    run(`7z x -y gropt.7z`)
end
