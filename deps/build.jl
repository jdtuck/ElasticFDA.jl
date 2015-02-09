@unix_only begin
    cd(joinpath(dirname(@__FILE__), "src", "fdasrsf"))
    suffix = @osx? "dylib" : "so"
    run(`make CC=icc SUFFIX=$suffix`)
end
