function plot_fdawarp(x::fdawarp)
    M = length(x.time)

    mean_f0 = mean(x.f,dims=2)
    std_f0 = std(x.f, dims=2)
    mean_fn = mean(x.fn, dims=2)
    std_fn = std(x.fn, dims=2)

    plot(x.time,x.f,title="Original Functions",legend=:none)

    plot(x.time,x.fn,title="Warped Functions",legend=:none)

    plot(LinRange(0,1,M),x.gam,title="Warping Functions",legend=:none,
         aspect_ratio=:equal)

    plot(x.time,[mean_f0, mean_f0+std_f0, mean_f0-std_f0], label=["Mean","Mean + STD","Mean - STD"],
         title="Original Data: Mean +/- STD")

    plot(x.time,[mean_fn, mean_fn+std_fn, mean_fn-std_fn], label=["Mean","Mean + STD","Mean - STD"],
         title="Warped Data: Mean +/- STD")

    if x.method=="mean"
        plot(x.time,x.fmean,title="fmean",legend=:none)
    else
        plot(x.time,x.fmean,title="fmedian",legend=:none)
    end
end
