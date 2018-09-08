@recipe function f(h::fdawarp)
    M = length(h.time)
    x = h.time

    mean_f0 = mean(h.f,dims=2)
    std_f0 = std(h.f, dims=2)
    mean_fn = mean(h.fn, dims=2)
    std_fn = std(h.fn, dims=2)

    stat0 = [mean_f0, mean_f0+std_f0, mean_f0-std_f0]
    statn = [mean_fn, mean_fn+std_fn, mean_fn-std_fn]

    if h.method=="mean"
        type_lbl = "fmean"
    else
        type_lbl = "median"
    end

    # set up the subplots
    legend := false
    layout := (2,3)

    @series begin
        subplot := 1
        title := "Original Functions"
        x, h.f
    end

    @series begin
        subplot := 2
        title := "Aligned Functions"
        x, h.fn
    end

    @series begin
        subplot := 3
        title := "Warping Functions"
        aspect_ratio := :equal
        LinRange(0,1,length(x)), h.gam
    end

    @series begin
        subplot := 4
        legend := :bottom
        title := "Original Data: Mean +/- STD"
        label := ["Mean","Mean + STD","Mean - STD"]
        x, stat0
    end

    @series begin
        subplot := 5
        legend := :bottom
        title := "Warped Data: Mean +/- STD"
        label := ["Mean","Mean + STD","Mean - STD"]
        x, statn
    end

    @series begin
        subplot := 6
        title := type_lbl
        x, h.fmean
    end

end
