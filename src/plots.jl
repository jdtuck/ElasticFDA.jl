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

@recipe function f(h::vfpca; no=1:3)
    x = h.time
    M, N = size(h.q_pca)

    # set up the subplots
    legend := false
    layout := (2,3)

    @series begin
        subplot := 1
        title := @sprintf "f domain: PD %d" no[1]
        x, h.f_pca[:,:,no[1]]
    end

    @series begin
        subplot := 2
        title := @sprintf "f domain: PD %d" no[2]
        x, h.f_pca[:,:,no[2]]
    end

    @series begin
        subplot := 3
        title := @sprintf "f domain: PD %d" no[3]
        x, h.f_pca[:,:,no[3]]
    end

    @series begin
        subplot := 4
        title := @sprintf "q domain: PD %d" no[1]
        x, h.q_pca[1:(M-1),:,no[1]]
    end

    @series begin
        subplot := 5
        title := @sprintf "q domain: PD %d" no[2]
        x, h.q_pca[1:(M-1),:,no[2]]
    end

    @series begin
        subplot := 6
        title := @sprintf "q domain: PD %d" no[3]
        x, h.q_pca[1:(M-1),:,no[3]]
    end

end

@recipe function f(h::hfpca, no=1:3)
    M, N, no1 = size(h.gam_pca)
    x = LinRange(0,1,M)

    # set up the subplots
    legend := false
    layout := (1,3)
    aspect_ratio := :equal

    @series begin
        subplot := 1
        title := @sprintf "gam: PD %d" no[1]
        x, h.gam_pca[:,:,no[1]]
    end

    @series begin
        subplot := 2
        title := @sprintf "gam: PD %d" no[2]
        x, h.gam_pca[:,:,no[2]]
    end

    @series begin
        subplot := 3
        title := @sprintf "gam: PD %d" no[3]
        aspect_ratio := :equal
        x, h.gam_pca[:,:,no[3]]
    end

end
