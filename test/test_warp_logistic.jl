timet = LinRange(0,1,101);
M = length(timet);
N = 5;
lam = 0.001;
center = [.3; .35];
center2 = [4; 4];
sd1 = .1;
gam_sd = 1;
num_comp = 5;
f_orig = zeros(M,N*length(center));
omega = 2*pi;
cnt = 1;
for ii = 1:length(center)
    global cnt
    tmp = Normal(center[ii], .1);
    tmp2 = Normal(center2[ii], sd1);
    for jj = 1:N
        f_orig[:, cnt] = rand(tmp2) * pdf.(tmp, timet);
        cnt+= 1;
    end
end

q_orig = f_to_srsf(f_orig,collect(timet));
y_orig = ones(size(q_orig,2));
y_orig[N+1:end] .= -1;

f = zeros(M, size(f_orig,2));
q = zeros(M, size(f_orig,2));
gam_orig = rgam(M, gam_sd, 2*N);
cnt = 1;
for ii = 1:length(center)
    global cnt
    for jj = 1:N
        f[:,cnt] = warp_f_gamma(collect(timet), f_orig[:,cnt], gam_orig[:,jj]);
        q[:,cnt] = f_to_srsf(f[:,cnt], collect(timet));
        cnt += 1;
    end
end
