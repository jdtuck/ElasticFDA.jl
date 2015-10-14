timet = linspace(0,1,101);
M = length(timet);
N = 30;
lam = 0.001;
center = [.35; .5; .65];
center2 = [4; 3.7; 4];
sd1 = .05;
gam_sd = 2;
num_comp = 5;
f_orig = zeros(M,N*length(center));
omega = 2*pi;
cnt = 1;
for ii = 1:length(center)
    tmp = Normal(center[ii], .075);
    tmp2 = Normal(center2[ii], sd1);
    for jj = 1:N
        f_orig[:, cnt] = rand(tmp2) * pdf(tmp, timet);
        cnt += 1;
    end
end

q_orig = f_to_srsf(f_orig,collect(timet));
y_orig = zeros(size(q_orig,2));
alpha_t = 0;
b1 = sin(omega*timet);
b2 = cos(omega*timet);
B = hcat(b1,b2);
bt = .5 * b1 + 0.9 * b2;
cnt = 1;
for ii = 1:length(center)
    tmp = Normal(center[ii], .075);
    d = Normal(center2[ii], sd1);
    d1 = Normal(0, 0.01);
    tmp2 = f_to_srsf(rand(d) * pdf(tmp, timet), collect(timet));
    for jj = 1:N
        y_orig[cnt] = alpha_t + trapz(collect(timet), tmp2 .* bt) + rand(d);
        cnt += 1;
    end
end

f = zeros(M, size(f_orig,2));
q = zeros(M, size(f_orig,2));
cnt = 1;
for ii = 1:length(center)
    gam_orig = rgam(M, gam_sd, N);
    for jj = 1:N
        f[:,cnt] = warp_f_gamma(collect(timet), f_orig[:,cnt], gam_orig[:,jj]);
        q[:,cnt] = f_to_srsf(f[:,cnt], collect(timet));
        cnt += 1;
    end
end




