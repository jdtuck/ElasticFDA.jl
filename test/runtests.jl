using ElasticFDA
using Distributions
using JLD
using Base.Test

d = load("simu_data.jld");
f = d["f"];
timet = collect(linspace(0,1,size(f,1)));
f1 = f[:,1];
f2 = f[:,2];

# test smooth_data
g = smooth_data(f1,1);
@test sum(g)-sum(f1)<1e-6
smooth_data!(f1,1);
@test sum(f1)-sum(f[:,1])<1e-6

# test trapz
d = trapz(timet, timet);
@test_approx_eq(d,0.5)

# test f_to_srsf
q1 = f_to_srsf(f1, timet);
f1a = srsf_to_f(q1,timet,f1[1]);

# test optimum reparam
gam = optimum_reparam(q1,timet,q1,f1o=f1[1],f2o=f2[1]);
@test norm(gam-linspace(0,1,101)) < 1e-10

# test warping functions
qw = warp_q_gamma(timet, q1, gam);
fw = warp_f_gamma(timet, f1, gam);
@test norm(q1-qw) < 1e-10
@test norm(f1-fw) < 1e-10

# test srsf_align
out = srsf_align(f, timet);

# test align_fPCA
out1 = align_fPCA(f, timet, MaxItr=2);

# test vert_fPCA
out1 = vert_fPCA(out["qn"],timet,out["qn"]);

# test horiz_fPCA
out1 = horiz_fPCA(out["gam"], timet);

# test gauss_model
out1 = gauss_model(out["fn"], timet, out["qn"], out["gam"]);

# test pair_warping_bayes
out1 = pair_warping_bayes(f1[1:100], f2[1:100], iter=2);

# test group_warping_bayes
out1 = group_warping_bayes(f, iter=2);

# test elastic_distance
da, dp = elastic_distance(f1, f1, timet);
@test da < 1e-10
@test dp < 1e-6

# test rgam
gam = rgam(101,.1,10);

# test elastic_regression
include("test_warp_regress.jl")
timet = collect(timet);
out = elastic_regression(f, y_orig, timet, max_itr=1);
out1 = elastic_prediction(f, timet, out, y=y_orig);

# test elastic logistic regression
include("test_warp_logistic.jl")
timet = collect(timet);
out = elastic_logistic(f, y_orig, timet, max_itr=1);
out1 = elastic_prediction(f, timet, out, y=y_orig);

# test elastic m-logistic regression
include("test_warp_mlogistic.jl")
timet = collect(timet);
out = elastic_mlogistic(f, y_orig, timet, max_itr=1);
out1 = elastic_prediction(f, timet, out, y=y_orig);

d = load("beta.jld");
beta = d["beta"];

# test curve_to_q
q1 = curve_to_q(beta[:,:,1]);
beta1 = q_to_curve(q1);

# test calc_shape_dist
beta1 = beta[:,:,1];
d1 = calc_shape_dist(beta1,beta1);

# test curve_pair_align
beta2n, q2n, gam, q1 = curve_pair_align(beta1,beta1);

# test curve_geodesic
geod, geod_q = curve_geodesic(beta1,beta[:,:,2]);

# test curve functions 'O'
mu, betamean, v, q = curve_karcher_mean(beta, maxit=1);
betan, qn, betamean1, q_mu = curve_srvf_align(beta, maxit=1);
K = curve_karcher_cov(betamean, beta);
pd = curve_principal_directions(betamean, mu, K);
s = sample_shapes(mu, K);

# test curve functions 'C'
mu, betamean, v, q = curve_karcher_mean(beta, mode='C', maxit=1);
betan, qn, betamean1, q_mu = curve_srvf_align(beta, mode='C', maxit=1);
K = curve_karcher_cov(betamean, beta, mode='C');
pd = curve_principal_directions(betamean, mu, K, mode='C');
s = sample_shapes(mu, K, mode='C');

# test image functions
vars = load("image.jld");
I1 = vars["I1"];
I2 = vars["I2"];
Inew, gam = pair_align_image(I1,I2,itermax=1);
