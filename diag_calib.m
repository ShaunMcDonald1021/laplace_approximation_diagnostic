function [post_mean, w, wce, alph] =...
    diag_calib(gam, lambda, d, v, Us, s_star)

% Calibrate the LA diagnostic in d dimensions, using a multivariate T
% density as the test function. Make sure you download FSKQ prior to
% running this (see https://github.com/tskarvone/fskq).

arguments
    % Hyperparameters gamma and lambda (see Sections 4-5 of manuscript)
    gam (1,1) {mustBePositive}
    lambda (1,1) {mustBePositive}
    % Dimensionality of domain
    d (1,1) {mustBePositive, mustBeInteger}
    % Degrees of freedom for calibration function
    v (1,1) {mustBePositive}
    % The preliminary interrogation grid as a cellular array. Will
    % usually be obtained by running fss_gen (from FSKQ) on a matrix of
    % generator vectors (see demo.m in FSKQ)
    Us cell
    % Preliminary grid in the form of n*d array. Not necessary, but
    % saves time if precomputed
    s_star double = cell2mat(Us)' 
end

% OUTPUTS:
% post_mean: posterior integral mean for diagnostic applied to test fn
% w: BQ weights. Since these only depend on hyperparameters and s_star,
% they can be used for any function to save time
% wce: worst-case error, or unscaled posterior integral standard deviation.
% wce^2*[function-specific factors]/alph^d = posterior integral variance.
% alph: the precision hyperparameter, calculated to put the test fn on the
% boundary of rejection

logf_inter = log(mvtpdf(s_star.*sqrt(v/(v+d)), eye(d), v));
logf_at_mode = (gammaln((v+d)/2) - gammaln(v/2) - d*log(v*pi)/2);
logdet_T = d*log(v/(v+d))/2;
lap_app = exp(logdet_T + logf_at_mode + d*log(2*pi)/2);

% Re-weight interrogations and subtract prior mean interrogation values
Y = exp(d*log(gam) + logdet_T + logf_inter - log(mvnpdf(s_star/gam)))-...
    exp(d*log(2*pi*gam) + logdet_T + logf_at_mode +...
    log(mvnpdf(sqrt(gam^2-1)*s_star/gam)));

% Kernel stuff
k = @(r)exp(-r.^2/(2*lambda^2));
kmean = @(x)(lambda^2/(gam^2+lambda^2))^(d/2)*...
    exp(-norm(x)^2/(2*(gam^2+lambda^2)));
Ikmean = (lambda^2/(2*gam^2 + lambda^2))^(d/2);

[Q, wce, w] = kq_fss(Y, Us, k, kmean, Ikmean);
post_mean = Q + lap_app;

alph = (abs(Q)/(1.96*wce*exp(logdet_T + logf_at_mode)))^(-2/d);
