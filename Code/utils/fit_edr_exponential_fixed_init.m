function [A, lambda] = fit_edr_exponential_fixed_init(bin_centers, bin_means, fit_start_bin, A0, bounds)
%FIT_EDR_EXPONENTIAL_FIXED_INIT Fit y = A*exp(-lambda*x) with an explicit
%   initial guess and bounds, following the exact optimization setup
%   found in Vohryzek et al.'s longrange_derivation_project_laplacian_v2.m
%   (lsqcurvefit, A0=[0.15 0.18], bounds [-100 100]).
%
%   NEW FILE -- does not modify FIT_EDR_EXPONENTIAL.M (Stage 1 original,
%   which used a log-linear initial guess instead of a fixed one). Both
%   are preserved for comparison, since nonlinear exponential fits on
%   noisy real data are known to be sensitive to starting point.

if nargin < 3 || isempty(fit_start_bin)
    fit_start_bin = 25;
end
if nargin < 4 || isempty(A0)
    A0 = [0.15, 0.18];
end
if nargin < 5 || isempty(bounds)
    bounds = [-100, 100];
end

bin_centers = double(bin_centers(:));
bin_means = double(bin_means(:));
x = bin_centers(fit_start_bin:end);
y = bin_means(fit_start_bin:end);
valid = ~isnan(y);
x = x(valid);
y = y(valid);

model = @(A, x) A(1) * exp(-A(2) * x);
opts = optimset('MaxFunEvals', 10000, 'MaxIter', 1000, 'Display', 'off');

Afit = lsqcurvefit(model, A0, x, y, [bounds(1) bounds(1)], [bounds(2) bounds(2)], opts);

A = Afit(1);
lambda = Afit(2);
end
