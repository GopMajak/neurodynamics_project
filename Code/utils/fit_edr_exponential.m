function [A, lambda] = fit_edr_exponential(bin_centers, bin_means, exclude_bins)
%FIT_EDR_EXPONENTIAL Fit y = A*exp(-lambda*x) to binned connectome weight
%   vs. Euclidean distance, following Vohryzek et al. 2025 Methods
%   ("EDR" section): 400 equal-width bins spanning 10-170 mm, excluding
%   the first 25 (shortest-distance) bins from the fit.
%
%   BIN_CENTERS, BIN_MEANS : from COMPUTE_DISTANCE_BIN_STATS
%   EXCLUDE_BINS            : number of low-distance bins to drop (25)
%
%   Returns fitted amplitude A and decay LAMBDA (mm^-1). Paper reports
%   A = 0.066, lambda = 0.162 mm^-1 for the structural connectome fit.

n_bins = numel(bin_centers);
fit_idx = (exclude_bins + 1):n_bins;
x = bin_centers(fit_idx)';
y = bin_means(fit_idx)';

valid = ~isnan(y) & y > 0;
x = x(valid);
y = y(valid);

% Initial guess via log-linear least squares: log(y) = log(A) - lambda*x
p_lin = polyfit(x, log(y), 1);
lambda0 = -p_lin(1);
A0 = exp(p_lin(2));

model = @(p, x) p(1) * exp(-p(2) * x);

try
    opts = optimoptions('lsqcurvefit', 'Display', 'off');
    p_fit = lsqcurvefit(model, [A0, lambda0], x, y, [0, 0], [Inf, Inf], opts);
    A = p_fit(1);
    lambda = p_fit(2);
catch
    % Optimization Toolbox unavailable: fall back to the log-linear fit.
    A = A0;
    lambda = lambda0;
end
end
