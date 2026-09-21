function r = corr_columns(Y, Yhat)
%CORR_COLUMNS Pearson correlation between corresponding columns of two
%   matrices, vectorized (no per-column corrcoef calls).
%   Y, Yhat: [n_features x n_columns] (e.g. parcels x subjects)
%   r: [1 x n_columns] correlation per column

Yc = Y - mean(Y, 1);
Yhatc = Yhat - mean(Yhat, 1);
numer = sum(Yc .* Yhatc, 1);
denom = sqrt(sum(Yc.^2, 1) .* sum(Yhatc.^2, 1));
r = numer ./ denom;
end
