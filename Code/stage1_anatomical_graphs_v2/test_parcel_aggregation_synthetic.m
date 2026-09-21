%% Smoke test for the parcel-pair aggregation logic (small synthetic data)
clear; clc;

rng(9);
n = 200;
num_parcels = 10;
parc = randi(num_parcels, n, 1); % random parcel assignment, 1..10

% synthetic "exception" pairs: concentrate them in just 2 parcel pairs
% to check the aggregation correctly reports a SMALL number of hot pairs
is_exception = false(n, n);
target_pairs = [1 2; 3 4]; % force exceptions only between parcel(1,2) and (3,4)
for tp = 1:size(target_pairs, 1)
    idx_a = find(parc == target_pairs(tp,1));
    idx_b = find(parc == target_pairs(tp,2));
    for a = idx_a'
        for b = idx_b'
            if rand < 0.3
                is_exception(a,b) = true;
            end
        end
    end
end

label_to_idx = 1:num_parcels; % identity here since labels are already 1..10
[rr, cc] = find(is_exception);
pi_idx = label_to_idx(parc(rr))';
pj_idx = label_to_idx(parc(cc))';
pair_idx = sub2ind([num_parcels, num_parcels], min(pi_idx,pj_idx), max(pi_idx,pj_idx));
counts = accumarray(pair_idx, 1, [num_parcels*num_parcels, 1]);
parcel_exception_count = reshape(counts, num_parcels, num_parcels);

n_possible = num_parcels*(num_parcels-1)/2;
n_with_exception = nnz(triu(parcel_exception_count,1) > 0);

fprintf('Injected exceptions concentrated in parcel pairs (1,2) and (3,4) only.\n');
fprintf('Possible parcel pairs: %d\n', n_possible);
fprintf('Parcel pairs with >=1 exception found: %d (expected: 2)\n', n_with_exception);
assert(n_with_exception == 2, 'expected exactly 2 hot parcel pairs');
assert(parcel_exception_count(1,2) > 0 || parcel_exception_count(2,1) > 0, 'pair (1,2) should be flagged');
assert(parcel_exception_count(3,4) > 0 || parcel_exception_count(4,3) > 0, 'pair (3,4) should be flagged');
fprintf('PASS: aggregation logic correctly isolates the 2 hot parcel pairs.\n');
