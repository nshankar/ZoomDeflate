

% dropoff_rate ratio of values that gets set to 0
function normalized = normalize_reconstruction(observed, reconstructed, dropoff_rate)
    [r, c] = size(observed);
    normalized = observed;
    for i = 1:r
        row = observed(i, 1:c);
        num_zeros = sum(row == 0); % num zeros in row i
        num_true_zeros = max(num_zeros - c * dropoff_rate, 0) / (1 - dropoff_rate);
        observed_dist = row(row ~= 0);
        % add in zeros so that we observed_dist samples (1 - dropoff_rate)
        % of each element from ground truth (including the true zeros)
        % observed_dist = [observed_dist zeros(1, floor(num_true_zeros * (1 - dropoff_rate)))];
        normalized(i, row==0) = match_dist(observed_dist, reconstructed(i, row ==0));
    end
end