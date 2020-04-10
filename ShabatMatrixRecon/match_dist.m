% transforms the values so that distribution of the value is similar to
% that of the desired distribution
function normalized = match_dist(desired_dist, values)
    desired_dist = sort(desired_dist);
    num_desired = length(desired_dist);
    [v, index] = sort(values);
    num_values = length(values);
    normalized = zeros(1, num_values);
    for i = 1:num_values
        desired_idx = (i - 1) / (num_values - 1) * (num_desired - 1);
        desired_idx = desired_idx + 1;
        
        idx = floor(desired_idx);
        offset = desired_idx - idx;
        if offset == 0
            normalized(i) = desired_dist(idx);
        else
            normalized(i) = desired_dist(idx) * (1 - offset) + desired_dist(idx + 1) * offset;
        end
    end
end