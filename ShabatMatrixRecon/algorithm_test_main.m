%% Seed RNG
rng(41);
%% Make ground truth
num_cells = 400;
num_genes = 800;
true_zero_prob = 0; 

ground_truth = randi([0,10],num_genes,10)*randi([0,10],10,num_cells); % 400x400 of rank 10
true_zeros = rand(size(ground_truth)) < true_zero_prob; % 0 percent true zeros
ground_truth(true_zeros) = 0;

%% Make zero-inflated "Observation" matrix.
zero_inflate_prob = 0.1;
zero_inflation_mask = rand(size(ground_truth)) < zero_inflate_prob; % 0%  of entries are 0
observed = ground_truth;
observed(zero_inflation_mask) = 0;

%% Reconstruct from "Observation" matrix
num_mask = 2;

lambda_tol = 10;
tol = 1e-8;
N = 100;
    
% reconstructed = masked_reconstruction(observed, num_mask, lambda_tol, tol, N);
% or 
% reconstructed = all_zero_reconstruction(observed, lambda_tol, tol, N);

reconstructed = masked_reconstruction(observed, num_mask, lambda_tol, tol, N);

%% Print performance results
zero_entries = (observed == 0);
recon_mask = logical(zero_entries); 
known_mask = 1 - recon_mask; 

fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(observed)));
fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(putative_array_matrix)));
Diff_sq = abs(reconstructed - ground_truth).^2;

fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq .*known_mask) / sum(known_mask(:)) ));
fprintf('MSE on unknown entries: %g \n',sqrt(sum2(Diff_sq.*recon_mask)/sum(recon_mask(:)) ));

num_clusters = 10;
fprintf('clustering perf of ground truth: %g \n ', clustering_perf(ground_truth, num_clusters))
fprintf('clustering perf of observed: %g \n ', clustering_perf(observed, num_clusters))
fprintf('clustering perf of reconstructed: %g \n ', clustering_perf(reconstructed, num_clusters))
%% Plot results
do_plot = false;

if do_plot
    figure(1)
    hold on
    histogram(Diff_sq .* recon_mask,1000); 

    figure(2)
    hold on
    histogram(putative_array_matrix .* known_mask,100);


    figure(3)
    hold on
    histogram(Diff_sq .* known_mask,100);
end

%% Helper Functions


% Computes the average Silhoutte coeficient on k clustering (clustering
% found by k means)
% X (n x p matrix) each column represents a single point
% k number of clusters
% Closer to 1 the better.
% Output seems to be randomized based on how kmeans operates
% The variance in the output should decrease for clustered data
function perf = clustering_perf(X, k)
    num_tries = 3; % Takes best out of 5 tries
    perf = 0;
    for i = 1:num_tries
        [idx, C, sumd] = kmeans(X',k);
        S = silhouette(X',idx,'Euclidean');
        trial_perf = mean(S);
        if trial_perf > perf;
            perf = trial_perf;
        end
    end
end

% Given an observed matrix, each zero entry has a chance of being set to
% "unknown" element. These unknown elements computed by minimized nuclear
% norm * const + l2 norm
function reconstructed = masked_reconstruction(observed, num_masks, lambda_tol, tol, N)
    %% Collect the set of zero entries
    zero_entries = (observed == 0);
    zero_inds = find(zero_entries);
    num_zeros = length(zero_inds);
    [num_genes, num_cells] = size(observed);

    %% Make masks
    theoretical_limit = 0.75; % Fix this to be something from CBT paper
    num_non_observed = min(floor(num_zeros * 0.9),floor(theoretical_limit*numel(observed)));

    masks = cell(num_masks,1);

    % Tried using all the zeros as the "unobserved population"
    for i = 1:num_masks
        masks{i} = zero_inds(randperm(num_zeros,num_non_observed));
    %     masks{i} = find(zero_inflation_mask); 
    end

    %% Reconstruct matrix with each mask (parallel step)
    recon_array = cell(num_masks,1);

    parfor j = 1:num_masks
        mask_matrix = zeros(num_genes,num_cells);
        mask_matrix(masks{j}) = 1;

        observed_entries = logical(1 - mask_matrix);

        [CompletedMat, ier] = MatrixCompletion(...
            observed, observed_entries,N, 'nuclear', lambda_tol, tol, 0);

        recon_array{j} = CompletedMat;
    end

    %% Arrange reconstructed matrices into an array for easy median-computation
    recon_array_arrangement = zeros(num_masks, num_genes * num_cells); 
    for i = 1:num_masks
        recon_array_arrangement(i,:) = recon_array{i}(:);
    end 

    %% Form putative matrix 
    reconstructed = zeros(num_genes,num_cells); 
    reconstructed(:) = median(recon_array_arrangement,1); 

end

% Given an observed matrix, every zero entry is set as an "unknown" 
% element. These unknown elements computed by minimized nuclear
% norm * const + l2 norm
function reconstructed = all_zero_reconstruction(observed, lambda_tol, tol, N)
    [reconstructed, ier] = MatrixCompletion(...
            observed, observed ~= 0,N, 'nuclear', lambda_tol, tol, 0);
end