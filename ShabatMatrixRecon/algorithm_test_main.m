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
zero_inflate_prob = 0.4;
zero_inflation_mask = rand(size(ground_truth)) < zero_inflate_prob; % 0%  of entries are 0
observed = ground_truth;
observed(zero_inflation_mask) = 0;

%% Collect the set of zero entries
zero_entries = (observed == 0);
zero_inds = find(zero_entries);
num_zeros = length(zero_inds);

%% Make masks
num_masks = 2;

theoretical_limit = 0.75; % Fix this to be something from CBT paper
num_non_observed = min(floor(num_zeros * 0.9),floor(theoretical_limit*numel(ground_truth)));

masks = cell(num_masks,1);

% Tried using all the zeros as the "unobserved population"
for i = 1:num_masks
    masks{i} = zero_inds(randi(num_zeros,num_non_observed,1));
%     masks{i} = find(zero_inflation_mask); 
end

%% Reconstruct matrix with each mask (parallel step)
recon_array = cell(num_masks,1);

lamnbda_tol = 10;
tol = 1e-8;
N = 100;
parfor j = 1:num_masks
    mask_matrix = zeros(num_genes,num_cells);
    mask_matrix(masks{j}) = 1;
    
    observed_entries = logical(1 - mask_matrix);
    
    [CompletedMat, ier] = MatrixCompletion(...
        observed, observed_entries,N, 'nuclear', lamnbda_tol, tol, 0);
    
    recon_array{j} = CompletedMat;
end

%% Arrange reconstructed matrices into an array for easy median-computation
recon_array_arrangement = zeros(num_masks, num_genes * num_cells); 
for i = 1:num_masks
    recon_array_arrangement(i,:) = recon_array{i}(:);
end 

%% Form putative matrix 
putative_array_matrix = zeros(num_genes,num_cells); 
putative_array_matrix(:) = median(recon_array_arrangement,1); 

%% 
recon_mask = logical(zero_entries); 
known_mask = 1 - recon_mask; 

fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(observed)));
fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(putative_array_matrix)));
Diff_sq = abs(putative_array_matrix - ground_truth).^2;

fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq .*known_mask) / sum(known_mask(:)) ));
fprintf('MSE on unknown entries: %g \n',sqrt(sum2(Diff_sq.*recon_mask)/sum(recon_mask(:)) ));
%% 
figure(1)
hold on
histogram(Diff_sq .* recon_mask,1000); 

figure(2)
hold on
histogram(putative_array_matrix .* known_mask,100);

%%

figure(3)
hold on
histogram(Diff_sq .* known_mask,100);