%% Seed RNG
rng(41);
%% Make ground truth
% num_cells = 24;
% num_genes = 8;
% true_zero_prob = 0; 
% 
% ground_truth = randi([0,10],num_genes,10)*randi([0,10],10,num_cells); % 400x400 of rank 10
% true_zeros = rand(size(ground_truth)) < true_zero_prob; % 0 percent true zeros
% ground_truth(true_zeros) = 0;
% 
%% Make zero-inflated "Observation" matrix.
% zero_inflate_prob = 0.1;
% zero_inflation_mask = rand(size(ground_truth)) < zero_inflate_prob; % 0%  of entries are 0
% observed = ground_truth;
% observed(zero_inflation_mask) = 0;

%% Setup for plots
rose_gold = [246/255,170/255,180/255];
sea_foam = [0 0.925 0.7];
burgundy = [0.502,0,0.125];
gold = [255/255 215/255 0];
make_figs = true;
save_figs = true;
formats = {'png','epsc','fig'};
% 
%% messing around with filepaths 
% put the path of your github repo here. End it with a slash
filepath = '../'; 

filepath_in = [filepath,'SplatGenData/one_cell_types_50_sparse/'];
filepath_in = PathSlashCorrector(filepath_in);

filepath_out = SubfolderMaker(filepath,'Figures/one_cell_types_50_sparse/'); 
filepath_out = PathSlashCorrector(filepath_out); 

%% Read in file 
% (1,1) allows csvread to skip the 0th row and 0th column (labels),
% start with numerical values
ground_truth = csvread([filepath_in,' true_counts.csv'],1,1);
observed = csvread([filepath_in,' counts.csv'],1,1);

[num_genes, num_cells] = size(observed);

%% look at some gene distributions
% figure(3)
% hold on
% for i=1:5:26
%     histogram(observed(i,:),'LineStyle','none','FaceAlpha',0.5)
% end 
%% Preprocessing step 
scaling_flag = 'std';   % allowable values are 'none','std','noise'
[observed,scalar_multiples] = preprocess(observed,scaling_flag);
ground_truth = scalar_multiples * ground_truth; 

%% look at some gene distributions
% figure(4)
% hold on
% for i=1:5:26
%     histogram(observed(i,:),'LineStyle','none','FaceAlpha',0.5)
% end 


%% Reconstruct from "Observation" matrix
num_mask = 20;

lambda_tol = 10;
tol = 1e-8;
N = 100;
    
% reconstructed = masked_reconstruction(observed, num_mask, lambda_tol, tol, N);
% or 
% reconstructed = all_zero_reconstruction(observed, lambda_tol, tol, N);

recon_array_arrangement = masked_reconstructions(observed, num_mask, lambda_tol, tol, N);

%% 
zero_entries = (observed == 0);
recon_mask = logical(zero_entries); 
known_mask = logical(1 - recon_mask); 

max_proj = max(recon_array_arrangement,1);
min_proj = min(recon_array_arrangement,1);

est_ranges = max_proj - min_proj; 

% figure,histogram(est_ranges(recon_mask),20000,'LineStyle','none');
%% Compute RMSE for different quantiles
num_quantiles = 20;
quantile_array = linspace(0.05,1,20);
RMSE_array = zeros(1,num_quantiles);

for i = 1:num_quantiles
    test_X = zeros(num_genes,num_cells);
    test_X(:) = quantile(recon_array_arrangement,quantile_array(i),1); 
    RMSE_array(i) = RMSE(test_X,ground_truth,recon_mask); 
end
    
%% 
reconstructed = zeros(num_genes,num_cells);
reconstructed(:) = median(recon_array_arrangement,1);
%% Print performance results

% 
% fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(observed)));
% fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(reconstructed)));
Diff_sq = abs(reconstructed - ground_truth).^2;
% 
% fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq .*known_mask) / sum(known_mask(:)) ));
% fprintf('MSE on unknown entries: %g \n',sqrt(sum2(Diff_sq.*recon_mask)/sum(recon_mask(:)) ));
% 
% num_clusters = 10;
% fprintf('clustering perf of ground truth: %g \n ', clustering_perf(ground_truth, num_clusters))
% fprintf('clustering perf of observed: %g \n ', clustering_perf(observed, num_clusters))
% fprintf('clustering perf of reconstructed: %g \n ', clustering_perf(max_proj, num_clusters))

%% Plot results
num_bins = 1000;

if make_figs
    % This plot is derpy/dumb but I wanted to give an example of plotting
    % two histograms on the same plot 
    h= figure(1);
    hold on
    set(gca,'FontName','Helvetica Neue','FontSize',20,'FontWeight','Bold')
    set(gca,'TickDir','out');
    histogram(Diff_sq(recon_mask),num_bins,'FaceColor',rose_gold,'LineStyle','none','FaceAlpha',0.7);
%     histogram(Diff_sq(known_mask),num_bins,'FaceColor',sea_foam,'LineStyle','none','FaceAlpha',0.7);
    xlabel('squared diff (ts^2)')
    ylabel('frequency')
    legend('Unknown Values')%,'Known values')
    title('Error of reconstructed matrix')
    
    % save plots in various figure formats
    if save_figs
        for k = 1:length(formats)
            save_name = [filepath_out,'hist_of_recon_diffs'];
            saveas(h,save_name,char(formats(k)));
        end
%         close all
    end
    
    
    h2 = figure(2); 
    hold on
    set(gca,'FontName','Helvetica Neue','FontSize',20,'FontWeight','Bold')
    set(gca,'TickDir','out');
    plot(quantile_array,RMSE_array,'o','MarkerSize',12,'LineWidth',2) 
    xlabel('Quantile for matrix entries')
    ylabel('RMSE')
    
    
    % save plots in various figure formats
    if save_figs
        for k = 1:length(formats)
            save_name = [filepath_out,'threshold_vs_RMSE'];
            saveas(h2,save_name,char(formats(k)));
        end
%         close all
    end
%     ylim([0,2000])
        %
    %     figure(2)
    %     hold on
    %     histogram(putative_array_matrix .* known_mask,100);
    %
    %
    %     figure(3)
    %     hold on
    %     histogram(Diff_sq .* known_mask,100);
    
end

%% Helper Functions

% computes RMSE of entries included in the mask. mask should be logical.
function [RMSE_value,scalar_multiples] = RMSE(X,truth,mask)
    squared_diff = (X - truth).^2;
    RMSE_value = sqrt(sum2(squared_diff(mask)) / sum(mask(:)));
end

function X_proc = preprocess(X,scaling_flag)
%Want to take stds only over the nonzero entries (NOT  row_stds =
%std(X,0,2);)  
scalar_multiples = eye(num_genes);
if strcmp(scaling_flag,'std') || strcmp(scaling_flag,'noise')
    num_genes = size(X,1); 
    row_stds = zeros(num_genes,1);
    for i = 1:num_genes
        nonzero_entries = double(~(X(i,:)==0));
        row_stds(i) = std(X(i,:),nonzero_entries,2); 
    end 
end
if strcmp(scaling_flag,'none')
    X_proc = X; 
    return
elseif strcmp(scaling_flag,'std')
    scalar_multiples = diag( 1 ./ row_stds); 
    X_proc = scalar_multiples * X; 
    return
elseif strcmp(scaling_flag,'noise')
    row_means = zeros(num_genes,1);
    for i = 1:num_genes
        row_means(i) = mean(X(i,(X(i,:)~=0)));
    end 
    row_noise = (row_stds ./ row_means).^2; 
    scalar_multiples = diag(1 ./ row_noise);
    X_proc =  scalar_multiples * X; 
    return
else
    fprintf('ERROR: INCORRECT FLAG IN PREPROCESSING. NO PROCESSING FOR YOU.')
    X_proc = X; 
end
end


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
function recon_array_arrangement = masked_reconstructions(observed, num_masks, lambda_tol, tol, N)
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
%         masks{i} = zero_inds; 
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
%     reconstructed = zeros(num_genes,num_cells); 
%     reconstructed(:) = median(recon_array_arrangement,1); 

end

% Given an observed matrix, every zero entry is set as an "unknown" 
% element. These unknown elements computed by minimized nuclear
% norm * const + l2 norm
function reconstructed = all_zero_reconstruction(observed, lambda_tol, tol, N)
    [reconstructed, ier] = MatrixCompletion(...
            observed, observed ~= 0,N, 'nuclear', lambda_tol, tol, 0);
end