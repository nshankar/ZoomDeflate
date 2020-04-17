function M_approx = kluger_low_rank_approx(M)
%% TEST
num_cells = 200;
num_genes = 200;

M = randi([0,10],num_genes,10)*randi([0,10],10,num_cells); % 400x400 of rank 10
%% Normalize
col_sums = sum(M,1); 
M_approx = 1e4 .* M ./ col_sums;

[u,s,v] = rsvd(M,100); 
% s = s(:); 

diffs = s(1:99) - s(2:100); 

last_21 = diffs(end-20:end); 
mean_diff = mean(last_21); 
std_diff = std(last_21); 

%%
% y = normpdf(diffs,