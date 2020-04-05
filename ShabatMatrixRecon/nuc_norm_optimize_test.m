%% Create synthetic data 
% clc

%seed rng
rng(37); 

frac_hidden = 1/10; 

ground_truth = randn(200,10)*randn(10,200); % 200x200 of rank 10
mask = rand(size(ground_truth)) < frac_hidden; % frac_hidden of values are 1 in this matrix 

measured = ground_truth; 
measured(mask) = 0; % frac_hidden values sent to 0. 

% number of hidden entries  
number_hidden = sum(mask(:))
%% Use spy to see which values are hidden
% figure(1)
% hold on
% spy(ground_truth)
% 
% figure(2)
% hold on
% spy(mask)
% 
% figure(3)
% hold on
% spy(measured) 

%% Fit imputed values to synthetic data by minimizing nuclear norm 
first_guesses = rand(number_hidden,1); 

options = optimset('MaxFunEvals',10000*number_hidden,'TolFun',1e-12); 

[imputed_vals, GoF] = fminsearch(@(vals_fitted) nuc_norm(...
    measured,mask,vals_fitted),first_guesses,...
   options);

%% 
recon_matrix = measured;
recon_matrix(mask) = imputed_vals; 

%%
fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(measured)));
fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(recon_matrix)));
Diff_sq = abs(recon_matrix - ground_truth).^2;
fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq.*mask)/sum(mask(:)) ));
fprintf('Sum of squared diff: %g \n',sum(Diff_sq(:)))