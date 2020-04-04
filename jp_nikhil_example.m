% Demo for matrix completion
% Create a random matrix of low rank, remove values and complete.
%
% clc
ground_truth = randn(200,10)*randn(10,200); % 2000x2000 of rank 20
mask = rand(size(ground_truth)) < 0.4; 

num_samples = sum(mask(:))
frac_samples = num_samples / numel(mask)
%%


lamnbda_tol = 10;
tol = 1e-8;
N = 100;
fprintf('Completion using nuclear norm minimization... \n');
[CompletedMat, ier] = MatrixCompletion(...
    ground_truth.*mask, mask,N, 'nuclear', lamnbda_tol, tol, 0);


%% 
recon_mask = logical((1 - mask)); 


fprintf('\n Corrupted matrix nuclear norm (initial): %g \n',sum(svd(ground_truth.*mask)));
fprintf('Restored matrix nuclear norm (final): %g \n',sum(svd(CompletedMat)));
Diff_sq = abs(CompletedMat-ground_truth).^2;
fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq.*mask)/sum(mask(:)) ));
fprintf('MSE on unknown entries: %g \n',sqrt(sum2(Diff_sq.*recon_mask)/sum(recon_mask(:)) ));

figure(1)
hold on
histogram(Diff_sq .* recon_mask,100); 

%%
figure
hold on
surf(ground_truth,'LineStyle','none')
view([0 90])
colorbar
%% 
figure, hold on
surf(recon_mask.*((ground_truth - CompletedMat).^2),'LineStyle','none')
view([0 90])
colorbar
%%
figure
histogram(((ground_truth(recon_mask) - CompletedMat(recon_mask)).^2))

% %% 
% fprintf('\n Completion using spectral norm minimization... \n');
% [CompletedMat, ier] = MatrixCompletion(ground_truth.*measure_mask, measure_mask,N, 'spectral', lamnbda_tol, tol, 0);
% 
% fprintf('\n Corrupted matrix spectral norm (initial): %g \n',norm(ground_truth.*measure_mask));
% fprintf('Restored matrix spectral norm (final): %g \n',norm(CompletedMat));
% Diff_sq = abs(CompletedMat-ground_truth).^2;
% fprintf('MSE on known entries: %g \n',sqrt(sum2(Diff_sq.*measure_mask)/sum(measure_mask(:))));
% 
% fprintf('\n Completion using weighted norm minimization (pushing to low rank, no global convergence)... \n');
% [CompletedMat, ier] = MatrixCompletion(ground_truth.*measure_mask, measure_mask,N, 'NuclearWeighted', lamnbda_tol, tol, 0, [ones(1,10) ones(1,190)*10000] ); % big penalty on small singular values
% 
% fprintf('Corrupted matrix rank (initial): %g \n',rank(ground_truth.*measure_mask));
% fprintf('Restored matrix rank (final): %g \n',rank(CompletedMat));
% Diff_sq = abs(CompletedMat-ground_truth).^2;
% fprintf('MSE on known entries: %g \n \n',sqrt(sum2(Diff_sq.*measure_mask)/sum(measure_mask(:))));
% 
