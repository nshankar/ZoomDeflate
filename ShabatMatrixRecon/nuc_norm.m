function norm = nuc_norm(measured,mask,vals_fitted)
A = measured;
A(mask) = vals_fitted; 

sigmas = svd(A); 
norm = sum(sigmas(:));  