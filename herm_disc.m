function out = herm_disc(A)
% Computing the discrepancy values of Hermitian matrix A

lambda = sort(eig(A));
out = sort(abs(lambda - flipud(lambda))/2, 'descend');

end