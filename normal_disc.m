function out = normal_disc(X)
% Computing the discrepancy values of Normal matrix A using SOCP

dim = length(X);
lambda = eig(X);

disc_norm = 0;
out = zeros(dim,1);

for k = 1:dim
    cvx_begin quiet
        variable al complex
        variable q
        variable u(dim,1)
        variable x(dim,1)
        minimize(ones(dim, 1)'*u + k*q);
        subject to
           u >= x - q*ones(dim,1);
           x >= abs(lambda - al*ones(dim,1));
           u >= 0;
    cvx_end
    
    out(k) = cvx_optval - disc_norm;
    disc_norm = disc_norm + out(k);
end

end