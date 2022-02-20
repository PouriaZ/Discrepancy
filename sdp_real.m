function out = sdp_real(X)
% Computing the discrepancy values of real square matrix A

dim = length(X);
C = [zeros(dim,dim), X; X', zeros(dim,dim)];

out = zeros(dim,1);
disc_norm = 0;

for k = 1:dim
    cvx_begin quiet
        variable al
        variable t
        variable s
        variable Z(2*dim, 2*dim);
        minimize(t);
        subject to
           Z == semidefinite(2*dim);
           t >= trace(Z) + k*s;
           Z + s*eye(2*dim,2*dim) - C + [zeros(dim,dim), al*eye(dim,dim); ...
           al*eye(dim,dim), zeros(dim,dim)] == semidefinite(2*dim);
    cvx_end

    out(k) = cvx_optval - disc_norm;
    disc_norm = disc_norm + out(k);
end

end