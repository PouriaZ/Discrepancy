function out = sdp_complex(X)
% Computing the discrepancy values of complex square matrix A

dim = length(X);
C = [zeros(dim,dim), X; X', zeros(dim,dim)];

out = zeros(dim,1);
disc_norm = 0;

for k = 1:dim
    cvx_begin quiet
        variable t
        variable s complex
        variable al complex
        variable Z(2*dim, 2*dim) complex;
        minimize(t);
        subject to
           Z == hermitian_semidefinite(2*dim);
           t >= real(trace(Z) + k*s);
           Z + s*eye(2*dim,2*dim) - C + [zeros(dim,dim), al*eye(dim,dim); ...
           conj(al)*eye(dim,dim), zeros(dim,dim)] == hermitian_semidefinite(2*dim);
    cvx_end

    out(k) = cvx_optval - disc_norm;
    disc_norm = disc_norm + out(k);
end

end