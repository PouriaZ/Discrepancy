function out = discrepancy(A)
% Computing the discrepancy values of matrix A

s = size(A);
if s(1) ~= s(2)
    error('Input matrix must be square.')
end

if ishermitian(A) % Hermitian matrices
    out = herm_disc(A);
    
elseif norm(A'*A-A*A') < 1e-10 % Normal matrices
    out = normal_disc(A);
    
elseif isreal(A) % Real matrices
    out = sdp_real(A);
    
else % Complex matrices
    out = sdp_complex(A);
end

end