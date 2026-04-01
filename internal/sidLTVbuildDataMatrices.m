function [D, Xl] = sidLTVbuildDataMatrices(X, U, N, p, q, L)
% SIDLTVBUILDDATAMATRICES Construct D(k) and X'(k) for all k.
%
%   D(k)  = [X(k)' U(k)'] / sqrt(N)    size (L x p+q)
%   Xl(k) = X(k+1)' / sqrt(N)           size (L x p)
%
%   Stored as 3D arrays: D is (L x p+q x N), Xl is (L x p x N).

    sqrtN = sqrt(N);
    D  = zeros(L, p + q, N);
    Xl = zeros(L, p, N);

    for k = 0:N-1
        D(:, :, k+1)  = [reshape(X(k+1, :, :), p, L)', ...
                          reshape(U(k+1, :, :), q, L)'] / sqrtN;
        Xl(:, :, k+1) = reshape(X(k+2, :, :), p, L)' / sqrtN;
    end
end
