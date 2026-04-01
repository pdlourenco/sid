function [S, T] = sidLTVbuildBlockTerms(D, Xl, lambda, N, p, q)
% SIDLTVBUILDBLOCKTERMS Compute the block diagonal S_kk and right-hand side T_k.
%   Handles both 3D array D (uniform) and cell array D (variable-length).

    d = p + q;
    S = zeros(d, d, N);
    T = zeros(d, p, N);
    useCell = iscell(D);

    for k = 1:N
        if useCell
            Dk  = D{k};
            Xlk = Xl{k};
        else
            Dk  = D(:, :, k);
            Xlk = Xl(:, :, k);
        end
        S(:, :, k) = Dk' * Dk;
        T(:, :, k) = Dk' * Xlk;
    end

    % Add regularization to diagonal
    I = eye(d);
    S(:, :, 1)   = S(:, :, 1)   + lambda(1) * I;
    S(:, :, N)   = S(:, :, N)   + lambda(N-1) * I;
    for k = 2:N-1
        S(:, :, k) = S(:, :, k) + (lambda(k-1) + lambda(k)) * I;
    end
end
