function [C, Lbd] = sidLTVcosmicSolve(S, T, lambda, N, p, q)
% SIDLTVCOSMSICSOLVE COSMIC forward-backward pass.
%
%   Solves the block tridiagonal system arising from the regularized
%   least squares formulation.
%
%   Forward pass: compute Lbd_k and Y_k for k = 0..N-1
%   Backward pass: recover C(k) for k = N-2..0
%
%   Returns Lbd (the forward Schur complements) for use in uncertainty
%   computation.

    d = p + q;
    Lbd = zeros(d, d, N);
    Y   = zeros(d, p, N);
    C   = zeros(d, p, N);
    I   = eye(d);

    % Forward pass
    Lbd(:, :, 1) = S(:, :, 1);
    Y(:, :, 1)   = Lbd(:, :, 1) \ T(:, :, 1);

    for k = 2:N
        Lbd(:, :, k) = S(:, :, k) - lambda(k-1)^2 * (Lbd(:, :, k-1) \ I);
        Y(:, :, k)   = Lbd(:, :, k) \ (T(:, :, k) + lambda(k-1) * Y(:, :, k-1));
    end

    % Backward pass
    C(:, :, N) = Y(:, :, N);

    for k = N-1:-1:1
        C(:, :, k) = Y(:, :, k) + lambda(k) * (Lbd(:, :, k) \ C(:, :, k+1));
    end
end
