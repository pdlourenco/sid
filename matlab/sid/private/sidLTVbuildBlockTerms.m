function [S, T] = sidLTVbuildBlockTerms(D, Xl, lambda, N, p, q)
% SIDLTVBUILDBLOCKTERMS Compute block diagonal and right-hand side terms.
%
%   [S, T] = sidLTVbuildBlockTerms(D, Xl, lambda, N, p, q)
%
%   Builds the block diagonal terms S(k) = D(k)'D(k) + reg(k)*I and the
%   right-hand side T(k) = D(k)'Xl(k) for the COSMIC block tridiagonal
%   system. Handles both 3D array D (uniform trajectories) and cell array
%   D (variable-length trajectories).
%
%   INPUTS:
%     D      - Data matrices. Either (L x d x N) 3D array or {N x 1} cell
%              array where D{k} is (L_k x d).
%     Xl     - Next-state matrices. Either (L x p x N) or {N x 1} cell
%              array where Xl{k} is (L_k x p).
%     lambda - (N-1 x 1) regularization weights.
%     N      - Number of time steps.
%     p      - State dimension.
%     q      - Input dimension.
%
%   OUTPUTS:
%     S - (d x d x N) block diagonal terms with regularization added.
%     T - (d x p x N) right-hand side terms.
%
%   EXAMPLES:
%     [D, Xl] = sidLTVbuildDataMatrices(X, U, N, p, q, L);
%     [S, T] = sidLTVbuildBlockTerms(D, Xl, lambda, N, p, q);
%
%   ALGORITHM:
%     1. S(k) = D(k)' * D(k),  T(k) = D(k)' * Xl(k)  for k = 1..N
%     2. Add regularization: S(1) += lambda(1)*I, S(N) += lambda(N-1)*I,
%        S(k) += (lambda(k-1) + lambda(k))*I  for k = 2..N-1
%
%   SPECIFICATION:
%     SPEC.md §8.3 — COSMIC Algorithm
%
%   See also: sidLTVcosmicSolve, sidLTVbuildDataMatrices, sidLTVdisc
%
%   Changelog:
%   2026-04-01: First version by Pedro Lourenço.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid
%  -----------------------------------------------------------------------

    d = p + q;
    S = zeros(d, d, N);  % (d x d x N) diagonal blocks of normal equations
    T = zeros(d, p, N);  % (d x p x N) right-hand side
    useCell = iscell(D);

    % S(k) = D(k)'D(k), T(k) = D(k)'X'(k) (SPEC.md §8.3.3)
    % Regularization lambda*I is added later in sidLTVcosmicSolve.
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
