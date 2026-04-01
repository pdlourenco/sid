function [C, Lbd] = sidLTVcosmicSolve(S, T, lambda, N, p, q)
% SIDLTVCOSMSICSOLVE COSMIC forward-backward block tridiagonal solver.
%
%   [C, Lbd] = sidLTVcosmicSolve(S, T, lambda, N, p, q)
%
%   Solves the block tridiagonal system arising from the regularized
%   least squares formulation in the COSMIC algorithm. Returns the
%   forward Schur complements Lbd for reuse in uncertainty computation.
%
%   INPUTS:
%     S      - (d x d x N) block diagonal terms from sidLTVbuildBlockTerms.
%     T      - (d x p x N) right-hand side terms.
%     lambda - (N-1 x 1) regularization weights.
%     N      - Number of time steps.
%     p      - State dimension.
%     q      - Input dimension.
%
%   OUTPUTS:
%     C   - (d x p x N) solution [A(k)'; B(k)'] for each time step.
%     Lbd - (d x d x N) forward Schur complements.
%
%   EXAMPLES:
%     [S, T] = sidLTVbuildBlockTerms(D, Xl, lambda, N, p, q);
%     [C, Lbd] = sidLTVcosmicSolve(S, T, lambda, N, p, q);
%
%   ALGORITHM:
%     Forward pass:  Lbd(k) = S(k) - lambda(k-1)^2 * Lbd(k-1)^{-1}
%                    Y(k) = Lbd(k) \ (T(k) + lambda(k-1) * Y(k-1))
%     Backward pass: C(k) = Y(k) + lambda(k) * Lbd(k) \ C(k+1)
%     Complexity: O(N * d^3), d = p + q.
%
%   REFERENCES:
%     Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
%     identification from large-scale data for LTV systems."
%     arXiv:2112.04355, 2022.
%
%   SPECIFICATION:
%     SPEC.md §8.3 — COSMIC Algorithm
%
%   See also: sidLTVbuildBlockTerms, sidLTVuncertaintyBackwardPass, sidLTVdisc
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
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

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
