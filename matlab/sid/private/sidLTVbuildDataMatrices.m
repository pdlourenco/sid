function [D, Xl] = sidLTVbuildDataMatrices(X, U, N, p, q, L)
% SIDLTVBUILDDATAMATRICES Construct COSMIC data matrices for uniform trajectories.
%
%   [D, Xl] = sidLTVbuildDataMatrices(X, U, N, p, q, L)
%
%   Builds the data matrices D(k) and next-state matrices Xl(k) for the
%   COSMIC algorithm. Each is normalized by 1/sqrt(N).
%
%     D(k)  = [X(k)' U(k)'] / sqrt(N)    size (L x p+q)
%     Xl(k) = X(k+1)' / sqrt(N)           size (L x p)
%
%   INPUTS:
%     X - (N+1 x p x L) state data.
%     U - (N x q x L) input data.
%     N - Number of time steps.
%     p - State dimension.
%     q - Input dimension.
%     L - Number of trajectories.
%
%   OUTPUTS:
%     D  - (L x p+q x N) data matrices.
%     Xl - (L x p x N) next-state matrices.
%
%   EXAMPLES:
%     [D, Xl] = sidLTVbuildDataMatrices(X, U, N, p, q, L);
%
%   SPECIFICATION:
%     SPEC.md §8.2 — Inputs
%
%   See also: sidLTVbuildDataMatricesVarLen, sidLTVbuildBlockTerms, sidLTVdisc
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

    % D(k) = [X(k); U(k)]' / sqrt(N), Xl(k) = X(k+1)' / sqrt(N) (SPEC.md §8.3.2)
    % The 1/sqrt(N) scaling makes the normal equations independent of N.
    sqrtN = sqrt(N);
    D  = zeros(L, p + q, N);  % (L x d x N) data matrices
    Xl = zeros(L, p, N);      % (L x p x N) next-state matrices

    for k = 0:N-1
        D(:, :, k+1)  = [reshape(X(k+1, :, :), p, L)', ...
                          reshape(U(k+1, :, :), q, L)'] / sqrtN;
        Xl(:, :, k+1) = reshape(X(k+2, :, :), p, L)' / sqrtN;
    end
end
