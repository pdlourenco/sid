function [cost, fidelity, reg] = sidLTVevaluateCost(A, B, D, Xl, lambda, N, p, q)
% SIDLTVEVALUATECOST Compute COSMIC cost function value.
%
%   [cost, fidelity, reg] = sidLTVevaluateCost(A, B, D, Xl, lambda, N, p, q)
%
%   Evaluates the total COSMIC cost as the sum of data fidelity and
%   temporal regularization terms. Handles both 3D array D (uniform
%   trajectories) and cell array D (variable-length trajectories).
%
%   INPUTS:
%     A      - (p x p x N) time-varying dynamics matrices.
%     B      - (p x q x N) time-varying input matrices.
%     D      - Data matrices. (L x d x N) 3D array or {N x 1} cell array.
%     Xl     - Next-state matrices. (L x p x N) or {N x 1} cell array.
%     lambda - (N-1 x 1) regularization weights.
%     N      - Number of time steps.
%     p      - State dimension.
%     q      - Input dimension.
%
%   OUTPUTS:
%     cost     - Total cost = fidelity + reg.
%     fidelity - (1/2) sum_k ||D(k) C(k) - Xl(k)||^2_F.
%     reg      - (1/2) sum_k lambda_k ||C(k) - C(k+1)||^2_F.
%
%   EXAMPLES:
%     [cost, fid, reg] = sidLTVevaluateCost(A, B, D, Xl, lambda, N, p, q);
%
%   ALGORITHM:
%     1. For each k: fidelity += ||D(k) * [A(k)'; B(k)'] - Xl(k)||^2_F
%     2. For each k < N: reg += lambda(k) * ||C(k) - C(k+1)||^2_F
%     3. cost = 0.5 * fidelity + 0.5 * reg
%
%   SPECIFICATION:
%     SPEC.md §8.3 — COSMIC Algorithm
%
%   See also: sidLTVcosmicSolve, sidLTVdisc
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

    % J = (1/2) sum_k ||D(k)C(k) - X'(k)||^2 + (1/2) sum_k lambda(k)||C(k+1)-C(k)||^2
    % (SPEC.md §8.3.3)
    fidelity = 0;
    priorVec = zeros(N - 1, 1);
    useCell = iscell(D);

    for k = 1:N
        % Data fidelity: ||D(k)*C(k) - X'(k)||_F^2
        % C(k) = [A(k)'; B(k)'] (SPEC.md §8.3.1)
        Ck = [A(:, :, k)'; B(:, :, k)'];
        if useCell
            residual = D{k} * Ck - Xl{k};
        else
            residual = D(:, :, k) * Ck - Xl(:, :, k);
        end
        fidelity = fidelity + norm(residual, 'fro')^2;

        % Regularization: ||C(k) - C(k-1)||^2_F
        if k < N
            Ck1 = [A(:, :, k+1)'; B(:, :, k+1)'];
            priorVec(k) = norm(Ck - Ck1, 'fro')^2;
        end
    end

    fidelity = 0.5 * fidelity;
    reg      = 0.5 * lambda' * priorVec;
    cost     = fidelity + reg;
end
