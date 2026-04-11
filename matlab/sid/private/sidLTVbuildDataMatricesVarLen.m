function [D, Xl] = sidLTVbuildDataMatricesVarLen(X, U, N, p, q, L, horizons)
% SIDLTVBUILDDATAMATRICESVARLEN Construct COSMIC data matrices for variable-length trajectories.
%
%   [D, Xl] = sidLTVbuildDataMatricesVarLen(X, U, N, p, q, L, horizons)
%
%   Like sidLTVbuildDataMatrices, but handles trajectories with different
%   lengths. At each time step k, only trajectories with horizon > k
%   contribute. Returns cell arrays instead of 3D arrays.
%
%     D{k}  has size (|L(k)| x p+q), normalized by 1/sqrt(N)
%     Xl{k} has size (|L(k)| x p),   normalized by 1/sqrt(N)
%
%   where N = max(horizons) is the longest horizon. This matches the
%   uniform-trajectory normalization in sidLTVbuildDataMatrices, ensuring
%   that the effective regularization strength lambda is consistent
%   regardless of how many trajectories are active at each time step.
%
%   INPUTS:
%     X        - Cell array {X1, X2, ...} of L trajectories, (N_l+1 x p).
%     U        - Cell array {U1, U2, ...} of L trajectories, (N_l x q).
%     N        - Maximum horizon across all trajectories.
%     p        - State dimension.
%     q        - Input dimension.
%     L        - Number of trajectories.
%     horizons - (L x 1) horizon length of each trajectory.
%
%   OUTPUTS:
%     D  - {N x 1} cell array, D{k} is (L_k x p+q).
%     Xl - {N x 1} cell array, Xl{k} is (L_k x p).
%
%   EXAMPLES:
%     [D, Xl] = sidLTVbuildDataMatricesVarLen(X, U, N, p, q, L, horizons);
%
%   SPECIFICATION:
%     SPEC.md §8.8 — Variable-Length Trajectories
%
%   See also: sidLTVbuildDataMatrices, sidLTVbuildBlockTerms, sidLTVdisc
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

    D  = cell(N, 1);
    Xl = cell(N, 1);

    for k = 0:N-1
        % Active trajectories at step k: those with horizon > k
        active = find(horizons > k);
        Lk = length(active);

        if Lk == 0
            % No trajectories active — empty matrices
            D{k+1}  = zeros(0, p + q);
            Xl{k+1} = zeros(0, p);
            continue;
        end

        sqrtN = sqrt(N);
        Dk  = zeros(Lk, p + q);
        Xlk = zeros(Lk, p);

        for ii = 1:Lk
            l = active(ii);
            Dk(ii, :)  = [X{l}(k+1, :), U{l}(k+1, :)] / sqrtN;
            Xlk(ii, :) = X{l}(k+2, :) / sqrtN;
        end

        D{k+1}  = Dk;
        Xl{k+1} = Xlk;
    end
end
