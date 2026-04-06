function [Sigma, dof] = sidEstimateNoiseCov(C, D, Xl, P, covMode, N, p, q)
% SIDESTIMATENOISECOV Estimate noise covariance from COSMIC residuals.
%
%   [Sigma, dof] = sidEstimateNoiseCov(C, D, Xl, P, covMode, N, p, q)
%
%   Estimates the (p x p) noise covariance matrix from the residuals of a
%   COSMIC solution. The data D and Xl are scaled by 1/sqrt(L) (COSMIC
%   convention). The scaled residuals E_s(k) have noise covariance Sigma/L.
%   This function returns the UNSCALED noise covariance Sigma.
%
%   INPUTS:
%     C       - (d x p x N) COSMIC solution matrices.
%     D       - Data matrices. Either (L x d x N) 3D array or {N x 1}
%               cell array where D{k} is (L_k x d).
%     Xl      - Next-state matrices. Either (L x p x N) or {N x 1} cell
%               array where Xl{k} is (L_k x p).
%     P       - (d x d x N) posterior covariance diagonal blocks from
%               sidLTVuncertaintyBackwardPass.
%     covMode - Covariance structure: 'diagonal', 'full', or 'isotropic'.
%     N       - Number of time steps.
%     p       - State dimension.
%     q       - Input dimension.
%
%   OUTPUTS:
%     Sigma - (p x p) estimated noise covariance matrix.
%     dof   - Effective degrees of freedom used.
%
%   EXAMPLES:
%     [Sigma, dof] = sidEstimateNoiseCov(C, D, Xl, P, 'diagonal', N, p, q);
%
%   SPECIFICATION:
%     SPEC.md §8.9.3 — Noise Covariance Estimation
%
%   See also: sidLTVdisc, sidLTVdiscIO, sidLTVuncertaintyBackwardPass
%
%   Changelog:
%   2026-04-05: Extracted from sidLTVdisc by Pedro Lourenco.
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
    useCell = iscell(D);

    % Accumulate scaled residual scatter matrix and count observations
    SSR_scaled = zeros(p, p);
    totalObs = 0;

    for k = 1:N
        Ck = C(:, :, k);
        if useCell
            Dk  = D{k};
            Xlk = Xl{k};
            Lk  = size(Dk, 1);
        else
            Dk  = D(:, :, k);
            Xlk = Xl(:, :, k);
            Lk  = size(Dk, 1);
        end

        if Lk == 0
            continue;
        end

        Ek = Xlk - Dk * Ck;  % (Lk x p), scaled residuals
        SSR_scaled = SSR_scaled + Ek' * Ek;
        totalObs = totalObs + Lk;
    end

    % Exact degrees of freedom using unscaled hat-matrix trace
    % (SPEC.md §8.9.3)
    traceSum = 0;
    for k = 1:N
        if useCell
            Dk = D{k};
        else
            Dk = D(:, :, k);
        end
        if size(Dk, 1) > 0
            DtD = Dk' * Dk;  % (d x d), scaled
            traceSum = traceSum + sum(sum(DtD .* P(:, :, k)));
        end
    end

    dof = totalObs - N * traceSum;

    % Conservative fallback if exact dof is non-positive
    if dof <= 0
        dof = totalObs - N * d;
        if dof <= 0
            dof = max(totalObs, 1);
        end
    end

    % Unscaled noise covariance: scaled residuals have variance Sigma/L,
    % so Sigma = L * SSR_scaled / dof. Since D is already scaled by
    % 1/sqrt(L), we use N (the scale factor) here.
    Sigma = N * SSR_scaled / dof;

    % Apply covariance mode restriction (SPEC.md §8.9.3)
    switch covMode
        case 'diagonal'
            Sigma = diag(diag(Sigma));
        case 'isotropic'
            Sigma = (trace(Sigma) / p) * eye(p);
        case 'full'
            % Keep as is
    end
end
