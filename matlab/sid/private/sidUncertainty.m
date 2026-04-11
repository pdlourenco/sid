function [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W, nTraj, PhiU)
% SIDUNCERTAINTY Asymptotic standard deviations for spectral estimates.
%
%   [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W)
%   [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W, nTraj)
%   [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W, nTraj, PhiU)
%
%   Computes the asymptotic standard deviations of the frequency response
%   and noise spectrum estimates, based on Ljung (1999), pp. 184 and 188.
%
%   For multi-trajectory ensemble averaging, pass nTraj > 1. The variance
%   scales as 1/(N*nTraj) instead of 1/N, reflecting the L-fold reduction
%   from ensemble averaging independent trajectories.
%
%   For MIMO systems, pass PhiU (the input auto-spectrum) to enable a
%   diagonal approximation of the per-element variance:
%     Var{G_{ij}(w)} ≈ (C_W / Neff) * Phi_v_{ii}(w) / Phi_u_{jj}(w)
%   This treats each (i,j) element independently and ignores cross-channel
%   correlations, but gives meaningful error bars. Without PhiU, MIMO
%   uncertainty returns NaN.
%
%   INPUTS:
%     G     - Complex frequency response estimate, (n_f x n_y x n_u) or [].
%     PhiV  - Noise spectrum estimate, (n_f x n_y x n_y) or (n_f x 1).
%     Coh   - Squared coherence (n_f x 1), or [] for time series / MIMO.
%     N     - Number of data samples per trajectory.
%     W     - Hann window values for lags 0..M, (M+1 x 1) vector.
%     nTraj - (optional) Number of trajectories, default 1.
%     PhiU  - (optional) Input auto-spectrum, (n_f x n_u x n_u). Required
%             for MIMO uncertainty; ignored for SISO.
%
%   OUTPUTS:
%     GStd    - Standard deviation of G, same size as G, or [].
%     PhiVStd - Standard deviation of PhiV, same size as PhiV.
%
%   EXAMPLES:
%     W = sidHannWin(30);
%     [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W);
%
%   SPECIFICATION:
%     SPEC.md §3 — Uncertainty Estimation
%
%   See also: sidFreqBT, sidHannWin, sidCov
%
%   Changelog:
%   2026-03-24: First version by Pedro Lourenço.
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

    if nargin < 6 || isempty(nTraj)
        nTraj = 1;
    end
    if nargin < 7
        PhiU = [];
    end

    % Effective sample size: ensemble averaging L trajectories of length N
    % reduces variance by factor L (see SPEC.md §2)
    Neff = N * nTraj;

    % ---- Window norm: C_W = sum_{tau=-M}^{M} W(tau)^2 ----
    % W contains values for tau = 0..M. The full sum is W(0)^2 + 2*sum(W(1..M)^2).
    CW = W(1)^2 + 2 * sum(W(2:end).^2);

    % ---- Noise spectrum variance (SPEC.md §3.1) ----
    % Var{Phi_v(w)} = (2 * C_W / Neff) * Phi_v(w)^2
    PhiVStd = sqrt(2 * CW / Neff) * abs(PhiV);

    % ---- Transfer function variance (SPEC.md §3.2) ----
    if isempty(G)
        GStd = [];
        return;
    end

    if ~isempty(Coh)
        % SISO case
        % Var{G(w)} = (C_W / Neff) * |G(w)|^2 * (1 - gamma^2(w)) / gamma^2(w)
        eps_floor = 1e-10;
        cohSafe = max(Coh, eps_floor);
        GVar = (CW / Neff) .* abs(G).^2 .* (1 - cohSafe) ./ cohSafe;
        GStd = sqrt(GVar);
    elseif ~isempty(PhiU)
        % MIMO case: diagonal approximation using noise and input spectra.
        % Var{G_{ij}(w)} ≈ (C_W / Neff) * Phi_v_{ii}(w) / Phi_u_{jj}(w)
        % This treats each element independently (Ljung 1999, §9.4).
        nf = size(G, 1);
        ny = size(G, 2);
        nu = size(G, 3);
        GStd = zeros(nf, ny, nu);
        eps_floor = 1e-10;
        for k = 1:nf
            for ii = 1:ny
                % Diagonal of noise spectrum at this frequency
                if ndims(PhiV) == 3 || (ndims(PhiV) == 2 && ny > 1) %#ok<ISMAT>
                    phiV_ii = real(PhiV(k, ii, ii));
                else
                    phiV_ii = real(PhiV(k));
                end
                for jj = 1:nu
                    phiU_jj = real(PhiU(k, jj, jj));
                    if phiU_jj > eps_floor
                        GStd(k, ii, jj) = sqrt(CW / Neff * phiV_ii / phiU_jj);
                    else
                        GStd(k, ii, jj) = Inf;
                    end
                end
            end
        end
    else
        % MIMO case without PhiU: cannot compute uncertainty.
        GStd = nan(size(G));
    end
end
