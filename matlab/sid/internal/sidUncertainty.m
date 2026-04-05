function [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W, nTraj)
% SIDUNCERTAINTY Asymptotic standard deviations for spectral estimates.
%
%   [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W)
%   [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W, nTraj)
%
%   Computes the asymptotic standard deviations of the frequency response
%   and noise spectrum estimates, based on Ljung (1999), pp. 184 and 188.
%
%   For multi-trajectory ensemble averaging, pass nTraj > 1. The variance
%   scales as 1/(N*nTraj) instead of 1/N, reflecting the L-fold reduction
%   from ensemble averaging independent trajectories.
%
%   INPUTS:
%     G     - Complex frequency response estimate (n_f x 1), or [].
%     PhiV  - Noise spectrum estimate (n_f x 1), real, non-negative.
%     Coh   - Squared coherence (n_f x 1), or [] for time series / MIMO.
%     N     - Number of data samples per trajectory.
%     W     - Hann window values for lags 0..M, (M+1 x 1) vector.
%     nTraj - (optional) Number of trajectories, default 1.
%
%   OUTPUTS:
%     GStd    - Standard deviation of G (n_f x 1), or [].
%     PhiVStd - Standard deviation of PhiV (n_f x 1).
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
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    if nargin < 6 || isempty(nTraj)
        nTraj = 1;
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
    else
        % MIMO case: per-element uncertainty is not computed in v1.0.
        % Return NaN array of the same size as G.
        GStd = nan(size(G));
    end
end
