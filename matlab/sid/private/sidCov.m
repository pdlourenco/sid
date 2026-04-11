function R = sidCov(x, z, maxLag)
% SIDCOV Biased cross-covariance estimate for lags 0..maxLag.
%
%   R = sidCov(x, z, maxLag)
%
%   Computes the biased sample cross-covariance:
%
%     R(tau) = (1/N) * sum_{t=1}^{N-tau} x(t+tau) * z(t)'
%
%   for tau = 0, 1, ..., maxLag.
%
%   Multi-trajectory support: if x and z are 3D arrays (N x p x L) and
%   (N x q x L), the covariance is ensemble-averaged across L trajectories:
%
%     R(tau) = (1/L) * sum_{l=1}^{L} R_l(tau)
%
%   where R_l is the per-trajectory biased covariance. This reduces
%   variance by a factor of L without affecting frequency resolution.
%
%   INPUTS:
%     x      - (N x p) or (N x p x L) matrix, first signal
%     z      - (N x q) or (N x q x L) matrix, second signal (can equal x)
%     maxLag - maximum lag M (non-negative integer)
%
%   OUTPUTS:
%     R      - (maxLag+1 x p x q) array of covariance estimates.
%              R(tau+1, :, :) is the covariance at lag tau.
%              For scalar signals, R is (maxLag+1 x 1).
%
%   EXAMPLES:
%     x = randn(200, 1);
%     R = sidCov(x, x, 30);  % auto-covariance for lags 0..30
%
%   SPECIFICATION:
%     SPEC.md §2.3 — Covariance Estimation
%
%   See also: sidFreqBT, sidHannWin, sidWindowedDFT
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

    N = size(x, 1);
    p = size(x, 2);
    q = size(z, 2);

    % Detect multi-trajectory (3D arrays)
    if ndims(x) == 3 %#ok<ISMAT>
        L = size(x, 3);
    else
        L = 1;
    end

    R = zeros(maxLag + 1, p, q);  % (maxLag+1 x p x q)

    if L == 1
        % R-hat_xz(tau) = (1/N) * sum_{t=tau+1}^{N} x(t) z(t-tau)' (SPEC.md §2.3)
        for tau = 0:maxLag
            R(tau + 1, :, :) = (x(tau + 1:N, :)' * z(1:N - tau, :)) / N;
        end
    else
        % Ensemble average: R-hat = (1/L) sum_l R-hat^{(l)} (SPEC.md §2.3)
        for tau = 0:maxLag
            Rsum = zeros(p, q);
            for l = 1:L
                Rsum = Rsum + x(tau + 1:N, :, l)' * z(1:N - tau, :, l);
            end
            R(tau + 1, :, :) = Rsum / (L * N);
        end
    end

    % Squeeze trailing singleton dimensions for scalar signals
    if p == 1 && q == 1
        R = R(:);
    end
end
