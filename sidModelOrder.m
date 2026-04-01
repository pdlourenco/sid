function [n, sv] = sidModelOrder(result, varargin)
% SIDMODELORDER Estimate model order from frequency response via Hankel SVD.
%
%   [n, sv] = sidModelOrder(result)
%   [n, sv] = sidModelOrder(result, 'Horizon', r)
%   [n, sv] = sidModelOrder(result, 'Threshold', tau)
%   sidModelOrder(result, 'Plot', true)
%
%   Estimates the state dimension n of a linear system from the singular
%   value decomposition of a block Hankel matrix built from the impulse
%   response coefficients. The impulse response is obtained via IFFT of
%   the frequency response estimate.
%
%   INPUTS:
%     result - Struct from any sidFreq* function, with at least:
%                .Response   (n_f x n_y x n_u) complex frequency response
%                .Frequency  (n_f x 1) frequency vector, rad/sample
%              For time-series (no input), uses .NoiseSpectrum instead.
%
%   NAME-VALUE OPTIONS:
%     'Horizon'    - Block Hankel prediction horizon r.
%                    Default: min(floor(N_imp/3), 50).
%     'Threshold'  - If specified, count singular values with sigma_k/sigma_1
%                    above this threshold instead of gap detection.
%                    Default: [] (use gap method).
%     'Plot'       - If true, display bar chart of singular values with
%                    the detected model order marked. Default: false.
%
%   OUTPUTS:
%     n  - Estimated model order (state dimension).
%     sv - Struct with fields:
%            .SingularValues  (m x 1) singular values of the Hankel matrix
%            .Horizon         scalar, prediction horizon used
%
%   ALGORITHM:
%     1. Compute impulse response g(k) via IFFT of the frequency response.
%     2. Build block Hankel matrix H with r block-rows and r block-cols.
%        For MIMO systems, each entry g(k) is an n_y x n_u block.
%     3. Compute SVD of H: [U, S, V] = svd(H).
%     4. Detect model order as argmax_k (sigma_k / sigma_{k+1}).
%        With threshold tau: n = number of sigma_k / sigma_1 > tau.
%
%   EXAMPLES:
%     % Automated model order detection
%     G = sidFreqBT(y, u);
%     [n, sv] = sidModelOrder(G);
%
%     % Visual inspection
%     sidModelOrder(G, 'Plot', true);
%
%     % Use with output-COSMIC
%     p_y = size(y, 2);
%     H = [eye(p_y), zeros(p_y, n - p_y)];
%     res = sidLTVdiscIO(y, u, H, 'Lambda', 1e5);
%
%   REFERENCES:
%     Kung, S.Y. "A new identification and model reduction algorithm via
%     singular value decomposition." Proc. 12th Asilomar Conference, 1978.
%
%   SPECIFICATION:
%     SPEC.md §8.12 — Output-COSMIC: Partial State Observation
%
%   See also: sidFreqBT, sidFreqETFE, sidLTVdiscIO
%
%   Changelog:
%   2026-04-01: First version by Pedro Lourenco.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenco, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    % ---- Validate input struct ----
    if ~isstruct(result)
        error('sid:badInput', 'First argument must be a sid result struct.');
    end

    isTimeSeries = ~isfield(result, 'Response') || isempty(result.Response);
    if isTimeSeries
        if ~isfield(result, 'NoiseSpectrum') || isempty(result.NoiseSpectrum)
            error('sid:badInput', ...
                'Result struct must have a Response or NoiseSpectrum field.');
        end
    end

    if ~isfield(result, 'Frequency') || isempty(result.Frequency)
        error('sid:badInput', 'Result struct must have a Frequency field.');
    end

    % ---- Parse name-value options ----
    horizon = [];
    threshold = [];
    doPlot = false;

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'horizon'
                    horizon = varargin{k+1};
                    if ~isscalar(horizon) || horizon < 1 || floor(horizon) ~= horizon
                        error('sid:badHorizon', ...
                            'Horizon must be a positive integer.');
                    end
                    k = k + 2;
                case 'threshold'
                    threshold = varargin{k+1};
                    if ~isscalar(threshold) || threshold <= 0
                        error('sid:badThreshold', ...
                            'Threshold must be a positive scalar.');
                    end
                    k = k + 2;
                case 'plot'
                    doPlot = varargin{k+1};
                    k = k + 2;
                otherwise
                    error('sid:unknownOption', 'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badInput', 'Expected option name at position %d.', k);
        end
    end

    % ---- Extract frequency response ----
    if isTimeSeries
        % Time-series: use noise spectrum (square root for "transfer function")
        G = result.NoiseSpectrum;
        ny = size(G, 2);
        nu = size(G, 3);
    else
        G = result.Response;
        ny = size(G, 2);
        nu = size(G, 3);
    end

    freqs = result.Frequency;
    nf = size(G, 1);

    % ---- Compute impulse response via IFFT ----
    % The Response is evaluated at nf frequencies in (0, pi].
    % Build conjugate-symmetric full-circle spectrum for real impulse response.
    % Assume linearly spaced frequencies: w = (1:nf) * pi/nf.

    % Construct full-circle spectrum: [G(0); G(w1); ...; G(pi); conj(G(pi-1)); ...; conj(G(w1))]
    % Approximate G(0) = real(G(w1)) (DC gain approximation)
    Nfft = 2 * nf;
    g_all = zeros(Nfft, ny, nu);

    for iy = 1:ny
        for iu = 1:nu
            Gvec = squeeze(G(:, iy, iu));  % (nf x 1) complex

            % Build full-circle: DC, positive freqs, mirror of negative freqs
            Gfull = zeros(Nfft, 1);
            Gfull(1) = real(Gvec(1));                    % DC approximation
            Gfull(2:nf) = Gvec(1:nf-1);                  % w1 to w_{nf-1}
            Gfull(nf+1) = real(Gvec(nf));                 % Nyquist (real)
            Gfull(nf+2:Nfft) = conj(Gvec(nf-1:-1:1));    % mirror

            g_all(:, iy, iu) = real(ifft(Gfull));
        end
    end

    % Use causal part (first half) as impulse response coefficients
    N_imp = nf;
    g = g_all(1:N_imp, :, :);

    % ---- Determine horizon ----
    if isempty(horizon)
        horizon = min(floor(N_imp / 3), 50);
    end

    if horizon < 2
        error('sid:tooShort', ...
            'Impulse response too short for Hankel matrix (need N_imp >= 6, got %d).', N_imp);
    end

    % Need at least 2*horizon - 1 impulse response coefficients
    if 2 * horizon - 1 > N_imp
        horizon = floor((N_imp + 1) / 2);
        if horizon < 2
            error('sid:tooShort', ...
                'Impulse response too short for Hankel matrix.');
        end
    end

    r = horizon;

    % ---- Build block Hankel matrix ----
    % H_hankel is (r*ny) x (r*nu)
    % H_hankel(i,j) block = g(i+j-1), for block indices i,j = 1..r
    H_hankel = zeros(r * ny, r * nu);

    for bi = 1:r
        for bj = 1:r
            idx = bi + bj - 1;
            if idx <= N_imp
                H_hankel((bi-1)*ny+1:bi*ny, (bj-1)*nu+1:bj*nu) = ...
                    reshape(g(idx, :, :), ny, nu);
            end
        end
    end

    % ---- SVD ----
    [~, Sigma, ~] = svd(H_hankel, 0);
    sigmas = diag(Sigma);

    % Remove near-zero singular values that are just numerical noise
    nSigma = length(sigmas);

    % ---- Detect model order ----
    if ~isempty(threshold)
        % Threshold method: count sigma_k / sigma_1 > tau
        if sigmas(1) < eps
            n = 1;
            warning('sid:allZeroSV', ...
                'All singular values near zero. Returning n = 1.');
        else
            n = sum(sigmas / sigmas(1) > threshold);
            n = max(n, 1);
        end
    else
        % Gap method: argmax_k (sigma_k / sigma_{k+1})
        if nSigma < 2 || sigmas(1) < eps
            n = 1;
            if sigmas(1) < eps
                warning('sid:allZeroSV', ...
                    'All singular values near zero. Returning n = 1.');
            end
        else
            % Only consider ratios among singular values above a noise
            % floor to avoid spurious gaps in the numerical tail.
            % The floor scales with sigma_1 and the matrix dimension.
            noiseFloor = sigmas(1) * sqrt(nSigma) * eps;
            lastSig = find(sigmas > noiseFloor, 1, 'last');
            maxK = min(lastSig, floor(nSigma / 2));
            maxK = max(maxK, 1);

            ratios = sigmas(1:maxK) ./ sigmas(2:maxK+1);

            % Find largest gap
            [~, n] = max(ratios);
        end
    end

    % ---- Optional plot ----
    if doPlot
        figure;
        bar(sigmas, 'FaceColor', [0.3 0.5 0.8]);
        hold on;
        % Mark the detected order with a vertical line
        if n < nSigma
            plot([n + 0.5, n + 0.5], ylim, 'r--', 'LineWidth', 1.5);
        end
        hold off;
        xlabel('Index');
        ylabel('Singular value');
        title(sprintf('Hankel Singular Values (detected n = %d)', n));
        grid on;
    end

    % ---- Pack output struct ----
    sv.SingularValues = sigmas;
    sv.Horizon = r;
end
