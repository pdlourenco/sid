function result = sidFreqBTFDR(y, u, varargin)
% SIDFREQBTFDR Blackman-Tukey spectral analysis with frequency-dependent resolution.
%
%   result = sidFreqBTFDR(y, u)
%   result = sidFreqBTFDR(y, [])
%   result = sidFreqBTFDR(y, u, 'Resolution', R)
%   result = sidFreqBTFDR(y, u, 'Resolution', R, 'Frequencies', w, 'SampleTime', Ts)
%
%   Like sidFreqBT, but the window size varies across frequencies.
%   The user specifies a resolution parameter (in rad/sample) instead
%   of a fixed window size. Finer resolution (smaller R) uses a larger
%   window and gives lower variance but coarser frequency detail, while
%   coarser resolution (larger R) uses a smaller window.
%
%   This is an open-source replacement for the System Identification
%   Toolbox function 'spafdr'.
%
%   INPUTS:
%     y    - Output data, (N x n_y) matrix. Column vector for SISO.
%            For multiple trajectories: (N x n_y x L) array or cell array
%            {y1, y2, ...} for variable-length data. Spectral estimates
%            are ensemble-averaged across trajectories.
%     u    - Input data, (N x n_u) matrix. Column vector for SISO.
%            For multiple trajectories: (N x n_u x L) or cell array.
%            Use [] for time series (output spectrum only).
%
%   NAME-VALUE OPTIONS:
%     'Resolution'    - Frequency resolution in rad/sample. Scalar (uniform)
%                       or vector of same length as frequency grid (per-freq).
%                       Default: 2*pi / min(floor(N/10), 30).
%     'Frequencies'   - Frequency vector in rad/sample, in (0, pi].
%                       Default: 128 linearly spaced values.
%     'SampleTime'    - Sample time in seconds. Default: 1.0.
%
%   OUTPUTS:
%     result - Struct with fields:
%       .Frequency        - (n_f x 1) frequency vector, rad/sample
%       .FrequencyHz      - (n_f x 1) frequency vector, Hz
%       .Response         - (n_f x n_y x n_u) complex frequency response
%       .ResponseStd      - (n_f x n_y x n_u) standard deviation of Response
%       .NoiseSpectrum    - (n_f x n_y x n_y) noise spectrum
%       .NoiseSpectrumStd - (n_f x n_y x n_y) standard deviation
%       .Coherence        - (n_f x 1) squared coherence (SISO only)
%       .SampleTime       - sample time in seconds
%       .WindowSize       - (n_f x 1) vector of window sizes M_k
%       .DataLength       - number of samples N
%       .NumTrajectories  - number of trajectories L
%       .Method           - 'sidFreqBTFDR'
%
%   EXAMPLES:
%     N = 1000; u = randn(N, 1);
%     y = filter([1 0.5], [1 -0.8], u) + 0.1*randn(N, 1);
%     result = sidFreqBTFDR(y, u, 'Resolution', 0.3);
%     sidBodePlot(result);
%
%   ALGORITHM:
%     For each frequency w_k:
%       1. Determine local window size M_k = ceil(2*pi / R_k).
%       2. Compute Hann window W_{M_k} and biased covariances up to lag M_k.
%       3. Compute windowed spectral estimates via direct DFT.
%       4. Form G(w_k) and Phi_v(w_k) as in sidFreqBT.
%       5. Compute asymptotic uncertainty using local window norm.
%
%   REFERENCES:
%     Ljung, L. "System Identification: Theory for the User", 2nd ed.,
%     Prentice Hall, 1999. Sections 6.3-6.4.
%
%   SPECIFICATION:
%     SPEC.md §5 — Frequency-Dependent Resolution
%
%   See also: sidFreqBT, sidFreqETFE, sidBodePlot, sidSpectrumPlot
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

    % ---- Parse inputs ----
    [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u);
    Neff = N * nTraj;  % effective sample size for variance scaling

    defs.Resolution = [];
    defs.Frequencies = [];
    defs.SampleTime = 1.0;
    opts = sidParseOptions(defs, varargin);
    R = opts.Resolution;
    freqs = opts.Frequencies;
    Ts = opts.SampleTime;

    % ---- Defaults ----
    if isempty(freqs)
        freqs = (1:128)' * pi / 128;
    else
        freqs = freqs(:);
    end

    nf = length(freqs);

    if isempty(R)
        Mdefault = min(floor(N / 10), 30);
        if Mdefault < 2
            Mdefault = 2;
        end
        R = 2 * pi / Mdefault;
    end

    % ---- Validate parameters ----
    if Ts <= 0
        error('sid:badTs', 'Sample time must be positive.');
    end

    if any(freqs <= 0) || any(freqs > pi)
        error('sid:badFreqs', 'Frequencies must be in the range (0, pi] rad/sample.');
    end

    if any(R <= 0)
        error('sid:badResolution', 'Resolution must be positive.');
    end

    % Expand scalar R to vector
    if isscalar(R)
        R = R * ones(nf, 1);
    else
        R = R(:);
        if length(R) ~= nf
            error('sid:badResolution', ...
                'Resolution vector length (%d) must match frequency vector length (%d).', ...
                length(R), nf);
        end
    end

    % ---- Resolution to window size (SPEC.md §5.2) ----
    % M_k = ceil(2*pi / R_k). Reaching a bound is reported, not silently
    % clamped (SPEC.md §5.2): an error when the requested resolution implies
    % M_k < 2, and a windowReduced warning when any M_k is reduced to
    % floor(N/2). Mirrors the fixed-window sidFreqBT handling.
    MkRaw = ceil(2 * pi ./ R);
    if any(MkRaw < 2)
        error('sid:badResolution', ...
            ['Resolution too coarse: it implies a lag-window size M_k < 2. ' ...
             'Decrease the resolution value.']);
    end
    Nhalf = floor(N / 2);
    if any(MkRaw > Nhalf)
        warning('sid:windowReduced', ...
            ['Resolution finer than the data supports at some frequencies; ' ...
             'window size reduced to N/2 = %d.'], Nhalf);
    end
    Mk = min(MkRaw, Nhalf);            % M <= N/2

    % ---- Pre-compute biased covariances up to max(M_k) (SPEC.md §2.3) ----
    Mmax = max(Mk);
    Ryy = sidCov(y, y, Mmax);          % (Mmax+1 x ny x ny)

    if ~isTimeSeries
        Ruu = sidCov(u, u, Mmax);      % (Mmax+1 x nu x nu)
        Ryu = sidCov(y, u, Mmax);      % (Mmax+1 x ny x nu)
        Ruy = sidCov(u, y, Mmax);      % (Mmax+1 x nu x ny) negative lags
    end

    % ---- Per-frequency spectral estimation (SPEC.md §5.2) ----
    % Determine signal dimensions for output pre-allocation
    isSISO = (ny == 1 && nu == 1 && ~isTimeSeries);

    if isTimeSeries
        G = [];
        PhiV = zeros(nf, ny, ny);
        GStd = [];
        PhiVStd = zeros(nf, ny, ny);
        Coh = [];

        for kk = 1:nf
            Mk_k = Mk(kk);

            % Truncate covariance and compute window for this frequency
            W = sidHannWin(Mk_k);
            Ryy_k = truncateCov(Ryy, Mk_k, ny, ny);

            % Direct DFT at single frequency
            PhiY_k = singleFreqDFT(Ryy_k, W, freqs(kk), ny, ny);
            PhiV(kk, :, :) = real(PhiY_k);

            % Var{Phi_y} = (2*C_W/N) * Phi_y^2 (SPEC.md §5.3)
            CW = W(1)^2 + 2 * sum(W(2:end).^2);
            PhiVStd(kk, :, :) = sqrt(2 * CW / Neff) * abs(PhiY_k);
        end

        % Squeeze if scalar
        if ny == 1
            PhiV = PhiV(:);
            PhiVStd = PhiVStd(:);
        end

    else
        % ---- Input/output path ----
        % Assemble the per-frequency windowed spectra into stacked arrays, then
        % apply the shared degenerate handling (SPEC.md §2.6/§2.7/§10.3) so BT
        % and BTFDR cannot drift. #141: a singular Phi_u now yields NaN + Inf
        % rather than the MATLAB "singular to working precision" Inf garbage.
        Wstore = cell(nf, 1);
        CWall = zeros(nf, 1);
        for kk = 1:nf
            W = sidHannWin(Mk(kk));
            Wstore{kk} = W;
            CWall(kk) = W(1)^2 + 2 * sum(W(2:end).^2);
        end

        forceDegenerate = sidInputExcitationDegenerate(u);
        if ~forceDegenerate
            dead = sidDeadInputChannels(u);
            if any(dead)
                warning('sid:deadInputChannel', ...
                    ['Input channel(s) %s are (near-)constant; their ' ...
                     'frequency-response columns are unreliable ' ...
                     '(SPEC.md 10.3).'], mat2str(find(dead)));
            end
        end

        if isSISO
            PhiY = zeros(nf, 1);
            PhiU = zeros(nf, 1);
            PhiYU = zeros(nf, 1);
            for kk = 1:nf
                Mk_k = Mk(kk);
                W = Wstore{kk};
                PhiY(kk)  = scalarSingleFreqDFT(Ryy(1:Mk_k+1), W, freqs(kk));
                PhiU(kk)  = scalarSingleFreqDFT(Ruu(1:Mk_k+1), W, freqs(kk));
                PhiYU(kk) = scalarSingleFreqDFT(Ryu(1:Mk_k+1), W, freqs(kk), Ruy(1:Mk_k+1));
            end
        else
            PhiY = zeros(nf, ny, ny);
            PhiU = zeros(nf, nu, nu);
            PhiYU = zeros(nf, ny, nu);
            for kk = 1:nf
                Mk_k = Mk(kk);
                W = Wstore{kk};
                PhiY(kk, :, :)  = singleFreqDFT( ...
                    truncateCov(Ryy, Mk_k, ny, ny), W, freqs(kk), ny, ny);
                PhiU(kk, :, :)  = singleFreqDFT( ...
                    truncateCov(Ruu, Mk_k, nu, nu), W, freqs(kk), nu, nu);
                PhiYU(kk, :, :) = singleFreqDFT( ...
                    truncateCov(Ryu, Mk_k, ny, nu), W, freqs(kk), ny, nu, ...
                    truncateCov(Ruy, Mk_k, nu, ny));
            end
        end

        [G, PhiV, Coh] = sidRegularizeResponse(PhiYU, PhiU, PhiY, ny, nu, forceDegenerate);

        % ---- Per-frequency asymptotic uncertainty (SPEC.md §5.3, §3.3) ----
        eps_floor = 1e-10;
        if isSISO
            GVar = (CWall / Neff) .* abs(G).^2 .* (1 - Coh) ./ Coh;
            GStd = sqrt(GVar);
            % §3.3: sigma_G = Inf where the input carries no usable information.
            GStd(Coh < eps_floor) = Inf;
            PhiVStd = sqrt(2 * CWall / Neff) .* abs(PhiV);
        else
            GStd = zeros(nf, ny, nu);
            PhiVStd = zeros(nf, ny, ny);
            for kk = 1:nf
                PhiV_k = reshape(PhiV(kk, :, :), ny, ny);
                PhiU_k = reshape(PhiU(kk, :, :), nu, nu);
                PhiVStd(kk, :, :) = sqrt(2 * CWall(kk) / Neff) * abs(PhiV_k);
                for ii = 1:ny
                    phiV_ii = max(real(PhiV_k(ii, ii)), 0);
                    for jj = 1:nu
                        phiU_jj = real(PhiU_k(jj, jj));
                        if phiU_jj > eps_floor
                            GStd(kk, ii, jj) = sqrt(CWall(kk) / Neff * phiV_ii / phiU_jj);
                        else
                            GStd(kk, ii, jj) = Inf;
                        end
                    end
                end
            end
            % §3.3 (equivalently): Inf where a singular Phi_u NaN'd the G slice.
            GStd(isnan(G)) = Inf;
        end
    end

    % ---- Pack result ----
    result.Frequency        = freqs(:);
    result.FrequencyHz      = freqs(:) / (2 * pi * Ts);
    result.Response         = G;
    result.ResponseStd      = GStd;
    result.NoiseSpectrum    = PhiV;
    result.NoiseSpectrumStd = PhiVStd;
    result.Coherence        = Coh;
    result.SampleTime       = Ts;
    result.WindowSize       = Mk(:);
    result.DataLength       = N;
    result.NumTrajectories  = nTraj;
    result.Method           = 'sidFreqBTFDR';
end

function Rk = truncateCov(R, Mk, p, q)
% TRUNCATECOV Extract covariances for lags 0..Mk from a larger array.
    if p == 1 && q == 1
        Rk = R(1:Mk+1);
    else
        Rk = R(1:Mk+1, :, :);
    end
end

function Phi = singleFreqDFT(R, W, w, p, q, Rneg)
% SINGLEFREQDFT Windowed DFT at a single frequency for matrix signals.
%   R:     (Mk+1 x p x q) covariance array
%   W:     (Mk+1 x 1) window values
%   w:     scalar frequency (rad/sample)
%   Rneg:  (optional) (Mk+1 x q x p) reverse covariance for negative lags
%   Returns: (p x q) complex spectral matrix
    if nargin < 6
        Rneg = [];
    end
    Mk = length(W) - 1;
    Phi = zeros(p, q);

    for ii = 1:p
        for jj = 1:q
            if p == 1 && q == 1
                Rvec = R(:);
                if isempty(Rneg)
                    Rneg_vec = [];
                else
                    Rneg_vec = Rneg(:);
                end
            else
                Rvec = R(:, ii, jj);
                if isempty(Rneg)
                    Rneg_vec = [];
                else
                    Rneg_vec = Rneg(:, jj, ii);
                end
            end
            Phi(ii, jj) = scalarSingleFreqDFT(Rvec, W, w, Rneg_vec);
        end
    end
end

function val = scalarSingleFreqDFT(R, W, w, Rneg)
% SCALARSINGLEFREQDFT Windowed DFT at one frequency for scalar covariance.
%   R:    (M+1 x 1) covariance for lags 0..M
%   W:    (M+1 x 1) window values
%   w:    scalar frequency
%   Rneg: (optional) (M+1 x 1) reverse covariance for negative lags
%   Returns: scalar complex spectral estimate
    if nargin < 4
        Rneg = [];
    end
    M = length(W) - 1;

    % Lag 0 contribution
    val = W(1) * R(1);

    % Lags 1..M: combine positive and negative lag contributions
    % For auto-covariance: R(-tau) = conj(R(tau))
    % For cross-covariance: R_xy(-tau) = Rneg(tau) = R_yx(tau)
    for tau = 1:M
        e = exp(-1j * w * tau);
        if isempty(Rneg)
            val = val + W(tau + 1) * (R(tau + 1) * e + conj(R(tau + 1)) * conj(e));
        else
            val = val + W(tau + 1) * (R(tau + 1) * e + Rneg(tau + 1) * conj(e));
        end
    end
end
