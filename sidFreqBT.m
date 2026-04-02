function result = sidFreqBT(y, u, varargin)
% SIDFREQBT Estimate frequency response via Blackman-Tukey spectral analysis.
%
%   result = sidFreqBT(y, u)
%   result = sidFreqBT(y, [])
%   result = sidFreqBT(y, u, 'WindowSize', M, 'Frequencies', w)
%   result = sidFreqBT(y, u, M)
%   result = sidFreqBT(y, u, M, w)
%
%   Estimates the frequency response G(e^{jw}) and noise spectrum from
%   time-domain input/output data using the Blackman-Tukey method.
%
%   This is an open-source replacement for the System Identification
%   Toolbox function 'spa'.
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
%     'WindowSize'    - Hann window lag size M. Default: min(floor(N/10), 30).
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
%       .WindowSize       - window size M used
%       .DataLength       - number of samples N
%       .NumTrajectories  - number of trajectories L
%       .Method           - 'sidFreqBT'
%
%   EXAMPLES:
%     % SISO system identification
%     N = 1000; u = randn(N,1);
%     y = filter([1], [1 -0.9], u) + 0.1*randn(N,1);
%     result = sidFreqBT(y, u);
%     sidBodePlot(result);
%
%     % Time series spectrum estimation
%     y = randn(500,1);
%     result = sidFreqBT(y, []);
%     sidSpectrumPlot(result);
%
%     % Custom window size and frequencies
%     w = linspace(0.01, pi, 256)';
%     result = sidFreqBT(y, u, 'WindowSize', 50, 'Frequencies', w);
%
%     % Multi-trajectory: average 5 independent experiments
%     L = 5; N = 1000;
%     y3d = zeros(N, 1, L); u3d = zeros(N, 1, L);
%     for l = 1:L
%         u3d(:,1,l) = randn(N, 1);
%         y3d(:,1,l) = filter([1], [1 -0.9], u3d(:,1,l)) + 0.1*randn(N,1);
%     end
%     result = sidFreqBT(y3d, u3d);
%
%   ALGORITHM:
%     1. Compute biased sample covariances R_y, R_u, R_yu for lags 0..M
%     2. Apply Hann window and Fourier transform to obtain spectral estimates
%     3. Form G = Phi_yu / Phi_u and Phi_v = Phi_y - |Phi_yu|^2 / Phi_u
%     4. Compute asymptotic standard deviations
%
%   REFERENCES:
%     Ljung, L. "System Identification: Theory for the User", 2nd ed.,
%     Prentice Hall, 1999. Sections 2.3, 6.3-6.4.
%
%   SPECIFICATION:
%     SPEC.md §2 — Blackman-Tukey Spectral Analysis
%
%   See also: sidFreqBTFDR, sidFreqETFE, sidBodePlot, sidSpectrumPlot
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

    % ---- Parse inputs ----
    [y, u, M, freqs, Ts, isTimeSeries, nTraj] = sidValidate(y, u, varargin{:});

    N = size(y, 1);
    ny = size(y, 2);
    nf = length(freqs);
    useFFT = isempty(varargin) || isDefaultFreqs(freqs, nf);

    % ---- Compute covariances ----
    Ryy = sidCov(y, y, M);                       % (M+1) x ny x ny

    if ~isTimeSeries
        nu = size(u, 2);
        Ruu = sidCov(u, u, M);                   % (M+1) x nu x nu
        Ryu = sidCov(y, u, M);                   % (M+1) x ny x nu
        Ruy = sidCov(u, y, M);                   % (M+1) x nu x ny — for negative lags
    end

    % ---- Window ----
    W = sidHannWin(M);                            % (M+1) x 1, for lags 0..M

    % ---- Spectral estimates ----
    PhiY = sidWindowedDFT(Ryy, W, freqs, useFFT, Ryy); % (nf x ny x ny)

    if ~isTimeSeries
        PhiU  = sidWindowedDFT(Ruu, W, freqs, useFFT, Ruu); % (nf x nu x nu)
        PhiYU = sidWindowedDFT(Ryu, W, freqs, useFFT, Ruy); % (nf x ny x nu)
    end

    % ---- Form transfer function and noise spectrum ----
    if isTimeSeries
        G = [];
        PhiV = real(PhiY);
        Coh = [];
    elseif ny == 1 && nu == 1
        % SISO
        G = PhiYU ./ PhiU;
        PhiV = real(PhiY) - abs(PhiYU).^2 ./ real(PhiU);
        PhiV = max(PhiV, 0);
        Coh = abs(PhiYU).^2 ./ (real(PhiY) .* real(PhiU));
        Coh = min(max(Coh, 0), 1);
    else
        % MIMO
        G = zeros(nf, ny, nu);
        PhiV = zeros(nf, ny, ny);
        for k = 1:nf
            PhiU_k = reshape(PhiU(k, :, :), nu, nu);
            PhiYU_k = reshape(PhiYU(k, :, :), ny, nu);
            PhiY_k = reshape(PhiY(k, :, :), ny, ny);
            G(k, :, :) = PhiYU_k / PhiU_k;
            PhiV(k, :, :) = PhiY_k - PhiYU_k / PhiU_k * PhiYU_k';
        end
        PhiV = real(PhiV);
        Coh = [];
    end

    % ---- Uncertainty ----
    if isTimeSeries
        [~, PhiVStd] = sidUncertainty([], PhiV, [], N, W, nTraj);
        GStd = [];
    else
        [GStd, PhiVStd] = sidUncertainty(G, PhiV, Coh, N, W, nTraj);
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
    result.WindowSize       = M;
    result.DataLength       = N;
    result.NumTrajectories  = nTraj;
    result.Method           = 'sidFreqBT';
end

function tf = isDefaultFreqs(freqs, nf)
% ISDEFAULTFREQS Check if frequency vector matches the default linear grid.
    if nf ~= 128
        tf = false;
        return;
    end
    defaultFreqs = (1:128)' * pi / 128;
    tf = max(abs(freqs(:) - defaultFreqs)) < 1e-12;
end
