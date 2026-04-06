function result = sidFreqETFE(y, u, varargin)
% SIDFREQETFE Empirical transfer function estimate.
%
%   result = sidFreqETFE(y, u)
%   result = sidFreqETFE(y, [])
%   result = sidFreqETFE(y, u, 'Smoothing', S)
%   result = sidFreqETFE(y, u, 'Smoothing', S, 'Frequencies', w, 'SampleTime', Ts)
%
%   Estimates the frequency response as the ratio of output and input
%   discrete Fourier transforms. Provides maximum frequency resolution
%   but high variance. Optional smoothing reduces variance.
%
%   This is an open-source replacement for the System Identification
%   Toolbox function 'etfe'.
%
%   INPUTS:
%     y    - Output data, (N x n_y) matrix. Column vector for SISO.
%            For multiple trajectories: (N x n_y x L) array or cell array
%            {y1, y2, ...} for variable-length data. Cross-periodograms
%            are ensemble-averaged across trajectories.
%     u    - Input data, (N x n_u) matrix. Column vector for SISO.
%            For multiple trajectories: (N x n_u x L) or cell array.
%            Use [] for time series (periodogram).
%
%   NAME-VALUE OPTIONS:
%     'Smoothing'     - Smoothing window length S (positive odd integer).
%                       Default: 1 (no smoothing).
%     'Frequencies'   - Frequency vector in rad/sample, in (0, pi].
%                       Default: 128 linearly spaced values.
%     'SampleTime'    - Sample time in seconds. Default: 1.0.
%
%   OUTPUTS:
%     result - Struct with fields:
%       .Frequency        - (n_f x 1) frequency vector, rad/sample
%       .FrequencyHz      - (n_f x 1) frequency vector, Hz
%       .Response         - (n_f x n_y x n_u) complex frequency response
%       .ResponseStd      - (n_f x n_y x n_u) standard deviation (NaN)
%       .NoiseSpectrum    - (n_f x n_y x n_y) noise spectrum
%       .NoiseSpectrumStd - (n_f x n_y x n_y) standard deviation (NaN)
%       .Coherence        - [] (not applicable for ETFE)
%       .SampleTime       - sample time in seconds
%       .WindowSize       - N (data length)
%       .DataLength       - number of samples N
%       .NumTrajectories  - number of trajectories L
%       .Method           - 'sidFreqETFE'
%
%   EXAMPLES:
%     N = 1000; u = randn(N, 1);
%     y = filter([1], [1 -0.9], u) + 0.1*randn(N, 1);
%     result = sidFreqETFE(y, u, 'Smoothing', 5);
%     sidBodePlot(result);
%
%   ALGORITHM:
%     1. Compute DFTs Y(w) and U(w) of the data signals.
%     2. Form raw ETFE: G(w) = Y(w) / U(w) (SISO) or matrix division (MIMO).
%     3. Optionally smooth G with a length-S boxcar window.
%     4. Noise spectrum: Phi_v(w) = (1/N) * |Y(w) - G(w) * U(w)|^2.
%     5. Time series mode: periodogram Phi_y(w) = (1/N) * |Y(w)|^2.
%
%   REFERENCES:
%     Ljung, L. "System Identification: Theory for the User", 2nd ed.,
%     Prentice Hall, 1999. Sections 2.3, 6.3.
%
%   SPECIFICATION:
%     SPEC.md §4 — Empirical Transfer Function Estimate
%
%   See also: sidFreqBT, sidFreqBTFDR, sidBodePlot, sidSpectrumPlot
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
    [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u);

    defs.Smoothing = 1;
    defs.Frequencies = [];
    defs.SampleTime = 1.0;
    opts = sidParseOptions(defs, varargin);
    S = opts.Smoothing;
    freqs = opts.Frequencies;
    Ts = opts.SampleTime;

    % ---- Validate parameters ----
    if Ts <= 0
        error('sid:badTs', 'Sample time must be positive.');
    end

    if S < 1 || S ~= round(S)
        error('sid:badSmoothing', 'Smoothing parameter S must be a positive integer.');
    end
    if S > 1 && mod(S, 2) == 0
        error('sid:badSmoothing', 'Smoothing parameter S must be odd.');
    end

    if isempty(freqs)
        freqs = (1:128)' * pi / 128;
        useFFT = true;
    else
        freqs = freqs(:);
        useFFT = sidIsDefaultFreqs(freqs, length(freqs));
    end

    nf = length(freqs);

    if any(freqs <= 0) || any(freqs > pi)
        error('sid:badFreqs', 'Frequencies must be in the range (0, pi] rad/sample.');
    end

    % ---- Compute DFTs (SPEC.md §4.1) ----
    % Y(w) = sum_{n=1}^{N} y(n) e^{-jwn}, U(w) analogously
    if nTraj == 1
        Ydft = sidDFT(y, freqs, useFFT);    % (nf x ny)
        if ~isTimeSeries
            Udft = sidDFT(u, freqs, useFFT); % (nf x nu)
        end
    else
        % Multi-trajectory: compute DFTs per trajectory. For the
        % input-output case, cross-periodograms are averaged before
        % forming the ratio (H1 estimator, SPEC.md §4.1):
        %   G = (1/L) sum_l Y_l conj(U_l) / ((1/L) sum_l |U_l|^2)
        % Single-trajectory DFTs are stored only for the time-series path;
        % for I/O, cross- and auto-periodograms are accumulated in the
        % transfer function formation section below.
        Ydft = [];
        Udft = [];
    end

    % ---- Form transfer function and noise spectrum (SPEC.md §4.2-4.3) ----
    if isTimeSeries
        % Periodogram: Phi_y(w) = (1/L) * sum_l (1/N) |Y_l(w)|^2
        G = [];
        if nTraj == 1
            if ny == 1
                PhiV = (1/N) * abs(Ydft).^2;
            else
                PhiV = zeros(nf, ny, ny);
                for kk = 1:nf
                    Yk = Ydft(kk, :).';
                    PhiV(kk, :, :) = (1/N) * (Yk * Yk');
                end
                PhiV = real(PhiV);
            end
        else
            % Ensemble average of per-trajectory periodograms
            if ny == 1
                PhiV = zeros(nf, 1);
                for l = 1:nTraj
                    Yl = sidDFT(y(:, :, l), freqs, useFFT);
                    PhiV = PhiV + (1/N) * abs(Yl).^2;
                end
                PhiV = PhiV / nTraj;
            else
                PhiV = zeros(nf, ny, ny);
                for l = 1:nTraj
                    Yl = sidDFT(y(:, :, l), freqs, useFFT);
                    for kk = 1:nf
                        Yk = Yl(kk, :).';
                        PhiV(kk, :, :) = PhiV(kk, :, :) + reshape((1/N) * (Yk * Yk'), [1 ny ny]);
                    end
                end
                PhiV = real(PhiV) / nTraj;
            end
        end
        Coh = [];

    elseif ny == 1 && nu == 1
        % SISO: G(w) = Phi_yu(w) / Phi_u(w) via H1 estimator (SPEC.md §4.2)
        epsReg = 1e-10;

        if nTraj == 1
            % Single trajectory: G = Y / U directly
            Uabs = abs(Udft);
            Umax = max(Uabs);

            G = zeros(nf, 1);
            for kk = 1:nf
                if Uabs(kk) < epsReg * Umax
                    G(kk) = NaN + 1j*NaN;
                else
                    G(kk) = Ydft(kk) / Udft(kk);
                end
            end
        else
            % Multi-trajectory H1: average cross-periodograms (SPEC.md §4.1)
            % Phi_yu = (1/L) sum_l Y_l conj(U_l)
            % Phi_u  = (1/L) sum_l |U_l|^2
            PhiYU = zeros(nf, 1);
            PhiU  = zeros(nf, 1);
            for l = 1:nTraj
                Yl = sidDFT(y(:, :, l), freqs, useFFT);
                Ul = sidDFT(u(:, :, l), freqs, useFFT);
                PhiYU = PhiYU + Yl .* conj(Ul);
                PhiU  = PhiU  + abs(Ul).^2;
            end
            PhiYU = PhiYU / nTraj;
            PhiU  = PhiU  / nTraj;

            Umax = max(PhiU);
            G = zeros(nf, 1);
            for kk = 1:nf
                if PhiU(kk) < epsReg * Umax
                    G(kk) = NaN + 1j*NaN;
                else
                    G(kk) = PhiYU(kk) / PhiU(kk);
                end
            end
        end

        % Optional smoothing
        if S > 1
            G = boxcarSmooth(G, S);
        end

        % Noise spectrum: Phi_v(w) = (1/N) * |Y(w) - G(w) * U(w)|^2
        % For multi-trajectory: average per-trajectory noise periodograms
        if nTraj == 1
            residual = Ydft - G .* Udft;
            PhiV = (1/N) * abs(residual).^2;
        else
            PhiV = zeros(nf, 1);
            for l = 1:nTraj
                Yl = sidDFT(y(:, :, l), freqs, useFFT);
                Ul = sidDFT(u(:, :, l), freqs, useFFT);
                res_l = Yl - G .* Ul;
                PhiV = PhiV + (1/N) * abs(res_l).^2;
            end
            PhiV = PhiV / nTraj;
        end
        PhiV = max(PhiV, 0);
        Coh = [];

    else
        % MIMO: G(w) = Phi_yu(w) * Phi_u(w)^{-1} via H1 estimator (SPEC.md §4.2)
        epsReg = 1e-10;

        if nTraj == 1
            % Single trajectory: rank-1 cross/auto-periodograms
            PhiYU = zeros(nf, ny, nu);
            PhiU  = zeros(nf, nu, nu);
            for kk = 1:nf
                Yk = reshape(Ydft(kk, :), ny, 1);
                Uk = reshape(Udft(kk, :), nu, 1);
                PhiYU(kk, :, :) = Yk * Uk';
                PhiU(kk, :, :)  = Uk * Uk';
            end
        else
            % Multi-trajectory: average cross-periodograms (SPEC.md §4.1)
            PhiYU = zeros(nf, ny, nu);
            PhiU  = zeros(nf, nu, nu);
            for l = 1:nTraj
                Yl = sidDFT(y(:, :, l), freqs, useFFT);
                Ul = sidDFT(u(:, :, l), freqs, useFFT);
                for kk = 1:nf
                    Yk = reshape(Yl(kk, :), ny, 1);
                    Uk = reshape(Ul(kk, :), nu, 1);
                    PhiYU(kk, :, :) = PhiYU(kk, :, :) + reshape(Yk * Uk', [1 ny nu]);
                    PhiU(kk, :, :)  = PhiU(kk, :, :)  + reshape(Uk * Uk', [1 nu nu]);
                end
            end
            PhiYU = PhiYU / nTraj;
            PhiU  = PhiU  / nTraj;
        end

        G = zeros(nf, ny, nu);
        for kk = 1:nf
            PhiU_k  = reshape(PhiU(kk, :, :), nu, nu);
            PhiYU_k = reshape(PhiYU(kk, :, :), ny, nu);

            if nu == 1
                if abs(PhiU_k) < epsReg * max(abs(PhiU(:)))
                    G(kk, :, :) = NaN;
                else
                    G(kk, :, :) = PhiYU_k / PhiU_k;
                end
            else
                rc = rcond(PhiU_k);
                if rc < epsReg
                    G(kk, :, :) = NaN;
                else
                    G(kk, :, :) = PhiYU_k / PhiU_k;
                end
            end
        end

        % Optional smoothing (element-wise)
        if S > 1
            for ii = 1:ny
                for jj = 1:nu
                    G(:, ii, jj) = boxcarSmooth(G(:, ii, jj), S);
                end
            end
        end

        % Noise spectrum: average per-trajectory residual periodograms
        PhiV = zeros(nf, ny, ny);
        if nTraj == 1
            for kk = 1:nf
                Yk = reshape(Ydft(kk, :), ny, 1);
                Uk = reshape(Udft(kk, :), nu, 1);
                Gk = reshape(G(kk, :, :), ny, nu);
                res = Yk - Gk * Uk;
                PhiV(kk, :, :) = (1/N) * (res * res');
            end
        else
            for l = 1:nTraj
                Yl = sidDFT(y(:, :, l), freqs, useFFT);
                Ul = sidDFT(u(:, :, l), freqs, useFFT);
                for kk = 1:nf
                    Yk = reshape(Yl(kk, :), ny, 1);
                    Uk = reshape(Ul(kk, :), nu, 1);
                    Gk = reshape(G(kk, :, :), ny, nu);
                    res = Yk - Gk * Uk;
                    PhiV(kk, :, :) = PhiV(kk, :, :) + ...
                        reshape((1/N) * (res * res'), [1 ny ny]);
                end
            end
            PhiV = PhiV / nTraj;
        end
        PhiV = real(PhiV);
        Coh = [];
    end

    % ---- Uncertainty (SPEC.md §4.5) ----
    % ETFE has no closed-form asymptotic variance — return NaN.
    if isTimeSeries
        GStd = [];
        PhiVStd = nan(size(PhiV));
    else
        GStd = nan(size(G));
        PhiVStd = nan(size(PhiV));
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
    result.WindowSize       = N;
    result.DataLength       = N;
    result.NumTrajectories  = nTraj;
    result.Method           = 'sidFreqETFE';
end

function xSmooth = boxcarSmooth(x, S)
% BOXCARSMOOTH Apply length-S boxcar (moving average) smoothing.
%   At boundaries, uses available neighbors (shrinking window).
    nf = length(x);
    xSmooth = zeros(size(x));
    halfS = (S - 1) / 2;

    for k = 1:nf
        lo = max(1, k - halfS);
        hi = min(nf, k + halfS);
        xSmooth(k) = mean(x(lo:hi));
    end
end
