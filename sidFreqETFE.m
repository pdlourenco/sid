function result = sidFreqETFE(y, u, varargin)
%SIDFREQETFE Empirical transfer function estimate.
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
%     u    - Input data, (N x n_u) matrix. Column vector for SISO.
%            Use [] for time series (periodogram).
%
%   NAME-VALUE OPTIONS:
%     'Smoothing'     - Smoothing window length S (positive odd integer).
%                       Default: 1 (no smoothing).
%     'Frequencies'   - Frequency vector in rad/sample, in (0, pi].
%                       Default: 128 linearly spaced values.
%     'SampleTime'    - Sample time in seconds. Default: 1.0.
%
%   See also: sidFreqBT, sidFreqBTFDR
%
%   Example:
%   TODO add example code here
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
%   OUTPUT:
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
%       .Method           - 'sidFreqETFE'
%
%   ALGORITHM:
%     1. Compute DFTs Y(w) and U(w) of the data signals.
%     2. Form raw ETFE: G(w) = Y(w) / U(w) (SISO) or matrix division (MIMO).
%     3. Optionally smooth G with a length-S boxcar window.
%     4. Noise spectrum: Phi_v(w) = (1/N) * |Y(w) - G(w) * U(w)|^2.
%     5. Time series mode: periodogram Phi_y(w) = (1/N) * |Y(w)|^2.
%
%   REFERENCE:
%     Ljung, L. "System Identification: Theory for the User", 2nd ed.,
%     Prentice Hall, 1999. Sections 2.3, 6.3.
%
%   See also: sidFreqBT, sidFreqBTFDR, sidBodePlot, sidSpectrumPlot

    % ---- Parse inputs ----
    [y, u, N, ny, nu, isTimeSeries] = sidValidateData(y, u);

    S = 1;
    freqs = [];
    Ts = 1.0;

    args = varargin;
    k = 1;
    while k <= length(args)
        if ischar(args{k})
            switch lower(args{k})
                case 'smoothing'
                    S = args{k+1};
                    k = k + 2;
                case 'frequencies'
                    freqs = args{k+1};
                    k = k + 2;
                case 'sampletime'
                    Ts = args{k+1};
                    k = k + 2;
                otherwise
                    error('sid:unknownOption', 'Unknown option: %s', args{k});
            end
        else
            error('sid:badInput', 'Expected a string option name at position %d.', k);
        end
    end

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
        useFFT = isDefaultFreqs(freqs, length(freqs));
    end

    nf = length(freqs);

    if any(freqs <= 0) || any(freqs > pi)
        error('sid:badFreqs', 'Frequencies must be in the range (0, pi] rad/sample.');
    end

    % ---- Compute DFTs ----
    Ydft = sidDFT(y, freqs, useFFT);    % (nf x ny)

    if ~isTimeSeries
        Udft = sidDFT(u, freqs, useFFT);  % (nf x nu)
    end

    % ---- Form transfer function and noise spectrum ----
    if isTimeSeries
        % Periodogram: Phi_y(w) = (1/N) * |Y(w)|^2
        G = [];
        if ny == 1
            PhiV = (1/N) * abs(Ydft).^2;
        else
            PhiV = zeros(nf, ny, ny);
            for kk = 1:nf
                Yk = Ydft(kk, :).';  % (ny x 1)
                PhiV(kk, :, :) = (1/N) * (Yk * Yk');
            end
            PhiV = real(PhiV);
        end
        Coh = [];

    elseif ny == 1 && nu == 1
        % SISO ETFE
        epsReg = 1e-10;
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

        % Optional smoothing
        if S > 1
            G = boxcarSmooth(G, S);
        end

        % Noise spectrum: Phi_v(w) = (1/N) * |Y(w) - G(w) * U(w)|^2
        residual = Ydft - G .* Udft;
        PhiV = (1/N) * abs(residual).^2;
        PhiV = max(PhiV, 0);
        Coh = [];

    else
        % MIMO ETFE
        G = zeros(nf, ny, nu);
        PhiV = zeros(nf, ny, ny);
        epsReg = 1e-10;

        for kk = 1:nf
            Uk = reshape(Udft(kk, :), nu, 1);     % (nu x 1) — but we need (nu x 1) for each freq
            Yk = reshape(Ydft(kk, :), ny, 1);

            % For MIMO: G = Y * U' * inv(U * U') — but we only have single-snapshot DFT
            % With a single snapshot, U is (nu x 1), so U*U' is rank-1 (nu x nu).
            % For nu > 1 this is singular. MIMO ETFE requires multiple snapshots or smoothing.
            % For single snapshot with nu=1: works fine.
            % For nu > 1: use pseudoinverse.
            if nu == 1
                if abs(Uk) < epsReg * max(abs(Udft(:)))
                    G(kk, :, :) = NaN;
                else
                    G(kk, :, :) = Yk / Uk;
                end
            else
                % Use left division: G = Yk * pinv(Uk') where Uk' is (1 x nu)
                % Actually for single snapshot: G(w) = Y(w) * U(w)^H / (U(w)^H * U(w))
                UkH = Uk';
                denom = UkH * Uk;  % scalar
                if abs(denom) < epsReg * max(abs(Udft(:)))^2
                    G(kk, :, :) = NaN;
                else
                    G(kk, :, :) = (Yk * UkH) / denom;
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

        % Noise spectrum
        for kk = 1:nf
            Uk = reshape(Udft(kk, :), nu, 1);
            Yk = reshape(Ydft(kk, :), ny, 1);
            Gk = reshape(G(kk, :, :), ny, nu);
            res = Yk - Gk * Uk;
            PhiV(kk, :, :) = (1/N) * (res * res');
        end
        PhiV = real(PhiV);
        Coh = [];
    end

    % ---- Uncertainty ----
    % ETFE has no simple asymptotic uncertainty formula.
    % Return NaN arrays to indicate uncertainty is not available.
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
    result.Method           = 'sidFreqETFE';
end


function xSmooth = boxcarSmooth(x, S)
%BOXCARSMOOTH Apply length-S boxcar (moving average) smoothing.
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


function tf = isDefaultFreqs(freqs, nf)
%ISDEFAULTFREQS Check if frequency vector matches the default linear grid.
    if nf ~= 128
        tf = false;
        return;
    end
    defaultFreqs = (1:128)' * pi / 128;
    tf = max(abs(freqs(:) - defaultFreqs)) < 1e-12;
end
