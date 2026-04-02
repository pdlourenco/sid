function result = sidSpectrogram(x, varargin)
% SIDSPECTROGRAM Short-time FFT spectrogram.
%
%   result = sidSpectrogram(x)
%   result = sidSpectrogram(x, 'WindowLength', 256, 'Overlap', 128)
%   result = sidSpectrogram(x, 'Window', 'hamming', 'NFFT', 512)
%
%   Computes the short-time Fourier transform (STFT) spectrogram of one
%   or more signals. Replaces the Signal Processing Toolbox spectrogram
%   function with no toolbox dependencies.
%
%   INPUTS:
%     x - Signal data, (N x n_ch) real matrix. Column vector for single
%         channel. Each column is treated as a separate channel.
%         For multiple trajectories: (N x n_ch x L) array. Power spectral
%         density is ensemble-averaged across trajectories within each segment.
%
%   NAME-VALUE OPTIONS:
%     'WindowLength'  - Segment length L. Default: 256.
%     'Overlap'       - Overlap P between segments, 0 <= P < L.
%                       Default: floor(L/2).
%     'NFFT'          - FFT length. Default: max(256, 2^nextpow2(L)).
%     'Window'        - Window type: 'hann' (default), 'hamming', 'rect',
%                       or a numeric vector of length L.
%     'SampleTime'    - Sample time in seconds. Default: 1.0.
%
%   OUTPUTS:
%     result - Struct with fields:
%       .Time          - (K x 1) center time of each segment (seconds)
%       .Frequency     - (n_bins x 1) frequency vector (Hz)
%       .FrequencyRad  - (n_bins x 1) frequency vector (rad/s)
%       .Power         - (n_bins x K x n_ch) power spectral density
%       .PowerDB       - (n_bins x K x n_ch) power in dB
%       .Complex       - (n_bins x K x n_ch) complex STFT coefficients
%       .SampleTime    - sample time in seconds
%       .WindowLength  - segment length L
%       .Overlap       - overlap P
%       .NFFT          - FFT length
%       .Method        - 'sidSpectrogram'
%
%   EXAMPLES:
%     % Spectrogram of a chirp signal
%     Fs = 1000; Ts = 1/Fs; N = 5000;
%     t = (0:N-1)' * Ts;
%     x = cos(2*pi * (50 + 100*t/max(t)) .* t);
%     result = sidSpectrogram(x, 'WindowLength', 256, 'SampleTime', Ts);
%     sidSpectrogramPlot(result);
%
%   ALGORITHM:
%     1. Divide signal into overlapping segments of length L
%     2. Apply time-domain window to each segment
%     3. Compute FFT of each windowed segment
%     4. Compute one-sided power spectral density
%     5. If L trajectories: ensemble-average PSD across realizations
%
%   REFERENCES:
%     Oppenheim, A.V. and Schafer, R.W. "Discrete-Time Signal Processing",
%     3rd ed., Prentice Hall, 2010.
%
%   SPECIFICATION:
%     SPEC.md §7 — Short-Time Spectral Analysis
%
%   See also: sidFreqMap, sidSpectrogramPlot
%
%   Changelog:
%   2026-03-28: First version by Pedro Lourenço.
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
    p = inputParser;
    p.addParameter('WindowLength', 256);
    p.addParameter('Overlap', []);
    p.addParameter('NFFT', []);
    p.addParameter('Window', 'hann');
    p.addParameter('SampleTime', 1.0);
    p.parse(varargin{:});
    opts = p.Results;

    L  = opts.WindowLength;
    Ts = opts.SampleTime;

    % ---- Validate signal ----
    if isvector(x)
        x = x(:);
    end
    N = size(x, 1);
    nCh = size(x, 2);
    if ndims(x) == 3 %#ok<ISMAT>
        nTrajS = size(x, 3);
    else
        nTrajS = 1;
    end

    if ~isreal(x)
        error('sid:complexData', 'Complex data is not supported. Signal x must be real.');
    end
    if any(~isfinite(x(:)))
        error('sid:nonFinite', 'Signal x contains NaN or Inf values.');
    end
    if ~isscalar(L) || L < 1 || L ~= round(L)
        error('sid:invalidWindowLength', 'WindowLength must be a positive integer.');
    end
    if N < L
        error('sid:tooShort', ...
            'Signal length (%d) is shorter than WindowLength (%d).', N, L);
    end
    if ~isscalar(Ts) || Ts <= 0
        error('sid:invalidSampleTime', 'SampleTime must be a positive scalar.');
    end

    % ---- Defaults that depend on L ----
    if isempty(opts.Overlap)
        P = floor(L / 2);
    else
        P = opts.Overlap;
    end
    if isempty(opts.NFFT)
        nfft = max(256, 2^nextpow2(L));
    else
        nfft = opts.NFFT;
    end

    % ---- Validate derived parameters ----
    if ~isscalar(P) || P < 0 || P >= L || P ~= round(P)
        error('sid:invalidOverlap', 'Overlap must be an integer in [0, L-1].');
    end
    if ~isscalar(nfft) || nfft < L || nfft ~= round(nfft)
        error('sid:invalidNFFT', 'NFFT must be an integer >= WindowLength.');
    end

    % ---- Build window vector ----
    w = buildWindow(opts.Window, L);

    % ---- Segmentation ----
    step = L - P;
    K = floor((N - L) / step) + 1;
    if K < 1
        error('sid:tooFewSegments', ...
            'Data too short for even one segment with L=%d, P=%d.', L, P);
    end

    Fs = 1 / Ts;
    nBins = floor(nfft / 2) + 1;
    S1 = sum(w .^ 2);  % window power normalization

    % ---- Pre-allocate ----
    stftCoeffs = zeros(nBins, K, nCh);
    Pxx = zeros(nBins, K, nCh);

    % ---- Compute STFT ----
    for ch = 1:nCh
        for k = 1:K
            startIdx = (k - 1) * step + 1;
            PkSum = zeros(nBins, 1);
            XSum  = zeros(nBins, 1);

            for lt = 1:nTrajS
                if nTrajS > 1
                    seg = x(startIdx:startIdx + L - 1, ch, lt) .* w;
                else
                    seg = x(startIdx:startIdx + L - 1, ch) .* w;
                end
                X = fft(seg, nfft);
                X = X(1:nBins);
                XSum = XSum + X;

                % One-sided PSD for this trajectory
                Pk = (1 / (Fs * S1)) * abs(X) .^ 2;
                % Double positive-frequency bins (not DC or Nyquist)
                if mod(nfft, 2) == 0
                    Pk(2:end-1) = 2 * Pk(2:end-1);
                else
                    Pk(2:end) = 2 * Pk(2:end);
                end
                PkSum = PkSum + Pk;
            end

            stftCoeffs(:, k, ch) = XSum / nTrajS;
            Pxx(:, k, ch) = PkSum / nTrajS;
        end
    end

    % ---- Time vector (center of each segment) ----
    timeVec = ((0:K-1)' * step + L / 2) * Ts;

    % ---- Frequency vector ----
    freqHz  = (0:nBins-1)' * Fs / nfft;
    freqRad = 2 * pi * freqHz;

    % ---- Power in dB ----
    PxxDB = 10 * log10(max(Pxx, eps));

    % ---- Pack result ----
    result.Time          = timeVec;
    result.Frequency     = freqHz;
    result.FrequencyRad  = freqRad;
    result.Power         = Pxx;
    result.PowerDB       = PxxDB;
    result.Complex       = stftCoeffs;
    result.SampleTime    = Ts;
    result.WindowLength  = L;
    result.Overlap       = P;
    result.NFFT          = nfft;
    result.NumTrajectories = nTrajS;
    result.Method        = 'sidSpectrogram';
end

function w = buildWindow(winSpec, L)
% BUILDWINDOW Create a window vector of length L from a specification.
    if isnumeric(winSpec)
        w = winSpec(:);
        if length(w) ~= L
            error('sid:windowSizeMismatch', ...
                'Custom window vector length (%d) must equal WindowLength (%d).', ...
                length(w), L);
        end
        return;
    end

    if ~ischar(winSpec)
        error('sid:invalidWindow', ...
            'Window must be ''hann'', ''hamming'', ''rect'', or a numeric vector.');
    end

    n = (0:L-1)';
    switch lower(winSpec)
        case 'hann'
            w = 0.5 * (1 - cos(2 * pi * n / (L - 1)));
        case 'hamming'
            w = 0.54 - 0.46 * cos(2 * pi * n / (L - 1));
        case {'rect', 'rectangular'}
            w = ones(L, 1);
        otherwise
            error('sid:invalidWindow', ...
                'Unknown window type ''%s''. Use ''hann'', ''hamming'', or ''rect''.', ...
                winSpec);
    end
end
