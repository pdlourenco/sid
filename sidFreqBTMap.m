function result = sidFreqBTMap(y, u, varargin)
%SIDFREQBTMAP Time-varying frequency response map via Blackman-Tukey.
%
%   result = sidFreqBTMap(y, u)
%   result = sidFreqBTMap(y, [], 'SegmentLength', 256)
%   result = sidFreqBTMap(y, u, 'SegmentLength', L, 'Overlap', P)
%
%   Applies sidFreqBT to overlapping segments of data, producing a
%   time-varying frequency response estimate G(w,t). For an LTI system
%   the map is constant along time; for an LTV system it reveals how the
%   transfer function, noise spectrum, and coherence evolve.
%
%   INPUTS:
%     y - Output data, (N x n_y) matrix. Column vector for SISO.
%     u - Input data, (N x n_u) matrix, or [] for time series mode.
%
%   NAME-VALUE OPTIONS:
%     'SegmentLength'  - Number of samples per segment L.
%                        Default: min(floor(N/4), 256).
%     'Overlap'        - Overlap P between segments, 0 <= P < L.
%                        Default: floor(L/2).
%     'WindowSize'     - Hann window lag size M for each segment.
%                        Default: min(floor(L/10), 30).
%     'Frequencies'    - Frequency vector in rad/sample, in (0, pi].
%                        Default: [] (128-point linear grid via sidFreqBT).
%     'SampleTime'     - Sample time in seconds. Default: 1.0.
%
%   OUTPUT:
%     result - Struct with fields:
%       .Time             - (K x 1) center time of each segment (seconds)
%       .Frequency        - (n_f x 1) frequency vector, rad/sample
%       .FrequencyHz      - (n_f x 1) frequency vector, Hz
%       .Response         - (n_f x K [x n_y x n_u]) complex, [] in time series
%       .ResponseStd      - (n_f x K [x n_y x n_u]) real, [] in time series
%       .NoiseSpectrum    - (n_f x K [x n_y x n_y]) real
%       .NoiseSpectrumStd - (n_f x K [x n_y x n_y]) real
%       .Coherence        - (n_f x K) real (SISO), [] for MIMO/time series
%       .SampleTime       - sample time
%       .SegmentLength    - L
%       .Overlap          - P
%       .WindowSize       - M
%       .Method           - 'sidFreqBTMap'
%
%   ALGORITHM:
%     1. Divide data into K overlapping segments of length L
%     2. Run sidFreqBT on each segment
%     3. Collect per-segment results into time-frequency arrays
%
%   See also: sidFreqBT, sidSpectrogram, sidMapPlot
%
%   Example:
%     % Time-varying analysis of a SISO system
%     N = 4000; u = randn(N, 1);
%     y = filter([1], [1 -0.9], u) + 0.1*randn(N, 1);
%     result = sidFreqBTMap(y, u, 'SegmentLength', 512);
%     sidMapPlot(result);
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

    % ---- Validate data ----
    [y, u, N, ny, nu, isTimeSeries] = sidValidateData(y, u);

    % ---- Parse options ----
    p = inputParser;
    p.addParameter('SegmentLength', min(floor(N / 4), 256));
    p.addParameter('Overlap', []);
    p.addParameter('WindowSize', []);
    p.addParameter('Frequencies', []);
    p.addParameter('SampleTime', 1.0);
    p.parse(varargin{:});
    opts = p.Results;

    L  = opts.SegmentLength;
    Ts = opts.SampleTime;

    % ---- Validate segment length ----
    if ~isscalar(L) || L < 4 || L ~= round(L)
        error('sid:invalidSegmentLength', 'SegmentLength must be an integer >= 4.');
    end
    if L > N
        error('sid:segmentTooLong', ...
            'SegmentLength (%d) exceeds data length (%d).', L, N);
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
    if isempty(opts.WindowSize)
        M = min(floor(L / 10), 30);
    else
        M = opts.WindowSize;
    end
    freqs = opts.Frequencies;

    % ---- Validate derived parameters ----
    if ~isscalar(P) || P < 0 || P >= L || P ~= round(P)
        error('sid:invalidOverlap', 'Overlap must be an integer in [0, L-1].');
    end
    if ~isscalar(M) || M < 2 || M ~= round(M)
        error('sid:invalidWindowSize', 'WindowSize must be an integer >= 2.');
    end
    if L <= 2 * M
        error('sid:segmentTooShort', ...
            'SegmentLength (%d) must be greater than 2*WindowSize (%d).', L, 2*M);
    end

    % ---- Segmentation ----
    step = L - P;
    K = floor((N - L) / step) + 1;
    if K < 1
        error('sid:tooFewSegments', ...
            'Data too short for even one segment with L=%d, P=%d.', L, P);
    end

    % ---- Build sidFreqBT argument list ----
    btArgs = {'WindowSize', M, 'SampleTime', Ts};
    if ~isempty(freqs)
        btArgs = [btArgs, {'Frequencies', freqs}];
    end

    % ---- Run sidFreqBT on first segment to get dimensions ----
    s1 = (0) * step + 1;
    e1 = s1 + L - 1;
    yk = y(s1:e1, :);
    if isTimeSeries
        uk = [];
    else
        uk = u(s1:e1, :);
    end
    r1 = sidFreqBT(yk, uk, btArgs{:});
    nf = length(r1.Frequency);

    % ---- Pre-allocate based on first segment result ----
    if isTimeSeries
        % NoiseSpectrum holds Phi_y
        nsSize = size(r1.NoiseSpectrum);
        NS    = zeros([nf, K, nsSize(2:end)]);
        NSStd = zeros([nf, K, nsSize(2:end)]);
        G     = [];
        GStd  = [];
        Coh   = [];
    else
        rSize = size(r1.Response);
        nsSize = size(r1.NoiseSpectrum);
        G     = zeros([nf, K, rSize(2:end)], 'like', 1i);
        GStd  = zeros([nf, K, rSize(2:end)]);
        NS    = zeros([nf, K, nsSize(2:end)]);
        NSStd = zeros([nf, K, nsSize(2:end)]);
        if ~isempty(r1.Coherence)
            Coh = zeros(nf, K);
        else
            Coh = [];
        end
    end

    % ---- Store first segment ----
    storeSegment(1, r1);

    % ---- Loop over remaining segments ----
    for k = 2:K
        startIdx = (k - 1) * step + 1;
        endIdx   = startIdx + L - 1;
        yk = y(startIdx:endIdx, :);
        if isTimeSeries
            uk = [];
        else
            uk = u(startIdx:endIdx, :);
        end
        rk = sidFreqBT(yk, uk, btArgs{:});
        storeSegment(k, rk);
    end

    % ---- Time vector ----
    timeVec = ((0:K-1)' * step + L / 2) * Ts;

    % ---- Pack result ----
    result.Time             = timeVec;
    result.Frequency        = r1.Frequency;
    result.FrequencyHz      = r1.FrequencyHz;
    result.Response         = G;
    result.ResponseStd      = GStd;
    result.NoiseSpectrum    = NS;
    result.NoiseSpectrumStd = NSStd;
    result.Coherence        = Coh;
    result.SampleTime       = Ts;
    result.SegmentLength    = L;
    result.Overlap          = P;
    result.WindowSize       = M;
    result.Method           = 'sidFreqBTMap';

    % ---- Nested helper to store segment results ----
    function storeSegment(idx, rk)
        if ~isTimeSeries
            if ndims(rk.Response) <= 2 %#ok<ISMAT>
                G(:, idx) = rk.Response;
                GStd(:, idx) = rk.ResponseStd;
            else
                G(:, idx, :, :) = rk.Response;
                GStd(:, idx, :, :) = rk.ResponseStd;
            end
            if ~isempty(Coh)
                Coh(:, idx) = rk.Coherence;
            end
        end
        if ndims(rk.NoiseSpectrum) <= 2 %#ok<ISMAT>
            NS(:, idx) = rk.NoiseSpectrum;
            NSStd(:, idx) = rk.NoiseSpectrumStd;
        else
            NS(:, idx, :, :) = rk.NoiseSpectrum;
            NSStd(:, idx, :, :) = rk.NoiseSpectrumStd;
        end
    end
end
