function result = sidFreqMap(y, u, varargin)
% SIDFREQMAP Time-varying frequency response map.
%
%   result = sidFreqMap(y, u)
%   result = sidFreqMap(y, [], 'SegmentLength', 256)
%   result = sidFreqMap(y, u, 'Algorithm', 'welch')
%
%   Estimates a time-varying frequency response G(w,t) by applying spectral
%   analysis to overlapping segments of input-output data. For an LTI system
%   the map is constant along time; for an LTV system it reveals how the
%   transfer function, noise spectrum, and coherence evolve.
%
%   Two algorithms are supported:
%     'bt'    (default) - Blackman-Tukey correlogram via sidFreqBT
%     'welch' - Welch averaged periodogram (tfestimate compatible)
%
%   INPUTS:
%     y - Output data, (N x n_y) matrix. Column vector for SISO.
%         For multiple trajectories: (N x n_y x L) array or cell array.
%         Spectral estimates within each segment are ensemble-averaged.
%     u - Input data, (N x n_u) matrix, or [] for time series mode.
%         For multiple trajectories: (N x n_u x L) or cell array.
%
%   NAME-VALUE OPTIONS (common):
%     'SegmentLength'  - Number of samples per segment L.
%                        Default: min(floor(N/4), 256).
%     'Overlap'        - Overlap P between segments, 0 <= P < L.
%                        Default: floor(L/2).
%     'Algorithm'      - 'bt' (default) or 'welch'.
%     'SampleTime'     - Sample time in seconds. Default: 1.0.
%
%   BT-specific options:
%     'WindowSize'     - Hann lag window size M. Default: min(floor(L/10), 30).
%     'Frequencies'    - Frequency vector in rad/sample, in (0, pi].
%                        Default: 128-point linear grid.
%
%   Welch-specific options:
%     'SubSegmentLength' - Sub-segment length within each segment.
%                          Default: floor(L/4.5).
%     'SubOverlap'       - Sub-segment overlap. Default: floor(SubSegmentLength/2).
%     'Window'           - 'hann' (default), 'hamming', 'rect', or numeric vector.
%     'NFFT'             - FFT length. Default: max(256, 2^nextpow2(SubSegmentLength)).
%
%   OUTPUTS:
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
%       .WindowSize       - M (BT) or [] (Welch)
%       .Algorithm        - 'bt' or 'welch'
%       .NumTrajectories  - number of trajectories L
%       .Method           - 'sidFreqMap'
%
%   EXAMPLES:
%     % Time-varying frequency map (Blackman-Tukey)
%     N = 4000; u = randn(N, 1);
%     y = filter([1], [1 -0.9], u) + 0.1*randn(N, 1);
%     result = sidFreqMap(y, u, 'SegmentLength', 512);
%     sidMapPlot(result);
%
%   SPECIFICATION:
%     SPEC.md §6 — Time-Varying Frequency Response Map
%
%   See also: sidFreqBT, sidSpectrogram, sidMapPlot
%
%   Changelog:
%   2026-03-29: Refactored from sidFreqBTMap, added Welch support.
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

    % ---- Validate data (preserve variable-length cells for SPEC.md §6.2) ----
    [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u, true);

    % If the input was a variable-length cell array, y and u remain cell
    % arrays here.  We detect that and fall into the per-segment filtering
    % path (SPEC.md §6.2): at each segment k, only trajectories that span
    % segment k contribute to the ensemble.
    varLen = iscell(y);
    if varLen
        horizons = zeros(nTraj, 1);
        for ll = 1:nTraj
            horizons(ll) = size(y{ll}, 1);
        end
    else
        horizons = [];
    end

    % ---- Parse options ----
    defs.SegmentLength = min(floor(N / 4), 256);
    defs.Overlap = [];
    defs.Algorithm = 'bt';
    defs.SampleTime = 1.0;
    % BT options
    defs.WindowSize = [];
    defs.Frequencies = [];
    % Welch options
    defs.SubSegmentLength = [];
    defs.SubOverlap = [];
    defs.Window = 'hann';
    defs.NFFT = [];
    opts = sidParseOptions(defs, varargin);

    L  = opts.SegmentLength;
    Ts = opts.SampleTime;
    algorithm = lower(opts.Algorithm);

    % ---- Validate common parameters ----
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
    if ~ismember(algorithm, {'bt', 'welch'})
        error('sid:invalidAlgorithm', ...
            'Algorithm must be ''bt'' or ''welch''. Got ''%s''.', algorithm);
    end

    % ---- Defaults that depend on L ----
    if isempty(opts.Overlap)
        P = floor(L / 2);
    else
        P = opts.Overlap;
    end

    % ---- Validate overlap ----
    if ~isscalar(P) || P < 0 || P >= L || P ~= round(P)
        error('sid:invalidOverlap', 'Overlap must be an integer in [0, L-1].');
    end

    % ---- Algorithm-specific parameter setup ----
    if strcmp(algorithm, 'bt')
        if isempty(opts.WindowSize)
            M = min(floor(L / 10), 30);
        else
            M = opts.WindowSize;
        end
        freqs = opts.Frequencies;

        if ~isscalar(M) || M < 2 || M ~= round(M)
            error('sid:invalidWindowSize', 'WindowSize must be an integer >= 2.');
        end
        if L <= 2 * M
            error('sid:segmentTooShort', ...
                'SegmentLength (%d) must be greater than 2*WindowSize (%d).', L, 2*M);
        end

        btArgs = {'WindowSize', M, 'SampleTime', Ts};
        if ~isempty(freqs)
            btArgs = [btArgs, {'Frequencies', freqs}];
        end
    else
        % Welch parameters
        if isempty(opts.SubSegmentLength)
            Lsub = floor(L / 4.5);
        else
            Lsub = opts.SubSegmentLength;
        end
        if isempty(opts.SubOverlap)
            Psub = floor(Lsub / 2);
        else
            Psub = opts.SubOverlap;
        end
        if isempty(opts.NFFT)
            nfft = max(256, 2^nextpow2(Lsub));
        else
            nfft = opts.NFFT;
        end
        winType = opts.Window;

        if ~isscalar(Lsub) || Lsub < 2 || Lsub ~= round(Lsub)
            error('sid:invalidSubSegmentLength', 'SubSegmentLength must be an integer >= 2.');
        end
        if Lsub > L
            error('sid:subSegmentTooLong', ...
                'SubSegmentLength (%d) exceeds SegmentLength (%d).', Lsub, L);
        end
        if ~isscalar(Psub) || Psub < 0 || Psub >= Lsub || Psub ~= round(Psub)
            error('sid:invalidSubOverlap', 'SubOverlap must be in [0, SubSegmentLength-1].');
        end
        if ~isscalar(nfft) || nfft < Lsub || nfft ~= round(nfft)
            error('sid:invalidNFFT', ...
                'NFFT (%d) must be an integer >= SubSegmentLength (%d).', ...
                nfft, Lsub);
        end

        M = [];  % not applicable for Welch
    end

    % ---- Outer segmentation (SPEC.md §6.2) ----
    % K = floor((N - L) / step) + 1 segments, stride = L - P
    step = L - P;
    K = floor((N - L) / step) + 1;
    if K < 1
        error('sid:tooFewSegments', ...
            'Data too short for even one segment with L=%d, P=%d.', L, P);
    end

    % ---- Track number of trajectories used per segment (SPEC.md §6.8) ----
    nTrajPerSeg = zeros(K, 1);

    % ---- Run inner estimator on first segment to get dimensions ----
    [yk, uk, nAct] = extractSegment(1);
    nTrajPerSeg(1) = nAct;

    if strcmp(algorithm, 'bt')
        r1 = sidFreqBT(yk, uk, btArgs{:});
    else
        r1 = welchEstimate(yk, uk, isTimeSeries, ny, nu, Lsub, Psub, nfft, winType, Ts);
    end
    nf = length(r1.Frequency);

    % ---- Pre-allocate based on first segment result ----
    if isTimeSeries
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

    % ---- Loop over remaining segments (SPEC.md §6.3) ----
    for k = 2:K
        [yk, uk, nAct] = extractSegment(k);
        nTrajPerSeg(k) = nAct;
        if strcmp(algorithm, 'bt')
            rk = sidFreqBT(yk, uk, btArgs{:});
        else
            rk = welchEstimate(yk, uk, isTimeSeries, ny, nu, Lsub, Psub, nfft, winType, Ts);
        end
        storeSegment(k, rk);
    end

    % ---- Time vector (SPEC.md §6.2) ----
    % t_k = (k * step + L/2) * Ts — center of each segment
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
    result.Algorithm        = algorithm;
    % NumTrajectories: scalar when every segment uses the same count,
    % (K x 1) vector when variable-length trajectories drop in/out
    % (SPEC.md §6.8).
    if varLen && any(nTrajPerSeg ~= nTrajPerSeg(1))
        result.NumTrajectories = nTrajPerSeg;
    elseif varLen
        result.NumTrajectories = nTrajPerSeg(1);
    else
        result.NumTrajectories = nTraj;
    end
    result.Method           = 'sidFreqMap';

    % ---- Nested helper to extract one segment (SPEC.md §6.2 / §6.3) ----
    % For uniform-length inputs, slice the 3-D / 2-D arrays directly.
    % For variable-length cell inputs, collect only those trajectories that
    % fully span the segment [startIdx, endIdx], matching the spec:
    %     "At each segment k, only trajectories that span segment k
    %      contribute to the ensemble."
    function [yk_, uk_, nAct_] = extractSegment(k_)
        startIdx_ = (k_ - 1) * step + 1;
        endIdx_   = startIdx_ + L - 1;
        if varLen
            active = find(horizons >= endIdx_);
            nAct_ = numel(active);
            if nAct_ == 0
                error('sid:noActiveTraj', ...
                    ['Segment starting at sample %d has no trajectories ' ...
                     'spanning it. Check segment length / overlap vs horizons.'], ...
                    startIdx_);
            end
            if nAct_ == 1
                % 2-D form so the inner estimator stays on its
                % single-trajectory code path (sidCov expects 2-D).
                li = active(1);
                yk_ = y{li}(startIdx_:endIdx_, :);
                if isTimeSeries
                    uk_ = [];
                else
                    uk_ = u{li}(startIdx_:endIdx_, :);
                end
            else
                yk_ = zeros(L, ny, nAct_);
                if isTimeSeries
                    uk_ = [];
                else
                    uk_ = zeros(L, nu, nAct_);
                end
                for ii_ = 1:nAct_
                    li = active(ii_);
                    yk_(:, :, ii_) = y{li}(startIdx_:endIdx_, :);
                    if ~isTimeSeries
                        uk_(:, :, ii_) = u{li}(startIdx_:endIdx_, :);
                    end
                end
            end
            return;
        end

        if nTraj > 1
            yk_ = y(startIdx_:endIdx_, :, :);
            if isTimeSeries
                uk_ = [];
            else
                uk_ = u(startIdx_:endIdx_, :, :);
            end
        else
            yk_ = y(startIdx_:endIdx_, :);
            if isTimeSeries
                uk_ = [];
            else
                uk_ = u(startIdx_:endIdx_, :);
            end
        end
        nAct_ = nTraj;
    end

    % ---- Nested helper to store segment results ----
    function storeSegment(idx, rk)
        if ~isTimeSeries
            if ny == 1 && nu == 1
                G(:, idx) = rk.Response(:);
                GStd(:, idx) = rk.ResponseStd(:);
            else
                G(:, idx, :, :) = reshape(rk.Response, [], 1, ny, nu);
                GStd(:, idx, :, :) = reshape(rk.ResponseStd, [], 1, ny, nu);
            end
            if ~isempty(Coh)
                Coh(:, idx) = rk.Coherence;
            end
        end
        if ny == 1 || isTimeSeries
            NS(:, idx) = rk.NoiseSpectrum(:);
            NSStd(:, idx) = rk.NoiseSpectrumStd(:);
        else
            NS(:, idx, :, :) = reshape(rk.NoiseSpectrum, [], 1, ny, ny);
            NSStd(:, idx, :, :) = reshape(rk.NoiseSpectrumStd, [], 1, ny, ny);
        end
    end

end

% ========================================================================
%  WELCH INNER ESTIMATOR
% ========================================================================

function result = welchEstimate(y, u, isTimeSeries, ny, nu, Lsub, Psub, nfft, winType, Ts)
% WELCHESTIMATE Welch averaged periodogram estimate within one segment.
%
%   Returns a result struct with the same fields as sidFreqBT output.
%   Supports multi-trajectory input: y is (Lseg x ny x L), u is (Lseg x nu x L).

    Lseg = size(y, 1);
    if ndims(y) == 3 %#ok<ISMAT>
        nTrajW = size(y, 3);
    else
        nTrajW = 1;
    end

    % ---- Build window ----
    if isnumeric(winType)
        w = winType(:);
        if length(w) ~= Lsub
            error('sid:invalidWindow', ...
                'Window vector length (%d) must match SubSegmentLength (%d).', length(w), Lsub);
        end
    else
        switch lower(winType)
            case 'hann'
                n = (0:Lsub-1)';
                w = 0.5 * (1 - cos(2 * pi * n / (Lsub - 1)));
            case 'hamming'
                n = (0:Lsub-1)';
                w = 0.54 - 0.46 * cos(2 * pi * n / (Lsub - 1));
            case 'rect'
                w = ones(Lsub, 1);
            otherwise
                error('sid:invalidWindow', ...
                    'Window must be ''hann'', ''hamming'', ''rect'', or a numeric vector.');
        end
    end

    S1 = sum(w .^ 2);  % window power normalization

    % ---- Sub-segmentation ----
    subStep = Lsub - Psub;
    J = floor((Lseg - Lsub) / subStep) + 1;
    if J < 1
        error('sid:tooFewSubSegments', ...
            'Segment too short for sub-segmentation with Lsub=%d, Psub=%d.', Lsub, Psub);
    end

    % ---- One-sided frequency grid (skip DC, up to Nyquist) ----
    nBins = floor(nfft / 2);  % number of one-sided bins (excluding DC)
    freqs = (1:nBins)' * (2 * pi / nfft);  % rad/sample, in (0, pi]

    % ---- Accumulate averaged periodograms ----
    PhiY  = zeros(nBins, ny, ny);
    if ~isTimeSeries
        PhiU  = zeros(nBins, nu, nu);
        PhiYU = zeros(nBins, ny, nu);
    end

    for lt = 1:nTrajW
        for j = 1:J
            s = (j - 1) * subStep + 1;
            e = s + Lsub - 1;

            % Windowed FFT of output
            if nTrajW > 1
                Yj = fft(bsxfun(@times, y(s:e, :, lt), w), nfft, 1);
            else
                Yj = fft(bsxfun(@times, y(s:e, :), w), nfft, 1);
            end
            Yj = Yj(2:nBins+1, :);  % one-sided, skip DC

            % Auto-spectrum of y
            for a = 1:ny
                for b = 1:ny
                    PhiY(:, a, b) = PhiY(:, a, b) + Yj(:, a) .* conj(Yj(:, b));
                end
            end

            if ~isTimeSeries
                % Windowed FFT of input
                if nTrajW > 1
                    Uj = fft(bsxfun(@times, u(s:e, :, lt), w), nfft, 1);
                else
                    Uj = fft(bsxfun(@times, u(s:e, :), w), nfft, 1);
                end
                Uj = Uj(2:nBins+1, :);

                % Auto-spectrum of u
                for a = 1:nu
                    for b = 1:nu
                        PhiU(:, a, b) = PhiU(:, a, b) + Uj(:, a) .* conj(Uj(:, b));
                    end
                end

                % Cross-spectrum y*u'
                for a = 1:ny
                    for b = 1:nu
                        PhiYU(:, a, b) = PhiYU(:, a, b) + Yj(:, a) .* conj(Uj(:, b));
                    end
                end
            end
        end
    end

    % ---- Average and normalize ----
    % Factor of 2 for one-sided spectrum (we exclude DC and use positive
    % frequencies only). This cancels in G = Pyu/Puu and coherence, but is
    % needed for correct PSD magnitude.
    % Total number of periodogram segments: J sub-segments × nTrajW trajectories
    Jtotal = J * nTrajW;
    PhiY = 2 * PhiY / (Jtotal * S1);
    if ~isTimeSeries
        PhiU  = 2 * PhiU  / (Jtotal * S1);
        PhiYU = 2 * PhiYU / (Jtotal * S1);
    end

    % ---- Degrees of freedom for uncertainty ----
    % For Hann window at 50% overlap: nu_dof ≈ 1.8 * J per trajectory
    % Multi-trajectory multiplies by nTrajW (independent realizations)
    overlapRatio = Psub / Lsub;
    if overlapRatio <= 0
        nuDof = 2 * J * nTrajW;
    else
        nuDof = max(2, 1.8 * J * nTrajW);  % conservative for Hann at ~50% overlap
    end

    % ---- Form transfer function, noise spectrum, coherence ----
    nf = nBins;
    if isTimeSeries
        G = [];
        GStd = [];
        PhiV = real(squeeze(PhiY));
        if ny == 1
            PhiV = PhiV(:);
        end
        PhiVStd = PhiV * sqrt(2 / nuDof);
        Coh = [];
    elseif ny == 1 && nu == 1
        % SISO
        PhiY_s  = PhiY(:);
        PhiU_s  = PhiU(:);
        PhiYU_s = PhiYU(:);
        G = PhiYU_s ./ PhiU_s;
        PhiV = real(PhiY_s) - abs(PhiYU_s).^2 ./ real(PhiU_s);
        PhiV = max(PhiV, 0);
        Coh = abs(PhiYU_s).^2 ./ (real(PhiY_s) .* real(PhiU_s));
        Coh = min(max(Coh, 0), 1);
        GStd = abs(G) .* sqrt((1 - Coh) ./ (Coh * nuDof));
        GStd(Coh < 1e-10) = NaN;
        PhiVStd = PhiV * sqrt(2 / nuDof);
    else
        % MIMO
        G = zeros(nf, ny, nu);
        PhiV = zeros(nf, ny, ny);
        for k = 1:nf
            PhiU_k  = reshape(PhiU(k, :, :), nu, nu);
            PhiYU_k = reshape(PhiYU(k, :, :), ny, nu);
            PhiY_k  = reshape(PhiY(k, :, :), ny, ny);
            G(k, :, :) = PhiYU_k / PhiU_k;
            PhiV(k, :, :) = PhiY_k - PhiYU_k / PhiU_k * PhiYU_k';
        end
        PhiV = real(PhiV);
        GStd = NaN(size(G));  % MIMO uncertainty not supported for Welch
        PhiVStd = abs(PhiV) * sqrt(2 / nuDof);
        Coh = [];
    end

    % ---- Pack result (matching sidFreqBT output structure) ----
    result.Frequency        = freqs;
    result.FrequencyHz      = freqs / (2 * pi * Ts);
    result.Response         = G;
    result.ResponseStd      = GStd;
    result.NoiseSpectrum    = PhiV;
    result.NoiseSpectrumStd = PhiVStd;
    result.Coherence        = Coh;
    result.SampleTime       = Ts;
    result.WindowSize       = [];
    result.DataLength       = Lseg;
    result.Method           = 'welch';
end
