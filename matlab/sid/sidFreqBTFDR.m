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
    % M_k = ceil(2*pi / R_k) — local window size at frequency k
    Mk = ceil(2 * pi ./ R);
    Mk = max(Mk, 2);                   % M >= 2
    Mk = min(Mk, floor(N / 2));        % M <= N/2

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

    elseif isSISO
        G = zeros(nf, 1);
        PhiV = zeros(nf, 1);
        GStd = zeros(nf, 1);
        PhiVStd = zeros(nf, 1);
        Coh = zeros(nf, 1);
        PhiUall = zeros(nf, 1);
        epsReg = 1e-10;

        % First pass: compute all spectral estimates
        PhiYall = zeros(nf, 1);
        PhiYUall = zeros(nf, 1);
        Wstore = cell(nf, 1);
        for kk = 1:nf
            Mk_k = Mk(kk);
            W = sidHannWin(Mk_k);
            Wstore{kk} = W;

            PhiYall(kk) = scalarSingleFreqDFT(Ryy(1:Mk_k+1), W, freqs(kk));
            PhiUall(kk) = scalarSingleFreqDFT(Ruu(1:Mk_k+1), W, freqs(kk));
            PhiYUall(kk) = scalarSingleFreqDFT(Ryu(1:Mk_k+1), W, freqs(kk), Ruy(1:Mk_k+1));
        end

        PhiUmax = max(abs(PhiUall));

        % Second pass: form G, PhiV, uncertainty
        for kk = 1:nf
            W = Wstore{kk};
            PhiY_k = PhiYall(kk);
            PhiU_k = PhiUall(kk);
            PhiYU_k = PhiYUall(kk);

            % G(w_k) = Phi_yu(w_k) / Phi_u(w_k) (SPEC.md §5.2)
            if abs(PhiU_k) < epsReg * PhiUmax
                G(kk) = NaN + 1j*NaN;
                PhiV(kk) = real(PhiY_k);
                Coh(kk) = 0;
                GStd(kk) = Inf;
            else
                G(kk) = PhiYU_k / PhiU_k;
                PhiV(kk) = max(real(PhiY_k) - abs(PhiYU_k)^2 / real(PhiU_k), 0);
                Coh(kk) = min(max(abs(PhiYU_k)^2 / (real(PhiY_k) * real(PhiU_k)), 0), 1);

                % Uncertainty with local window norm (SPEC.md §5.3)
                CW = W(1)^2 + 2 * sum(W(2:end).^2);
                cohSafe = max(Coh(kk), epsReg);
                GStd(kk) = sqrt((CW / Neff) * abs(G(kk))^2 * (1 - cohSafe) / cohSafe);
            end

            % Noise uncertainty
            CW = W(1)^2 + 2 * sum(W(2:end).^2);
            PhiVStd(kk) = sqrt(2 * CW / Neff) * abs(PhiV(kk));
        end

    else
        % MIMO: G(w_k) = Phi_yu(w_k) * Phi_u(w_k)^{-1} (SPEC.md §5.2)
        G = zeros(nf, ny, nu);          % (nf x ny x nu) complex
        PhiV = zeros(nf, ny, ny);       % (nf x ny x ny) noise spectrum
        GStd = zeros(nf, ny, nu);
        PhiVStd = zeros(nf, ny, ny);
        Coh = [];
        eps_floor = 1e-10;

        for kk = 1:nf
            Mk_k = Mk(kk);
            W = sidHannWin(Mk_k);

            Ryy_k = truncateCov(Ryy, Mk_k, ny, ny);
            Ruu_k = truncateCov(Ruu, Mk_k, nu, nu);
            Ryu_k = truncateCov(Ryu, Mk_k, ny, nu);
            Ruy_k = truncateCov(Ruy, Mk_k, nu, ny);

            PhiY_k = singleFreqDFT(Ryy_k, W, freqs(kk), ny, ny);
            PhiU_k = singleFreqDFT(Ruu_k, W, freqs(kk), nu, nu);
            PhiYU_k = singleFreqDFT(Ryu_k, W, freqs(kk), ny, nu, Ruy_k);

            G(kk, :, :) = PhiYU_k / PhiU_k;
            % Phi_v = Phi_y - Phi_yu * Phi_u^{-1} * Phi_uy
            PhiV_k = PhiY_k - PhiYU_k / PhiU_k * PhiYU_k';
            PhiV(kk, :, :) = real(PhiV_k);

            % Noise uncertainty (SPEC.md §5.3)
            CW = W(1)^2 + 2 * sum(W(2:end).^2);
            PhiVStd(kk, :, :) = sqrt(2 * CW / Neff) * abs(PhiV_k);

            % Diagonal MIMO G uncertainty: Var{G_{ij}} ≈ C_W/Neff * Phi_v_{ii} / Phi_u_{jj}
            for ii = 1:ny
                for jj = 1:nu
                    phiU_jj = real(PhiU_k(jj, jj));
                    phiV_ii = max(real(PhiV_k(ii, ii)), 0);
                    if phiU_jj > eps_floor
                        GStd(kk, ii, jj) = sqrt(CW / Neff * phiV_ii / phiU_jj);
                    else
                        GStd(kk, ii, jj) = Inf;
                    end
                end
            end
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
