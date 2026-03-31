function result = sidResidual(model, y, u, varargin)
%SIDRESIDUAL Compute model residuals and perform diagnostic tests.
%
%   result = sidResidual(model, y, u)
%   result = sidResidual(model, y, u, 'MaxLag', M)
%   sidResidual(model, y, u, 'Plot', true)
%
%   Computes residuals from an estimated model and performs whiteness and
%   independence tests to assess model quality.
%
%   INPUTS:
%     model - Result struct from any sid estimator (sidFreqBT, sidLTVdisc, etc.)
%     y     - (N x ny) measured output (or (N+1 x p x L) state data for COSMIC)
%     u     - (N x nu) input, or [] for time-series models
%
%   NAME-VALUE OPTIONS:
%     'MaxLag' - Maximum lag for correlation tests (default: min(25, floor(N/5)))
%     'Plot'   - Display diagnostic plot (default: true if nargout==0)
%
%   OUTPUT:
%     result.Residual         - (N x ny) residual time series e(t)
%     result.AutoCorr         - (M+1 x 1) normalised autocorrelation r_ee(tau)
%     result.CrossCorr        - (2M+1 x 1) normalised cross-corr r_eu(tau), or []
%     result.ConfidenceBound  - scalar, 99% bound 2.58/sqrt(N)
%     result.WhitenessPass    - logical, true if autocorrelation test passes
%     result.IndependencePass - logical, true if cross-correlation test passes
%     result.DataLength       - N
%
%   See also: sidCompare, sidFreqBT, sidLTVdisc

    % ---- Parse options ----
    maxLag = [];
    doPlot = (nargout == 0);

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'maxlag'
                    maxLag = varargin{k+1};
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

    isTimeSeries = isempty(u);

    % ---- Dispatch on model type ----
    if isfield(model, 'A') && isfield(model, 'B')
        % State-space model (sidLTVdisc)
        [e, N_eff] = computeResidualSS(model, y, u);
    elseif isfield(model, 'Response')
        % Frequency-domain model (sidFreqBT, sidFreqETFE, etc.)
        % Response may be empty for time-series models — that's valid
        [e, N_eff] = computeResidualFreq(model, y, u);
    else
        error('sid:badModel', ...
            'Model struct must have Response field (freq-domain) or A,B fields (state-space).');
    end

    % ---- Default MaxLag ----
    if isempty(maxLag)
        maxLag = min(25, floor(N_eff / 5));
    end

    % ---- Whiteness test: normalised autocorrelation ----
    ny = size(e, 2);
    % Use first channel for SISO tests; average across channels for MIMO
    if ny == 1
        e_test = e;
    else
        % Average autocorrelation across channels
        e_test = e(:, 1);  % use first channel for scalar diagnostics
    end

    Ree = sidCov(e_test, e_test, maxLag);  % (M+1 x 1)
    Ree0 = Ree(1);
    if Ree0 > 0
        autoCorr = Ree / Ree0;
    else
        autoCorr = zeros(maxLag + 1, 1);
    end

    confBound = 2.58 / sqrt(N_eff);
    whitenessPass = all(abs(autoCorr(2:end)) < confBound);

    % ---- Independence test: normalised cross-correlation ----
    if ~isTimeSeries
        % Extract first channel, first trajectory as column vector
        if ndims(u) == 3 %#ok<ISMAT>
            u_test = u(:, 1, 1);
        elseif size(u, 2) > 1
            u_test = u(:, 1);
        else
            u_test = u(:);
        end
        % Ensure same length as residual
        Nu = min(size(e_test, 1), size(u_test, 1));
        e_test_trim = e_test(1:Nu, 1);
        u_test_trim = u_test(1:Nu, 1);

        Reu_pos = sidCov(e_test_trim, u_test_trim, maxLag);  % positive lags
        Rue_pos = sidCov(u_test_trim, e_test_trim, maxLag);  % for negative lags

        Ruu0 = u_test_trim' * u_test_trim / Nu;
        denom = sqrt(Ree0 * Ruu0);
        if denom > 0
            crossCorr_pos = Reu_pos / denom;
            crossCorr_neg = Rue_pos / denom;
        else
            crossCorr_pos = zeros(maxLag + 1, 1);
            crossCorr_neg = zeros(maxLag + 1, 1);
        end

        % Assemble: tau = -M, ..., -1, 0, 1, ..., M
        crossCorr = [flipud(crossCorr_neg(2:end)); crossCorr_pos];
        independencePass = all(abs(crossCorr) < confBound);
    else
        crossCorr = [];
        independencePass = true;  % vacuously true for time-series
    end

    % ---- Pack result ----
    result.Residual         = e;
    result.AutoCorr         = autoCorr;
    result.CrossCorr        = crossCorr;
    result.ConfidenceBound  = confBound;
    result.WhitenessPass    = whitenessPass;
    result.IndependencePass = independencePass;
    result.DataLength       = N_eff;

    % ---- Plot ----
    if doPlot
        plotResidualDiagnostics(autoCorr, crossCorr, confBound, maxLag, isTimeSeries);
    end
end


% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [e, N] = computeResidualSS(model, X, U)
%COMPUTERESIDUALSS Residuals from state-space model (COSMIC).
%   e(k) = x(k+1) - A(k)*x(k) - B(k)*u(k)

    Nm = model.DataLength;  % number of time steps
    p = model.StateDim;

    % Handle multi-trajectory
    if ndims(X) == 3 %#ok<ISMAT>
        L = size(X, 3);
    else
        L = 1;
    end

    e_all = zeros(Nm, p);
    for l = 1:L
        if L > 1
            Xl = X(:, :, l);   % (N+1 x p)
            Ul = U(:, :, l);   % (N x q)
        else
            Xl = X;
            Ul = U;
        end

        for k = 1:Nm
            x_pred = model.A(:, :, k) * Xl(k, :)' + model.B(:, :, k) * Ul(k, :)';
            e_all(k, :) = e_all(k, :) + (Xl(k+1, :) - x_pred');
        end
    end
    e = e_all / L;
    N = Nm;
end


function [e, N] = computeResidualFreq(model, y, u)
%COMPUTERESIDUALFREQ Residuals from frequency-domain model via IFFT.
%   Filter input through G(w) in frequency domain, IFFT, subtract from y.

    if isvector(y)
        y = y(:);
    end
    if ~isempty(u) && isvector(u)
        u = u(:);
    end

    N = size(y, 1);
    ny = size(y, 2);

    if isempty(u)
        % Time-series: residual = y - estimated noise spectrum contribution
        % For time-series, there's no G to filter through. Return y as residual
        % (residual analysis tests whether y itself is white).
        e = y;
        return;
    end

    nu = size(u, 2);

    % Interpolate G onto full FFT grid
    freqs_model = model.Frequency;  % (nf x 1) rad/sample
    G_model = model.Response;       % (nf x ny x nu) complex

    nfft = N;
    freqs_fft = (1:floor(nfft/2))' * (2 * pi / nfft);  % FFT freq grid, rad/sample

    % FFT of input
    U_fft = fft(u, nfft, 1);  % (nfft x nu)

    % Predicted output in frequency domain
    Y_pred_fft = zeros(nfft, ny);

    if ny == 1 && nu == 1
        % SISO: interpolate G onto FFT grid
        G_interp = interpG(freqs_model, G_model(:), freqs_fft);
        % Apply to positive frequencies
        npos = length(freqs_fft);
        Y_pred_fft(2:npos+1, 1) = G_interp .* U_fft(2:npos+1, 1);
        % Mirror for negative frequencies (conjugate symmetry)
        if mod(nfft, 2) == 0
            Y_pred_fft(npos+2:end, 1) = conj(Y_pred_fft(npos:-1:2, 1));
        else
            Y_pred_fft(npos+2:end, 1) = conj(Y_pred_fft(npos+1:-1:2, 1));
        end
    else
        % MIMO: interpolate each channel
        for iy = 1:ny
            for iu = 1:nu
                Gij = reshape(G_model(:, iy, iu), [], 1);
                G_interp = interpG(freqs_model, Gij, freqs_fft);
                npos = length(freqs_fft);
                Y_pred_fft(2:npos+1, iy) = Y_pred_fft(2:npos+1, iy) + ...
                    G_interp .* U_fft(2:npos+1, iu);
            end
        end
        % Mirror
        npos = length(freqs_fft);
        for iy = 1:ny
            if mod(nfft, 2) == 0
                Y_pred_fft(npos+2:end, iy) = conj(Y_pred_fft(npos:-1:2, iy));
            else
                Y_pred_fft(npos+2:end, iy) = conj(Y_pred_fft(npos+1:-1:2, iy));
            end
        end
    end

    y_pred = real(ifft(Y_pred_fft, nfft, 1));
    e = y - y_pred;
end


function G_interp = interpG(freqs_model, G, freqs_target)
%INTERPG Interpolate complex transfer function onto target frequency grid.
%   Uses linear interpolation of real and imaginary parts.
%   Always returns a column vector.

    G_interp = interp1(freqs_model(:), real(G(:)), freqs_target(:), 'linear', 'extrap') + ...
        1i * interp1(freqs_model(:), imag(G(:)), freqs_target(:), 'linear', 'extrap');
    G_interp = G_interp(:);  % ensure column
end


function plotResidualDiagnostics(autoCorr, crossCorr, confBound, maxLag, isTimeSeries)
%PLOTRESIDUALDIAGNOSTICS Two-panel diagnostic plot.

    if isTimeSeries
        nPanels = 1;
    else
        nPanels = 2;
    end

    figure;

    % Top panel: autocorrelation
    subplot(nPanels, 1, 1);
    lags_auto = (0:maxLag)';
    bar(lags_auto, autoCorr, 0.5);
    hold on;
    plot([0 maxLag], [confBound confBound], 'r--', 'LineWidth', 1);
    plot([0 maxLag], [-confBound -confBound], 'r--', 'LineWidth', 1);
    hold off;
    xlabel('Lag');
    ylabel('r_{ee}(\tau)');
    title('Residual Autocorrelation (Whiteness Test)');

    % Bottom panel: cross-correlation
    if ~isTimeSeries && ~isempty(crossCorr)
        subplot(nPanels, 1, 2);
        lags_cross = (-maxLag:maxLag)';
        bar(lags_cross, crossCorr, 0.5);
        hold on;
        plot([-maxLag maxLag], [confBound confBound], 'r--', 'LineWidth', 1);
        plot([-maxLag maxLag], [-confBound -confBound], 'r--', 'LineWidth', 1);
        hold off;
        xlabel('Lag');
        ylabel('r_{eu}(\tau)');
        title('Residual-Input Cross-Correlation (Independence Test)');
    end
end
