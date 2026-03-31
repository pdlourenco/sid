function result = sidCompare(model, y, u, varargin)
%SIDCOMPARE Compare model predicted output to measured data.
%
%   result = sidCompare(model, y, u)
%   result = sidCompare(model, y, u, 'InitialState', x0)
%   sidCompare(model, y, u, 'Plot', true)
%
%   Simulates the model's predicted output given the input signal and
%   compares it to the measured output using the NRMSE fit metric.
%
%   INPUTS:
%     model - Result struct from any sid estimator
%     y     - (N x ny) measured output, or (N+1 x p x L) state data for COSMIC
%     u     - (N x nu) input, or (N x q x L) for multi-trajectory COSMIC
%
%   NAME-VALUE OPTIONS:
%     'InitialState' - (p x 1) initial state for state-space simulation
%                      (default: first row of y)
%     'Plot'         - Display comparison plot (default: true if nargout==0)
%
%   OUTPUT:
%     result.Predicted - (N x ny) model-predicted output
%     result.Measured  - (N x ny) measured output (copy)
%     result.Fit       - (1 x ny) NRMSE fit percentage per channel
%     result.Residual  - (N x ny) residual (y - y_pred)
%     result.Method    - char, method of the source model
%
%   FIT METRIC:
%     fit = 100 * (1 - ||y - y_pred|| / ||y - mean(y)||)
%     100% = perfect, 0% = no better than mean, negative = worse than mean.
%
%   See also: sidResidual, sidFreqBT, sidLTVdisc

    % ---- Parse options ----
    x0 = [];
    doPlot = (nargout == 0);

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'initialstate'
                    x0 = varargin{k+1};
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

    % ---- Dispatch on model type ----
    if isfield(model, 'A') && isfield(model, 'B')
        [y_pred, y_meas] = simulateSS(model, y, u, x0);
    elseif isfield(model, 'Response') && ~isempty(model.Response)
        [y_pred, y_meas] = simulateFreq(model, y, u);
    else
        error('sid:badModel', ...
            'Model struct must have Response field (freq-domain) or A,B fields (state-space).');
    end

    % ---- Compute NRMSE fit ----
    ny = size(y_meas, 2);
    fitVec = zeros(1, ny);
    for ch = 1:ny
        ym = y_meas(:, ch);
        yp = y_pred(:, ch);
        denom = norm(ym - mean(ym));
        if denom > 0
            fitVec(ch) = 100 * (1 - norm(ym - yp) / denom);
        else
            fitVec(ch) = NaN;
        end
    end

    % ---- Pack result ----
    result.Predicted = y_pred;
    result.Measured  = y_meas;
    result.Fit       = fitVec;
    result.Residual  = y_meas - y_pred;
    if isfield(model, 'Method')
        result.Method = model.Method;
    else
        result.Method = 'unknown';
    end

    % ---- Plot ----
    if doPlot
        plotComparison(y_meas, y_pred, fitVec, result.Method);
    end
end


% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [y_pred, y_meas] = simulateSS(model, X, U, x0)
%SIMULATESS Simulate state-space model and return predicted vs measured.

    Nm = model.DataLength;
    p = model.StateDim;

    % Handle multi-trajectory: average fit across trajectories
    if ndims(X) == 3 %#ok<ISMAT>
        L = size(X, 3);
    else
        L = 1;
    end

    y_pred_sum = zeros(Nm, p);
    y_meas_sum = zeros(Nm, p);

    for l = 1:L
        if L > 1
            Xl = X(:, :, l);
            Ul = U(:, :, l);
        else
            Xl = X;
            Ul = U;
        end

        % Initial state
        if ~isempty(x0)
            xk = x0(:);
        else
            xk = Xl(1, :)';
        end

        x_hat = zeros(Nm, p);
        for k = 1:Nm
            x_next = model.A(:, :, k) * xk + model.B(:, :, k) * Ul(k, :)';
            x_hat(k, :) = x_next';
            xk = x_next;
        end

        % Measured: x(2:N+1), Predicted: x_hat(1:N)
        y_meas_sum = y_meas_sum + Xl(2:end, :);
        y_pred_sum = y_pred_sum + x_hat;
    end

    y_meas = y_meas_sum / L;
    y_pred = y_pred_sum / L;
end


function [y_pred, y_meas] = simulateFreq(model, y, u)
%SIMULATEFREQ Simulate frequency-domain model output via IFFT.

    if isvector(y)
        y = y(:);
    end
    if ~isempty(u) && isvector(u)
        u = u(:);
    end

    N = size(y, 1);
    ny = size(y, 2);

    if isempty(u)
        % Time-series: no input to filter, predicted = 0
        y_pred = zeros(N, ny);
        y_meas = y;
        return;
    end

    nu = size(u, 2);

    % Interpolate G onto FFT grid and filter
    freqs_model = model.Frequency;
    G_model = model.Response;

    nfft = N;
    freqs_fft = (1:floor(nfft/2))' * (2 * pi / nfft);

    U_fft = fft(u, nfft, 1);
    Y_pred_fft = zeros(nfft, ny);

    npos = length(freqs_fft);
    for iy = 1:ny
        for iu = 1:nu
            if ny == 1 && nu == 1
                Gij = G_model(:);
            else
                Gij = reshape(G_model(:, iy, iu), [], 1);
            end
            G_interp = interp1(freqs_model, real(Gij), freqs_fft, 'linear', 'extrap') + ...
                1i * interp1(freqs_model, imag(Gij), freqs_fft, 'linear', 'extrap');
            Y_pred_fft(2:npos+1, iy) = Y_pred_fft(2:npos+1, iy) + ...
                G_interp .* U_fft(2:npos+1, iu);
        end
    end

    % Conjugate symmetry for real output
    for iy = 1:ny
        if mod(nfft, 2) == 0
            Y_pred_fft(npos+2:end, iy) = conj(Y_pred_fft(npos:-1:2, iy));
        else
            Y_pred_fft(npos+2:end, iy) = conj(Y_pred_fft(npos+1:-1:2, iy));
        end
    end

    y_pred = real(ifft(Y_pred_fft, nfft, 1));
    y_meas = y;
end


function plotComparison(y_meas, y_pred, fitVec, methodName)
%PLOTCOMPARISON Overlay measured and predicted outputs.

    ny = size(y_meas, 2);
    N = size(y_meas, 1);
    t = (1:N)';

    figure;
    for ch = 1:ny
        subplot(ny, 1, ch);
        plot(t, y_meas(:, ch), 'b-', 'LineWidth', 1); hold on;
        plot(t, y_pred(:, ch), 'r--', 'LineWidth', 1);
        hold off;
        if ny > 1
            title(sprintf('Channel %d — Fit: %.1f%%', ch, fitVec(ch)));
        else
            title(sprintf('Model: %s — Fit: %.1f%%', methodName, fitVec(ch)));
        end
        xlabel('Sample');
        ylabel('Output');
        legend('Measured', 'Predicted', 'Location', 'best');
    end
end
