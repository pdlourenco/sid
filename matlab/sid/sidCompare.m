function result = sidCompare(model, y, u, varargin)
% SIDCOMPARE Compare model predicted output to measured data.
%
%   result = sidCompare(model, y, u)
%   result = sidCompare(model, y, u, 'InitialState', x0)
%   sidCompare(model, y, u, 'Plot', true)
%
%   Simulates the model's predicted output given the input signal and
%   compares it to the measured output using the NRMSE fit metric.
%
%   INPUTS:
%     model - Result struct from any sid estimator (see sidResultTypes).
%             Freq-domain (§1): requires .Response, .Frequency, .SampleTime
%             State-space (§4/§5): requires .A, .B, .StateDim, .InputDim,
%               .DataLength
%     y     - (N x ny) measured output, or (N+1 x p x L) state data for COSMIC
%     u     - (N x nu) input, or (N x q x L) for multi-trajectory COSMIC
%
%   NAME-VALUE OPTIONS:
%     'InitialState' - (p x 1) initial state for state-space simulation
%                      (default: first row of y)
%     'Plot'         - Display comparison plot (default: true if nargout==0)
%
%   OUTPUTS:
%     result.Predicted - (N x ny) model-predicted output
%     result.Measured  - (N x ny) measured output (copy)
%     result.Fit       - (1 x ny) NRMSE fit percentage per channel
%     result.Residual  - (N x ny) residual (y - y_pred)
%     result.Method    - char, method of the source model
%
%   EXAMPLES:
%     % Compare frequency-domain model to data
%     G = sidFreqBT(y, u);
%     result = sidCompare(G, y, u);
%
%     % Visual comparison with plot
%     sidCompare(G, y, u, 'Plot', true);
%
%   FIT METRIC:
%     fit = 100 * (1 - ||y - y_pred|| / ||y - mean(y)||)
%     100% = perfect, 0% = no better than mean, negative = worse than mean.
%
%   SPECIFICATION:
%     (Model output comparison — not yet in SPEC.md)
%
%   See also: sidResidual, sidFreqBT, sidLTVdisc
%
%   Changelog:
%   2026-03-29: First version by Pedro Lourenço.
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

    % ---- Parse options ----
    defs.InitialState = [];
    defs.Plot = (nargout == 0);
    opts = sidParseOptions(defs, varargin);
    x0 = opts.InitialState;
    doPlot = opts.Plot;

    % ---- Dispatch on model type ----
    if isfield(model, 'A') && isfield(model, 'B')
        [y_pred, y_meas] = simulateSS(model, y, u, x0);
    elseif isfield(model, 'Response')
        [y_pred, y_meas] = simulateFreq(model, y, u);
    else
        error('sid:badModel', ...
            'Model struct must have Response field (freq-domain) or A,B fields (state-space).');
    end

    % ---- Compute NRMSE fit ----
    % fit = 100 * (1 - ||y - y_pred|| / ||y - mean(y)||)
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
% SIMULATESS Simulate state-space model and return predicted vs measured.

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
% SIMULATEFREQ Simulate frequency-domain model output via IFFT.

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
    G_model = model.Response;

    % Ensure G_model is 3D (Octave may drop trailing singleton dims)
    if ndims(G_model) == 2 && (ny > 1 || nu > 1) %#ok<ISMAT>
        G_model = reshape(G_model, size(G_model, 1), ny, nu);
    end

    y_pred = sidFreqDomainSim(G_model, model.Frequency, u, N);
    y_meas = y;
end

function plotComparison(y_meas, y_pred, fitVec, methodName)
% PLOTCOMPARISON Overlay measured and predicted outputs.

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
            title(sprintf('Channel %d - Fit: %.1f%%', ch, fitVec(ch)));
        else
            title(sprintf('Model: %s - Fit: %.1f%%', methodName, fitVec(ch)));
        end
        xlabel('Sample');
        ylabel('Output');
        legend('Measured', 'Predicted', 'Location', 'best');
    end
end
