function result = sidResidual(model, y, u, varargin)
% SIDRESIDUAL Compute model residuals and perform diagnostic tests.
%
%   result = sidResidual(model, y, u)
%   result = sidResidual(model, y, u, 'MaxLag', M)
%   sidResidual(model, y, u, 'Plot', true)
%
%   Computes residuals from an estimated model and performs whiteness and
%   independence tests to assess model quality.
%
%   INPUTS:
%     model - Result struct from any sid estimator (see sidResultTypes).
%             Freq-domain (§1): requires .Response, .Frequency, .SampleTime
%             State-space (§4/§5): requires .A, .B, .StateDim, .InputDim,
%               .DataLength
%     y     - (N x ny) measured output (or (N+1 x p x L) state data for COSMIC)
%     u     - (N x nu) input, or [] for time-series models
%
%   NAME-VALUE OPTIONS:
%     'MaxLag' - Maximum lag for correlation tests (default: min(25, floor(N/5)))
%     'Plot'   - Display diagnostic plot (default: true if nargout==0)
%
%   OUTPUTS:
%     result.Residual         - (N x ny) residual time series e(t)
%     result.AutoCorr         - (M+1 x 1) normalised autocorrelation r_ee(tau)
%     result.CrossCorr        - (2M+1 x 1) normalised cross-corr r_eu(tau), or []
%     result.ConfidenceBound  - scalar, 99% bound 2.58/sqrt(N)
%     result.WhitenessPass    - logical, true if autocorrelation test passes
%     result.IndependencePass - logical, true if cross-correlation test passes
%     result.DataLength       - N
%
%   EXAMPLES:
%     % Residual analysis for a frequency-domain model
%     G = sidFreqBT(y, u);
%     result = sidResidual(G, y, u);
%
%     % Visual diagnostic plot
%     sidResidual(G, y, u, 'Plot', true);
%
%   SPECIFICATION:
%     (Model residual analysis — not yet in SPEC.md)
%
%   See also: sidCompare, sidFreqBT, sidLTVdisc
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
    defs.MaxLag = [];
    defs.Plot = (nargout == 0);
    opts = sidParseOptions(defs, varargin);
    maxLag = opts.MaxLag;
    doPlot = opts.Plot;

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

    % ---- Per-channel whiteness and independence tests ----
    % Run tests on EVERY output channel, not just the first.
    % The overall pass/fail requires ALL channels to pass.
    ny = size(e, 2);
    confBound = 2.58 / sqrt(N_eff);

    % Determine number of input channels for independence tests
    if ~isTimeSeries
        if ndims(u) == 3 %#ok<ISMAT>
            nu = size(u, 2);
            u_2d = u(:, :, 1);  % use first trajectory
        else
            nu = size(u, 2);
            u_2d = u;
        end
        Nu = min(size(e, 1), size(u_2d, 1));
    else
        nu = 0;
        u_2d = [];
        Nu = size(e, 1);
    end

    % Storage for per-channel results
    autoCorr_all   = zeros(maxLag + 1, ny);
    whitenessPass_all = true(ny, 1);

    if ~isTimeSeries
        % Cross-correlation: test each (output, input) pair
        nPairs = ny * nu;
        crossCorr_all   = zeros(2 * maxLag + 1, nPairs);
        indepPass_all   = true(nPairs, 1);
    end

    for ch = 1:ny
        e_ch = e(1:Nu, ch);

        % ---- Whiteness: normalised autocorrelation ----
        Ree = sidCov(e_ch, e_ch, maxLag);
        Ree0 = Ree(1);
        if Ree0 > 0
            autoCorr_all(:, ch) = Ree / Ree0;
        end
        whitenessPass_all(ch) = all(abs(autoCorr_all(2:end, ch)) < confBound);

        % ---- Independence: normalised cross-correlation per input ----
        if ~isTimeSeries
            for iu = 1:nu
                pairIdx = (ch - 1) * nu + iu;
                u_ch = u_2d(1:Nu, iu);

                Reu_pos = sidCov(e_ch, u_ch, maxLag);
                Rue_pos = sidCov(u_ch, e_ch, maxLag);

                Ruu0 = u_ch' * u_ch / Nu;
                denom = sqrt(Ree0 * Ruu0);
                if denom > 0
                    cc_pos = Reu_pos / denom;
                    cc_neg = Rue_pos / denom;
                else
                    cc_pos = zeros(maxLag + 1, 1);
                    cc_neg = zeros(maxLag + 1, 1);
                end

                crossCorr_all(:, pairIdx) = [flipud(cc_neg(2:end)); cc_pos];
                indepPass_all(pairIdx) = all(abs(crossCorr_all(:, pairIdx)) < confBound);
            end
        end
    end

    % Aggregate: pass only if ALL channels/pairs pass
    whitenessPass    = all(whitenessPass_all);
    if ~isTimeSeries
        independencePass = all(indepPass_all);
        crossCorr = crossCorr_all;
    else
        independencePass = true;
        crossCorr = [];
    end

    % For backward compatibility and plotting, provide first-channel summary
    autoCorr = autoCorr_all(:, 1);

    % ---- Pack result ----
    result.Residual           = e;
    result.AutoCorr           = autoCorr;
    result.AutoCorrAll        = autoCorr_all;
    result.CrossCorr          = crossCorr;
    result.ConfidenceBound    = confBound;
    result.WhitenessPass      = whitenessPass;
    result.WhitenessPassAll   = whitenessPass_all;
    result.IndependencePass   = independencePass;
    if ~isTimeSeries
        result.IndependencePassAll = indepPass_all;
    end
    result.DataLength         = N_eff;

    % ---- Plot ----
    if doPlot
        % Plot first channel for backward-compatible display
        if ~isTimeSeries
            crossCorr_plot = crossCorr(:, 1);
        else
            crossCorr_plot = [];
        end
        plotResidualDiagnostics(autoCorr, crossCorr_plot, confBound, maxLag, isTimeSeries);
    end
end

% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [e, N] = computeResidualSS(model, X, U)
% COMPUTERESIDUALSS Residuals from state-space model (COSMIC).
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
% COMPUTERESIDUALFREQ Residuals from frequency-domain model via IFFT.
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
    G_model = model.Response;

    % Ensure G_model is 3D (Octave may drop trailing singleton dims)
    if ndims(G_model) == 2 %#ok<ISMAT>
        G_model = reshape(G_model, size(G_model, 1), ny, nu);
    end

    y_pred = sidFreqDomainSim(G_model, model.Frequency, u, N);
    e = y - y_pred;
end

function plotResidualDiagnostics(autoCorr, crossCorr, confBound, maxLag, isTimeSeries)
% PLOTRESIDUALDIAGNOSTICS Two-panel diagnostic plot.

    if isTimeSeries
        nPanels = 1;
    else
        nPanels = 2;
    end

    figure;

    % Top panel: autocorrelation
    subplot(nPanels, 1, 1);
    lags_auto = (0:maxLag)';
    bar(lags_auto, autoCorr, 0.5, 'FaceColor', [0.3 0.5 0.8]);
    hold on;
    plot([0 maxLag], [confBound confBound], 'r--', 'LineWidth', 1);
    plot([0 maxLag], [-confBound -confBound], 'r--', 'LineWidth', 1);
    % Highlight violations
    violations = abs(autoCorr(2:end)) >= confBound;
    if any(violations)
        vIdx = find(violations);
        bar(vIdx, autoCorr(vIdx + 1), 0.5, 'FaceColor', [0.9 0.2 0.2]);
    end
    xlabel('Lag');
    ylabel('r_{ee}(\tau)');
    title('Residual Autocorrelation (Whiteness Test)');
    hold off;

    % Bottom panel: cross-correlation
    if ~isTimeSeries && ~isempty(crossCorr)
        subplot(nPanels, 1, 2);
        lags_cross = (-maxLag:maxLag)';
        bar(lags_cross, crossCorr, 0.5, 'FaceColor', [0.3 0.5 0.8]);
        hold on;
        plot([-maxLag maxLag], [confBound confBound], 'r--', 'LineWidth', 1);
        plot([-maxLag maxLag], [-confBound -confBound], 'r--', 'LineWidth', 1);
        violations_c = abs(crossCorr) >= confBound;
        if any(violations_c)
            vIdx_c = find(violations_c);
            bar(lags_cross(vIdx_c), crossCorr(vIdx_c), 0.5, 'FaceColor', [0.9 0.2 0.2]);
        end
        xlabel('Lag');
        ylabel('r_{eu}(\tau)');
        title('Residual-Input Cross-Correlation (Independence Test)');
        hold off;
    end
end
