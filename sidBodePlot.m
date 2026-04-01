function varargout = sidBodePlot(result, varargin)
% SIDBODEPLOT Bode diagram with confidence bands.
%
%   sidBodePlot(result)
%   sidBodePlot(result, 'Confidence', 3)
%   h = sidBodePlot(...)
%
%   Plots the magnitude (dB) and phase (degrees) of the estimated
%   frequency response from a sidFreq* result struct.
%
%   INPUTS:
%     result     - Struct returned by sidFreqBT, sidFreqBTFDR, or sidFreqETFE.
%
%   NAME-VALUE OPTIONS:
%     'Confidence'    - Number of standard deviations for the shaded
%                       confidence band. Default: 3. Set to 0 to hide.
%     'FrequencyUnit' - 'rad/s' (default) or 'Hz'.
%     'Color'         - Line color. Default: [0 0.447 0.741] (MATLAB blue).
%     'LineWidth'     - Line width. Default: 1.5.
%     'Axes'          - [ax_mag, ax_phase] axes handles. Creates new
%                       figure if empty.
%
%   OUTPUTS:
%     h - Struct with fields .fig, .axMag, .axPhase, .lineMag, .linePhase.
%
%   EXAMPLES:
%     N = 1000; u = randn(N, 1);
%     y = filter([1], [1 -0.9], u) + 0.1*randn(N, 1);
%     result = sidFreqBT(y, u);
%     sidBodePlot(result, 'Confidence', 3);
%
%   SPECIFICATION:
%     SPEC.md §11.1 — sidBodePlot
%
%   See also: sidFreqBT, sidSpectrumPlot
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

    % ---- Parse options ----
    p = inputParser;
    p.addParameter('Confidence', 3);
    p.addParameter('FrequencyUnit', 'rad/s');
    p.addParameter('Color', [0 0.447 0.741]);
    p.addParameter('LineWidth', 1.5);
    p.addParameter('Axes', []);
    p.parse(varargin{:});
    opts = p.Results;

    if isempty(result.Response)
        error('sid:noResponse', ...
            'Result contains no frequency response (time series mode). Use sidSpectrumPlot.');
    end

    % ---- Frequency axis ----
    if strcmpi(opts.FrequencyUnit, 'Hz')
        freq = result.FrequencyHz;
        freqLabel = 'Frequency (Hz)';
    else
        freq = result.Frequency / result.SampleTime;
        freqLabel = 'Frequency (rad/s)';
    end

    % ---- Extract SISO response (first channel pair for MIMO) ----
    G = result.Response(:, 1);
    if ~isempty(result.ResponseStd)
        GStd = result.ResponseStd(:, 1);
    else
        GStd = [];
    end

    mag = abs(G);
    magdB = 20 * log10(mag);
    phase = angle(G) * 180 / pi;

    % ---- Create axes ----
    if isempty(opts.Axes)
        fig = figure;
        axMag = subplot(2, 1, 1);
        axPhase = subplot(2, 1, 2);
    else
        axMag = opts.Axes(1);
        axPhase = opts.Axes(2);
        fig = get(axMag, 'Parent');
    end

    % ---- Magnitude plot ----
    axes(axMag);  %#ok<LAXES>
    lineMag = semilogx(freq, magdB, 'Color', opts.Color, 'LineWidth', opts.LineWidth);
    hold on;

    if opts.Confidence > 0 && ~isempty(GStd)
        magUpper = 20 * log10(mag + opts.Confidence * GStd);
        magLower = 20 * log10(max(mag - opts.Confidence * GStd, 1e-20));
        fillX = [freq; flipud(freq)];
        fillY = [magUpper; flipud(magLower)];
        fill(fillX, fillY, opts.Color, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    ylabel('Magnitude (dB)');
    title(sprintf('Bode Diagram (%s)', result.Method));
    grid on;
    set(axMag, 'XScale', 'log');
    hold off;

    % ---- Phase plot ----
    axes(axPhase);  %#ok<LAXES>
    linePhase = semilogx(freq, phase, 'Color', opts.Color, 'LineWidth', opts.LineWidth);
    hold on;

    if opts.Confidence > 0 && ~isempty(GStd)
        phaseStd = opts.Confidence * GStd ./ max(mag, 1e-20) * 180 / pi;
        fillY = [phase + phaseStd; flipud(phase - phaseStd)];
        fill(fillX, fillY, opts.Color, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    xlabel(freqLabel);
    ylabel('Phase (deg)');
    grid on;
    set(axPhase, 'XScale', 'log');
    hold off;

    % ---- Output ----
    if nargout > 0
        h.fig = fig;
        h.axMag = axMag;
        h.axPhase = axPhase;
        h.lineMag = lineMag;
        h.linePhase = linePhase;
        varargout{1} = h;
    end
end
