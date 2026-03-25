function varargout = sidSpectrumPlot(result, varargin)
%SIDSPECTRUMPLOT Power spectrum plot with confidence bands.
%
%   sidSpectrumPlot(result)
%   sidSpectrumPlot(result, 'Confidence', 3)
%   h = sidSpectrumPlot(...)
%
%   Plots the noise spectrum (or output spectrum for time series) in dB.
%
%   INPUTS:
%     result - Struct returned by sidFreqBT, sidFreqBTFDR, or sidFreqETFE.
%
%   NAME-VALUE OPTIONS:
%     'Confidence'    - Number of standard deviations. Default: 3.
%     'FrequencyUnit' - 'rad/s' (default) or 'Hz'.
%     'Color'         - Line color. Default: [0.850 0.325 0.098].
%     'LineWidth'     - Line width. Default: 1.5.
%     'Axes'          - Axes handle. Creates new figure if empty.
%
%   OUTPUT:
%     h - Struct with fields .fig, .ax, .line.
%
%   See also: sidFreqBT, sidBodePlot
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

    % ---- Parse options ----
    p = inputParser;
    p.addParameter('Confidence', 3);
    p.addParameter('FrequencyUnit', 'rad/s');
    p.addParameter('Color', [0.850 0.325 0.098]);
    p.addParameter('LineWidth', 1.5);
    p.addParameter('Axes', []);
    p.parse(varargin{:});
    opts = p.Results;

    % ---- Frequency axis ----
    if strcmpi(opts.FrequencyUnit, 'Hz')
        freq = result.FrequencyHz;
        freqLabel = 'Frequency (Hz)';
    else
        freq = result.Frequency / result.SampleTime;
        freqLabel = 'Frequency (rad/s)';
    end

    % ---- Extract spectrum (first output channel) ----
    PhiV = result.NoiseSpectrum(:, 1);
    PhiVStd = result.NoiseSpectrumStd(:, 1);

    specDB = 10 * log10(max(PhiV, 1e-20));

    % ---- Create axes ----
    if isempty(opts.Axes)
        fig = figure;
        ax = axes;
    else
        ax = opts.Axes;
        fig = get(ax, 'Parent');
    end

    % ---- Plot ----
    axes(ax);  %#ok<LAXES>
    lineH = semilogx(freq, specDB, 'Color', opts.Color, 'LineWidth', opts.LineWidth);
    hold on;

    if opts.Confidence > 0 && ~isempty(PhiVStd)
        upper = 10 * log10(max(PhiV + opts.Confidence * PhiVStd, 1e-20));
        lower = 10 * log10(max(PhiV - opts.Confidence * PhiVStd, 1e-20));
        fillX = [freq; flipud(freq)];
        fillY = [upper; flipud(lower)];
        fill(fillX, fillY, opts.Color, ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
    end

    xlabel(freqLabel);
    if isempty(result.Response)
        ylabel('Output Spectrum (dB)');
        titleStr = 'Output Power Spectrum';
    else
        ylabel('Noise Spectrum (dB)');
        titleStr = 'Noise Spectrum';
    end
    title(sprintf('%s (%s)', titleStr, result.Method));
    grid on;
    set(ax, 'XScale', 'log');
    hold off;

    % ---- Output ----
    if nargout > 0
        h.fig = fig;
        h.ax = ax;
        h.line = lineH;
        varargout{1} = h;
    end
end
