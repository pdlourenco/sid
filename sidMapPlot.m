function varargout = sidMapPlot(result, varargin)
%SIDMAPPLOT Time-frequency color map for sidFreqBTMap results.
%
%   sidMapPlot(result)
%   sidMapPlot(result, 'PlotType', 'coherence')
%   h = sidMapPlot(...)
%
%   Plots the time-varying frequency response as a color map, with time
%   on the x-axis and frequency on the y-axis.
%
%   INPUTS:
%     result - Struct returned by sidFreqBTMap.
%
%   NAME-VALUE OPTIONS:
%     'PlotType'      - What to plot:
%                        'magnitude' (default) - 20*log10(|G(w,t)|) in dB
%                        'phase'     - angle(G(w,t)) in degrees
%                        'noise'     - 10*log10(Phi_v(w,t)) in dB
%                        'coherence' - gamma^2(w,t) on [0, 1]
%                        'spectrum'  - 10*log10(Phi_y(w,t)) for time series
%     'FrequencyUnit' - 'rad/s' (default) or 'Hz'.
%     'CLim'          - Color axis limits [cmin cmax]. Default: [] (auto).
%     'Axes'          - Axes handle. Creates new figure if empty.
%
%   OUTPUT:
%     h - Struct with fields .fig, .ax, .surf.
%
%   See also: sidFreqBTMap, sidSpectrogramPlot, sidBodePlot
%
%   Example:
%     % Time-varying magnitude map
%     N = 4000; u = randn(N, 1);
%     y = filter([1], [1 -0.9], u) + 0.1*randn(N, 1);
%     result = sidFreqBTMap(y, u, 'SegmentLength', 512);
%     sidMapPlot(result, 'PlotType', 'magnitude');
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

    % ---- Parse options ----
    p = inputParser;
    p.addParameter('PlotType', 'magnitude');
    p.addParameter('FrequencyUnit', 'rad/s');
    p.addParameter('CLim', []);
    p.addParameter('Axes', []);
    p.parse(varargin{:});
    opts = p.Results;

    if ~isfield(result, 'Method') || ~strcmp(result.Method, 'sidFreqBTMap')
        error('sid:invalidResult', ...
            'Input must be a result struct from sidFreqBTMap.');
    end

    % ---- Frequency axis ----
    if strcmpi(opts.FrequencyUnit, 'Hz')
        freq = result.FrequencyHz;
        freqLabel = 'Frequency (Hz)';
    else
        freq = result.Frequency / result.SampleTime;
        freqLabel = 'Frequency (rad/s)';
    end

    T = result.Time(:)';  % row vector for pcolor
    F = freq(:);           % column vector

    % ---- Select data based on PlotType ----
    plotType = lower(opts.PlotType);
    isTimeSeriesResult = isempty(result.Response);

    switch plotType
        case 'magnitude'
            if isTimeSeriesResult
                error('sid:noResponse', ...
                    'PlotType ''magnitude'' requires input-output data (not time series).');
            end
            Z = 20 * log10(max(abs(result.Response(:, :, 1)), eps));
            colorLabel = 'Magnitude (dB)';
            titleStr = 'Time-Varying Magnitude';

        case 'phase'
            if isTimeSeriesResult
                error('sid:noResponse', ...
                    'PlotType ''phase'' requires input-output data (not time series).');
            end
            Z = angle(result.Response(:, :, 1)) * 180 / pi;
            colorLabel = 'Phase (deg)';
            titleStr = 'Time-Varying Phase';

        case 'noise'
            Z = 10 * log10(max(result.NoiseSpectrum(:, :, 1), eps));
            colorLabel = 'Noise PSD (dB)';
            titleStr = 'Time-Varying Noise Spectrum';

        case 'coherence'
            if isempty(result.Coherence)
                error('sid:noCoherence', ...
                    'Coherence is only available for SISO input-output data.');
            end
            Z = result.Coherence;
            colorLabel = 'Coherence';
            titleStr = 'Time-Varying Coherence';

        case 'spectrum'
            Z = 10 * log10(max(result.NoiseSpectrum(:, :, 1), eps));
            colorLabel = 'PSD (dB)';
            titleStr = 'Time-Varying Power Spectrum';

        otherwise
            error('sid:invalidPlotType', ...
                'Unknown PlotType ''%s''. Use magnitude, phase, noise, coherence, or spectrum.', ...
                opts.PlotType);
    end

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
    surfH = pcolor(ax, T, F, Z);
    shading(ax, 'flat');
    set(ax, 'YScale', 'log');
    cb = colorbar(ax);
    ylabel(cb, colorLabel);

    xlabel('Time (s)');
    ylabel(freqLabel);
    title(sprintf('%s (L=%d, M=%d)', titleStr, result.SegmentLength, result.WindowSize));

    if ~isempty(opts.CLim)
        caxis(ax, opts.CLim);
    end

    set(ax, 'Layer', 'top');

    % ---- Output ----
    if nargout > 0
        h.fig = fig;
        h.ax = ax;
        h.surf = surfH;
        varargout{1} = h;
    end
end
