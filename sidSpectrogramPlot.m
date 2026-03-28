function varargout = sidSpectrogramPlot(result, varargin)
%SIDSPECTROGRAMPLOT Spectrogram color map plot.
%
%   sidSpectrogramPlot(result)
%   sidSpectrogramPlot(result, 'FrequencyScale', 'log')
%   h = sidSpectrogramPlot(...)
%
%   Plots the spectrogram as a time-frequency color map with power in dB.
%
%   INPUTS:
%     result - Struct returned by sidSpectrogram.
%
%   NAME-VALUE OPTIONS:
%     'FrequencyScale' - 'linear' (default) or 'log'.
%     'Channel'        - Channel index to plot for multi-channel data.
%                        Default: 1.
%     'CLim'           - Color axis limits [cmin cmax] in dB.
%                        Default: [] (auto).
%     'Axes'           - Axes handle. Creates new figure if empty.
%
%   OUTPUT:
%     h - Struct with fields .fig, .ax, .surf.
%
%   See also: sidSpectrogram, sidMapPlot
%
%   Example:
%     % Plot spectrogram of a chirp signal
%     Fs = 1000; Ts = 1/Fs; N = 5000;
%     t = (0:N-1)' * Ts;
%     x = cos(2*pi * (50 + 100*t/max(t)) .* t);
%     result = sidSpectrogram(x, 'WindowLength', 256, 'SampleTime', Ts);
%     sidSpectrogramPlot(result);
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
    p.addParameter('FrequencyScale', 'linear');
    p.addParameter('Channel', 1);
    p.addParameter('CLim', []);
    p.addParameter('Axes', []);
    p.parse(varargin{:});
    opts = p.Results;

    if ~isfield(result, 'Method') || ~strcmp(result.Method, 'sidSpectrogram')
        error('sid:invalidResult', ...
            'Input must be a result struct from sidSpectrogram.');
    end

    ch = opts.Channel;
    nCh = size(result.Power, 3);
    if ch < 1 || ch > nCh
        error('sid:invalidChannel', ...
            'Channel %d out of range (data has %d channels).', ch, nCh);
    end

    % ---- Extract data ----
    Z = result.PowerDB(:, :, ch);  % (n_bins x K)
    T = result.Time(:)';           % row for meshgrid
    F = result.Frequency(:);       % column

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
    cb = colorbar(ax);
    ylabel(cb, 'Power (dB)');

    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title(sprintf('Spectrogram (L=%d, P=%d, NFFT=%d)', ...
        result.WindowLength, result.Overlap, result.NFFT));

    if strcmpi(opts.FrequencyScale, 'log')
        set(ax, 'YScale', 'log');
    end

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
