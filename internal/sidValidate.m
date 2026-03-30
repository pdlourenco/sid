function [y, u, M, freqs, Ts, isTimeSeries, nTraj] = sidValidate(y, u, varargin)
%SIDVALIDATE Parse and validate inputs for sidFreq* functions.
%
%   [y, u, M, freqs, Ts, isTimeSeries, nTraj] = sidValidate(y, u, ...)
%
%   Supports both positional and name-value calling conventions:
%     sidValidate(y, u, M)
%     sidValidate(y, u, M, freqs)
%     sidValidate(y, u, 'WindowSize', M, 'Frequencies', freqs)
%
%   Example:
%     y = randn(500, 1); u = randn(500, 1);
%     [y, u, M, freqs, Ts, isTS] = sidValidate(y, u, 'WindowSize', 20);
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

    % ---- Ensure column orientation ----
    if isvector(y)
        y = y(:);
    end

    isTimeSeries = isempty(u);
    if ~isTimeSeries
        if isvector(u)
            u = u(:);
        end
    end

    % ---- Detect multi-trajectory (3D arrays) ----
    if ndims(y) == 3 %#ok<ISMAT>
        nTraj = size(y, 3);
    else
        nTraj = 1;
    end

    N = size(y, 1);

    % ---- Validate data ----
    if N < 2
        error('sid:tooShort', 'Data must have at least 2 samples.');
    end
    if ~isreal(y)
        error('sid:complexData', 'Complex data is not supported in v1.0. Input y must be real.');
    end
    if any(~isfinite(y(:)))
        error('sid:nonFinite', 'Data y contains NaN or Inf values.');
    end

    if ~isTimeSeries
        if size(u, 1) ~= N
            error('sid:sizeMismatch', ...
                'Input u (%d samples) and output y (%d samples) must have the same length.', ...
                size(u, 1), N);
        end
        if ~isreal(u)
            error('sid:complexData', 'Complex data is not supported in v1.0. Input u must be real.');
        end
        if any(~isfinite(u(:)))
            error('sid:nonFinite', 'Data u contains NaN or Inf values.');
        end
        % Multi-trajectory: u must have same number of trajectories
        if nTraj > 1
            if ndims(u) ~= 3 || size(u, 3) ~= nTraj %#ok<ISMAT>
                error('sid:trajMismatch', ...
                    'y has %d trajectories but u does not match.', nTraj);
            end
        end
    end

    % ---- Parse optional arguments ----
    % Detect whether positional or name-value style is used.
    M = [];
    freqs = [];
    Ts = 1.0;

    args = varargin;
    if ~isempty(args) && isnumeric(args{1})
        % Positional style: sidFreqBT(y, u, M) or sidFreqBT(y, u, M, freqs)
        M = args{1};
        if length(args) >= 2 && isnumeric(args{2})
            freqs = args{2};
            args = args(3:end);
        else
            args = args(2:end);
        end
    end

    % Name-value pairs (remaining args)
    k = 1;
    while k <= length(args)
        if ischar(args{k})
            switch lower(args{k})
                case 'windowsize'
                    M = args{k+1};
                    k = k + 2;
                case 'frequencies'
                    freqs = args{k+1};
                    k = k + 2;
                case 'sampletime'
                    Ts = args{k+1};
                    k = k + 2;
                otherwise
                    error('sid:unknownOption', 'Unknown option: %s', args{k});
            end
        else
            error('sid:badInput', 'Expected a string option name at position %d.', k);
        end
    end

    % ---- Defaults ----
    if isempty(M)
        M = min(floor(N / 10), 30);
    end

    if isempty(freqs)
        freqs = (1:128)' * pi / 128;
    else
        freqs = freqs(:);
    end

    % ---- Validate parameters ----
    if Ts <= 0
        error('sid:badTs', 'Sample time must be positive.');
    end

    if M < 2
        error('sid:badWindowSize', 'Window size M must be at least 2.');
    end

    if M > floor(N / 2)
        Morig = M;
        M = floor(N / 2);
        warning('sid:windowReduced', ...
            'Window size %d exceeds N/2 = %d. Reduced to %d.', Morig, floor(N/2), M);
    end

    if N < 10
        warning('sid:shortData', ...
            'Very short data (N = %d). Estimates will be unreliable.', N);
    end

    if any(freqs <= 0) || any(freqs > pi)
        error('sid:badFreqs', 'Frequencies must be in the range (0, pi] rad/sample.');
    end
end
