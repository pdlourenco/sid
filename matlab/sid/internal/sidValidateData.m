function [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u)
% SIDVALIDATEDATA Validate and orient data for sidFreq* functions.
%
%   [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u)
%
%   Shared data validation used by sidFreqBT, sidFreqETFE, sidFreqBTFDR.
%   Ensures column orientation, checks for NaN/Inf, complex data, and
%   size consistency.
%
%   Supports multi-trajectory input:
%     - 3D arrays: y is (N x n_y x L), u is (N x n_u x L)
%     - 2D arrays: y is (N x n_y), treated as L=1
%
%   INPUTS:
%     y - Output data, (N x n_y), (N x n_y x L), or vector
%     u - Input data, (N x n_u), (N x n_u x L), vector, or [] for time series
%
%   OUTPUTS:
%     y            - (N x ny x nTraj) oriented output data
%     u            - (N x nu x nTraj) oriented input data, or []
%     N            - Number of samples per trajectory
%     ny           - Number of output channels
%     nu           - Number of input channels (0 for time series)
%     isTimeSeries - Logical, true when u is empty
%     nTraj        - Number of trajectories (1 for single-trajectory)
%
%   EXAMPLES:
%     y = randn(500, 1); u = randn(500, 1);
%     [y, u, N, ny, nu, isTS] = sidValidateData(y, u);
%
%   SPECIFICATION:
%     (Input validation — not yet in SPEC.md)
%
%   See also: sidParseOptions, sidFreqBT
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

    % ---- Handle cell array input (variable-length trajectories) ----
    % Cell arrays are trimmed to the shortest trajectory length and
    % stacked into a 3D array. This follows SPEC.md §2.1.
    if iscell(y)
        isTimeSeries = isempty(u) || (iscell(u) && isempty(u));
        L = numel(y);
        if L == 0
            error('sid:badInput', 'Cell arrays must not be empty.');
        end
        if ~isTimeSeries
            if ~iscell(u)
                error('sid:badInput', ...
                    'When y is a cell array, u must also be a cell array or [].');
            end
            if numel(u) ~= L
                error('sid:dimMismatch', ...
                    'y has %d trajectories but u has %d.', L, numel(u));
            end
        end

        % Determine common length (trim to shortest)
        lengths = zeros(L, 1);
        for l = 1:L
            if isvector(y{l})
                y{l} = y{l}(:);
            end
            lengths(l) = size(y{l}, 1);
        end
        N_common = min(lengths);
        ny = size(y{1}, 2);

        % Stack into 3D array
        y_3d = zeros(N_common, ny, L);
        for l = 1:L
            if size(y{l}, 2) ~= ny
                error('sid:dimMismatch', ...
                    'y{%d} has %d columns, expected %d.', l, size(y{l}, 2), ny);
            end
            y_3d(:, :, l) = y{l}(1:N_common, :);
        end
        y = y_3d;

        if ~isTimeSeries
            nu = size(u{1}, 2);
            u_3d = zeros(N_common, nu, L);
            for l = 1:L
                if isvector(u{l})
                    u{l} = u{l}(:);
                end
                if size(u{l}, 2) ~= nu
                    error('sid:dimMismatch', ...
                        'u{%d} has %d columns, expected %d.', l, size(u{l}, 2), nu);
                end
                if size(u{l}, 1) < N_common
                    error('sid:sizeMismatch', ...
                        'u{%d} has %d samples but y requires at least %d.', ...
                        l, size(u{l}, 1), N_common);
                end
                u_3d(:, :, l) = u{l}(1:N_common, :);
            end
            u = u_3d;
        end

        if any(lengths ~= N_common)
            warning('sid:trimmedTrajectories', ...
                'Variable-length trajectories trimmed to shortest length N = %d.', ...
                N_common);
        end
    end

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
    ny = size(y, 2);

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
        nu = size(u, 2);
        if size(u, 1) ~= N
            error('sid:sizeMismatch', ...
                'Input u (%d samples) and output y (%d samples) must have the same length.', ...
                size(u, 1), N);
        end
        if ~isreal(u)
            error('sid:complexData', ...
                'Complex data is not supported in v1.0. Input u must be real.');
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
    else
        nu = 0;
    end

    if N < 10
        warning('sid:shortData', ...
            'Very short data (N = %d). Estimates will be unreliable.', N);
    end
end
