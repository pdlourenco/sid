function [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u, preserveLengths)
% SIDVALIDATEDATA Validate and orient data for sidFreq* functions.
%
%   [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u)
%   [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u, preserveLengths)
%
%   Shared data validation used by sidFreqBT, sidFreqETFE, sidFreqBTFDR.
%   Ensures column orientation, checks for NaN/Inf, complex data, and
%   size consistency.
%
%   Supports multi-trajectory input:
%     - 3D arrays: y is (N x n_y x L), u is (N x n_u x L)
%     - 2D arrays: y is (N x n_y), treated as L=1
%     - Cell arrays: {y1, y2, ...} for variable-length trajectories
%
%   By default, variable-length cell input is trimmed to the shortest
%   trajectory length and a warning is emitted.  When preserveLengths is
%   true, variable-length cell input is returned as a cell array of
%   per-trajectory matrices instead; callers (e.g. sidFreqMap) are then
%   responsible for per-segment filtering (SPEC.md §6.2).
%
%   INPUTS:
%     y                - Output data, (N x n_y), (N x n_y x L), cell, or vector
%     u                - Input data, (N x n_u), (N x n_u x L), cell, vector,
%                        or [] for time series
%     preserveLengths  - (optional) logical, default false. When true and y
%                        is a variable-length cell array, the returned y/u
%                        are cell arrays of per-trajectory data and N is
%                        max(N_l).
%
%   OUTPUTS:
%     y            - (N x ny x nTraj) oriented output data, or cell array
%                    of per-trajectory (N_l x ny) when preserving lengths
%     u            - (N x nu x nTraj) oriented input data, cell array, or []
%     N            - Number of samples per trajectory; max(N_l) in preserved
%                    variable-length mode, otherwise the common trimmed length
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
%   https://github.com/pdlourenco/sid
%  -----------------------------------------------------------------------

    if nargin < 3 || isempty(preserveLengths)
        preserveLengths = false;
    end

    % ---- Handle cell array input (variable-length trajectories) ----
    % Default: trimmed to the shortest trajectory length and stacked into a
    % 3D array (SPEC.md §2.1).  When preserveLengths is true, variable-length
    % cell input is returned as a cell array of per-trajectory matrices so
    % callers can do per-segment filtering (SPEC.md §6.2 for sidFreqMap).
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

        % Orient every trajectory and validate contents before we decide
        % whether to trim or preserve lengths.
        lengths = zeros(L, 1);
        for l = 1:L
            if isvector(y{l})
                y{l} = y{l}(:);
            end
            if ~isreal(y{l})
                error('sid:complexData', ...
                    'Complex data is not supported in v1.0. y{%d} must be real.', l);
            end
            if any(~isfinite(y{l}(:)))
                error('sid:nonFinite', 'Data y{%d} contains NaN or Inf values.', l);
            end
            lengths(l) = size(y{l}, 1);
        end
        ny = size(y{1}, 2);
        for l = 2:L
            if size(y{l}, 2) ~= ny
                error('sid:dimMismatch', ...
                    'y{%d} has %d columns, expected %d.', l, size(y{l}, 2), ny);
            end
        end

        nu_raw = 0;
        if ~isTimeSeries
            nu_raw = size(u{1}, 2);
            for l = 1:L
                if isvector(u{l})
                    u{l} = u{l}(:);
                end
                if ~isreal(u{l})
                    error('sid:complexData', ...
                        'Complex data is not supported in v1.0. u{%d} must be real.', l);
                end
                if any(~isfinite(u{l}(:)))
                    error('sid:nonFinite', 'Data u{%d} contains NaN or Inf values.', l);
                end
                if size(u{l}, 2) ~= nu_raw
                    error('sid:dimMismatch', ...
                        'u{%d} has %d columns, expected %d.', l, size(u{l}, 2), nu_raw);
                end
                if size(u{l}, 1) ~= lengths(l)
                    error('sid:sizeMismatch', ...
                        ['u{%d} has %d samples but y{%d} has %d. ' ...
                         'Trajectory-wise lengths must match.'], ...
                        l, size(u{l}, 1), l, lengths(l));
                end
            end
        end

        variableLength = any(lengths ~= lengths(1));

        if preserveLengths && variableLength
            % Return cell arrays as-is so callers can do per-segment filtering.
            N = max(lengths);
            if N < 2
                error('sid:tooShort', 'Data must have at least 2 samples.');
            end
            if ~isTimeSeries
                nu = nu_raw;
            else
                nu = 0;
            end
            nTraj = L;
            if N < 10
                warning('sid:shortData', ...
                    'Very short data (N_max = %d). Estimates will be unreliable.', N);
            end
            return;
        end

        % Default path: trim to shortest and stack into 3-D arrays.
        N_common = min(lengths);

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
