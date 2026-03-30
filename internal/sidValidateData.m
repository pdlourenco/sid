function [y, u, N, ny, nu, isTimeSeries, nTraj] = sidValidateData(y, u)
%SIDVALIDATEDATA Validate and orient data for sidFreq* functions.
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
%   Example:
%     y = randn(500, 1); u = randn(500, 1);
%     [y, u, N, ny, nu, isTS] = sidValidateData(y, u);
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
    else
        nu = 0;
    end

    if N < 10
        warning('sid:shortData', ...
            'Very short data (N = %d). Estimates will be unreliable.', N);
    end
end
