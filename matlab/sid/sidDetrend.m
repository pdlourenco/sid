function [x_detrended, trend] = sidDetrend(x, varargin)
% SIDDETREND Remove polynomial trend from time-domain data.
%
%   x_dt = sidDetrend(x)
%   x_dt = sidDetrend(x, 'Order', d)
%   [x_dt, trend] = sidDetrend(x, 'SegmentLength', L)
%
%   Removes a polynomial trend of degree d from each channel of x.
%   Detrending is standard preprocessing before spectral estimation:
%   unremoved trends bias low-frequency spectral estimates and violate
%   stationarity assumptions.
%
%   INPUTS:
%     x - (N x n_ch) or (N x n_ch x L) real data matrix. Vectors are
%         treated as single-channel column data.
%
%   NAME-VALUE OPTIONS:
%     'Order'         - Polynomial degree to remove (default: 1).
%                         0 = remove mean, 1 = linear, 2 = quadratic, etc.
%     'SegmentLength' - Detrend each non-overlapping segment independently
%                       (default: N, i.e., full record).
%
%   OUTPUTS:
%     x_detrended - (N x n_ch) or (N x n_ch x L), same size as input
%     trend       - (N x n_ch) or (N x n_ch x L), the removed trend
%                   (x = x_detrended + trend)
%
%   EXAMPLES:
%     y_dt = sidDetrend(y);                    % remove linear trend
%     y_dm = sidDetrend(y, 'Order', 0);        % remove mean only
%     [y_dt, trend] = sidDetrend(y);           % also get the trend
%     y_ds = sidDetrend(y, 'SegmentLength', 500); % segment-wise
%
%     % Typical workflow
%     y_dt = sidDetrend(y);
%     u_dt = sidDetrend(u);
%     result = sidFreqBT(y_dt, u_dt);
%
%   SPECIFICATION:
%     (Data preprocessing — not yet in SPEC.md)
%
%   See also: sidFreqBT, sidFreqETFE, sidFreqMap
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
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    % ---- Input validation ----
    if isvector(x)
        x = x(:);
    end
    if ~isreal(x)
        error('sid:complexData', 'Input x must be real.');
    end
    if any(~isfinite(x(:)))
        error('sid:nonFinite', 'Input x contains NaN or Inf values.');
    end

    N = size(x, 1);
    nCh = size(x, 2);
    if ndims(x) == 3 %#ok<ISMAT>
        nTraj = size(x, 3);
    else
        nTraj = 1;
    end

    % ---- Parse options ----
    defs.Order = 1;
    defs.SegmentLength = N;
    opts = sidParseOptions(defs, varargin);
    order = opts.Order;
    segLen = opts.SegmentLength;

    % ---- Validate parameters ----
    if ~isscalar(order) || order < 0 || order ~= round(order)
        error('sid:badOrder', 'Order must be a non-negative integer.');
    end
    if ~isscalar(segLen) || segLen < 1 || segLen ~= round(segLen)
        error('sid:badSegmentLength', 'SegmentLength must be a positive integer.');
    end
    if segLen > N
        segLen = N;
    end

    % ---- Detrend ----
    trend = zeros(size(x));

    for lt = 1:nTraj
        for ch = 1:nCh
            % Extract channel data
            if nTraj > 1
                col = x(:, ch, lt);
            else
                col = x(:, ch);
            end

            trendCol = zeros(N, 1);

            % Process each segment
            idx = 1;
            while idx <= N
                segEnd = min(idx + segLen - 1, N);
                segN = segEnd - idx + 1;
                t = (0:segN-1)';

                seg = col(idx:segEnd);
                actualOrder = min(order, segN - 1);
                if actualOrder < order
                    warning('sid:detrendOrderReduced', ...
                        ['Segment of length %d is too short for polynomial order %d. ' ...
                         'Reduced to order %d.'], segN, order, actualOrder);
                end
                coeffs = polyfit(t, seg, actualOrder);
                trendCol(idx:segEnd) = polyval(coeffs, t);

                idx = segEnd + 1;
            end

            if nTraj > 1
                trend(:, ch, lt) = trendCol;
            else
                trend(:, ch) = trendCol;
            end
        end
    end

    x_detrended = x - trend;
end
