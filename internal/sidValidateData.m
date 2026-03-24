function [y, u, N, ny, nu, isTimeSeries] = sidValidateData(y, u)
%SIDVALIDATEDATA Validate and orient data for sidFreq* functions.
%
%   [y, u, N, ny, nu, isTimeSeries] = sidValidateData(y, u)
%
%   Shared data validation used by sidFreqBT, sidFreqETFE, sidFreqBTFDR.
%   Ensures column orientation, checks for NaN/Inf, complex data, and
%   size consistency.
%
%   INPUTS:
%     y - Output data, (N x n_y) or vector
%     u - Input data, (N x n_u) or vector, or [] for time series
%
%   OUTPUTS:
%     y            - (N x ny) oriented output data
%     u            - (N x nu) oriented input data, or []
%     N            - Number of samples
%     ny           - Number of output channels
%     nu           - Number of input channels (0 for time series)
%     isTimeSeries - Logical, true when u is empty

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
    else
        nu = 0;
    end

    if N < 10
        warning('sid:shortData', ...
            'Very short data (N = %d). Estimates will be unreliable.', N);
    end
end
