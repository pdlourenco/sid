function R = sidCov(x, z, maxLag)
%SIDCOV Biased cross-covariance estimate for lags 0..maxLag.
%
%   R = sidCov(x, z, maxLag)
%
%   Computes the biased sample cross-covariance:
%
%     R(tau) = (1/N) * sum_{t=1}^{N-tau} x(t+tau) * z(t)'
%
%   for tau = 0, 1, ..., maxLag.
%
%   INPUTS:
%     x      - (N x p) matrix, first signal
%     z      - (N x q) matrix, second signal (can equal x for auto-cov)
%     maxLag - maximum lag M (non-negative integer)
%
%   OUTPUT:
%     R      - (maxLag+1 x p x q) array of covariance estimates.
%              R(tau+1, :, :) is the covariance at lag tau.
%              For scalar signals, R is (maxLag+1 x 1).

    N = size(x, 1);
    p = size(x, 2);
    q = size(z, 2);

    R = zeros(maxLag + 1, p, q);

    for tau = 0:maxLag
        % x(t+tau) for t = 1..N-tau  =>  x(tau+1:N, :)
        % z(t)     for t = 1..N-tau  =>  z(1:N-tau, :)
        R(tau + 1, :, :) = (x(tau + 1 : N, :)' * z(1 : N - tau, :)) / N;
    end

    % Squeeze trailing singleton dimensions for scalar signals
    if p == 1 && q == 1
        R = R(:);
    end
end
