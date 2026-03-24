function W = sidHannWin(M)
%SIDHANNWIN Hann (Hanning) lag window for spectral analysis.
%
%   W = sidHannWin(M)
%
%   Returns the Hann window values for lags 0, 1, ..., M:
%
%     W(tau) = 0.5 * (1 + cos(pi * tau / M))
%
%   W(0) = 1 and W(M) = 0.
%
%   INPUT:
%     M - Window size (positive integer, M >= 2)
%
%   OUTPUT:
%     W - (M+1 x 1) vector of window values for lags 0..M

    tau = (0:M)';
    W = 0.5 * (1 + cos(pi * tau / M));
end
