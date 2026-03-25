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
%
%   Example:
%   TODO add example code here
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

    tau = (0:M)';
    W = 0.5 * (1 + cos(pi * tau / M));
end
