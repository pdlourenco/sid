function Phi = sidWindowedDFT(R, W, freqs, useFFT, Rneg)
%SIDWINDOWEDDFT Windowed Fourier transform of covariance estimates.
%
%   Phi = sidWindowedDFT(R, W, freqs, useFFT)
%   Phi = sidWindowedDFT(R, W, freqs, useFFT, Rneg)
%
%   Computes the spectral estimate:
%
%     Phi(w) = sum_{tau=-M}^{M} R(tau) * W(|tau|) * exp(-j*w*tau)
%
%   Using either the FFT (when useFFT is true and the frequency grid is
%   the default linear grid) or direct summation.
%
%   For auto-covariance, R(-tau) = conj(R(tau)) is used automatically.
%   For cross-covariance, pass Rneg = R_zy (the reverse covariance) so
%   that R(-tau) = Rneg(tau) is used for negative lags.
%
%   INPUTS:
%     R      - Covariance estimates for lags 0..M.
%              Scalar signals: (M+1 x 1) vector.
%              Matrix signals: (M+1 x p x q) array.
%     W      - Hann window values for lags 0..M. (M+1 x 1) vector.
%     freqs  - (n_f x 1) frequency vector in rad/sample.
%     useFFT - Logical. If true, use FFT fast path.
%     Rneg   - (optional) Covariance for negative lags. Same size as R.
%              When provided, Rneg(tau) is used for R(-tau) instead of
%              conj(R(tau)). Required for cross-covariance estimation.
%
%   OUTPUT:
%     Phi    - Spectral estimate at each frequency.
%              Scalar signals: (n_f x 1) complex vector.
%              Matrix signals: (n_f x p x q) complex array.
%
%   Example:
%     R = sidCov(x, x, 30); W = sidHannWin(30);
%     freqs = (1:128)' * pi / 128;
%     Phi = sidWindowedDFT(R, W, freqs, true, R);
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

    if nargin < 5
        Rneg = [];
    end

    M = length(W) - 1;
    nf = length(freqs);

    % Determine signal dimensions
    if isvector(R)
        p = 1; q = 1;
        isScalar = true;
    else
        p = size(R, 2); q = size(R, 3);
        isScalar = false;
    end

    if useFFT && isScalar
        Rneg_s = Rneg;
        Phi = windowedDFT_FFT(R, W, nf, Rneg_s);
    elseif useFFT && ~isScalar
        Phi = zeros(nf, p, q);
        for ii = 1:p
            for jj = 1:q
                if isempty(Rneg)
                    Rneg_ij = [];
                else
                    Rneg_ij = Rneg(:, jj, ii);
                end
                Phi(:, ii, jj) = windowedDFT_FFT(R(:, ii, jj), W, nf, Rneg_ij);
            end
        end
    elseif ~useFFT && isScalar
        Rneg_s = Rneg;
        Phi = windowedDFT_direct(R, W, freqs, Rneg_s);
    else
        Phi = zeros(nf, p, q);
        for ii = 1:p
            for jj = 1:q
                if isempty(Rneg)
                    Rneg_ij = [];
                else
                    Rneg_ij = Rneg(:, jj, ii);
                end
                Phi(:, ii, jj) = windowedDFT_direct(R(:, ii, jj), W, freqs, Rneg_ij);
            end
        end
    end
end


function Phi = windowedDFT_FFT(R, W, nf, Rneg)
%WINDOWEDDFT_FFT FFT fast path for default linear frequency grid.
%
%   Assumes frequencies are (1:nf) * pi / nf, i.e. the bins produced by
%   a length-2*nf FFT.

    M = length(W) - 1;
    L = 2 * nf;   % FFT length: 256 for the default 128 frequencies

    % Build the windowed covariance sequence in FFT order.
    % Lags 0..M go into positions 1..M+1.
    % Lags -M..-1 go into positions L-M+1..L (wrapped negative lags).
    s = zeros(L, 1);

    % Lag 0
    s(1) = W(1) * R(1);

    % Positive lags 1..M  (using R(tau) for tau >= 0)
    for tau = 1:M
        s(tau + 1) = W(tau + 1) * R(tau + 1);
    end

    % Negative lags -1..-M
    % For auto-covariance: R(-tau) = conj(R(tau))
    % For cross-covariance: R_xy(-tau) = R_yx(tau), passed as Rneg
    for tau = 1:M
        if isempty(Rneg)
            s(L - tau + 1) = W(tau + 1) * conj(R(tau + 1));
        else
            s(L - tau + 1) = W(tau + 1) * Rneg(tau + 1);
        end
    end

    % FFT
    S = fft(s, L);

    % Extract bins 2..nf+1 (skipping DC at bin 1).
    % Bin k+1 corresponds to frequency k * 2*pi / L = k * pi / nf.
    Phi = S(2 : nf + 1);
end


function Phi = windowedDFT_direct(R, W, freqs, Rneg)
%WINDOWEDDFT_DIRECT Direct DFT at arbitrary frequencies.

    M = length(W) - 1;
    nf = length(freqs);
    Phi = zeros(nf, 1);

    for k = 1:nf
        w = freqs(k);

        % Lag 0 contribution
        val = W(1) * R(1);

        % Lags 1..M: positive and negative combined
        % For auto-covariance: R(-tau) = conj(R(tau))
        % For cross-covariance: R_xy(-tau) = Rneg(tau) = R_yx(tau)
        for tau = 1:M
            e = exp(-1j * w * tau);
            if isempty(Rneg)
                val = val + W(tau + 1) * (R(tau + 1) * e + conj(R(tau + 1)) * conj(e));
            else
                val = val + W(tau + 1) * (R(tau + 1) * e + Rneg(tau + 1) * conj(e));
            end
        end

        Phi(k) = val;
    end
end
