function X = sidDFT(x, freqs, useFFT)
% SIDDFT Compute DFT of time-domain signal at specified frequencies.
%
%   X = sidDFT(x, freqs, useFFT)
%
%   Computes the discrete Fourier transform:
%
%     X(w_k) = sum_{t=1}^{N} x(t) * exp(-j * w_k * t)
%
%   INPUTS:
%     x      - (N x p) real data matrix (p channels)
%     freqs  - (nf x 1) frequency vector in rad/sample, in (0, pi]
%     useFFT - Logical. If true, use FFT fast path (requires default
%              linear grid of 128 points).
%
%   OUTPUTS:
%     X      - (nf x p) complex DFT values
%
%   EXAMPLES:
%     x = randn(500, 1);
%     freqs = (1:128)' * pi / 128;
%     X = sidDFT(x, freqs, true);
%
%   SPECIFICATION:
%     SPEC.md §4.1 — ETFE (DFT computation)
%
%   See also: sidWindowedDFT, sidFreqETFE
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

    [N, p] = size(x);
    nf = length(freqs);

    if useFFT
        % FFT fast path for default 128-point linear grid.
        % Frequencies are k * pi / 128 for k = 1..128, which correspond
        % to bins k+1 of a 256-point FFT (since bin k+1 = freq k*2*pi/256
        % = k*pi/128).
        L = 2 * nf;   % 256

        % Zero-pad x to length L if needed, or just compute FFT with L pts
        % Note: fft(x, L) zero-pads if N < L, or truncates if N > L.
        % For N > L, we need the full-length FFT and extract the right bins.
        % The DFT at w_k = k*pi/nf for a length-N signal can be obtained
        % from a length-N FFT, but the bins won't align with the 128-grid
        % unless N is a multiple of 256.
        %
        % Instead, use the standard approach: zero-pad to at least N and
        % pick the closest bins, OR use the direct path.
        %
        % For consistency with sidWindowedDFT which uses L=256 FFT on
        % covariance sequences (length 2M+1 << 256), the ETFE uses a
        % length-N FFT and extracts at the desired frequencies.

        Nfft = max(N, L);
        % Round up to next power of 2 for efficiency
        Nfft = 2^nextpow2(Nfft);

        Xfft = fft(x, Nfft);  % (Nfft x p)

        % The frequency of bin k (0-indexed) is k * 2*pi / Nfft.
        % We want frequencies freqs(i) = i * pi / nf for i = 1..nf.
        % So the bin index (0-based) = freqs(i) * Nfft / (2*pi).
        binIdx = round(freqs(:) * Nfft / (2 * pi));
        % Clamp to valid range [1, Nfft/2]
        binIdx = max(min(binIdx, Nfft/2), 1);

        % MATLAB 1-indexed: bin 0 → index 1, bin k → index k+1
        X = Xfft(binIdx + 1, :);

        % Phase correction: MATLAB's fft computes sum_{t=0}^{N-1} x(t+1)*exp(-j*2*pi*k*t/N)
        % which equals sum_{t=1}^{N} x(t)*exp(-j*2*pi*k*(t-1)/N).
        % Our convention is X(w) = sum_{t=1}^{N} x(t)*exp(-j*w*t).
        % Difference: exp(-j*w*1) vs exp(-j*w*0) → multiply by exp(-j*w).
        % Actually: fft gives sum x(t)*exp(-j*w*(t-1)), we want sum x(t)*exp(-j*w*t)
        % = exp(-j*w) * sum x(t)*exp(-j*w*(t-1))
        w_actual = binIdx * 2 * pi / Nfft;
        correction = exp(-1j * w_actual(:));  % (nf x 1)
        X = X .* repmat(correction, 1, p);
    else
        % Direct DFT at arbitrary frequencies
        X = zeros(nf, p);
        t = (1:N)';  % (N x 1) time indices

        for k = 1:nf
            e = exp(-1j * freqs(k) * t);  % (N x 1)
            X(k, :) = e.' * x;            % (1 x p)
        end
    end
end
