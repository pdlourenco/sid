function Y_pred = sidFreqDomainSim(G_model, freqs_model, u, N)
% SIDFREQDOMAINSIM Simulate frequency-domain model output via IFFT.
%
%   Y_pred = sidFreqDomainSim(G_model, freqs_model, u, N)
%
%   Filters the input signal u through the frequency response G_model by:
%     1. Computing FFT of u
%     2. Interpolating G onto the FFT frequency grid using magnitude/phase
%     3. Multiplying and taking IFFT
%
%   Frequencies outside the model's grid are set to zero (no extrapolation).
%   Interpolation is performed on log-magnitude and unwrapped phase to
%   avoid artifacts near resonances.
%
%   INPUTS:
%     G_model     - (nf x ny x nu) complex frequency response.
%     freqs_model - (nf x 1) frequency vector in rad/sample, in (0, pi].
%     u           - (N x nu) real input signal.
%     N           - Number of samples (= size(u, 1)).
%
%   OUTPUTS:
%     Y_pred - (N x ny) predicted output signal.
%
%   EXAMPLES:
%     G = sidFreqBT(y, u);
%     Y_pred = sidFreqDomainSim(G.Response, G.Frequency, u, length(u));
%
%   See also: sidCompare, sidResidual
%
%   Changelog:
%   2026-04-05: First version by Pedro Lourenco.
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

    ny = size(G_model, 2);
    if ndims(G_model) < 3 %#ok<ISMAT>
        nu = size(G_model, 3);
        if isempty(nu) || nu == 0
            nu = 1;
        end
        G_model = reshape(G_model, size(G_model, 1), ny, nu);
    else
        nu = size(G_model, 3);
    end

    nfft = N;
    freqs_fft = (1:floor(nfft / 2))' * (2 * pi / nfft);
    npos = length(freqs_fft);

    U_fft = fft(u, nfft, 1);
    Y_pred_fft = zeros(nfft, ny);

    for iy = 1:ny
        for iu = 1:nu
            Gij = reshape(G_model(:, iy, iu), [], 1);
            G_interp = interpGsafe(freqs_model(:), Gij, freqs_fft);
            Y_pred_fft(2:npos + 1, iy) = Y_pred_fft(2:npos + 1, iy) + ...
                G_interp .* U_fft(2:npos + 1, iu);
        end
    end

    % Conjugate symmetry for real output
    for iy = 1:ny
        if mod(nfft, 2) == 0
            Y_pred_fft(npos + 2:end, iy) = conj(Y_pred_fft(npos:-1:2, iy));
        else
            Y_pred_fft(npos + 2:end, iy) = conj(Y_pred_fft(npos + 1:-1:2, iy));
        end
    end

    Y_pred = real(ifft(Y_pred_fft, nfft, 1));
end

function G_interp = interpGsafe(freqs_model, G, freqs_target)
% INTERPGSAFE Interpolate complex transfer function with safe boundaries.
%
%   Uses magnitude/phase interpolation (log-magnitude + unwrapped phase)
%   which is more numerically stable near resonances than real/imag.
%   Frequencies outside the model grid are set to zero (no extrapolation).

    freqs_model = freqs_model(:);
    G = G(:);
    freqs_target = freqs_target(:);

    % Identify target frequencies within the model range
    fmin = freqs_model(1);
    fmax = freqs_model(end);
    inRange = (freqs_target >= fmin) & (freqs_target <= fmax);

    G_interp = zeros(length(freqs_target), 1);

    if ~any(inRange)
        return;
    end

    % Interpolate on log-magnitude and unwrapped phase for stability
    mag = abs(G);
    ph  = unwrap(angle(G));

    % Avoid log(0): clamp magnitude floor
    mag = max(mag, eps);

    logmag_interp = interp1(freqs_model, log(mag), freqs_target(inRange), 'linear');
    ph_interp     = interp1(freqs_model, ph,       freqs_target(inRange), 'linear');

    G_interp(inRange) = exp(logmag_interp) .* exp(1i * ph_interp);
end
