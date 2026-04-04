function [A0, B0] = sidLTIfreqIO(Y, U, H, varargin)
% SIDLTIFREQIO Identify LTI state-space model from input-output data.
%
%   [A0, B0] = sidLTIfreqIO(Y, U, H)
%   [A0, B0] = sidLTIfreqIO(Y, U, H, 'Horizon', r)
%   [A0, B0] = sidLTIfreqIO(Y, U, H, 'MaxStabilize', s)
%
%   Estimates a constant (LTI) state-space realization from input-output
%   data, given a known observation matrix H:
%
%       x(k+1) = A0 x(k) + B0 u(k)
%       y(k)   = H  x(k)
%
%   Uses the Ho-Kalman realization algorithm applied to the frequency
%   response estimated via sidFreqBT. The realization is transformed
%   to the H-basis so that the observation equation y = H x holds
%   exactly with the returned (A0, B0).
%
%   INPUTS:
%     Y - Output data, (N+1 x py), (N+1 x py x L), or cell array {L x 1}
%         where Y{l} is (N_l+1 x py). For cells, trajectories shorter
%         than 2/3 of max horizon are discarded; the rest are trimmed
%         to a common length for spectral estimation.
%     U - Input data, (N x q), (N x q x L), or cell array {L x 1}
%         where U{l} is (N_l x q). Must match Y format.
%     H - Observation matrix, (py x n).
%
%   NAME-VALUE OPTIONS:
%     'Horizon'       - Hankel matrix depth r. Default: min(floor(nf/3), 50)
%                       where nf is the number of frequency bins.
%     'MaxStabilize'  - Maximum eigenvalue magnitude after stabilization.
%                       Unstable eigenvalues (|lambda| > threshold) are
%                       reflected inside the unit circle: lambda <- 1/conj(lambda),
%                       then clamped to this radius. Default: 0.999.
%
%   OUTPUTS:
%     A0 - (n x n) estimated LTI dynamics matrix.
%     B0 - (n x q) estimated LTI input matrix.
%
%   EXAMPLES:
%     % Estimate LTI dynamics from partial observations
%     [A0, B0] = sidLTIfreqIO(Y, U, H);
%
%     % With custom Hankel horizon
%     [A0, B0] = sidLTIfreqIO(Y, U, H, 'Horizon', 30);
%
%   ALGORITHM:
%     1. Estimate transfer function G(e^{jw}) = H(zI-A)^{-1}B via sidFreqBT
%     2. Compute Markov parameters g(k) = H A^{k-1} B via IFFT
%     3. Build block Hankel matrices H_0 and H_1 (shifted)
%     4. SVD of H_0, truncate to order n (Ho-Kalman realization)
%     5. Transform realization to H-basis: find T s.t. C_r T^{-1} = H
%     6. Stabilize eigenvalues if needed
%
%   REFERENCES:
%     Ho, B.L. and Kalman, R.E. "Effective construction of linear
%     state-variable models from input/output functions." Regelungstechnik,
%     14(12):545-548, 1966.
%
%   SPECIFICATION:
%     SPEC.md section 8.12 -- Output-COSMIC (LTI initialization)
%
%   See also: sidFreqBT, sidModelOrder, sidLTVdiscIO, sidLTVStateEst
%
%   Changelog:
%   2026-04-02: First version by Pedro Lourenco.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenco, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    % ---- Parse inputs ----
    [Y_trim, U, H, horizon, maxStab, py, n, q] = parseInputs( ...
        Y, U, H, varargin{:});

    % ---- Step 1: Frequency response estimation (SPEC.md §8.13) ----
    G_result = sidFreqBT(Y_trim, U);
    G = G_result.Response;  % (nf x py x q) complex
    nf = size(G, 1);

    % ---- Step 2: Impulse response via IFFT (SPEC.md §8.13) ----
    % Markov parameters: g(k) = H A^{k-1} B
    g = freqToImpulse(G, nf, py, q);  % (N_imp x py x q)
    N_imp = size(g, 1);

    % ---- Step 3: Determine horizon ----
    if isempty(horizon)
        horizon = min(floor(N_imp / 3), 50);
    end

    if horizon < 2
        error('sid:tooShort', ...
            'Data too short for Hankel matrix (need N_imp >= 6, got %d).', ...
            N_imp);
    end

    % Need at least 2*horizon impulse response coefficients for H_1
    if 2 * horizon > N_imp
        horizon = floor(N_imp / 2);
        if horizon < 2
            error('sid:tooShort', ...
                'Data too short for shifted Hankel matrix.');
        end
    end

    r = horizon;

    % Check that Hankel matrix can support order n
    if r * py < n || r * q < n
        error('sid:tooShort', ...
            ['Hankel size (%d*%d x %d*%d) too small for ' ...
             'order n=%d. Increase data length or Horizon.'], ...
            r, py, r, q, n);
    end

    % ---- Step 4: Build block Hankel matrices (SPEC.md §8.13) ----
    % H_0{i,j} = g(i+j-1), H_1{i,j} = g(i+j)
    [H0, H1] = buildHankel(g, r, py, q, N_imp);

    % ---- Step 5: Ho-Kalman SVD realization (SPEC.md §8.13) ----
    [A_r, B_r, C_r] = hoKalman(H0, H1, n, py, q);

    % ---- Step 6: Transform to H-basis ----
    [A0, B0] = transformToHBasis(A_r, B_r, C_r, H, n);

    % ---- Step 7: Stabilize ----
    A0 = stabilize(A0, maxStab);

end

% ========================================================================
%  LOCAL FUNCTIONS
% ========================================================================

function [Y_trim, U, H, horizon, maxStab, py, n, q] = parseInputs( ...
    Y, U, H, varargin)
% PARSEINPUTS Validate and parse inputs for sidLTIfreqIO.
%   For cell array inputs, trajectories shorter than 2/3 of the maximum
%   horizon are discarded and the rest are trimmed to a common length.

    py = size(H, 1);
    n  = size(H, 2);

    if iscell(Y)
        % ---- Variable-length: trim and stack into 3D ----
        if ~iscell(U)
            error('sid:badInput', ...
                'When Y is a cell array, U must also be a cell array.');
        end
        L = numel(Y);
        if numel(U) ~= L
            error('sid:dimMismatch', ...
                'Y has %d trajectories but U has %d.', L, numel(U));
        end
        if L == 0
            error('sid:badInput', 'Cell arrays must not be empty.');
        end

        q = size(U{1}, 2);
        horizons = zeros(L, 1);
        for l = 1:L
            if size(Y{l}, 2) ~= py
                error('sid:dimMismatch', ...
                    'Y{%d} has %d columns but H has %d rows.', ...
                    l, size(Y{l}, 2), py);
            end
            if size(U{l}, 2) ~= q
                error('sid:dimMismatch', ...
                    'U{%d} has %d columns, expected %d.', ...
                    l, size(U{l}, 2), q);
            end
            horizons(l) = size(U{l}, 1);
        end

        % Keep trajectories with horizon >= ceil(2/3 * max)
        maxH = max(horizons);
        threshold = ceil(2 * maxH / 3);
        keep = find(horizons >= threshold);
        if isempty(keep)
            error('sid:tooShort', ...
                'No trajectories long enough for spectral estimation.');
        end

        N_common = min(horizons(keep));
        L_kept = length(keep);

        % Stack into 3D arrays, trimmed to common length
        Y_trim = zeros(N_common, py, L_kept);
        U_out  = zeros(N_common, q, L_kept);
        for ii = 1:L_kept
            l = keep(ii);
            Y_trim(:, :, ii) = Y{l}(1:N_common, :);
            U_out(:, :, ii)  = U{l}(1:N_common, :);
        end
        U = U_out;
    else
        % ---- Uniform-horizon: existing path ----
        if ndims(Y) == 2  %#ok<ISMAT>
            Y = reshape(Y, size(Y, 1), size(Y, 2), 1);
        end
        if ndims(U) == 2  %#ok<ISMAT>
            U = reshape(U, size(U, 1), size(U, 2), 1);
        end

        N = size(U, 1);
        q = size(U, 2);

        if size(Y, 1) < N
            error('sid:dimMismatch', ...
                'Y must have at least N=%d rows, got %d.', ...
                N, size(Y, 1));
        end
        if size(Y, 2) ~= py
            error('sid:dimMismatch', ...
                'Y has %d columns but H has %d rows.', ...
                size(Y, 2), py);
        end
        if size(U, 3) ~= size(Y, 3)
            error('sid:dimMismatch', ...
                'U has %d trajectories but Y has %d.', ...
                size(U, 3), size(Y, 3));
        end

        % Trim Y to first N rows to match U
        Y_trim = Y(1:N, :, :);
    end

    % Defaults
    horizon = [];
    maxStab = 0.999;

    % Parse name-value options
    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'horizon'
                    horizon = varargin{k + 1};
                    k = k + 2;
                case {'maxstabilize', 'maxstability'}
                    maxStab = varargin{k + 1};
                    k = k + 2;
                otherwise
                    error('sid:badOption', ...
                        'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badOption', ...
                'Expected option name (string), got %s.', ...
                class(varargin{k}));
        end
    end
end

function g = freqToImpulse(G, nf, ny, nu)
% FREQTOIMPULSE Convert frequency response to impulse response via IFFT.
%
%   Reuses the conjugate-symmetric IFFT pattern from sidModelOrder.

    Nfft = 2 * nf;
    g_all = zeros(Nfft, ny, nu);

    for iy = 1:ny
        for iu = 1:nu
            Gvec = squeeze(G(:, iy, iu));  % (nf x 1) complex

            % Build full-circle: DC, positive freqs, mirror of negative
            Gfull = zeros(Nfft, 1);
            Gfull(1) = real(Gvec(1));                     % DC approximation
            Gfull(2:nf) = Gvec(1:nf - 1);                % w1 to w_{nf-1}
            Gfull(nf + 1) = real(Gvec(nf));               % Nyquist (real)
            Gfull(nf + 2:Nfft) = conj(Gvec(nf - 1:-1:1)); % mirror

            g_all(:, iy, iu) = real(ifft(Gfull));
        end
    end

    % Use causal part starting from lag 1 (skip direct feedthrough at lag 0).
    % For x(k+1) = A x(k) + B u(k), y(k) = H x(k), the Markov parameters
    % are h(k) = H A^{k-1} B for k >= 1, with h(0) = 0 (no feedthrough).
    % g_all(1) is lag 0, g_all(2) is lag 1 = H*B, etc.
    N_imp = nf - 1;
    g = g_all(2:nf, :, :);
end

function [H0, H1] = buildHankel(g, r, py, q, N_imp)
% BUILDHANKEL Build block Hankel matrix H_0 and shifted Hankel H_1.
%
%   H_0: block(i,j) = g(i+j-1), i,j = 1..r
%   H_1: block(i,j) = g(i+j),   i,j = 1..r

    H0 = zeros(r * py, r * q);
    H1 = zeros(r * py, r * q);

    for bi = 1:r
        for bj = 1:r
            idx0 = bi + bj - 1;
            idx1 = bi + bj;
            if idx0 <= N_imp
                H0((bi - 1) * py + 1:bi * py, ...
                   (bj - 1) * q + 1:bj * q) = ...
                    reshape(g(idx0, :, :), py, q);
            end
            if idx1 <= N_imp
                H1((bi - 1) * py + 1:bi * py, ...
                   (bj - 1) * q + 1:bj * q) = ...
                    reshape(g(idx1, :, :), py, q);
            end
        end
    end
end

function [A_r, B_r, C_r] = hoKalman(H0, H1, n, py, q)
% HOKALMAN Ho-Kalman realization from Hankel matrices.
%
%   Given H_0 and H_1, compute the minimal realization (A_r, B_r, C_r)
%   of order n such that G(z) = C_r (zI - A_r)^{-1} B_r.

    [U_svd, Sigma, V_svd] = svd(H0, 0);

    sigmas = diag(Sigma);
    if length(sigmas) < n
        error('sid:tooFewSV', ...
            'Hankel SVD has only %d singular values but n=%d requested.', ...
            length(sigmas), n);
    end

    % Truncate to order n
    U_n = U_svd(:, 1:n);
    S_n = diag(sigmas(1:n));
    V_n = V_svd(:, 1:n);

    S_n_sqrt  = diag(sqrt(sigmas(1:n)));
    S_n_isqrt = diag(1 ./ sqrt(sigmas(1:n)));

    % Realization in arbitrary basis
    A_r = S_n_isqrt * U_n' * H1 * V_n * S_n_isqrt;  % (n x n)
    C_r = U_n(1:py, :) * S_n_sqrt;                    % (py x n)
    B_r = S_n_sqrt * V_n(1:q, :)';                    % (n x q)
end

function [A0, B0] = transformToHBasis(A_r, B_r, C_r, H, n)
% TRANSFORMTOHBASIS Transform realization to the H-basis.
%
%   Find T such that C_r * T^{-1} = H, then A0 = T * A_r * T^{-1},
%   B0 = T * B_r.

    % Tinv = pinv(C_r) * H + (I - pinv(C_r) * C_r)
    Cr_pinv = pinv(C_r);
    Tinv = Cr_pinv * H + eye(n) - Cr_pinv * C_r;

    % Check conditioning
    rc = rcond(Tinv);
    if rc < eps * 1e3
        warning('sid:illConditioned', ...
            'Basis transform Tinv is near-singular (rcond=%.2e). Using raw realization.', ...
            rc);
        A0 = real(A_r);
        B0 = real(B_r);
        return;
    end

    T = Tinv \ eye(n);
    A0 = real(T * A_r * Tinv);
    B0 = real(T * B_r);
end

function A = stabilize(A, maxStab)
% STABILIZE Reflect unstable eigenvalues inside the unit circle.
%
%   Eigenvalues with |lambda| > 1 are reflected: lambda <- 1/conj(lambda)
%   (SPEC.md §8.13.1 step 6). The result is then clamped so that
%   |lambda| <= maxStab.

    [V, D] = eig(A);
    d = diag(D);
    magnitudes = abs(d);

    if max(magnitudes) <= maxStab
        return;
    end

    % Step 1: Reflect eigenvalues outside the unit circle (SPEC.md §8.13.1)
    outside = magnitudes > 1;
    d(outside) = 1 ./ conj(d(outside));
    magnitudes(outside) = abs(d(outside));

    % Step 2: Clamp to maxStab radius
    toolarge = magnitudes > maxStab;
    d(toolarge) = maxStab * d(toolarge) ./ magnitudes(toolarge);

    A = real(V * diag(d) / V);
end
