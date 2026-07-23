function [G, PhiV, Coh] = sidRegularizeResponse( ...
        PhiYU, PhiU, PhiY, ny, nu, forceDegenerate, epsReg, doWarn)
% SIDREGULARIZERESPONSE Form G, Phi_v and coherence with degenerate handling.
%
%   [G, PhiV, Coh] = sidRegularizeResponse(PhiYU, PhiU, PhiY, ny, nu)
%   [G, PhiV, Coh] = sidRegularizeResponse(..., forceDegenerate, epsReg, doWarn)
%
%   Shared frequency-response / noise-spectrum former for every Blackman-
%   Tukey-style estimator (sidFreqBT, sidFreqBTFDR, and the Welch path of
%   sidFreqMap), implementing the SPEC §2.6 near-singular Phi_u guard, the
%   §2.7 PSD clamp, and the §10.3 whole-signal degenerate override. Factored
%   out so the estimators cannot drift apart on these rules, and kept in lock
%   step with the Python regularize_response helper (issues #143, #141).
%
%   The sigma_G = Inf sentinel of §3.3 is applied downstream in
%   sidUncertainty (which sees Coh < eps or NaN in G); it is not formed here.
%
%   INPUTS:
%     PhiYU           - Cross-spectrum. SISO: (nf x 1). MIMO: (nf x ny x nu).
%     PhiU            - Input auto-spectrum. SISO: (nf x 1). MIMO: (nf x nu x nu).
%     PhiY            - Output auto-spectrum. SISO: (nf x 1). MIMO: (nf x ny x ny).
%     ny, nu          - Output / input channel counts (disambiguate SISO from
%                       a MIMO case that MATLAB would squeeze to 2-D).
%     forceDegenerate - When true (the §10.3 input-excitation check fired),
%                       treat every frequency as degenerate: G = NaN, Phi_v =
%                       Phi_y, coherence 0. Default false.
%     epsReg          - Relative floor (SISO) / rcond threshold (MIMO), §2.6.
%                       Default 1e-10.
%     doWarn          - Emit the near-singular warning when any frequency is
%                       degenerate. Default true.
%
%   OUTPUTS:
%     G    - Frequency response, (nf x 1) (SISO) or (nf x ny x nu) (MIMO).
%     PhiV - Noise spectrum, (nf x 1) (SISO) or (nf x ny x ny) (MIMO).
%     Coh  - Squared coherence (nf x 1) for SISO; [] for MIMO.
%
%   EXAMPLES:
%     [G, PhiV, Coh] = sidRegularizeResponse(PhiYU, PhiU, PhiY, 1, 1);
%
%   SPECIFICATION:
%     SPEC.md §2.6 — Regularization; §2.7 — Noise Spectrum PSD clamp
%
%   See also: sidInputExcitationDegenerate, sidDeadInputChannels, sidUncertainty
%
%   Changelog:
%   2026-07-23: First version by Pedro Lourenço.
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

    if nargin < 6 || isempty(forceDegenerate), forceDegenerate = false; end
    if nargin < 7 || isempty(epsReg),          epsReg = 1e-10;          end
    if nargin < 8 || isempty(doWarn),          doWarn = true;           end

    if ny == 1 && nu == 1
        % ---- SISO ----
        PhiYUc = PhiYU(:);
        PhiUr = real(PhiU(:));
        PhiYr = real(PhiY(:));

        PhiUabs = abs(PhiUr);
        PhiUmax = max(PhiUabs);
        % Defense-in-depth: an all-zero Phi_u makes the relative mask vacuous
        % (0 < 0 is false), so treat it as fully degenerate even if a caller
        % skipped the §10.3 input-excitation check -- otherwise the silent
        % NaN / no-warning shape of #143 returns by omission.
        if forceDegenerate || PhiUmax <= realmin
            singularMask = true(size(PhiUabs));
            forceDegenerate = true;
        else
            singularMask = PhiUabs < epsReg * PhiUmax;
        end

        G = PhiYUc ./ PhiU(:);   % divide by the complex Phi_u (parity w/ old)
        G(singularMask) = NaN + 1j * NaN;

        PhiV = PhiYr - abs(PhiYUc).^2 ./ PhiUr;
        PhiV(singularMask) = PhiYr(singularMask);
        PhiV = clampPSDscalar(PhiV);

        Coh = abs(PhiYUc).^2 ./ (PhiYr .* PhiUr);
        Coh(singularMask) = 0;
        Coh = min(max(Coh, 0), 1);

        if doWarn && any(singularMask)
            emitDegenerateWarning(forceDegenerate);
        end
        return;
    end

    % ---- MIMO ----
    nf = size(PhiYU, 1);
    G = zeros(nf, ny, nu);
    PhiV = zeros(nf, ny, ny);
    Coh = [];
    anySingular = false;

    for k = 1:nf
        PhiU_k = reshape(PhiU(k, :, :), nu, nu);
        PhiYU_k = reshape(PhiYU(k, :, :), ny, nu);
        PhiY_k = real(reshape(PhiY(k, :, :), ny, ny));

        if forceDegenerate || ~all(isfinite(PhiU_k(:))) || rcond(PhiU_k) < epsReg
            G(k, :, :) = NaN;
            PhiV(k, :, :) = clampPSDmatrix(PhiY_k);   % all output is noise
            anySingular = true;
        else
            Gk = PhiYU_k / PhiU_k;                     % PhiYU_k * inv(PhiU_k)
            G(k, :, :) = Gk;
            PhiV(k, :, :) = clampPSDmatrix(PhiY_k - Gk * PhiYU_k');
        end
    end

    if doWarn && anySingular
        emitDegenerateWarning(forceDegenerate);
    end
end

function emitDegenerateWarning(forceDegenerate)
% Single canonical wording shared with the Python helper.
    if forceDegenerate
        warning('sid:constantInput', ...
            ['Input u has negligible variation (constant or zero). The ' ...
             'frequency response is unidentifiable; G set to NaN and ' ...
             'sigma_G to Inf everywhere.']);
    else
        warning('sid:singularPhiU', ...
            ['Input spectrum Phi_u is near-singular at some frequencies. ' ...
             'G set to NaN at those points.']);
    end
end

function out = clampPSDscalar(phiV)
% Clamp a real scalar noise spectrum to be non-negative (§2.7). Only finite
% negatives are zeroed; NaN (degenerate frequency) is preserved.
    out = phiV;
    neg = isfinite(out) & (out < 0);
    out(neg) = 0;
end

function Vc = clampPSDmatrix(V)
% Clamp a single (ny x ny) noise-spectrum matrix to PSD (§2.7): take the real
% part FIRST, symmetrise, then zero any negative eigenvalues. A matrix with
% NaN (degenerate frequency) is returned unchanged so the NaN survives.
    V = real(V);
    if ~all(isfinite(V(:)))
        Vc = V;
        return;
    end
    Vk = (V + V.') / 2;
    [vecs, vals] = eig(Vk);
    d = diag(vals);
    if any(d < 0)
        d = max(d, 0);
        Vk = vecs * diag(d) * vecs.';
    end
    Vc = real(Vk);
end
