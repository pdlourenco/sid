function P = sidLTVuncertaintyBackwardPass(S_scaled, lambda, N, d)
% SIDLTVUNCERTAINTYBACKWARDPASS Bayesian posterior covariance diagonal blocks.
%
%   P = sidLTVuncertaintyBackwardPass(S_scaled, lambda, N, d)
%
%   Computes P(k) = [A_unscaled^{-1}]_{kk}, the diagonal blocks of the
%   inverse of the unscaled COSMIC Hessian. These are the row covariance
%   matrices needed for Bayesian uncertainty estimation.
%
%   The COSMIC algorithm normalizes data by 1/sqrt(N), so S_scaled contains
%   D_s'D_s + regularization. This function reconstructs the unscaled
%   diagonal blocks, then computes left and right Schur complements to
%   obtain the diagonal blocks of the inverse:
%
%     P(k) = (Lbd_k^L + Lbd_k^R - S_kk)^{-1}
%
%   INPUTS:
%     S_scaled - (d x d x N) scaled block diagonal terms (from
%                sidLTVbuildBlockTerms, with 1/sqrt(N) normalization).
%     lambda   - (N-1 x 1) regularization weights.
%     N        - Number of time steps.
%     d        - Combined dimension, d = p + q.
%
%   OUTPUTS:
%     P - (d x d x N) diagonal blocks of the inverse Hessian.
%         Cov(vec(C(k))) = Sigma ⊗ P(k).
%
%   EXAMPLES:
%     P = sidLTVuncertaintyBackwardPass(S, lambda, N, d);
%
%   ALGORITHM:
%     1. Reconstruct unscaled S_u(k) = N * DtD(k) + reg(k)
%     2. Forward pass: left Schur complements Lbd^L(k)
%     3. Backward pass: right Schur complements Lbd^R(k)
%     4. Combine: P(k) = (Lbd^L(k) + Lbd^R(k) - S_u(k))^{-1}
%     Complexity: O(N * d^3).
%
%   REFERENCES:
%     Carvalho, Soares, Lourenco, Ventura. "COSMIC: fast closed-form
%     identification from large-scale data for LTV systems."
%     arXiv:2112.04355, 2022.
%
%   SPECIFICATION:
%     SPEC.md §8.9 — Bayesian Uncertainty Estimation
%
%   See also: sidLTVcosmicSolve, sidLTVbuildBlockTerms, sidLTVdisc
%
%   Changelog:
%   2026-04-01: First version by Pedro Lourenço.
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

    I = eye(d);

    % ---- Reconstruct unscaled Hessian diagonal blocks (SPEC.md §8.9.2) ----
    % S_u(k) = N * D'D(k) + reg(k), undoing the 1/sqrt(N) data scaling
    S = zeros(d, d, N);
    for k = 1:N
        if k == 1
            reg = lambda(1) * I;
        elseif k == N
            reg = lambda(N-1) * I;
        else
            reg = (lambda(k-1) + lambda(k)) * I;
        end
        DtD_scaled = S_scaled(:, :, k) - reg;
        S(:, :, k) = N * DtD_scaled + reg;
    end

    % ---- Left Schur complements — forward pass (SPEC.md §8.9.3) ----
    LbdL = zeros(d, d, N);
    LbdL(:, :, 1) = S(:, :, 1);
    for k = 2:N
        LbdL(:, :, k) = S(:, :, k) - lambda(k-1)^2 * (LbdL(:, :, k-1) \ I);
    end

    % ---- Right Schur complements — backward pass (SPEC.md §8.9.3) ----
    LbdR = zeros(d, d, N);
    LbdR(:, :, N) = S(:, :, N);
    for k = N-1:-1:1
        LbdR(:, :, k) = S(:, :, k) - lambda(k)^2 * (LbdR(:, :, k+1) \ I);
    end

    % ---- Combine: P(k) = (LbdL(k) + LbdR(k) - S(k))^{-1} (SPEC.md §8.9.3) ----
    P = zeros(d, d, N);
    for k = 1:N
        M = LbdL(:, :, k) + LbdR(:, :, k) - S(:, :, k);
        P(:, :, k) = M \ I;
    end
end
