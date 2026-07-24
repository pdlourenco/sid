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

    % ---- Schur recursion on the SCALED Hessian blocks (SPEC.md §8.9.2) ----
    % The returned MAP minimises ||unscaled residual||^2 + N*lambda*||dC||^2
    % (effective prior weight N*lambda, from the 1/sqrt(N) scaling of §8.3.2), so
    % its Hessian is A_est = N * A_scaled and the reported posterior is
    % P(k) = [A_est^{-1}]_kk = P_scaled(k) / N. Run the left/right Schur
    % complements on the scaled blocks (un-inflated lambda_k couplings) to get
    % P_scaled, then divide by N -- exactly equivalent to inflating S -> N*S with
    % couplings (N*lambda_k)^2, but minimal. The previous reconstruction rebuilt
    % V'V + lambda F'F (unscaled data, un-inflated prior -- a different
    % estimator), overstating P by up to a factor N (issue #137).
    S = S_scaled;

    % ---- Left Schur complements — forward pass (SPEC.md §8.9.2) ----
    LbdL = zeros(d, d, N);
    LbdL(:, :, 1) = S(:, :, 1);
    for k = 2:N
        LbdL(:, :, k) = S(:, :, k) - lambda(k-1)^2 * (LbdL(:, :, k-1) \ I);
    end

    % ---- Right Schur complements — backward pass (SPEC.md §8.9.2) ----
    LbdR = zeros(d, d, N);
    LbdR(:, :, N) = S(:, :, N);
    for k = N-1:-1:1
        LbdR(:, :, k) = S(:, :, k) - lambda(k)^2 * (LbdR(:, :, k+1) \ I);
    end

    % ---- Combine: P_scaled(k) = (LbdL + LbdR - S)^{-1}, then P = P_scaled / N -
    P = zeros(d, d, N);
    for k = 1:N
        M = LbdL(:, :, k) + LbdR(:, :, k) - S(:, :, k);
        P(:, :, k) = (M \ I) / N;
    end
end
