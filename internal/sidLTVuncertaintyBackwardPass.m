function P = sidLTVuncertaintyBackwardPass(S_scaled, lambda, N, d)
% SIDLTVUNCERTAINTYBACKWARDPASS Compute P(k) = [A_unscaled^{-1}]_{kk}.
%
%   The COSMIC algorithm normalizes data by 1/sqrt(N), so S_scaled contains
%   D_s'D_s + regularization. The posterior covariance requires the UNSCALED
%   Hessian where D'D = N * D_s'D_s appears. This function reconstructs the
%   unscaled diagonal blocks, then computes left and right Schur complements
%   to obtain the diagonal blocks of the inverse.
%
%   P(k) = (Lbd_k^L + Lbd_k^R - S_kk)^{-1}
%
%   Complexity: O(N * d^3).

    I = eye(d);

    % ---- Reconstruct unscaled diagonal blocks S_u(k) = N*DtD(k) + reg(k) ----
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

    % ---- Left Schur complements (forward) ----
    LbdL = zeros(d, d, N);
    LbdL(:, :, 1) = S(:, :, 1);
    for k = 2:N
        LbdL(:, :, k) = S(:, :, k) - lambda(k-1)^2 * (LbdL(:, :, k-1) \ I);
    end

    % ---- Right Schur complements (backward) ----
    LbdR = zeros(d, d, N);
    LbdR(:, :, N) = S(:, :, N);
    for k = N-1:-1:1
        LbdR(:, :, k) = S(:, :, k) - lambda(k)^2 * (LbdR(:, :, k+1) \ I);
    end

    % ---- Combine: P(k) = (LbdL(k) + LbdR(k) - S(k))^{-1} ----
    P = zeros(d, d, N);
    for k = 1:N
        M = LbdL(:, :, k) + LbdR(:, :, k) - S(:, :, k);
        P(:, :, k) = M \ I;
    end
end
