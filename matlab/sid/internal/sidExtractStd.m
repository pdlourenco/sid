function [AStd, BStd] = sidExtractStd(P, Sigma, N, p, q)
% SIDEXTRACTSTD Standard deviations of A(k) and B(k) entries.
%
%   [AStd, BStd] = sidExtractStd(P, Sigma, N, p, q)
%
%   Computes element-wise standard deviations of the time-varying system
%   matrices from the Kronecker posterior structure Cov(vec(C(k))) = Sigma
%   kron P(k).
%
%   INPUTS:
%     P     - (d x d x N) posterior covariance diagonal blocks, d = p + q.
%     Sigma - (p x p) noise covariance matrix.
%     N     - Number of time steps.
%     p     - State dimension.
%     q     - Input dimension.
%
%   OUTPUTS:
%     AStd - (p x p x N) standard deviations of A(k) entries.
%     BStd - (p x q x N) standard deviations of B(k) entries.
%
%   EXAMPLES:
%     [AStd, BStd] = sidExtractStd(P, Sigma, N, p, q);
%
%   SPECIFICATION:
%     SPEC.md §8.9.4 — Standard Deviations
%
%   See also: sidLTVdisc, sidLTVdiscIO, sidEstimateNoiseCov
%
%   Changelog:
%   2026-04-05: Extracted from sidLTVdisc by Pedro Lourenco.
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

    % Var(A(k)_{b,a}) = Sigma_{bb} * P(k)_{aa}       for a = 1,...,p
    % Var(B(k)_{b,a}) = Sigma_{bb} * P(k)_{p+a,p+a}  for a = 1,...,q
    %
    % Note: C(k) = [A(k)'; B(k)'], so row a of C(k) corresponds to:
    %   a = 1..p   -> column a of A(k) -> A(k)_{:,a}
    %   a = p+1..d -> column (a-p) of B(k) -> B(k)_{:,a-p}

    AStd = zeros(p, p, N);
    BStd = zeros(p, q, N);

    sigDiag = diag(Sigma);  % (p x 1)

    for k = 1:N
        pDiag = diag(P(:, :, k));  % (d x 1)

        for a = 1:p
            AStd(:, a, k) = sqrt(sigDiag * pDiag(a));
        end

        for a = 1:q
            BStd(:, a, k) = sqrt(sigDiag * pDiag(p + a));
        end
    end
end
