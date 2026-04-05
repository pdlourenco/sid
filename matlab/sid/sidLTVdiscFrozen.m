function result = sidLTVdiscFrozen(ltvResult, varargin)
% SIDLTVDISCFROZEN Frozen transfer function from LTV state-space model.
%
%   result = sidLTVdiscFrozen(ltvResult)
%   result = sidLTVdiscFrozen(ltvResult, 'Frequencies', w)
%   result = sidLTVdiscFrozen(ltvResult, 'TimeSteps', kVec)
%   result = sidLTVdiscFrozen(ltvResult, 'SampleTime', Ts)
%
%   Computes the frozen (instantaneous) transfer function at each time
%   step k and frequency w:
%
%       G(w, k) = (e^{jw} I - A(k))^{-1} B(k)
%
%   If the ltvResult includes uncertainty (from sidLTVdisc with
%   'Uncertainty', true), the standard deviation of G is propagated via
%   first-order (Jacobian) linearization.
%
%   INPUTS:
%     ltvResult - Result struct from sidLTVdisc.
%
%   NAME-VALUE OPTIONS:
%     'Frequencies' - (nf x 1) frequency vector in rad/sample.
%                     Default: 128 linearly spaced in (0, pi].
%     'TimeSteps'   - (nk x 1) indices of time steps to evaluate
%                     (1-based). Default: all 1:N.
%     'SampleTime'  - Sample time in seconds. Default: 1.0.
%
%   OUTPUTS:
%     result - Struct with fields:
%       .Frequency      - (nf x 1) rad/sample
%       .FrequencyHz    - (nf x 1) Hz
%       .TimeSteps      - (nk x 1) selected time step indices
%       .Response       - (nf x p x q x nk) complex transfer function
%       .ResponseStd    - (nf x p x q x nk) std dev ([] if no uncertainty)
%       .SampleTime     - scalar
%       .Method         - 'sidLTVdiscFrozen'
%
%   EXAMPLES:
%     % Basic usage
%     ltv = sidLTVdisc(X, U, 'Lambda', 1e5, 'Uncertainty', true);
%     frz = sidLTVdiscFrozen(ltv);
%
%     % Custom frequencies and selected time steps
%     w = logspace(-2, log10(pi), 200)';
%     frz = sidLTVdiscFrozen(ltv, 'Frequencies', w, 'TimeSteps', [1 50 100]);
%
%   SPECIFICATION:
%     SPEC.md §8.9 — Bayesian Uncertainty Estimation
%
%   See also: sidLTVdisc, sidBodePlot, sidMapPlot
%
%   Changelog:
%   2026-03-29: First version by Pedro Lourenço.
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

    % ---- Parse inputs ----
    nf_default = 128;
    defs.Frequencies = ((1:nf_default)' * pi / nf_default);
    defs.TimeSteps = [];
    defs.SampleTime = 1.0;
    opts = sidParseOptions(defs, varargin);
    w = opts.Frequencies;
    w = w(:);
    kVec = opts.TimeSteps;
    if ~isempty(kVec), kVec = kVec(:); end
    Ts = opts.SampleTime;

    % ---- Extract from ltvResult ----
    A = ltvResult.A;   % (p x p x N)
    B = ltvResult.B;   % (p x q x N)
    p = ltvResult.StateDim;
    q = ltvResult.InputDim;
    N = ltvResult.DataLength;

    hasUncertainty = isfield(ltvResult, 'P') && ~isempty(ltvResult.P);

    if isempty(kVec)
        kVec = (1:N)';
    end
    nk = length(kVec);
    nf = length(w);

    % Validate time step indices
    if any(kVec < 1) || any(kVec > N)
        error('sid:badTimeSteps', ...
            'TimeSteps must be in range [1, %d].', N);
    end

    % ---- Compute frozen transfer function (SPEC.md §8.6) ----
    % G(w, k) = (e^{jw} I - A(k))^{-1} B(k)
    G    = zeros(nf, p, q, nk);
    GStd = [];
    if hasUncertainty
        GStd = zeros(nf, p, q, nk);
    end

    Ip = eye(p);

    for ik = 1:nk
        ki = kVec(ik);
        Ak = A(:, :, ki);    % (p x p)
        Bk = B(:, :, ki);    % (p x q)

        for iw = 1:nf
            z = exp(1i * w(iw));
            % R = (zI - A(k))^{-1}
            R = (z * Ip - Ak) \ Ip;   % (p x p)
            Gk = R * Bk;              % (p x q)
            G(iw, :, :, ik) = Gk;
        end

        if hasUncertainty
            Pk = ltvResult.P(:, :, ki);   % (d x d), d = p+q
            Sigma = ltvResult.NoiseCov;    % (p x p)
            sigDiag = diag(Sigma);         % (p x 1)
            pDiag = diag(Pk);              % (d x 1)

            for iw = 1:nf
                z = exp(1i * w(iw));
                R = (z * Ip - Ak) \ Ip;    % (p x p)
                Gk = R * Bk;               % (p x q)

                % First-order uncertainty propagation (SPEC.md §8.6):
                % Cov(vec(C(k))) = Sigma kron P(k), so
                % Var(G_{ab}) = sum_{r,j} |dG_{ab}/dC_{rj}|^2 * Sigma_{jj} * P_{rr}
                %
                % Jacobians: dG_{ab}/dA_{ji} = R_{aj} * [R*B]_{ib}
                %            dG_{ab}/dB_{ji} = R_{aj} * delta_{ib}
                varG = zeros(p, q);
                for b = 1:q
                    for a = 1:p
                        v = 0;
                        % Contribution from A entries: dG_{ab}/dA_{ji} = R_{aj} * Gk_{ib}
                        % where Gk = R*B (already computed).
                        % C(k) stores A' in rows 1:p, so A_{ji} = C_{i,j}, row index i, col j
                        % Var(C_{i,j}) = Sigma_{jj} * P_{ii}
                        for j = 1:p
                            for i = 1:p
                                dG = R(a, j) * Gk(i, b);
                                v = v + abs(dG)^2 * sigDiag(j) * pDiag(i);
                            end
                        end
                        % Contribution from B entries: dG_{ab}/dB_{ji} = R_{aj} * delta_{ib}
                        % C(k) stores B' in rows p+1:d, so B_{ji} = C_{p+i,j}, row p+i, col j
                        % Var(C_{p+i,j}) = Sigma_{jj} * P_{p+i,p+i}
                        for j = 1:p
                            dG = R(a, j);
                            v = v + abs(dG)^2 * sigDiag(j) * pDiag(p + b);
                        end
                        varG(a, b) = v;
                    end
                end
                GStd(iw, :, :, ik) = sqrt(varG);
            end
        end
    end

    % ---- Pack result ----
    result.Frequency    = w;
    result.FrequencyHz  = w / (2 * pi * Ts);
    result.TimeSteps    = kVec;
    result.Response     = G;
    result.ResponseStd  = GStd;
    result.SampleTime   = Ts;
    result.Method       = 'sidLTVdiscFrozen';
end
