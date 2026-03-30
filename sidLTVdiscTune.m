function [bestResult, bestLambda, thirdOut] = sidLTVdiscTune(varargin)
%SIDLTVDISCTUNE Lambda tuning for sidLTVdisc (validation or frequency-based).
%
%   VALIDATION METHOD (requires held-out state data):
%
%   [bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, ...
%       X_val, U_val, 'Method', 'validation', ...)
%
%   FREQUENCY METHOD (no validation data needed):
%
%   [bestResult, bestLambda, info] = sidLTVdiscTune(X, U, ...
%       'Method', 'frequency', ...)
%
%   The frequency method compares the COSMIC frozen transfer function
%   against a non-parametric sidFreqMap estimate using Mahalanobis-like
%   consistency scoring. It selects the largest lambda where the two
%   estimates agree statistically at most (omega, t) grid points.
%
%   VALIDATION METHOD INPUTS:
%     X_train - Training state data, (N+1 x p x L_train)
%     U_train - Training input data, (N x q x L_train)
%     X_val   - Validation state data, (N+1 x p x L_val)
%     U_val   - Validation input data, (N x q x L_val)
%
%   FREQUENCY METHOD INPUTS:
%     X       - State data, (N+1 x p x L) or cell array
%     U       - Input data, (N x q x L) or cell array
%
%   NAME-VALUE OPTIONS (both methods):
%     'Method'       - 'validation' (default) or 'frequency'.
%     'LambdaGrid'   - Vector of candidate lambda values.
%                      Default: logspace(0, 10, 25) (frequency),
%                               logspace(-3, 15, 50) (validation).
%     'Precondition'  - Passed through to sidLTVdisc. Default: false.
%     'Algorithm'     - Passed through to sidLTVdisc. Default: 'cosmic'.
%
%   FREQUENCY METHOD OPTIONS:
%     'SegmentLength'          - Outer segment length for sidFreqMap.
%                                Default: min(floor(N/4), 256).
%     'ConsistencyThreshold'   - Fraction of grid points required to be
%                                consistent. Default: 0.90.
%     'CoherenceThreshold'     - Minimum coherence for a grid point to be
%                                included. Default: 0.3.
%
%   OUTPUTS:
%     bestResult  - sidLTVdisc result struct at optimal lambda.
%     bestLambda  - Optimal scalar lambda value.
%     thirdOut    - For 'validation': (nGrid x 1) trajectory RMSE per lambda.
%                   For 'frequency': info struct with fields:
%                     .lambdaGrid, .fractions, .bestFraction, .freqMapResult
%
%   EXAMPLES:
%     % Validation-based
%     [best, lam, losses] = sidLTVdiscTune(Xtr, Utr, Xval, Uval);
%
%     % Frequency-based (no validation data)
%     [best, lam, info] = sidLTVdiscTune(X, U, 'Method', 'frequency');
%
%   REFERENCE:
%     Carvalho et al., "COSMIC", arXiv:2112.04355, 2022.
%     Ljung, L. "System Identification", 2nd ed., Prentice Hall, 1999.
%
%   See also: sidLTVdisc, sidLTVdiscFrozen, sidFreqMap

    % ---- Detect method from varargin ----
    method = 'validation';
    for k = 1:length(varargin)
        if ischar(varargin{k}) && strcmpi(varargin{k}, 'method')
            if k < length(varargin)
                method = lower(varargin{k+1});
            end
            break;
        end
    end

    switch method
        case 'validation'
            [bestResult, bestLambda, thirdOut] = validationTune(varargin{:});
        case 'frequency'
            [bestResult, bestLambda, thirdOut] = frequencyTune(varargin{:});
        otherwise
            error('sid:badMethod', ...
                'Method must be ''validation'' or ''frequency''. Got ''%s''.', method);
    end
end


% ========================================================================
%  VALIDATION-BASED TUNING
% ========================================================================

function [bestResult, bestLambda, allLosses] = validationTune(X_train, U_train, X_val, U_val, varargin)
%VALIDATIONTUNE Grid search over lambda, evaluated by trajectory RMSE.

    % ---- Parse options ----
    lambdaGrid = logspace(-3, 15, 50);
    extraArgs = {};

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'lambdagrid'
                    lambdaGrid = varargin{k+1};
                    k = k + 2;
                case 'method'
                    k = k + 2;  % skip, already handled
                case {'precondition', 'algorithm'}
                    extraArgs = [extraArgs, varargin(k), varargin(k+1)]; %#ok<AGROW>
                    k = k + 2;
                otherwise
                    error('sid:unknownOption', 'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badInput', 'Expected option name at position %d.', k);
        end
    end

    lambdaGrid = lambdaGrid(:)';
    nGrid = length(lambdaGrid);

    % Dimensions
    N = size(X_val, 1) - 1;
    p = size(X_val, 2);
    L_val = size(X_val, 3);

    % ---- Grid search ----
    allLosses = zeros(nGrid, 1);

    for j = 1:nGrid
        res = sidLTVdisc(X_train, U_train, 'Lambda', lambdaGrid(j), extraArgs{:});
        allLosses(j) = trajectoryRMSE(res.A, res.B, X_val, U_val, N, p, L_val);
    end

    % ---- Select best ----
    [~, bestIdx] = min(allLosses);
    bestLambda = lambdaGrid(bestIdx);
    bestResult = sidLTVdisc(X_train, U_train, 'Lambda', bestLambda, extraArgs{:});
end


% ========================================================================
%  FREQUENCY-RESPONSE CONSISTENCY TUNING
% ========================================================================

function [bestResult, bestLambda, info] = frequencyTune(X, U, varargin)
%FREQUENCYTUNE Select lambda via frequency-response consistency scoring.
%
%   Compares COSMIC frozen TF against non-parametric sidFreqMap estimate.
%   Selects the largest lambda where >=threshold fraction of (omega, t)
%   grid points are consistent at 95% confidence.

    % ---- Parse options ----
    lambdaGrid = logspace(0, 10, 25);
    segLen = [];
    consistThresh = 0.90;
    cohThresh = 0.3;
    extraArgs = {};

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'lambdagrid'
                    lambdaGrid = varargin{k+1};
                    k = k + 2;
                case 'method'
                    k = k + 2;  % skip, already handled
                case 'segmentlength'
                    segLen = varargin{k+1};
                    k = k + 2;
                case 'consistencythreshold'
                    consistThresh = varargin{k+1};
                    k = k + 2;
                case 'coherencethreshold'
                    cohThresh = varargin{k+1};
                    k = k + 2;
                case {'precondition', 'algorithm'}
                    extraArgs = [extraArgs, varargin(k), varargin(k+1)]; %#ok<AGROW>
                    k = k + 2;
                otherwise
                    error('sid:unknownOption', 'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badInput', 'Expected option name at position %d.', k);
        end
    end

    lambdaGrid = sort(lambdaGrid(:)');
    nGrid = length(lambdaGrid);

    % ---- Dimensions ----
    if iscell(X)
        N = size(X{1}, 1) - 1;
        p = size(X{1}, 2);
    else
        N = size(X, 1) - 1;
        p = size(X, 2);
    end
    if isempty(segLen)
        segLen = min(floor(N / 4), 256);
    end

    % ---- Step 1: Run sidFreqMap on training data (state as output) ----
    % Treat each state component x_i as an output channel, with u as input.
    % Use all trajectories (multi-trajectory support from Phase 9a).
    if iscell(X)
        % Variable-length cell arrays: use first trajectory for sidFreqMap
        y_freq = X{1}(1:end-1, :);  % (N x p) — drop last row (x has N+1 rows)
        u_freq = U{1};              % (N x q)
    else
        nTrajData = size(X, 3);
        y_freq = X(1:N, :, :);     % (N x p x L) — state at times 1..N
        u_freq = U;                 % (N x q x L)
        if nTrajData == 1
            y_freq = y_freq(:, :);
            u_freq = u_freq(:, :);
        end
    end

    fmapResult = sidFreqMap(y_freq, u_freq, 'SegmentLength', segLen);
    fmapFreqs = fmapResult.Frequency;

    % ---- Align time grids ----
    % sidFreqMap segment center times (in samples, 1-based)
    segCenterSamples = fmapResult.Time / fmapResult.SampleTime;
    kNearest = max(1, min(N, round(segCenterSamples)));
    nk = length(kNearest);
    nf = length(fmapFreqs);

    % ---- Extract sidFreqMap data ----
    % Response and ResponseStd are (nf x K) for SISO, (nf x K x ny x nu) for MIMO
    G_data = fmapResult.Response;
    GStd_data = fmapResult.ResponseStd;

    % Coherence mask: only score points with sufficient coherence
    if ~isempty(fmapResult.Coherence)
        cohMask = fmapResult.Coherence >= cohThresh;  % (nf x K)
    else
        % MIMO: no per-point coherence, use all points
        cohMask = true(nf, nk);
    end

    % Chi-square threshold for 95% confidence, 2 DOF (complex scalar SISO)
    chi2thresh = 5.991;
    % For MIMO p x q complex entries: DOF = 2*p*q
    if p > 1
        dof = 2 * p * size(U, 2);
        % Approximate chi2inv(0.95, dof) using Wilson-Hilferty
        z = 1.6449;  % norminv(0.95)
        chi2thresh = dof * (1 - 2/(9*dof) + z*sqrt(2/(9*dof)))^3;
    end

    % ---- Step 2: Grid search with Mahalanobis scoring ----
    fractions = zeros(nGrid, 1);

    for j = 1:nGrid
        % Run COSMIC with uncertainty
        res = sidLTVdisc(X, U, 'Lambda', lambdaGrid(j), 'Uncertainty', true, extraArgs{:});

        % Frozen transfer function at aligned time steps and matching frequencies
        frz = sidLTVdiscFrozen(res, 'Frequencies', fmapFreqs, 'TimeSteps', kNearest);

        % Compute consistency fraction
        fractions(j) = computeConsistency(G_data, GStd_data, ...
            frz.Response, frz.ResponseStd, cohMask, chi2thresh, p);
    end

    % ---- Step 3: Select largest lambda with sufficient consistency ----
    consistent = find(fractions >= consistThresh);
    if ~isempty(consistent)
        bestIdx = consistent(end);  % largest lambda (grid is sorted ascending)
        bestLambda = lambdaGrid(bestIdx);
    else
        % Fallback: no lambda meets threshold, use best available
        [~, bestIdx] = max(fractions);
        bestLambda = lambdaGrid(bestIdx);
        warning('sid:noConsistentLambda', ...
            'No lambda achieved %.0f%% consistency. Using best (%.1f%% at lambda=%.2e).', ...
            consistThresh * 100, fractions(bestIdx) * 100, bestLambda);
    end

    % ---- Re-run at optimal lambda with uncertainty ----
    bestResult = sidLTVdisc(X, U, 'Lambda', bestLambda, 'Uncertainty', true, extraArgs{:});

    % ---- Pack info struct ----
    info.lambdaGrid     = lambdaGrid(:);
    info.fractions      = fractions;
    info.bestFraction   = fractions(bestIdx);
    info.freqMapResult  = fmapResult;
    info.chi2Threshold  = chi2thresh;
end


% ========================================================================
%  LOCAL HELPER FUNCTIONS
% ========================================================================

function frac = computeConsistency(G_data, GStd_data, G_frozen, GStd_frozen, ...
                                    cohMask, chi2thresh, p)
%COMPUTECONSISTENCY Mahalanobis-like consistency between two TF estimates.
%
%   Returns the fraction of valid (ω, t) grid points where the parametric
%   (frozen) and non-parametric (sidFreqMap) estimates are consistent.

    nf = size(G_data, 1);
    nk = size(G_data, 2);

    if p == 1
        % SISO: G_data is (nf x K), G_frozen is (nf x 1 x 1 x nk)
        G_frz = reshape(G_frozen, nf, nk);
        GStd_frz = reshape(GStd_frozen, nf, nk);
        G_dat = G_data;
        GStd_dat = GStd_data;

        % Mahalanobis distance: d² = |G_frozen - G_data|² / (σ²_f + σ²_d)
        denominator = GStd_frz.^2 + GStd_dat.^2;
        denominator(denominator < eps) = eps;
        d2 = abs(G_frz - G_dat).^2 ./ denominator;

        % Consistent where d² < chi² threshold AND coherence is sufficient
        isConsistent = (d2 < chi2thresh) & cohMask;
        isValid = cohMask;

    else
        % MIMO: compare element-wise, average d² across channels
        % G_data: (nf x K x p x q), G_frozen: (nf x p x q x nk)
        % Reshape frozen to match: (nf x nk x p x q)
        G_frz = permute(G_frozen, [1 4 2 3]);
        GStd_frz = permute(GStd_frozen, [1 4 2 3]);
        G_dat = G_data;
        GStd_dat = GStd_data;

        denominator = GStd_frz.^2 + GStd_dat.^2;
        denominator(denominator < eps) = eps;
        d2_all = abs(G_frz - G_dat).^2 ./ denominator;

        % Sum d² across channels → (nf x nk)
        d2 = sum(sum(d2_all, 3), 4);

        isConsistent = (d2 < chi2thresh) & cohMask;
        isValid = cohMask;
    end

    nValid = sum(isValid(:));
    if nValid == 0
        frac = 0;
    else
        frac = sum(isConsistent(:)) / nValid;
    end
end


function rmse = trajectoryRMSE(A, B, X_val, U_val, N, p, L_val)
%TRAJECTORYRMSE Average trajectory prediction RMSE over validation set.

    totalRMSE = 0;

    for l = 1:L_val
        x_hat = zeros(N + 1, p);
        x_hat(1, :) = X_val(1, :, l);   % initial condition from data

        for kk = 1:N
            x_hat(kk+1, :) = (A(:, :, kk) * x_hat(kk, :)' + ...
                             B(:, :, kk) * reshape(U_val(kk, :, l), [], 1))';
        end

        x_true = X_val(:, :, l);
        err = x_hat - x_true;
        totalRMSE = totalRMSE + sqrt(mean(sum(err.^2, 2)));
    end

    rmse = totalRMSE / L_val;
end
