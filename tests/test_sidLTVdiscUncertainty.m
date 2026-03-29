%% test_sidLTVdiscUncertainty - Unit tests for COSMIC Bayesian uncertainty
%
% Tests the uncertainty backward pass, noise covariance estimation,
% user-provided noise covariance, covariance modes, and Monte Carlo
% validation of posterior standard deviations.

fprintf('Running test_sidLTVdiscUncertainty...\n');

%% Test 1: Uncertainty fields present and correct dimensions
rng(1001);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
sigma = 0.02;
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + sigma * randn(1, p);
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true);

d = p + q;
assert(isfield(result, 'AStd'), 'Missing field: AStd');
assert(isfield(result, 'BStd'), 'Missing field: BStd');
assert(isfield(result, 'P'), 'Missing field: P');
assert(isfield(result, 'NoiseCov'), 'Missing field: NoiseCov');
assert(isfield(result, 'NoiseCovEstimated'), 'Missing field: NoiseCovEstimated');
assert(isfield(result, 'NoiseVariance'), 'Missing field: NoiseVariance');
assert(isfield(result, 'DegreesOfFreedom'), 'Missing field: DegreesOfFreedom');

assert(isequal(size(result.AStd), [p, p, N]), 'AStd should be (p x p x N)');
assert(isequal(size(result.BStd), [p, q, N]), 'BStd should be (p x q x N)');
assert(isequal(size(result.P), [d, d, N]), 'P should be (d x d x N)');
assert(isequal(size(result.NoiseCov), [p, p]), 'NoiseCov should be (p x p)');
assert(result.NoiseCovEstimated == true, 'NoiseCovEstimated should be true');
assert(isscalar(result.NoiseVariance), 'NoiseVariance should be scalar');
assert(isscalar(result.DegreesOfFreedom), 'DegreesOfFreedom should be scalar');
fprintf('  Test 1 passed: uncertainty fields present with correct dimensions.\n');

%% Test 2: AStd, BStd are positive and finite
assert(all(result.AStd(:) > 0), 'AStd entries should be positive');
assert(all(result.BStd(:) > 0), 'BStd entries should be positive');
assert(all(isfinite(result.AStd(:))), 'AStd entries should be finite');
assert(all(isfinite(result.BStd(:))), 'BStd entries should be finite');
fprintf('  Test 2 passed: AStd, BStd are positive and finite.\n');

%% Test 3: P(k) matrices are symmetric positive definite
for k = 1:N
    Pk = result.P(:, :, k);
    asymm = norm(Pk - Pk', 'fro') / max(norm(Pk, 'fro'), eps);
    assert(asymm < 1e-10, 'P(%d) is not symmetric (asymm=%.2e)', k, asymm);
    eigvals = eig(Pk);
    assert(all(eigvals > 0), 'P(%d) is not positive definite', k);
end
fprintf('  Test 3 passed: all P(k) are symmetric positive definite.\n');

%% Test 4: No-regularization limit (lambda -> 0)
% With very small lambda, P(k) should approach (D(k)'*D(k))^{-1} (OLS).
rng(1002);
p = 2; q = 1; N = 10; L = 8;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.01 * randn(1, p);
    end
end
lambda_small = 1e-8;
result = sidLTVdisc(X, U, 'Lambda', lambda_small, 'Uncertainty', true);

% Reconstruct D(k)'*D(k) and check P(k) ~ (D'D)^{-1}
d = p + q;
sqrtN = sqrt(N);
for k = 1:N
    Dk = [reshape(X(k, :, :), p, L)', reshape(U(k, :, :), q, L)'] / sqrtN;
    DtD = Dk' * Dk;
    P_ols = inv(DtD);  %#ok<MINV>
    Pk = result.P(:, :, k);
    relErr = norm(Pk - P_ols, 'fro') / norm(P_ols, 'fro');
    assert(relErr < 0.1, 'No-reg limit: P(%d) should be ~(D''D)^{-1}, relErr=%.3f', k, relErr);
end
fprintf('  Test 4 passed: no-regularization limit matches OLS covariance.\n');

%% Test 5: Large-lambda limit: P(k) shrinks
rng(1003);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end

res_lo = sidLTVdisc(X, U, 'Lambda', 1e2, 'Uncertainty', true);
res_hi = sidLTVdisc(X, U, 'Lambda', 1e8, 'Uncertainty', true);

% Average trace of P should be smaller for larger lambda (more pooling)
trP_lo = 0; trP_hi = 0;
for k = 1:N
    trP_lo = trP_lo + trace(res_lo.P(:, :, k));
    trP_hi = trP_hi + trace(res_hi.P(:, :, k));
end
assert(trP_hi < trP_lo, ...
    'Large lambda should give smaller P (trP_hi=%.4e, trP_lo=%.4e)', trP_hi, trP_lo);
fprintf('  Test 5 passed: large lambda shrinks P(k).\n');

%% Test 6: Single time step (N=2 gives 2 steps, use known OLS)
rng(1004);
p = 1; q = 1; N = 2; L = 10;
A_true = 0.8; B_true = 1.0;
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn;
    for k = 1:N
        X(k+1, :, l) = A_true * X(k, :, l) + B_true * U(k, :, l) + 0.01 * randn;
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e-10, 'Uncertainty', true);
% With near-zero lambda and N=2, each P(k) should be close to (D'D)^{-1}
for k = 1:N
    eigP = eig(result.P(:, :, k));
    assert(all(eigP > 0), 'P(%d) should be PD for single-step limit', k);
end
fprintf('  Test 6 passed: N=2 case produces valid uncertainty.\n');

%% Test 7: Isotropic noise -> equal AStd columns
rng(1005);
p = 3; q = 1; N = 30; L = 10;
A_true = 0.8 * eye(p);
B_true = ones(p, q);
sigma = 0.05;
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + sigma * randn(1, p);
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true, 'CovarianceMode', 'isotropic');

% With isotropic noise, all columns of AStd at each time step should have
% the same pattern (same P diagonal, same sigma for all)
for k = 1:N
    AStd_k = result.AStd(:, :, k);
    % Each column a has AStd(b,a,k) = sqrt(sigma^2 * P(a,a))
    % So column ratios should reflect P diagonal ratios, same across b
    col1 = AStd_k(:, 1);
    for a = 2:p
        ratio = AStd_k(:, a) ./ col1;
        % All entries in ratio should be the same scalar (sqrt(P(a,a)/P(1,1)))
        assert(max(ratio) - min(ratio) < 1e-10, ...
            'Isotropic noise: AStd columns should have consistent ratios at k=%d', k);
    end
end
fprintf('  Test 7 passed: isotropic noise gives consistent AStd structure.\n');

%% Test 8: CovarianceMode diagonal vs full vs isotropic
rng(1006);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end

res_diag = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true, 'CovarianceMode', 'diagonal');
res_full = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true, 'CovarianceMode', 'full');
res_iso  = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true, 'CovarianceMode', 'isotropic');

% Diagonal: off-diagonal of NoiseCov should be zero
offdiag = res_diag.NoiseCov - diag(diag(res_diag.NoiseCov));
assert(norm(offdiag, 'fro') < 1e-15, 'Diagonal mode: NoiseCov should be diagonal');

% Isotropic: NoiseCov should be scalar * I
expected_iso = (trace(res_iso.NoiseCov) / p) * eye(p);
assert(norm(res_iso.NoiseCov - expected_iso, 'fro') < 1e-15, ...
    'Isotropic mode: NoiseCov should be sigma^2 * I');

% Full: NoiseCov may have off-diagonals
assert(isequal(size(res_full.NoiseCov), [p, p]), 'Full mode: NoiseCov should be p x p');

% All modes should produce positive AStd
assert(all(res_diag.AStd(:) > 0), 'Diagonal AStd should be positive');
assert(all(res_full.AStd(:) > 0), 'Full AStd should be positive');
assert(all(res_iso.AStd(:) > 0), 'Isotropic AStd should be positive');
fprintf('  Test 8 passed: covariance modes produce correct structure.\n');

%% Test 9: User-provided NoiseCov
rng(1007);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end

Sigma_known = diag([0.01, 0.05]);
result = sidLTVdisc(X, U, 'Lambda', 1e4, 'NoiseCov', Sigma_known);

% NoiseCov should be returned exactly as provided
assert(isequal(result.NoiseCov, Sigma_known), 'User-provided NoiseCov should be stored exactly');
assert(result.NoiseCovEstimated == false, 'NoiseCovEstimated should be false');
assert(isnan(result.DegreesOfFreedom), 'DegreesOfFreedom should be NaN for user-provided cov');

% AStd should reflect the provided Sigma
assert(all(result.AStd(:) > 0), 'AStd should be positive with user NoiseCov');
assert(all(isfinite(result.AStd(:))), 'AStd should be finite with user NoiseCov');

% Providing NoiseCov should auto-enable uncertainty
assert(isfield(result, 'P'), 'NoiseCov should auto-enable uncertainty');
fprintf('  Test 9 passed: user-provided NoiseCov works correctly.\n');

%% Test 10: NoiseCov validation errors
rng(1008);
p = 2; q = 1; N = 10; L = 3;
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (0.9*eye(p) * X(k, :, l)' + ones(p,q) * U(k, :, l)')';
    end
end

% Wrong size NoiseCov
try
    sidLTVdisc(X, U, 'Lambda', 1e3, 'NoiseCov', eye(3));
    error('sid:testFailed', 'Should have thrown error for wrong-size NoiseCov');
catch e
    assert(strcmp(e.identifier, 'sid:badNoiseCov'), ...
        'Expected sid:badNoiseCov, got %s', e.identifier);
end

% NoiseCov with NaN
try
    badCov = eye(p); badCov(1,1) = NaN;
    sidLTVdisc(X, U, 'Lambda', 1e3, 'NoiseCov', badCov);
    error('sid:testFailed', 'Should have thrown error for NaN in NoiseCov');
catch e
    assert(strcmp(e.identifier, 'sid:badNoiseCov'), ...
        'Expected sid:badNoiseCov, got %s', e.identifier);
end
fprintf('  Test 10 passed: NoiseCov validation errors correct.\n');

%% Test 11: Uncertainty=false gives no uncertainty fields (backward compat)
rng(1009);
p = 2; q = 1; N = 15; L = 3;
X = randn(N+1, p, L); U = randn(N, q, L);
result = sidLTVdisc(X, U, 'Lambda', 1e3);
assert(~isfield(result, 'AStd'), 'Default should not have AStd');
assert(~isfield(result, 'P'), 'Default should not have P');
assert(~isfield(result, 'NoiseCov'), 'Default should not have NoiseCov');
fprintf('  Test 11 passed: backward compatibility (no uncertainty fields by default).\n');

%% Test 12: Variable-length trajectories with uncertainty
rng(1010);
p = 2; q = 1; L = 4;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
Nvec = [20, 15, 25, 18];
Xc = cell(L, 1); Uc = cell(L, 1);
for l = 1:L
    Nl = Nvec(l);
    Xc{l} = zeros(Nl+1, p);
    Uc{l} = randn(Nl, q);
    Xc{l}(1, :) = randn(1, p);
    for k = 1:Nl
        Xc{l}(k+1, :) = (A_true * Xc{l}(k, :)' + B_true * Uc{l}(k, :)')' + 0.02 * randn(1, p);
    end
end
result = sidLTVdisc(Xc, Uc, 'Lambda', 1e4, 'Uncertainty', true);
Nmax = max(Nvec);
d = p + q;
assert(isequal(size(result.AStd), [p, p, Nmax]), 'VarLen AStd dims');
assert(isequal(size(result.P), [d, d, Nmax]), 'VarLen P dims');
assert(all(result.AStd(:) > 0), 'VarLen AStd should be positive');
fprintf('  Test 12 passed: variable-length trajectories with uncertainty.\n');

%% Test 13: Monte Carlo validation
% Generate many realizations, check empirical std matches AStd.
% Use small lambda so Bayesian posterior ~ frequentist variance (§4.4 of
% cosmic_uncertainty_derivation.md: they converge as lambda -> 0).
% Use large L to stabilize the data Gram matrices across realizations.
rng(1011);
p = 2; q = 1; N = 15; L = 30;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
sigma = 0.05;
Sigma_true = sigma^2 * eye(p);
lambda_val = 0.001;
nMC = 200;

A_samples = zeros(p, p, N, nMC);
for mc = 1:nMC
    X = zeros(N+1, p, L); U = randn(N, q, L);
    for l = 1:L
        X(1, :, l) = randn(1, p);
        for k = 1:N
            X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + ...
                sigma * randn(1, p);
        end
    end
    res = sidLTVdisc(X, U, 'Lambda', lambda_val);
    A_samples(:, :, :, mc) = res.A;
end

% Run one more with uncertainty to get AStd
rng(1012);
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + ...
            sigma * randn(1, p);
    end
end
res_unc = sidLTVdisc(X, U, 'Lambda', lambda_val, 'Uncertainty', true, ...
    'NoiseCov', Sigma_true);

% Empirical std across MC realizations
A_emp_std = std(A_samples, 0, 4);  % (p x p x N)

% Check at interior time steps (edges have boundary effects)
k_interior = 4:N-3;
for k = k_interior
    for b = 1:p
        for a = 1:p
            emp = A_emp_std(b, a, k);
            pred = res_unc.AStd(b, a, k);
            if emp > 1e-8
                ratio = pred / emp;
                assert(ratio > 0.3 && ratio < 3.0, ...
                    'MC validation: A(%d,%d,k=%d) ratio=%.2f out of [0.3, 3.0]', ...
                    b, a, k, ratio);
            end
        end
    end
end
fprintf('  Test 13 passed: Monte Carlo validates AStd (200 realizations).\n');

%% Test 14: DegreesOfFreedom is reasonable
rng(1013);
p = 2; q = 1; N = 20; L = 5;
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (0.9*eye(p) * X(k, :, l)' + ones(p,q) * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true);

% Total observations = N * L (for uniform trajectories)
totalObs = N * L;
% dof should be between 0 and totalObs
assert(result.DegreesOfFreedom > 0, 'DoF should be positive');
assert(result.DegreesOfFreedom < totalObs, 'DoF should be < total observations');
fprintf('  Test 14 passed: degrees of freedom is reasonable (dof=%.1f).\n', result.DegreesOfFreedom);

fprintf('test_sidLTVdiscUncertainty: ALL TESTS PASSED\n');
