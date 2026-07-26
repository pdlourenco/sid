%% test_sidLTVdiscTune - Unit tests for validation-based lambda tuning
%
% Tests sidLTVdiscTune for output shapes, lambda selection quality,
% custom grids, and consistency with sidLTVdisc.

fprintf('Running test_sidLTVdiscTune...\n');
runner__nPassed = 0;

%% Helper: generate LTV data split into train/val
% Use an LTV system so the bias-variance tradeoff in lambda is genuine:
% too small lambda overfits per-step, too large lambda over-smooths.
rng(1000);
p = 2; q = 1; N = 30;
A0 = [0.95 0.1; -0.1 0.85];
dA = [-0.3 0.05; 0.05 -0.25];
B_true = [0.5; 0.3];
sigma = 0.15;

L_train = 3; L_val = 4;
X_train = zeros(N+1, p, L_train); U_train = randn(N, q, L_train);
X_val   = zeros(N+1, p, L_val);   U_val   = randn(N, q, L_val);

for l = 1:L_train
    X_train(1, :, l) = randn(1, p);
    for k = 1:N
        Ak = A0 + (k / N) * dA;
        X_train(k+1, :, l) = (Ak * X_train(k, :, l)' + ...
            B_true * U_train(k, :, l)')' + ...
            sigma * randn(1, p);
    end
end
for l = 1:L_val
    X_val(1, :, l) = randn(1, p);
    for k = 1:N
        Ak = A0 + (k / N) * dA;
        X_val(k+1, :, l) = (Ak * X_val(k, :, l)' + B_true * U_val(k, :, l)')' + sigma * randn(1, p);
    end
end

%% Test 1: Output shapes
grid = logspace(-2, 6, 15);
[bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, X_val, U_val, ...
    'LambdaGrid', grid);

assert(isstruct(bestResult), 'bestResult should be a struct');
assert(isscalar(bestLambda), 'bestLambda should be scalar');
assert(bestLambda > 0, 'bestLambda should be positive');
assert(isequal(size(allLosses), [15, 1]), 'allLosses should be (nGrid x 1)');
assert(all(allLosses > 0), 'allLosses should be positive');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: output shapes correct.\n');

%% Test 2: bestResult has correct fields
assert(isfield(bestResult, 'A'), 'bestResult should have A');
assert(isfield(bestResult, 'B'), 'bestResult should have B');
assert(strcmp(bestResult.Method, 'sidLTVdisc'), 'bestResult.Method should be sidLTVdisc');
assert(isequal(size(bestResult.A), [p, p, N]), 'bestResult.A dimensions');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 2 passed: bestResult has correct fields.\n');

%% Test 3: Optimal lambda not at grid boundary
% For an LTV problem, the optimum should be interior (bias-variance tradeoff)
[~, bestIdx] = min(allLosses);
assert(bestIdx > 1 && bestIdx < length(grid), ...
    'Optimal lambda should not be at grid boundary (idx=%d/%d)', bestIdx, length(grid));
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 3 passed: optimal lambda is interior (idx=%d/%d).\n', bestIdx, length(grid));

%% Test 4: Custom lambda grid with few points
grid_small = [1, 100, 10000];
[~, bestLam2, losses2] = sidLTVdiscTune(X_train, U_train, X_val, U_val, ...
    'LambdaGrid', grid_small);
assert(length(losses2) == 3, 'allLosses should have 3 entries');
assert(any(abs(bestLam2 - grid_small) < 1e-12), 'bestLambda should be from grid');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 4 passed: custom grid works.\n');

%% Test 5: Consistency with sidLTVdisc
% Re-run sidLTVdisc at the best lambda and verify results match
check = sidLTVdisc(X_train, U_train, 'Lambda', bestLambda);
assert(max(abs(bestResult.A(:) - check.A(:))) < 1e-10, ...
    'bestResult.A should match direct sidLTVdisc call');
assert(max(abs(bestResult.B(:) - check.B(:))) < 1e-10, ...
    'bestResult.B should match direct sidLTVdisc call');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: bestResult consistent with sidLTVdisc.\n');

%% Test 6: Passthrough of Precondition option (not_implemented in v1.0)
warning('off', 'sid:preconditionUnsupported');
[bestRes_pc, ~, ~] = sidLTVdiscTune(X_train, U_train, X_val, U_val, ...
    'LambdaGrid', [100, 1000, 10000], 'Precondition', true);
warning('on', 'sid:preconditionUnsupported');
% Preconditioning requested but not implemented in v1.0
assert(ischar(bestRes_pc.Preconditioned) && ...
    strcmp(bestRes_pc.Preconditioned, 'not_implemented'), ...
    'Precondition should be ''not_implemented'' when requested');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: Precondition option forwarded, reports not_implemented.\n');

%% ====================================================================
%  FREQUENCY-RESPONSE METHOD TESTS (Phase 8d)
%  ====================================================================

%% Test 7: Frequency method produces valid output
rng(2001);
p = 2; q = 1; N = 60; L = 5; sigma = 0.05;
A0 = [0.9 0.1; -0.05 0.8]; B_true = [0.5; 0.3];
Xf = zeros(N+1, p, L); Uf = randn(N, q, L);
for l = 1:L
    Xf(1, :, l) = 0.1*randn(1, p);
    for k = 1:N
        Xf(k+1, :, l) = (A0 * Xf(k, :, l)' + B_true * Uf(k, :, l)')' + sigma*randn(1, p);
    end
end

grid_freq = logspace(1, 8, 10);
[bestRes_f, bestLam_f, info_f] = sidLTVdiscTune(Xf, Uf, ...
    'Method', 'frequency', 'LambdaGrid', grid_freq, 'SegmentLength', 20);

assert(isstruct(bestRes_f), 'bestResult should be a struct');
assert(isscalar(bestLam_f) && bestLam_f > 0, 'bestLambda should be positive scalar');
assert(isfield(info_f, 'lambdaGrid'), 'info should have lambdaGrid');
assert(isfield(info_f, 'fractions'), 'info should have fractions');
assert(isfield(info_f, 'bestFraction'), 'info should have bestFraction');
assert(isfield(info_f, 'freqMapResults'), 'info should have freqMapResults');
assert(length(info_f.fractions) == length(grid_freq), 'fractions length should match grid');
assert(all(info_f.fractions >= 0 & info_f.fractions <= 1), 'fractions should be in [0,1]');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: frequency method produces valid output.\n');

%% Test 8: LTI system -> large lambda selected
% For a constant system, regularization helps — lambda should be large.
rng(2002);
p = 2; q = 1; N = 80; L = 8; sigma = 0.03;
A_lti = [0.9 0.1; -0.05 0.8]; B_lti = [0.5; 0.3];
Xlti = zeros(N+1, p, L); Ulti = randn(N, q, L);
for l = 1:L
    Xlti(1, :, l) = 0.1*randn(1, p);
    for k = 1:N
        Xlti(k+1, :, l) = (A_lti * Xlti(k, :, l)' + B_lti * Ulti(k, :, l)')' + sigma*randn(1, p);
    end
end

grid_lti = logspace(1, 10, 15);
[~, bestLam_lti, info_lti] = sidLTVdiscTune(Xlti, Ulti, ...
    'Method', 'frequency', 'LambdaGrid', grid_lti, 'SegmentLength', 25);

% For LTI: large lambda should be preferred (system is constant)
midGrid = sqrt(grid_lti(1) * grid_lti(end));
assert(bestLam_lti >= midGrid, ...
    'LTI system: bestLambda=%.2e should be >= midGrid=%.2e', bestLam_lti, midGrid);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: LTI system selects large lambda (%.2e).\n', bestLam_lti);

%% Test 9: LTV system -> moderate lambda (not extreme)
rng(2003);
p = 2; q = 1; N = 80; L = 8; sigma = 0.03;
A0_ltv = [0.95 0.1; -0.1 0.85];
dA_ltv = [-0.4 0.05; 0.05 -0.3];
B_ltv = [0.5; 0.3];
Xltv = zeros(N+1, p, L); Ultv = randn(N, q, L);
for l = 1:L
    Xltv(1, :, l) = 0.1*randn(1, p);
    for k = 1:N
        Ak = A0_ltv + (k/N) * dA_ltv;
        Xltv(k+1, :, l) = (Ak * Xltv(k, :, l)' + B_ltv * Ultv(k, :, l)')' + sigma*randn(1, p);
    end
end

grid_ltv = logspace(0, 10, 20);
[~, bestLam_ltv, ~] = sidLTVdiscTune(Xltv, Ultv, ...
    'Method', 'frequency', 'LambdaGrid', grid_ltv, 'SegmentLength', 25);

% LTV: lambda should be moderate (not at extremes)
assert(bestLam_ltv > grid_ltv(1), 'LTV: lambda should be > smallest candidate');
assert(bestLam_ltv < grid_ltv(end), 'LTV: lambda should be < largest candidate');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: LTV system selects moderate lambda (%.2e).\n', bestLam_ltv);

%% Test 10: Fallback when threshold is very restrictive
% Use LTV data (from Test 9) where transitions make perfect consistency
% impossible. With threshold=0.9999, the fallback should trigger.
rng(2004);
lastwarn('');  % clear last warning
[~, ~, info_strict] = sidLTVdiscTune(Xltv, Ultv, ...
    'Method', 'frequency', 'LambdaGrid', logspace(1, 6, 5), ...
    'SegmentLength', 25, 'ConsistencyThreshold', 0.9999);

% Should still return a valid result even if threshold not met
assert(isfield(info_strict, 'bestFraction'), 'info should have bestFraction');
assert(info_strict.bestFraction >= 0 && info_strict.bestFraction <= 1, ...
    'bestFraction should be in [0, 1]');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 10 passed: fallback with strict threshold (frac=%.3f).\n', ...
    info_strict.bestFraction);

%% Test 11: Backward compatibility — validation method unchanged
rng(1000);
[~, bestLam_compat, losses_compat] = sidLTVdiscTune(X_train, U_train, X_val, U_val, ...
    'Method', 'validation', 'LambdaGrid', grid);
assert(abs(bestLam_compat - bestLambda) < 1e-12, ...
    'Explicit Method=validation should match default');
assert(max(abs(losses_compat - allLosses)) < 1e-12, ...
    'Validation losses should match exactly');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 11 passed: backward compatibility confirmed.\n');

%% Test 12: Fractions are in valid range
% LTI fractions can legitimately all be 1.0 (system matches spectral data
% at every lambda), so we only check the valid range [0, 1].
assert(any(info_f.fractions > 0), 'At least some fractions should be > 0');
assert(all(info_f.fractions >= 0 & info_f.fractions <= 1), 'All fractions should be in [0, 1]');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 12 passed: fractions are in valid range.\n');

%% Test 13: 2-D single-trajectory train/val accepted (SPEC §8.11, #189)
% Single-trajectory data may be passed as (N+1 x p) / (N x q); the port
% promotes it to L=1 rather than requiring 3-D. Port-symmetric with the Python
% 1-D single-channel test (MATLAB has no 1-D vectors).
rng(1189);
p13 = 2; q13 = 1; N13 = 30;
A13 = [0.9 0.1; -0.1 0.8]; B13 = [0.5; 0.3];
Xtr13 = zeros(N13 + 1, p13); Utr13 = randn(N13, q13); Xtr13(1, :) = randn(1, p13);
Xva13 = zeros(N13 + 1, p13); Uva13 = randn(N13, q13); Xva13(1, :) = randn(1, p13);
for k = 1:N13
    Xtr13(k+1, :) = (A13 * Xtr13(k, :)' + B13 * Utr13(k, :)')' + 0.02 * randn(1, p13);
    Xva13(k+1, :) = (A13 * Xva13(k, :)' + B13 * Uva13(k, :)')' + 0.02 * randn(1, p13);
end
[bestRes13, bestLam13, losses13] = sidLTVdiscTune( ...
    Xtr13, Utr13, Xva13, Uva13, 'LambdaGrid', logspace(-2, 6, 15));
assert(isstruct(bestRes13) && isscalar(bestLam13) && bestLam13 > 0, ...
    '2-D single-trajectory tune should return a valid result');
assert(numel(losses13) == 15 && all(losses13 > 0), ...
    '2-D single-trajectory tune should return per-lambda losses');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 13 passed: 2-D single-trajectory train/val accepted.\n');

fprintf('test_sidLTVdiscTune: %d/%d passed\n', runner__nPassed, runner__nPassed);
