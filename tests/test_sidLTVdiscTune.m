%% test_sidLTVdiscTune - Unit tests for validation-based lambda tuning
%
% Tests sidLTVdiscTune for output shapes, lambda selection quality,
% custom grids, and consistency with sidLTVdisc.

fprintf('Running test_sidLTVdiscTune...\n');

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
        X_train(k+1, :, l) = (Ak * X_train(k, :, l)' + B_true * U_train(k, :, l)')' + sigma * randn(1, p);
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
fprintf('  Test 1 passed: output shapes correct.\n');

%% Test 2: bestResult has correct fields
assert(isfield(bestResult, 'A'), 'bestResult should have A');
assert(isfield(bestResult, 'B'), 'bestResult should have B');
assert(strcmp(bestResult.Method, 'sidLTVdisc'), 'bestResult.Method should be sidLTVdisc');
assert(isequal(size(bestResult.A), [p, p, N]), 'bestResult.A dimensions');
fprintf('  Test 2 passed: bestResult has correct fields.\n');

%% Test 3: Optimal lambda not at grid boundary
% For an LTV problem, the optimum should be interior (bias-variance tradeoff)
[~, bestIdx] = min(allLosses);
assert(bestIdx > 1 && bestIdx < length(grid), ...
    'Optimal lambda should not be at grid boundary (idx=%d/%d)', bestIdx, length(grid));
fprintf('  Test 3 passed: optimal lambda is interior (idx=%d/%d).\n', bestIdx, length(grid));

%% Test 4: Custom lambda grid with few points
grid_small = [1, 100, 10000];
[~, bestLam2, losses2] = sidLTVdiscTune(X_train, U_train, X_val, U_val, ...
    'LambdaGrid', grid_small);
assert(length(losses2) == 3, 'allLosses should have 3 entries');
assert(any(abs(bestLam2 - grid_small) < 1e-12), 'bestLambda should be from grid');
fprintf('  Test 4 passed: custom grid works.\n');

%% Test 5: Consistency with sidLTVdisc
% Re-run sidLTVdisc at the best lambda and verify results match
check = sidLTVdisc(X_train, U_train, 'Lambda', bestLambda);
assert(max(abs(bestResult.A(:) - check.A(:))) < 1e-10, ...
    'bestResult.A should match direct sidLTVdisc call');
assert(max(abs(bestResult.B(:) - check.B(:))) < 1e-10, ...
    'bestResult.B should match direct sidLTVdisc call');
fprintf('  Test 5 passed: bestResult consistent with sidLTVdisc.\n');

%% Test 6: Passthrough of Precondition option
[bestRes_pc, ~, ~] = sidLTVdiscTune(X_train, U_train, X_val, U_val, ...
    'LambdaGrid', [100, 1000, 10000], 'Precondition', true);
assert(bestRes_pc.Preconditioned == true, 'Precondition should be forwarded');
fprintf('  Test 6 passed: Precondition option forwarded.\n');

fprintf('test_sidLTVdiscTune: ALL TESTS PASSED\n');
