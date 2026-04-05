%% test_sidLTVdisc - Unit tests for discrete LTV state-space identification
%
% Tests sidLTVdisc for result structure, known LTI/LTV recovery,
% multi-trajectory support, lambda options, preconditioning, and
% input validation.

fprintf('Running test_sidLTVdisc...\n');

%% Test 1: Result struct has all required fields
rng(100);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.01 * randn(1, p);
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e3);

requiredFields = {'A', 'B', 'Lambda', 'Cost', 'DataLength', 'StateDim', ...
    'InputDim', 'NumTrajectories', 'Algorithm', 'Preconditioned', 'Method'};
for i = 1:length(requiredFields)
    assert(isfield(result, requiredFields{i}), 'Missing field: %s', requiredFields{i});
end
fprintf('  Test 1 passed: all required fields present.\n');

%% Test 2: Correct metadata values
assert(result.DataLength == N, 'DataLength should be N=%d', N);
assert(result.StateDim == p, 'StateDim should be p=%d', p);
assert(result.InputDim == q, 'InputDim should be q=%d', q);
assert(result.NumTrajectories == L, 'NumTrajectories should be L=%d', L);
assert(strcmp(result.Algorithm, 'cosmic'), 'Algorithm should be cosmic');
assert(result.Preconditioned == false, 'Preconditioned should be false');
assert(strcmp(result.Method, 'sidLTVdisc'), 'Method should be sidLTVdisc');
fprintf('  Test 2 passed: metadata values correct.\n');

%% Test 3: Output dimensions
assert(isequal(size(result.A), [p, p, N]), 'A should be (p x p x N)');
assert(isequal(size(result.B), [p, q, N]), 'B should be (p x q x N)');
assert(isequal(size(result.Lambda), [N-1, 1]), 'Lambda should be (N-1 x 1)');
assert(isequal(size(result.Cost), [1, 3]), 'Cost should be (1 x 3)');
fprintf('  Test 3 passed: output dimensions correct.\n');

%% Test 4: Known LTI system recovery
rng(200);
p = 2; q = 1; N = 50; L = 10;
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
result = sidLTVdisc(X, U, 'Lambda', 1e4);

% Recovered A(k) should be approximately constant and close to A_true
A_mean = mean(result.A, 3);
B_mean = mean(result.B, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
errB = norm(B_mean - B_true, 'fro') / norm(B_true, 'fro');
assert(errA < 0.15, 'LTI A recovery error too large: %.3f', errA);
assert(errB < 0.15, 'LTI B recovery error too large: %.3f', errB);

% Check temporal constancy: std across time steps should be small
A_std = std(reshape(result.A, [], N), 0, 2);
assert(max(A_std) < 0.1, 'A(k) should be approximately constant for LTI system');
fprintf('  Test 4 passed: LTI system recovered (errA=%.4f, errB=%.4f).\n', errA, errB);

%% Test 5: Known LTV system with ramp
rng(300);
p = 1; q = 1; N = 40; L = 15;
A_true_seq = 0.5 + 0.4 * (0:N-1)' / (N-1);  % ramp from 0.5 to 0.9
B_true = 1.0;
sigma = 0.01;
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = A_true_seq(k) * X(k, :, l) + B_true * U(k, :, l) + sigma * randn(1, p);
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e2);

% Check correlation between recovered A(k) and true ramp
A_recovered = squeeze(result.A);
corrMat = corrcoef(A_recovered(:), A_true_seq(:));
assert(corrMat(1,2) > 0.9, 'LTV ramp correlation too low: %.3f', corrMat(1,2));
fprintf('  Test 5 passed: LTV ramp tracked (corr=%.4f).\n', corrMat(1,2));

%% Test 6: Multi-trajectory improves accuracy
rng(400);
p = 2; q = 1; N = 30;
A_true = [0.8 0.05; -0.05 0.7];
B_true = [1; 0.5];
sigma = 0.1;

% Single trajectory
X1 = zeros(N+1, p, 1); U1 = randn(N, q, 1);
X1(1, :, 1) = randn(1, p);
for k = 1:N
    X1(k+1, :, 1) = (A_true * X1(k, :, 1)' + B_true * U1(k, :, 1)')' + sigma * randn(1, p);
end
res1 = sidLTVdisc(X1, U1, 'Lambda', 1e3);

% 10 trajectories
L = 10;
XL = zeros(N+1, p, L); UL = randn(N, q, L);
XL(:, :, 1) = X1; UL(:, :, 1) = U1;
for l = 2:L
    XL(1, :, l) = randn(1, p);
    for k = 1:N
        XL(k+1, :, l) = (A_true * XL(k, :, l)' + B_true * UL(k, :, l)')' + sigma * randn(1, p);
    end
end
resL = sidLTVdisc(XL, UL, 'Lambda', 1e3);

err1 = norm(mean(res1.A, 3) - A_true, 'fro');
errL = norm(mean(resL.A, 3) - A_true, 'fro');
assert(errL < err1, 'Multi-trajectory should improve accuracy');
fprintf('  Test 6 passed: multi-trajectory better (err1=%.4f, errL=%.4f).\n', err1, errL);

%% Test 7: Manual scalar lambda
rng(500);
p = 2; q = 1; N = 15; L = 3;
X = randn(N+1, p, L); U = randn(N, q, L);
result = sidLTVdisc(X, U, 'Lambda', 42.0);
assert(length(result.Lambda) == N-1, 'Lambda should have N-1 elements');
assert(all(abs(result.Lambda - 42.0) < 1e-12), 'All lambda values should be 42.0');
fprintf('  Test 7 passed: scalar lambda expanded correctly.\n');

%% Test 8: Manual per-step lambda vector
rng(501);
N = 15; p = 2; q = 1; L = 3;
X = randn(N+1, p, L); U = randn(N, q, L);
lam = logspace(1, 5, N-1)';
result = sidLTVdisc(X, U, 'Lambda', lam);
assert(max(abs(result.Lambda - lam)) < 1e-10, 'Per-step lambda should be stored exactly');
fprintf('  Test 8 passed: per-step lambda stored correctly.\n');

%% Test 9: L-curve automatic lambda
rng(600);
p = 2; q = 1; N = 30; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end
result = sidLTVdisc(X, U, 'Lambda', 'auto');
assert(all(result.Lambda > 0), 'Auto lambda should be positive');
assert(length(result.Lambda) == N-1, 'Auto lambda should have N-1 elements');
fprintf('  Test 9 passed: L-curve auto lambda works (lambda=%.2e).\n', result.Lambda(1));

%% Test 10: Preconditioning request is gracefully disabled (v1.0)
rng(700);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end
warning('off', 'sid:preconditionUnsupported');
result = sidLTVdisc(X, U, 'Lambda', 1e4, 'Precondition', true);
warning('on', 'sid:preconditionUnsupported');
% Preconditioning is disabled in v1.0 due to off-diagonal block bug
assert(result.Preconditioned == false, 'Preconditioned should be false (disabled in v1.0)');
assert(isequal(size(result.A), [p, p, N]), 'A dimensions should be valid');
fprintf('  Test 10 passed: preconditioning gracefully disabled.\n');

%% Test 11: Cost decomposition: total = fidelity + reg
rng(800);
p = 2; q = 1; N = 20; L = 5;
X = randn(N+1, p, L); U = randn(N, q, L);
result = sidLTVdisc(X, U, 'Lambda', 1e3);
assert(abs(result.Cost(1) - result.Cost(2) - result.Cost(3)) < 1e-10 * abs(result.Cost(1)), ...
    'Cost(1) should equal Cost(2) + Cost(3)');
assert(result.Cost(2) >= 0, 'Data fidelity should be non-negative');
assert(result.Cost(3) >= 0, 'Regularization should be non-negative');
fprintf('  Test 11 passed: cost decomposition correct.\n');

%% Test 12: Input validation - dimension mismatch
try
    sidLTVdisc(randn(10, 2, 3), randn(8, 1, 3), 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for dimension mismatch');
catch e
    assert(strcmp(e.identifier, 'sid:dimMismatch'), ...
        'Expected sid:dimMismatch, got %s', e.identifier);
end

% Trajectory count mismatch
try
    sidLTVdisc(randn(10, 2, 3), randn(9, 1, 2), 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for trajectory mismatch');
catch e
    assert(strcmp(e.identifier, 'sid:dimMismatch'), ...
        'Expected sid:dimMismatch, got %s', e.identifier);
end

% NaN in data
try
    Xbad = randn(10, 2); Xbad(3, 1) = NaN;
    sidLTVdisc(Xbad, randn(9, 1), 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for NaN data');
catch e
    assert(strcmp(e.identifier, 'sid:nonFinite'), ...
        'Expected sid:nonFinite, got %s', e.identifier);
end

% Too short data
try
    sidLTVdisc(randn(2, 2), randn(1, 1), 'Lambda', 1);
    error('sid:testFailed', 'Should have thrown error for too short data');
catch e
    assert(strcmp(e.identifier, 'sid:tooShort'), ...
        'Expected sid:tooShort, got %s', e.identifier);
end

% Negative lambda
try
    sidLTVdisc(randn(10, 2), randn(9, 1), 'Lambda', -1);
    error('sid:testFailed', 'Should have thrown error for negative lambda');
catch e
    assert(strcmp(e.identifier, 'sid:badLambda'), ...
        'Expected sid:badLambda, got %s', e.identifier);
end

% Bad algorithm
try
    sidLTVdisc(randn(10, 2), randn(9, 1), 'Lambda', 1, 'Algorithm', 'foo');
    error('sid:testFailed', 'Should have thrown error for bad algorithm');
catch e
    assert(strcmp(e.identifier, 'sid:badAlgorithm'), ...
        'Expected sid:badAlgorithm, got %s', e.identifier);
end
fprintf('  Test 12 passed: input validation errors correct.\n');

%% Test 13: Noiseless LTI recovery (near-exact)
rng(900);
p = 2; q = 1; N = 30; L = 5;
A_true = [0.85 0.1; -0.05 0.75];
B_true = [0.6; 0.4];
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
    end
end
result = sidLTVdisc(X, U, 'Lambda', 1e6);
A_mean = mean(result.A, 3);
B_mean = mean(result.B, 3);
errA = norm(A_mean - A_true, 'fro') / norm(A_true, 'fro');
errB = norm(B_mean - B_true, 'fro') / norm(B_true, 'fro');
assert(errA < 0.01, 'Noiseless LTI A recovery should be near-exact: %.4f', errA);
assert(errB < 0.01, 'Noiseless LTI B recovery should be near-exact: %.4f', errB);
fprintf(['  Test 13 passed: noiseless LTI near-exact recovery' ...
    ' (errA=%.6f, errB=%.6f).\n'], errA, errB);

%% Test 14: Cost decomposition — fidelity + reg = total (SPEC §8.3.3)
% Verify Cost(1) = Cost(2) + Cost(3) across different lambda regimes
rng(1400);
p = 2; q = 1; N = 25; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' ...
            + 0.05 * randn(1, p);
    end
end

for lam = [1e0, 1e3, 1e6, 1e9]
    r = sidLTVdisc(X, U, 'Lambda', lam);
    total = r.Cost(1);
    fidelity = r.Cost(2);
    reg = r.Cost(3);
    assert(abs(total - (fidelity + reg)) < 1e-8 * max(total, 1), ...
        'Cost decomposition: %.6f != %.6f + %.6f at lambda=%g', ...
        total, fidelity, reg, lam);
end
fprintf('  Test 14 passed: cost decomposition across lambda regimes.\n');

%% Test 15: High lambda → nearly constant A, B (maximum regularisation)
r_high = sidLTVdisc(X, U, 'Lambda', 1e12);
A_std = std(r_high.A, 0, 3);
A_var = max(A_std(:));
B_std = std(r_high.B, 0, 3);
B_var = max(B_std(:));
assert(A_var < 1e-6, ...
    'Very high lambda: A should be near-constant, var=%.2e', A_var);
assert(B_var < 1e-6, ...
    'Very high lambda: B should be near-constant, var=%.2e', B_var);
fprintf('  Test 15 passed: high lambda gives constant A,B (var=%.2e).\n', ...
    A_var);

%% Test 16: Low lambda → more temporal variation (SPEC §8.3.3)
r_low = sidLTVdisc(X, U, 'Lambda', 1e-3);
% Low lambda should allow more temporal variation than high lambda
A_std_low = std(r_low.A, 0, 3);
A_std_hi  = std(r_high.A, 0, 3);
assert(max(A_std_low(:)) > max(A_std_hi(:)), ...
    'Low lambda A should vary more than high lambda');
% Low lambda: fidelity should be at least as good (lower or equal)
assert(r_low.Cost(2) <= r_high.Cost(2) * 1.01, ...
    'Low lambda fidelity %.6f should be <= high lambda %.6f', ...
    r_low.Cost(2), r_high.Cost(2));
fprintf('  Test 16 passed: low vs high lambda variation.\n');

%% Test 17: Preconditioning disabled in v1.0, results match unpreconditioned
rng(1700);
p = 2; q = 1; N = 30; L = 5;
A_true = [0.9 0.1; -0.1 0.8];
B_true = [0.5; 0.3];
X = zeros(N+1, p, L);
U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' ...
            + 0.05 * randn(1, p);
    end
end

r_nopre = sidLTVdisc(X, U, 'Lambda', 1e4, 'Precondition', false);
warning('off', 'sid:preconditionUnsupported');
r_pre   = sidLTVdisc(X, U, 'Lambda', 1e4, 'Precondition', true);
warning('on', 'sid:preconditionUnsupported');

% Preconditioning is disabled in v1.0; both should report false
assert(r_nopre.Preconditioned == false, 'Should report Preconditioned=false');
assert(r_pre.Preconditioned == false, 'Should report Preconditioned=false (disabled in v1.0)');

% Since preconditioning is disabled, results should be identical
assert(max(abs(r_pre.A(:) - r_nopre.A(:))) < 1e-14, ...
    'Disabled precondition should match unpreconditioned');
assert(all(isfinite(r_pre.A(:))), 'A should be finite');
assert(all(isfinite(r_pre.B(:))), 'B should be finite');
fprintf('  Test 17 passed: preconditioning disabled, results match.\n');

%% Test 18: Uncertainty fields present when requested (SPEC §8.9)
r_unc = sidLTVdisc(X, U, 'Lambda', 1e3, 'Uncertainty', true);
assert(isfield(r_unc, 'AStd'), 'Should have AStd field');
assert(isfield(r_unc, 'BStd'), 'Should have BStd field');
assert(isfield(r_unc, 'P'), 'Should have P field');
assert(isequal(size(r_unc.AStd), [p, p, N]), 'AStd should be (p x p x N)');
assert(isequal(size(r_unc.BStd), [p, q, N]), 'BStd should be (p x q x N)');
assert(all(r_unc.AStd(:) >= 0), 'AStd should be non-negative');
assert(all(r_unc.BStd(:) >= 0), 'BStd should be non-negative');
fprintf('  Test 18 passed: uncertainty fields present and valid.\n');

%% Test 19: NaN in input data rejected (SPEC §8.4)
rng(1900);
p = 2; q = 1; N = 20;
X_nan = randn(N+1, p);
U_nan = randn(N, q);
X_nan(5, 1) = NaN;
passed = false;
try
    sidLTVdisc(X_nan, U_nan, 'Lambda', 1e3);
catch e
    if ~isempty(strfind(e.identifier, 'sid:'))
        passed = true;
    end
end
assert(passed, 'NaN in X should raise error');

% Inf in U
U_inf = randn(N, q);
U_inf(3) = Inf;
X_ok = randn(N+1, p);
passed = false;
try
    sidLTVdisc(X_ok, U_inf, 'Lambda', 1e3);
catch e
    if ~isempty(strfind(e.identifier, 'sid:'))
        passed = true;
    end
end
assert(passed, 'Inf in U should raise error');
fprintf('  Test 19 passed: NaN/Inf in data rejected.\n');

%% Test 20: Very short data (N=2, minimum viable)
rng(2000);
p = 2; q = 1; N = 2; L = 3;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
    end
end

result_short = sidLTVdisc(X, U, 'Lambda', 1e3);
assert(result_short.DataLength == 2, 'DataLength should be 2');
assert(isequal(size(result_short.A), [p, p, N]), 'A should be (p x p x 2)');
fprintf('  Test 20 passed: N=2 minimum viable data.\n');

fprintf('test_sidLTVdisc: ALL TESTS PASSED\n');
