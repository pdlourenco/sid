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

%% Test 10: Preconditioning runs without error
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
result = sidLTVdisc(X, U, 'Lambda', 1e4, 'Precondition', true);
assert(result.Preconditioned == true, 'Preconditioned should be true');
assert(isequal(size(result.A), [p, p, N]), 'Preconditioned A dimensions');
fprintf('  Test 10 passed: preconditioning runs correctly.\n');

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
fprintf('  Test 13 passed: noiseless LTI near-exact recovery (errA=%.6f, errB=%.6f).\n', errA, errB);

fprintf('test_sidLTVdisc: ALL TESTS PASSED\n');
