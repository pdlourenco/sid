%% test_sidLTVdiscFrozen - Unit tests for frozen transfer function
%
% Tests sidLTVdiscFrozen for result structure, LTI consistency,
% time-varying behavior, and uncertainty propagation.

fprintf('Running test_sidLTVdiscFrozen...\n');

%% Test 1: Result struct fields and dimensions
rng(2001);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.02 * randn(1, p);
    end
end
ltv = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true);
frz = sidLTVdiscFrozen(ltv);

assert(isfield(frz, 'Frequency'), 'Missing Frequency');
assert(isfield(frz, 'FrequencyHz'), 'Missing FrequencyHz');
assert(isfield(frz, 'TimeSteps'), 'Missing TimeSteps');
assert(isfield(frz, 'Response'), 'Missing Response');
assert(isfield(frz, 'ResponseStd'), 'Missing ResponseStd');
assert(isfield(frz, 'SampleTime'), 'Missing SampleTime');
assert(isfield(frz, 'Method'), 'Missing Method');
assert(strcmp(frz.Method, 'sidLTVdiscFrozen'), 'Method should be sidLTVdiscFrozen');

nf = 128; nk = N;
assert(isequal(size(frz.Frequency), [nf, 1]), 'Frequency should be (128 x 1)');
assert(isequal(size(frz.Response), [nf, p, q, nk]), ...
    'Response should be (nf x p x q x nk)');
assert(isequal(size(frz.ResponseStd), [nf, p, q, nk]), ...
    'ResponseStd should be (nf x p x q x nk)');
fprintf('  Test 1 passed: result struct fields and dimensions.\n');

%% Test 2: Custom frequencies and time steps
w = logspace(-2, log10(pi), 50)';
kVec = [1 10 N];
frz2 = sidLTVdiscFrozen(ltv, 'Frequencies', w, 'TimeSteps', kVec);
assert(isequal(size(frz2.Response), [50, p, q, 3]), 'Custom dims');
assert(isequal(frz2.TimeSteps, kVec'), 'TimeSteps stored');
assert(max(abs(frz2.Frequency - w)) < 1e-15, 'Frequencies stored');
fprintf('  Test 2 passed: custom frequencies and time steps.\n');

%% Test 3: LTI frozen TF matches analytic
rng(2002);
p = 2; q = 1; N = 30; L = 10;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')';
    end
end
ltv = sidLTVdisc(X, U, 'Lambda', 1e8);
w = linspace(0.01, pi, 64)';
frz = sidLTVdiscFrozen(ltv, 'Frequencies', w);

% Compute analytic TF: G(w) = (e^{jw} I - A)^{-1} B
Ip = eye(p);
maxErr = 0;
for iw = 1:length(w)
    z = exp(1i * w(iw));
    G_analytic = (z * Ip - A_true) \ B_true;
    G_frozen_mid = squeeze(frz.Response(iw, :, :, round(N/2)));  % mid time step
    err = norm(G_frozen_mid - G_analytic) / max(norm(G_analytic), eps);
    maxErr = max(maxErr, err);
end
assert(maxErr < 0.02, 'LTI frozen TF should match analytic (maxErr=%.4f)', maxErr);
fprintf('  Test 3 passed: LTI frozen TF matches analytic (maxErr=%.4f).\n', maxErr);

%% Test 4: Time-varying response changes across time steps
rng(2003);
p = 1; q = 1; N = 40; L = 15;
A_seq = 0.5 + 0.4 * (0:N-1)' / (N-1);  % ramp 0.5 -> 0.9
B_true = 1.0;
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn;
    for k = 1:N
        X(k+1, :, l) = A_seq(k) * X(k, :, l) + B_true * U(k, :, l) + 0.01 * randn;
    end
end
ltv = sidLTVdisc(X, U, 'Lambda', 1e2);
w_dc = 0.01;
frz = sidLTVdiscFrozen(ltv, 'Frequencies', w_dc);

% DC gain ~ B/(1-A) for SISO, should increase as A(k) increases
dcGain = abs(squeeze(frz.Response(1, 1, 1, :)));
corrMat = corrcoef(dcGain, A_seq);
assert(corrMat(1, 2) > 0.8, ...
    'DC gain should correlate with A(k) ramp (corr=%.3f)', corrMat(1, 2));
fprintf('  Test 4 passed: time-varying DC gain tracks A(k) ramp (corr=%.3f).\n', corrMat(1, 2));

%% Test 5: Uncertainty propagation produces finite positive std
rng(2004);
p = 2; q = 1; N = 20; L = 5;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + 0.05 * randn(1, p);
    end
end
ltv = sidLTVdisc(X, U, 'Lambda', 1e4, 'Uncertainty', true);
frz = sidLTVdiscFrozen(ltv);

assert(all(frz.ResponseStd(:) > 0), 'ResponseStd should be positive');
assert(all(isfinite(frz.ResponseStd(:))), 'ResponseStd should be finite');
fprintf('  Test 5 passed: uncertainty propagation gives finite positive std.\n');

%% Test 6: No uncertainty -> ResponseStd is empty
ltv_no_unc = sidLTVdisc(X, U, 'Lambda', 1e4);
frz_no = sidLTVdiscFrozen(ltv_no_unc);
assert(isempty(frz_no.ResponseStd), 'ResponseStd should be empty without uncertainty');
fprintf('  Test 6 passed: no uncertainty -> empty ResponseStd.\n');

%% Test 7: SampleTime affects FrequencyHz
frz_ts = sidLTVdiscFrozen(ltv_no_unc, 'SampleTime', 0.01);
w_default = frz_ts.Frequency;
expected_hz = w_default / (2 * pi * 0.01);
assert(max(abs(frz_ts.FrequencyHz - expected_hz)) < 1e-12, ...
    'FrequencyHz should be w / (2*pi*Ts)');
fprintf('  Test 7 passed: SampleTime affects FrequencyHz correctly.\n');

fprintf('test_sidLTVdiscFrozen: ALL TESTS PASSED\n');
