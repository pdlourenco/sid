%% test_sidLTVdiscFrozen - Unit tests for frozen transfer function
%
% Tests sidLTVdiscFrozen for result structure, LTI consistency,
% time-varying behavior, and uncertainty propagation.

fprintf('Running test_sidLTVdiscFrozen...\n');
runner__nPassed = 0;

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
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 1 passed: result struct fields and dimensions.\n');

%% Test 2: Custom frequencies and time steps
w = logspace(-2, log10(pi), 50)';
kVec = [1 10 N];
frz2 = sidLTVdiscFrozen(ltv, 'Frequencies', w, 'TimeSteps', kVec);
assert(isequal(size(frz2.Response), [50, p, q, 3]), 'Custom dims');
assert(isequal(frz2.TimeSteps, kVec'), 'TimeSteps stored');
assert(max(abs(frz2.Frequency - w)) < 1e-15, 'Frequencies stored');
runner__nPassed = runner__nPassed + 1;
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
    G_frozen_mid = reshape(frz.Response(iw, :, :, round(N/2)), [p, q]);  % mid time step
    err = norm(G_frozen_mid - G_analytic) / max(norm(G_analytic), eps);
    maxErr = max(maxErr, err);
end
assert(maxErr < 0.02, 'LTI frozen TF should match analytic (maxErr=%.4f)', maxErr);
runner__nPassed = runner__nPassed + 1;
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
runner__nPassed = runner__nPassed + 1;
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
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 5 passed: uncertainty propagation gives finite positive std.\n');

%% Test 6: No uncertainty -> ResponseStd is empty
ltv_no_unc = sidLTVdisc(X, U, 'Lambda', 1e4);
frz_no = sidLTVdiscFrozen(ltv_no_unc);
assert(isempty(frz_no.ResponseStd), 'ResponseStd should be empty without uncertainty');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 6 passed: no uncertainty -> empty ResponseStd.\n');

%% Test 7: SampleTime affects FrequencyHz
frz_ts = sidLTVdiscFrozen(ltv_no_unc, 'SampleTime', 0.01);
w_default = frz_ts.Frequency;
expected_hz = w_default / (2 * pi * 0.01);
assert(max(abs(frz_ts.FrequencyHz - expected_hz)) < 1e-12, ...
    'FrequencyHz should be w / (2*pi*Ts)');
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 7 passed: SampleTime affects FrequencyHz correctly.\n');

%% Test 8: Exact variance formula matches brute-force Kronecker product
% Verify the rank-1 compact formula against a full double sum over
% Cov(vec(C(k))) = Sigma kron P(k) with non-diagonal Sigma and P.
rng(2008);
p = 2; q = 1; N = 15; L = 8;
A_true = [0.9 0.1; -0.1 0.8]; B_true = [0.5; 0.3];
X = zeros(N+1, p, L); U = randn(N, q, L);
for l = 1:L
    X(1, :, l) = randn(1, p);
    for k = 1:N
        noise = [0.02; 0.08] .* randn(p, 1);  % anisotropic noise
        X(k+1, :, l) = (A_true * X(k, :, l)' + B_true * U(k, :, l)')' + noise';
    end
end
ltv = sidLTVdisc(X, U, 'Lambda', 1e3, 'Uncertainty', true, ...
    'CovarianceMode', 'full');
frz = sidLTVdiscFrozen(ltv, 'Frequencies', [0.5; 1.5], 'TimeSteps', [5; 10]);

% Brute-force: Var(G_{ab}) = J^H (Sigma kron P) J for each (a,b,w,k)
d = p + q;
Sigma = ltv.NoiseCov;
maxRelErr = 0;
for ik = 1:2
    ki = [5; 10];
    Pk = ltv.P(:, :, ki(ik));
    Ak = ltv.A(:, :, ki(ik));
    Bk = ltv.B(:, :, ki(ik));
    Ip = eye(p);
    for iw = 1:2
        ww = [0.5; 1.5];
        z = exp(1i * ww(iw));
        R = (z * Ip - Ak) \ Ip;
        Gk = R * Bk;

        for a = 1:p
            for b = 1:q
                % Full Kronecker double sum
                v_bf = 0;
                for r = 1:d
                    for j = 1:p
                        % dG_{ab}/dC_{rj}
                        if r <= p
                            dG_rj = R(a, j) * Gk(r, b);
                        elseif r == p + b
                            dG_rj = R(a, j);
                        else
                            dG_rj = 0;
                        end
                        for s = 1:d
                            for l2 = 1:p
                                if s <= p
                                    dG_sl = R(a, l2) * Gk(s, b);
                                elseif s == p + b
                                    dG_sl = R(a, l2);
                                else
                                    dG_sl = 0;
                                end
                                v_bf = v_bf + conj(dG_rj) * dG_sl ...
                                    * Sigma(j, l2) * Pk(r, s);
                            end
                        end
                    end
                end
                var_bf = real(v_bf);
                var_impl = frz.ResponseStd(iw, a, b, ik)^2;
                relErr = abs(var_bf - var_impl) / max(var_bf, 1e-30);
                maxRelErr = max(maxRelErr, relErr);
            end
        end
    end
end
assert(maxRelErr < 1e-10, ...
    'Exact formula should match brute-force (maxRelErr=%.2e)', maxRelErr);
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 8 passed: exact variance matches brute-force Kronecker (maxRelErr=%.2e).\n', ...
    maxRelErr);

%% Test 9: Frozen TF of an Output-COSMIC (IO) result honors H (#144)
% For an IO result with H != I (py < n) the response must be (nf x py x q x nk)
% and equal H*(zI-A)^-1 B -- not the n x q state response the pre-#144 code gave.
n = 2; q = 1; py = 1; N = 40; L = 6;
A_t = [0.9 0.2; -0.15 0.85];
B_t = [1.0; 0.5];
H9 = [1 0];
rng(909);
U9 = randn(N, q, L); Y9 = zeros(N + 1, py, L); X9 = zeros(N + 1, n, L);
for l = 1:L
    X9(1, :, l) = randn(1, n);
    Y9(1, :, l) = (H9 * X9(1, :, l)')';
    for k = 1:N
        X9(k + 1, :, l) = (A_t * X9(k, :, l)' + B_t * U9(k, :, l)')';
        Y9(k + 1, :, l) = (H9 * X9(k + 1, :, l)')';
    end
end
resIO = sidLTVdiscIO(Y9, U9, H9, 'Lambda', 1e3);
w9 = [0.3; 1.0; 2.5];
frzIO = sidLTVdiscFrozen(resIO, 'TimeSteps', 6, 'Frequencies', w9);

% Per-dimension checks: MATLAB size() drops trailing singleton dims.
assert(size(frzIO.Response, 1) == 3 && size(frzIO.Response, 2) == py && ...
       size(frzIO.Response, 3) == q && size(frzIO.Response, 4) == 1, ...
    'Frozen-of-IO Response must be (nf x py x q x nk), got %s', ...
    mat2str(size(frzIO.Response)));
assert(all(isfinite(frzIO.ResponseStd(:))) && all(frzIO.ResponseStd(:) >= 0), ...
    'Frozen-of-IO ResponseStd must be finite and non-negative');

% Response must equal H*(zI-A)^-1 B (the H fold-in is real, not a no-op).
Ak9 = resIO.A(:, :, 6); Bk9 = resIO.B(:, :, 6);
for iw = 1:3
    z = exp(1i * w9(iw));
    stateResp = (z * eye(n) - Ak9) \ Bk9;      % (n x q)
    expected  = H9 * stateResp;                % (py x q)
    got = frzIO.Response(iw, 1, 1, 1);
    assert(abs(got - expected) < 1e-10, ...
        'Frozen-of-IO must equal H*(zI-A)^-1 B at frequency %d', iw);
end
runner__nPassed = runner__nPassed + 1;
fprintf('  Test 9 passed: frozen-of-IO honors H (py x q output response).\n');

fprintf('test_sidLTVdiscFrozen: %d/%d passed\n', runner__nPassed, runner__nPassed);
