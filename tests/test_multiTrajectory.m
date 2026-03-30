% test_multiTrajectory.m - Test multi-trajectory support for spectral functions
%
% Verifies that sidFreqBT, sidFreqETFE, sidFreqMap, and sidSpectrogram
% correctly handle 3D (N x n_ch x L) multi-trajectory input, performing
% ensemble averaging of covariances/periodograms across trajectories.

fprintf('Running test_multiTrajectory...\n');

%% Test 1: sidFreqBT multi-trajectory produces valid output
rng(3001);
N = 1000; L = 5; a = 0.85;
y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
for l = 1:L
    ul = randn(N, 1);
    yl = filter(1, [1 -a], ul) + 0.1 * randn(N, 1);
    y3(:, 1, l) = yl;
    u3(:, 1, l) = ul;
end

r = sidFreqBT(y3, u3);
assert(isfield(r, 'NumTrajectories'), 'Should have NumTrajectories field');
assert(r.NumTrajectories == L, 'NumTrajectories should be %d, got %d', L, r.NumTrajectories);
assert(all(isfinite(r.Response)), 'Response should be finite');
assert(all(isfinite(r.ResponseStd)), 'ResponseStd should be finite');
assert(all(r.ResponseStd > 0), 'ResponseStd should be positive');
fprintf('  Test 1 passed: sidFreqBT multi-trajectory produces valid output.\n');

%% Test 2: Variance reduction — multi-trajectory std ≈ single-trajectory std / sqrt(L)
rng(3002);
N = 2000; L = 10; a = 0.85;
y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
for l = 1:L
    ul = randn(N, 1);
    yl = filter(1, [1 -a], ul) + 0.2 * randn(N, 1);
    y3(:, 1, l) = yl;
    u3(:, 1, l) = ul;
end

r_multi = sidFreqBT(y3, u3);
r_single = sidFreqBT(y3(:, :, 1), u3(:, :, 1));

% Predicted ratio: std_multi / std_single ≈ 1/sqrt(L)
expected_ratio = 1 / sqrt(L);
actual_ratio = median(r_multi.ResponseStd ./ r_single.ResponseStd);
assert(abs(actual_ratio - expected_ratio) < 0.05, ...
    'Variance reduction: expected ratio %.3f, got %.3f', expected_ratio, actual_ratio);
fprintf('  Test 2 passed: variance reduction ratio = %.3f (expected %.3f).\n', ...
    actual_ratio, expected_ratio);

%% Test 3: sidFreqETFE multi-trajectory produces valid output
rng(3003);
N = 500; L = 3;
y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
for l = 1:L
    ul = randn(N, 1);
    yl = filter(1, [1 -0.7], ul) + 0.1 * randn(N, 1);
    y3(:, 1, l) = yl;
    u3(:, 1, l) = ul;
end

r = sidFreqETFE(y3, u3, 'Smoothing', 5);
assert(r.NumTrajectories == L, 'ETFE NumTrajectories should be %d', L);
assert(all(isfinite(r.Response)), 'ETFE Response should be finite');
fprintf('  Test 3 passed: sidFreqETFE multi-trajectory produces valid output.\n');

%% Test 4: sidFreqMap multi-trajectory produces valid output
rng(3004);
N = 2000; L = 4; a = 0.85;
y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
for l = 1:L
    ul = randn(N, 1);
    yl = filter(1, [1 -a], ul) + 0.1 * randn(N, 1);
    y3(:, 1, l) = yl;
    u3(:, 1, l) = ul;
end

r = sidFreqMap(y3, u3, 'SegmentLength', 500);
assert(r.NumTrajectories == L, 'FreqMap NumTrajectories should be %d', L);
assert(all(isfinite(r.Response(:))), 'FreqMap Response should be finite');
fprintf('  Test 4 passed: sidFreqMap multi-trajectory produces valid output.\n');

%% Test 5: sidSpectrogram multi-trajectory produces valid output
rng(3005);
N = 1000; L = 5;
x3 = randn(N, 1, L);

r = sidSpectrogram(x3, 'WindowLength', 128);
assert(r.NumTrajectories == L, 'Spectrogram NumTrajectories should be %d', L);
assert(all(isfinite(r.Power(:))), 'Spectrogram Power should be finite');
assert(all(r.Power(:) >= 0), 'Spectrogram Power should be non-negative');
fprintf('  Test 5 passed: sidSpectrogram multi-trajectory produces valid output.\n');

%% Test 6: Spectrogram ensemble PSD has lower variance than single trajectory
rng(3006);
N = 2000; L = 8;
x3 = randn(N, 1, L);  % white noise

r_multi = sidSpectrogram(x3, 'WindowLength', 256);
r_single = sidSpectrogram(x3(:, :, 1), 'WindowLength', 256);

% For white noise, PSD should be roughly constant across frequency.
% Multi-trajectory average should have lower variance across frequency.
psd_multi_var = var(r_multi.Power(:, round(end/2), 1));
psd_single_var = var(r_single.Power(:, round(end/2), 1));
assert(psd_multi_var < psd_single_var, ...
    'Ensemble spectrogram should have lower variance (multi=%.4e, single=%.4e)', ...
    psd_multi_var, psd_single_var);
fprintf('  Test 6 passed: ensemble spectrogram has lower PSD variance.\n');

%% Test 7: Backward compatibility — single trajectory gives same result as 2D input
rng(3007);
N = 500;
y = randn(N, 1);
u = randn(N, 1);

r_2d = sidFreqBT(y, u, 'WindowSize', 25);
y3 = reshape(y, N, 1, 1);
u3 = reshape(u, N, 1, 1);
r_3d = sidFreqBT(y3, u3, 'WindowSize', 25);

assert(max(abs(r_2d.Response - r_3d.Response)) < 1e-12, ...
    'Single trajectory 3D should match 2D exactly');
assert(max(abs(r_2d.ResponseStd - r_3d.ResponseStd)) < 1e-12, ...
    'Single trajectory 3D std should match 2D exactly');
assert(r_3d.NumTrajectories == 1, 'NumTrajectories should be 1 for single');
fprintf('  Test 7 passed: backward compatibility (2D = single 3D).\n');

%% Test 8: Time-series (output-only) multi-trajectory
rng(3008);
N = 1000; L = 6;
y3 = zeros(N, 1, L);
for l = 1:L
    y3(:, 1, l) = filter(1, [1 -0.6], randn(N, 1));
end

r = sidFreqBT(y3, []);
assert(r.NumTrajectories == L, 'Time-series NumTrajectories should be %d', L);
assert(isempty(r.Response), 'Time-series Response should be empty');
assert(all(isfinite(r.NoiseSpectrum)), 'Time-series NoiseSpectrum should be finite');
assert(all(r.NoiseSpectrum > 0), 'Time-series NoiseSpectrum should be positive');
fprintf('  Test 8 passed: time-series multi-trajectory works.\n');

fprintf('  test_multiTrajectory: ALL PASSED\n');
