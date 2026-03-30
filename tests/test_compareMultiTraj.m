% test_compareMultiTraj.m - Compare multi-trajectory results vs MathWorks toolbox
%
% Requires: System Identification Toolbox (spa, etfe, spafdr, iddata, merge)
%
% Tests that sid multi-trajectory (3D array) results match MATLAB's
% spa/etfe/spafdr with merged multi-experiment iddata objects.

fprintf('Running test_compareMultiTraj...\n');

%% Toolbox check
if ~exist('spa', 'file') || ~exist('iddata', 'file')
    fprintf('  test_compareMultiTraj: SKIPPED (requires System Identification Toolbox)\n');
    return;
end

%% Test 1: sidFreqBT vs spa with merged iddata (SISO, L=5)
rng(4001);
N = 2000; L = 5; Ts = 1; M = 30;
a = 0.85;
y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
dList = cell(1, L);
for l = 1:L
    ul = randn(N, 1);
    yl = filter(1, [1 -a], ul) + 0.3 * randn(N, 1);
    y3(:, 1, l) = yl;
    u3(:, 1, l) = ul;
    dList{l} = iddata(yl, ul, Ts);
end

% sid: multi-trajectory
r_sid = sidFreqBT(y3, u3, 'WindowSize', M);
w = r_sid.Frequency;

% MathWorks: merge experiments
dMerged = merge(dList{:});
G_spa = spa(dMerged, M, w);
resp_spa = squeeze(G_spa.ResponseData);

% Compare magnitude
relErr_mag = max(abs(abs(r_sid.Response) - abs(resp_spa)) ./ max(abs(resp_spa), 1e-10));
assert(relErr_mag < 0.02, ...
    'Test 1: SISO multi-traj magnitude relErr=%.4f should be <2%%', relErr_mag);

% Compare phase
phaseErr = max(abs(angle(r_sid.Response) - angle(resp_spa)));
assert(phaseErr < 0.02, ...
    'Test 1: SISO multi-traj phase error=%.4f rad should be <0.02', phaseErr);

fprintf('  Test 1 passed: sidFreqBT matches spa(merge(...)) for L=%d (mag err=%.4f).\n', L, relErr_mag);

%% Test 2: sidFreqBT time-series vs spa merged (L=5)
rng(4002);
N = 2000; L = 5; Ts = 1; M = 30;
y3 = zeros(N, 1, L);
dList = cell(1, L);
for l = 1:L
    yl = filter(1, [1 -0.6], randn(N, 1));
    y3(:, 1, l) = yl;
    dList{l} = iddata(yl, [], Ts);
end

r_sid = sidFreqBT(y3, []);
w = r_sid.Frequency;
dMerged = merge(dList{:});
G_spa = spa(dMerged, M, w);
spec_spa = real(squeeze(G_spa.SpectrumData));

relErr = max(abs(real(r_sid.NoiseSpectrum) - spec_spa) ./ max(abs(spec_spa), 1e-10));
assert(relErr < 0.05, ...
    'Test 2: time-series multi-traj spectrum relErr=%.4f should be <5%%', relErr);
fprintf('  Test 2 passed: time-series multi-traj spectrum matches spa(merge(...)).\n');

%% Test 3: sidFreqETFE vs etfe with merged iddata (L=3)
rng(4003);
N = 1000; L = 3; Ts = 1;
y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
dList = cell(1, L);
for l = 1:L
    ul = randn(N, 1);
    yl = filter([1 0.3], [1 -0.8], ul) + 0.1 * randn(N, 1);
    y3(:, 1, l) = yl;
    u3(:, 1, l) = ul;
    dList{l} = iddata(yl, ul, Ts);
end

r_sid = sidFreqETFE(y3, u3, 'Smoothing', 11);
dMerged = merge(dList{:});

% MathWorks etfe with smoothing (etfe does not accept a frequency vector;
% use default frequencies and compare at those)
G_etfe = etfe(dMerged, 11);
% G_etfe.Frequency is in rad/TimeUnit. With Ts=1, this equals rad/sample.
w_etfe = G_etfe.Frequency(:);
% Clip to (0, pi] (etfe may include DC or frequencies above pi)
valid = w_etfe > 0 & w_etfe <= pi;
w_etfe = w_etfe(valid);
resp_etfe = squeeze(G_etfe.ResponseData);
resp_etfe = resp_etfe(valid);

% Re-run sid at the same frequencies for a fair comparison
r_sid = sidFreqETFE(y3, u3, 'Smoothing', 11, 'Frequencies', w_etfe);

relErr = max(abs(abs(r_sid.Response) - abs(resp_etfe)) ./ max(abs(resp_etfe), 1e-10));
% ETFE has high variance and MATLAB's merge averaging may differ from our
% DFT-sum approach. Use generous tolerance (single-traj test_compareEtfe
% also uses loose tolerances for ETFE).
assert(relErr < 0.50, ...
    'Test 3: ETFE multi-traj relErr=%.4f should be <50%%', relErr);
fprintf('  Test 3 passed: sidFreqETFE matches etfe(merge(...)) (err=%.4f).\n', relErr);

%% Test 4: sidFreqBTFDR vs spafdr with merged iddata (L=3)
rng(4004);
if ~exist('spafdr', 'file')
    fprintf('  Test 4: SKIPPED (spafdr not available)\n');
else
    N = 2000; L = 3; Ts = 1;
    y3 = zeros(N, 1, L); u3 = zeros(N, 1, L);
    dList = cell(1, L);
    for l = 1:L
        ul = randn(N, 1);
        yl = filter(1, [1 -0.9], ul) + 0.1 * randn(N, 1);
        y3(:, 1, l) = yl;
        u3(:, 1, l) = ul;
        dList{l} = iddata(yl, ul, Ts);
    end

    r_sid = sidFreqBTFDR(y3, u3);
    w = r_sid.Frequency;
    dMerged = merge(dList{:});
    G_spafdr = spafdr(dMerged, [], w);
    resp_spafdr = squeeze(G_spafdr.ResponseData);

    % Use median (not max) like existing test_compareSpafdr — adaptive window
    % sizes may differ between implementations, especially with multi-traj data
    relErr = median(abs(abs(r_sid.Response) - abs(resp_spafdr)) ./ max(abs(resp_spafdr), 1e-10));
    assert(relErr < 0.10, ...
        'Test 4: BTFDR multi-traj median relErr=%.4f should be <10%%', relErr);
    fprintf('  Test 4 passed: sidFreqBTFDR matches spafdr(merge(...)) (median err=%.4f).\n', relErr);
end

%% Test 5: MIMO multi-trajectory (2 outputs, 1 input, L=4)
rng(4005);
N = 3000; L = 4; Ts = 1; M = 30;
y3 = zeros(N, 2, L); u3 = zeros(N, 1, L);
dList = cell(1, L);
for l = 1:L
    ul = randn(N, 1);
    y1l = filter(1, [1 -0.5], ul) + 0.1 * randn(N, 1);
    y2l = filter(0.3, [1 -0.7], ul) + 0.1 * randn(N, 1);
    y3(:, :, l) = [y1l, y2l];
    u3(:, 1, l) = ul;
    dList{l} = iddata([y1l, y2l], ul, Ts);
end

r_sid = sidFreqBT(y3, u3, 'WindowSize', M);
w = r_sid.Frequency;
dMerged = merge(dList{:});
G_spa = spa(dMerged, M, w);
resp_spa = permute(G_spa.ResponseData, [3, 1, 2]);

for ch = 1:2
    relErr = max(abs(abs(r_sid.Response(:, ch)) - abs(resp_spa(:, ch))) ./ max(abs(resp_spa(:, ch)), 1e-10));
    assert(relErr < 0.02, ...
        'Test 5: MIMO ch%d multi-traj relErr=%.4f should be <2%%', ch, relErr);
end
fprintf('  Test 5 passed: MIMO multi-trajectory matches spa(merge(...)).\n');

fprintf('  test_compareMultiTraj: ALL PASSED\n');
