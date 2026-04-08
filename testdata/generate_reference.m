function generate_reference()
%GENERATE_REFERENCE Generate cross-language reference data.
%
%   generate_reference
%
%   Produces JSON files with canonical test vectors for validating numerical
%   equivalence across MATLAB, Python, and Julia implementations.
%
%   Usage:
%     run('testdata/generate_reference.m')

fprintf('=== Generating cross-language reference data ===\n\n');

% Add paths
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
sidDir = fullfile(rootDir, 'matlab', 'sid');
addpath(sidDir);
% MATLAB ignores addpath on directories named 'private'. Copy to a
% temporary non-private-named directory so private helpers are accessible.
privateDir = fullfile(sidDir, 'private');
shimDir = fullfile(thisDir, 'private_shim');
if exist(shimDir, 'dir')
    rmdir(shimDir, 's');
end
mkdir(shimDir);
copyfile(fullfile(privateDir, '*.m'), shimDir);
addpath(shimDir);
cleanupObj = onCleanup(@() cleanupShim(shimDir));

% ---- Test case 1: SISO Blackman-Tukey ----
fprintf('Generating reference_siso_bt.json...\n');
rng(42);
N = 500;
u1 = randn(N, 1);
y1 = filter([1], [1 -0.9], u1) + 0.1 * randn(N, 1);
r1 = sidFreqBT(y1, u1, 'WindowSize', 30);

ref1 = struct();
ref1.function_name = 'sidFreqBT';
ref1.params = struct('WindowSize', 30, 'SampleTime', 1.0);
ref1.input = struct('y', y1, 'u', u1);
ref1.output = struct( ...
    'Frequency', r1.Frequency, ...
    'Response_real', real(r1.Response), ...
    'Response_imag', imag(r1.Response), ...
    'NoiseSpectrum', r1.NoiseSpectrum, ...
    'Coherence', r1.Coherence);
ref1.tolerance = struct('Response_rel', 1e-10, 'NoiseSpectrum_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_siso_bt.json'), ref1);

% ---- Test case 1b: SISO Blackman-Tukey (time series mode) ----
fprintf('Generating reference_timeseries_bt.json...\n');
rng(42);
N_ts = 500;
y_ts = filter([1], [1 -0.9 0.2], randn(N_ts, 1));
r1b = sidFreqBT(y_ts, [], 'WindowSize', 30);

ref1b = struct();
ref1b.function_name = 'sidFreqBT';
ref1b.params = struct('WindowSize', 30, 'SampleTime', 1.0, 'mode', 'timeseries');
ref1b.input = struct('y', y_ts);
ref1b.output = struct( ...
    'Frequency', r1b.Frequency, ...
    'NoiseSpectrum', r1b.NoiseSpectrum, ...
    'NoiseSpectrumStd', r1b.NoiseSpectrumStd);
ref1b.tolerance = struct('NoiseSpectrum_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_timeseries_bt.json'), ref1b);

% ---- Test case 1c: Internal helpers (cov, windowed DFT, uncertainty) ----
fprintf('Generating reference_internals.json...\n');
rng(42);
N_int = 100;
x_int = randn(N_int, 1);
z_int = randn(N_int, 1);
M_int = 20;

% Biased covariance
R_xx = sidCov(x_int, x_int, M_int);
R_xz = sidCov(x_int, z_int, M_int);

% Hann window
W_int = sidHannWin(M_int);

% Windowed DFT (auto-spectrum, default freqs)
freqs_int = (1:128)' * pi / 128;
Phi_xx = sidWindowedDFT(R_xx, W_int, freqs_int, true, R_xx);

% Windowed DFT (cross-spectrum)
R_zx = sidCov(z_int, x_int, M_int);
Phi_xz = sidWindowedDFT(R_xz, W_int, freqs_int, true, R_zx);

% DFT of time-domain signal
X_dft = sidDFT(x_int, freqs_int, true);

ref_int = struct();
ref_int.function_name = 'internals';
ref_int.params = struct();
ref_int.input = struct('x', x_int, 'z', z_int, 'M', M_int);
ref_int.output = struct( ...
    'R_xx', R_xx, ...
    'R_xz', R_xz, ...
    'W', W_int, ...
    'Phi_xx_real', real(Phi_xx), ...
    'Phi_xx_imag', imag(Phi_xx), ...
    'Phi_xz_real', real(Phi_xz), ...
    'Phi_xz_imag', imag(Phi_xz), ...
    'DFT_real', real(X_dft), ...
    'DFT_imag', imag(X_dft));
ref_int.tolerance = struct( ...
    'R_xx_rel', 1e-12, 'R_xz_rel', 1e-12, ...
    'W_rel', 1e-15, ...
    'Phi_xx_real_rel', 1e-10, ...
    'Phi_xx_imag_rel', 1e-10, 'Phi_xx_imag_atol', 1e-14, ...
    'Phi_xz_real_rel', 1e-10, ...
    'Phi_xz_imag_rel', 1e-10, ...
    'DFT_real_rel', 1e-10, ...
    'DFT_imag_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_internals.json'), ref_int);

% ---- Test case 2: MIMO Blackman-Tukey ----
fprintf('Generating reference_mimo_bt.json...\n');
rng(43);
N = 1000;
u2 = randn(N, 2);
y2 = [filter([1 0.5], [1 -0.8], u2(:,1)) + 0.05 * randn(N, 1), ...
      filter([0.3], [1 -0.7 0.2], u2(:,2)) + 0.05 * randn(N, 1)];
r2 = sidFreqBT(y2, u2);

ref2 = struct();
ref2.function_name = 'sidFreqBT';
ref2.params = struct('WindowSize', r2.WindowSize, 'SampleTime', 1.0);
ref2.input = struct('y', y2, 'u', u2);
ref2.output = struct( ...
    'Frequency', r2.Frequency, ...
    'Response_real', real(r2.Response), ...
    'Response_imag', imag(r2.Response), ...
    'NoiseSpectrum', r2.NoiseSpectrum);
ref2.tolerance = struct('Response_rel', 1e-10, 'Response_atol', 1e-14, ...
                        'NoiseSpectrum_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_mimo_bt.json'), ref2);

% ---- Test case 3: SISO ETFE ----
fprintf('Generating reference_siso_etfe.json...\n');
rng(44);
N = 500;
u3 = randn(N, 1);
y3 = filter([1], [1 -0.85], u3) + 0.1 * randn(N, 1);
r3 = sidFreqETFE(y3, u3, 'Smoothing', 5);

ref3 = struct();
ref3.function_name = 'sidFreqETFE';
ref3.params = struct('Smoothing', 5, 'SampleTime', 1.0);
ref3.input = struct('y', y3, 'u', u3);
ref3.output = struct( ...
    'Frequency', r3.Frequency, ...
    'Response_real', real(r3.Response), ...
    'Response_imag', imag(r3.Response), ...
    'NoiseSpectrum', r3.NoiseSpectrum);
ref3.tolerance = struct('Response_rel', 1e-10, 'NoiseSpectrum_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_siso_etfe.json'), ref3);

% ---- Test case 4: SISO BTFDR ----
fprintf('Generating reference_siso_btfdr.json...\n');
rng(45);
N = 500;
u4 = randn(N, 1);
y4 = filter([1], [1 -0.9], u4) + 0.1 * randn(N, 1);
r4 = sidFreqBTFDR(y4, u4);

ref4 = struct();
ref4.function_name = 'sidFreqBTFDR';
ref4.params = struct('SampleTime', 1.0);
ref4.input = struct('y', y4, 'u', u4);
ref4.output = struct( ...
    'Frequency', r4.Frequency, ...
    'Response_real', real(r4.Response), ...
    'Response_imag', imag(r4.Response), ...
    'NoiseSpectrum', r4.NoiseSpectrum);
ref4.tolerance = struct('Response_rel', 1e-10, 'NoiseSpectrum_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_siso_btfdr.json'), ref4);

% ---- Test case 4b: SISO Spectrogram ----
fprintf('Generating reference_spectrogram.json...\n');
rng(47);
N_sp = 2000;
Ts_sp = 0.001;
x_sp = randn(N_sp, 1);
r_sp = sidSpectrogram(x_sp, 'WindowLength', 128, 'Overlap', 64, 'SampleTime', Ts_sp);

ref_sp = struct();
ref_sp.function_name = 'sidSpectrogram';
ref_sp.params = struct('WindowLength', 128, 'Overlap', 64, 'SampleTime', Ts_sp);
ref_sp.input = struct('x', x_sp);
ref_sp.output = struct( ...
    'Time', r_sp.Time, ...
    'Frequency', r_sp.Frequency, ...
    'Power', r_sp.Power);
ref_sp.tolerance = struct('Time_rel', 1e-12, 'Power_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_spectrogram.json'), ref_sp);

% ---- Test case 4c: SISO FreqMap (BT) ----
fprintf('Generating reference_freqmap_bt.json...\n');
rng(48);
N_fm = 4000;
u_fm = randn(N_fm, 1);
y_fm = filter([1], [1 -0.9], u_fm) + 0.1 * randn(N_fm, 1);
r_fm = sidFreqMap(y_fm, u_fm, 'SegmentLength', 512, 'Overlap', 256, ...
                               'WindowSize', 25);

ref_fm = struct();
ref_fm.function_name = 'sidFreqMap';
ref_fm.params = struct('SegmentLength', 512, 'Overlap', 256, 'WindowSize', 25, ...
                       'SampleTime', 1.0, 'Algorithm', 'bt');
ref_fm.input = struct('y', y_fm, 'u', u_fm);
ref_fm.output = struct( ...
    'Time', r_fm.Time, ...
    'Frequency', r_fm.Frequency, ...
    'Response_real', real(r_fm.Response), ...
    'Response_imag', imag(r_fm.Response), ...
    'NoiseSpectrum', r_fm.NoiseSpectrum, ...
    'Coherence', r_fm.Coherence);
ref_fm.tolerance = struct('Response_rel', 1e-10, 'NoiseSpectrum_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_freqmap_bt.json'), ref_fm);

% ---- Test case 5: LTV COSMIC ----
fprintf('Generating reference_ltv_cosmic.json...\n');
rng(46);
N = 50; p = 2; q = 1;
A_true = 0.95 * eye(p);
B_true = [1; 0.5];
X = zeros(N + 1, p);
U = randn(N, q);
X(1, :) = randn(1, p);
for k = 1:N
    X(k+1, :) = (A_true * X(k, :)' + B_true * U(k, :)')' ...
                + 0.01 * randn(1, p);
end
r5 = sidLTVdisc(X, U, 'Lambda', 1e5);

ref5 = struct();
ref5.function_name = 'sidLTVdisc';
ref5.params = struct('Lambda', 1e5, 'Precondition', false);
ref5.input = struct('X', X, 'U', U);
ref5.output = struct( ...
    'A', r5.A, ...
    'B', r5.B, ...
    'Cost', r5.Cost);
ref5.tolerance = struct('A_rel', 1e-6, 'B_rel', 1e-6, 'Cost_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_ltv_cosmic.json'), ref5);

% ---- Test case 6: sidDetrend ----
fprintf('Generating reference_detrend.json...\n');
rng(50);
N_dt = 200;
t_dt = (0:N_dt-1)';
x_dt = 3.0 + 0.02 * t_dt + 0.5 * randn(N_dt, 1);
[x_detrended, trend] = sidDetrend(x_dt, 'Order', 1);

ref6 = struct();
ref6.function_name = 'sidDetrend';
ref6.params = struct('Order', 1);
ref6.input = struct('x', x_dt);
ref6.output = struct('x_detrended', x_detrended, 'trend', trend);
ref6.tolerance = struct('x_detrended_rel', 1e-10, 'trend_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_detrend.json'), ref6);

% ---- Test case 7: sidModelOrder ----
fprintf('Generating reference_model_order.json...\n');
rng(51);
N_mo = 500;
u_mo = randn(N_mo, 1);
y_mo = filter([1 0.5], [1 -0.8 0.2], u_mo) + 0.05 * randn(N_mo, 1);
r_mo_bt = sidFreqBT(y_mo, u_mo, 'WindowSize', 40);
[n_mo, sv_mo] = sidModelOrder(r_mo_bt, 'Horizon', 30);

ref7 = struct();
ref7.function_name = 'sidModelOrder';
ref7.params = struct('Horizon', 30, 'bt_WindowSize', 40);
ref7.input = struct('y', y_mo, 'u', u_mo);
ref7.output = struct('n', n_mo, 'SingularValues', sv_mo.SingularValues);
ref7.tolerance = struct('SingularValues_rel', 1e-8);

writeJSON(fullfile(thisDir, 'reference_model_order.json'), ref7);

% ---- Test case 8: sidCompare (state-space model) ----
fprintf('Generating reference_compare.json...\n');
% Reuse the LTV COSMIC data from test case 5
rng(46);
N_cmp = 50; p_cmp = 2; q_cmp = 1;
A_cmp = 0.95 * eye(p_cmp);
B_cmp = [1; 0.5];
X_cmp = zeros(N_cmp + 1, p_cmp);
U_cmp = randn(N_cmp, q_cmp);
X_cmp(1, :) = randn(1, p_cmp);
for k = 1:N_cmp
    X_cmp(k+1, :) = (A_cmp * X_cmp(k, :)' + B_cmp * U_cmp(k, :)')' ...
                + 0.01 * randn(1, p_cmp);
end
r_cmp = sidLTVdisc(X_cmp, U_cmp, 'Lambda', 1e5);
comp = sidCompare(r_cmp, X_cmp, U_cmp);

ref8 = struct();
ref8.function_name = 'sidCompare';
ref8.params = struct('Lambda', 1e5, 'Precondition', false);
ref8.input = struct('X', X_cmp, 'U', U_cmp);
ref8.output = struct('Predicted', comp.Predicted, 'Fit', comp.Fit);
ref8.tolerance = struct('Predicted_rel', 1e-6, 'Fit_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_compare.json'), ref8);

% ---- Test case 9: sidResidual (frequency-domain model) ----
fprintf('Generating reference_residual.json...\n');
rng(52);
N_res = 500;
u_res = randn(N_res, 1);
y_res = filter([1], [1 -0.85], u_res) + 0.1 * randn(N_res, 1);
r_res_bt = sidFreqBT(y_res, u_res, 'WindowSize', 30);
res = sidResidual(r_res_bt, y_res, u_res, 'MaxLag', 20, 'Plot', false);

ref9 = struct();
ref9.function_name = 'sidResidual';
ref9.params = struct('MaxLag', 20, 'bt_WindowSize', 30);
ref9.input = struct('y', y_res, 'u', u_res);
ref9.output = struct( ...
    'Residual', res.Residual, ...
    'AutoCorr', res.AutoCorr, ...
    'CrossCorr', res.CrossCorr);
ref9.tolerance = struct( ...
    'Residual_rel', 1e-6, ...
    'AutoCorr_rel', 1e-6, ...
    'CrossCorr_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_residual.json'), ref9);

% ---- Test case 10: sidFreqDomainSim ----
fprintf('Generating reference_freq_domain_sim.json...\n');
rng(53);
N_sim = 300;
u_sim = randn(N_sim, 1);
y_sim = filter([1], [1 -0.9], u_sim);
r_sim_bt = sidFreqBT(y_sim, u_sim, 'WindowSize', 25);
Y_pred = sidFreqDomainSim(r_sim_bt.Response, r_sim_bt.Frequency, u_sim, N_sim);

ref10 = struct();
ref10.function_name = 'sidFreqDomainSim';
ref10.params = struct('bt_WindowSize', 25);
ref10.input = struct('y_noiseless', y_sim, 'u', u_sim);
ref10.output = struct('Y_pred', Y_pred);
ref10.tolerance = struct('Y_pred_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_freq_domain_sim.json'), ref10);

% ---- Test case 11: sidUncertainty ----
fprintf('Generating reference_uncertainty.json...\n');
rng(54);
N_unc = 500;
M_unc = 30;
u_unc = randn(N_unc, 1);
y_unc = filter([1], [1 -0.9], u_unc) + 0.1 * randn(N_unc, 1);
r_unc = sidFreqBT(y_unc, u_unc, 'WindowSize', M_unc);
W_unc = sidHannWin(M_unc);
[GStd, PhiVStd] = sidUncertainty(r_unc.Response, r_unc.NoiseSpectrum, ...
                                  r_unc.Coherence, N_unc, W_unc, 1);

ref11 = struct();
ref11.function_name = 'sidUncertainty';
ref11.params = struct('bt_WindowSize', M_unc);
ref11.input = struct('y', y_unc, 'u', u_unc);
ref11.output = struct('GStd', GStd, 'PhiVStd', PhiVStd);
ref11.tolerance = struct('GStd_rel', 1e-8, 'PhiVStd_rel', 1e-8);

writeJSON(fullfile(thisDir, 'reference_uncertainty.json'), ref11);

% ---- Test case 12: sidLTVdiscFrozen ----
fprintf('Generating reference_ltv_frozen.json...\n');
rng(46);
N_fr = 50; p_fr = 2; q_fr = 1;
A_fr_true = 0.95 * eye(p_fr);
B_fr_true = [1; 0.5];
X_fr = zeros(N_fr + 1, p_fr);
U_fr = randn(N_fr, q_fr);
X_fr(1, :) = randn(1, p_fr);
for k = 1:N_fr
    X_fr(k+1, :) = (A_fr_true * X_fr(k, :)' + B_fr_true * U_fr(k, :)')' ...
                + 0.01 * randn(1, p_fr);
end
r_fr = sidLTVdisc(X_fr, U_fr, 'Lambda', 1e5);
frozen = sidLTVdiscFrozen(r_fr, 'TimeSteps', [1 25 50]);

ref12 = struct();
ref12.function_name = 'sidLTVdiscFrozen';
ref12.params = struct('Lambda', 1e5, 'Precondition', false, ...
                      'frozen_TimeSteps', [1 25 50]);
ref12.input = struct('X', X_fr, 'U', U_fr);
ref12.output = struct( ...
    'Frequency', frozen.Frequency, ...
    'Response_real', real(frozen.Response), ...
    'Response_imag', imag(frozen.Response));
ref12.tolerance = struct('Frequency_rel', 1e-12, 'Response_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_ltv_frozen.json'), ref12);

% ---- Test case 13: COSMIC pipeline internals ----
fprintf('Generating reference_cosmic_internals.json...\n');
rng(55);
N_ci = 30; p_ci = 2; q_ci = 1;
A_ci = 0.9 * eye(p_ci);
B_ci = [1; 0.3];
X_ci = zeros(N_ci + 1, p_ci);
U_ci = randn(N_ci, q_ci);
X_ci(1, :) = randn(1, p_ci);
for k = 1:N_ci
    X_ci(k+1, :) = (A_ci * X_ci(k, :)' + B_ci * U_ci(k, :)')' ...
                + 0.01 * randn(1, p_ci);
end
d_ci = p_ci + q_ci;
lambda_ci = 1e4 * ones(N_ci - 1, 1);

% Build data matrices
[D_ci, Xl_ci] = sidLTVbuildDataMatrices(X_ci, U_ci, N_ci, p_ci, q_ci, 1);

% Build block terms
[S_ci, T_ci] = sidLTVbuildBlockTerms(D_ci, Xl_ci, lambda_ci, N_ci, p_ci, q_ci);

% COSMIC solve
[C_ci, Lbd_ci] = sidLTVcosmicSolve(S_ci, T_ci, lambda_ci, N_ci, p_ci, q_ci);

% Extract A, B from C
A_est_ci = zeros(p_ci, p_ci, N_ci);
B_est_ci = zeros(p_ci, q_ci, N_ci);
for k = 1:N_ci
    A_est_ci(:, :, k) = C_ci(1:p_ci, :, k)';
    B_est_ci(:, :, k) = C_ci(p_ci+1:d_ci, :, k)';
end

% Evaluate cost
[cost_ci, fid_ci, reg_ci] = sidLTVevaluateCost( ...
    A_est_ci, B_est_ci, D_ci, Xl_ci, lambda_ci, N_ci, p_ci, q_ci);

% Uncertainty backward pass
S_scaled_ci = S_ci / N_ci;
P_ci = sidLTVuncertaintyBackwardPass(S_scaled_ci, lambda_ci, N_ci, d_ci);

ref13 = struct();
ref13.function_name = 'cosmic_internals';
ref13.params = struct();
ref13.input = struct('X', X_ci, 'U', U_ci, 'p', p_ci, 'q', q_ci, ...
                     'lambda', lambda_ci(1));
ref13.output = struct( ...
    'D', D_ci, ...
    'Xl', Xl_ci, ...
    'S', S_ci, ...
    'T', T_ci, ...
    'C', C_ci, ...
    'cost', cost_ci, ...
    'fidelity', fid_ci, ...
    'regularization', reg_ci, ...
    'P', P_ci);
ref13.tolerance = struct( ...
    'D_rel', 1e-12, 'Xl_rel', 1e-12, ...
    'S_rel', 1e-10, 'T_rel', 1e-10, ...
    'C_rel', 1e-8, 'C_atol', 1e-10, ...
    'cost_rel', 1e-8, 'fidelity_rel', 1e-8, 'regularization_rel', 1e-8, ...
    'P_rel', 1e-8);

writeJSON(fullfile(thisDir, 'reference_cosmic_internals.json'), ref13);

% ---- Test case 14: sidLTVdiscIO (Output-COSMIC) ----
fprintf('Generating reference_ltv_io.json...\n');
rng(56);
N_io = 40; p_io = 3; q_io = 1; py_io = 2;
A_io = 0.9 * eye(p_io);
B_io = [1; 0.5; 0.2];
H_io = [1 0 0; 0 1 0];
X_io = zeros(N_io + 1, p_io);
U_io = randn(N_io, q_io);
X_io(1, :) = randn(1, p_io);
for k = 1:N_io
    X_io(k+1, :) = (A_io * X_io(k, :)' + B_io * U_io(k, :)')' ...
                    + 0.01 * randn(1, p_io);
end
Y_io = (H_io * X_io')';
r_io = sidLTVdiscIO(Y_io, U_io, H_io, 'Lambda', 1e5);

ref14 = struct();
ref14.function_name = 'sidLTVdiscIO';
ref14.params = struct('Lambda', 1e5);
ref14.input = struct('Y', Y_io, 'U', U_io, 'H', H_io);
ref14.output = struct('A', r_io.A, 'B', r_io.B, 'Cost', r_io.Cost);
ref14.tolerance = struct('A_rel', 1e-2, 'B_rel', 1e-2, 'Cost_rel', 1e-2);

writeJSON(fullfile(thisDir, 'reference_ltv_io.json'), ref14);

% ---- Test case 15: sidLTVStateEst ----
fprintf('Generating reference_ltv_state_est.json...\n');
rng(57);
N_se = 30; p_se = 2; q_se = 1; py_se = 2;
A_se = repmat(0.9 * eye(p_se), [1 1 N_se]);
B_se = repmat([1; 0.5], [1 1 N_se]);
H_se = eye(p_se);
U_se = randn(N_se, q_se);
X_true_se = zeros(N_se + 1, p_se);
X_true_se(1, :) = randn(1, p_se);
for k = 1:N_se
    X_true_se(k+1, :) = (A_se(:,:,k) * X_true_se(k, :)' + ...
                          B_se(:,:,k) * U_se(k, :)')' + 0.01 * randn(1, p_se);
end
Y_se = (H_se * X_true_se')' + 0.05 * randn(N_se + 1, py_se);
X_hat_se = sidLTVStateEst(Y_se, U_se, A_se, B_se, H_se);

ref15 = struct();
ref15.function_name = 'sidLTVStateEst';
ref15.params = struct();
ref15.input = struct('Y', Y_se, 'U', U_se, 'A', A_se, 'B', B_se, 'H', H_se);
ref15.output = struct('X_hat', X_hat_se);
ref15.tolerance = struct('X_hat_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_ltv_state_est.json'), ref15);

% ---- Test case 16: sidLTIfreqIO ----
fprintf('Generating reference_lti_freq_io.json...\n');
rng(58);
N_lti = 200; p_lti = 2; q_lti = 1;
A0_true = [0.8 0.1; -0.1 0.7];
B0_true = [1; 0.3];
H_lti = eye(p_lti);
X_lti = zeros(N_lti + 1, p_lti);
U_lti = randn(N_lti, q_lti);
X_lti(1, :) = randn(1, p_lti);
for k = 1:N_lti
    X_lti(k+1, :) = (A0_true * X_lti(k, :)' + B0_true * U_lti(k, :)')' ...
                     + 0.01 * randn(1, p_lti);
end
Y_lti = (H_lti * X_lti')';
[A0_est, B0_est] = sidLTIfreqIO(Y_lti, U_lti, H_lti);

ref16 = struct();
ref16.function_name = 'sidLTIfreqIO';
ref16.params = struct();
ref16.input = struct('Y', Y_lti, 'U', U_lti, 'H', H_lti);
ref16.output = struct('A0', A0_est, 'B0', B0_est);
ref16.tolerance = struct('A0_rel', 1e-6, 'B0_rel', 1e-6);

writeJSON(fullfile(thisDir, 'reference_lti_freq_io.json'), ref16);

% ---- Test case 17: sidTestMSD ----
fprintf('Generating reference_test_msd.json...\n');
m_msd = [2; 1; 3];
k_msd = [100; 200; 150];
c_msd = [5; 3; 4];
F_msd = [1 0; 0 0; 0 1];
Ts_msd = 0.01;
[Ad_msd, Bd_msd] = sidTestMSD(m_msd, k_msd, c_msd, F_msd, Ts_msd);

ref17 = struct();
ref17.function_name = 'sidTestMSD';
ref17.params = struct();
ref17.input = struct('m', m_msd, 'k_spring', k_msd, 'c_damp', c_msd, ...
                     'F', F_msd, 'Ts', Ts_msd);
ref17.output = struct('Ad', Ad_msd, 'Bd', Bd_msd);
ref17.tolerance = struct('Ad_rel', 1e-10, 'Bd_rel', 1e-10, 'Bd_atol', 1e-14);

writeJSON(fullfile(thisDir, 'reference_test_msd.json'), ref17);

fprintf('\n=== All reference data generated ===\n');

end


function cleanupShim(shimDir)
%CLEANUPSHIM Remove the temporary shim directory from the path and disk.
    if exist(shimDir, 'dir')
        rmpath(shimDir);
        rmdir(shimDir, 's');
    end
end


function writeJSON(filepath, data)
%WRITEJSON Write struct to JSON file.
    json = jsonencode(data);
    fid = fopen(filepath, 'w');
    if fid == -1
        error('Could not open %s for writing.', filepath);
    end
    fwrite(fid, json);
    fclose(fid);
end
