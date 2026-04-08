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
ref2.tolerance = struct('Response_rel', 1e-10, 'NoiseSpectrum_rel', 1e-10);

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
ref5.tolerance = struct('A_rel', 1e-10, 'B_rel', 1e-10);

writeJSON(fullfile(thisDir, 'reference_ltv_cosmic.json'), ref5);

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
