%% runAllExamples - Run all examples and report pass/fail summary
%
% Validates that every example script runs without error. Designed for
% headless CI execution (all figures are closed after each example).
%
% Usage:
%   run('examples/runAllExamples.m')

fprintf('=== sid-matlab Examples ===\n\n');

% Add paths
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
addpath(rootDir);
addpath(fullfile(rootDir, 'internal'));

exampleFiles = {
    'exampleSISO'
    'exampleETFE'
    'exampleFreqDepRes'
    'exampleCoherence'
    'exampleMethodComparison'
    'exampleMIMO'
    'exampleFreqMap'
    'exampleSpectrogram'
    'exampleLTVdisc'
    'exampleMultiTrajectory'
    'exampleOutputCOSMIC'
};

nExamples = length(exampleFiles);
passed = 0;
failed = 0;
failedNames = {};

for i = 1:nExamples
    fprintf('Running %s...\n', exampleFiles{i});
    try
        run(fullfile(thisDir, [exampleFiles{i} '.m']));
        close all;
        passed = passed + 1;
        fprintf('  %s: OK\n', exampleFiles{i});
    catch e
        close all;
        failed = failed + 1;
        failedNames{end+1} = exampleFiles{i}; %#ok<SAGROW>
        fprintf('  *** %s: FAILED ***\n', exampleFiles{i});
        fprintf('      Error: %s\n', e.message);
        % Emit GitHub Actions annotation
        fprintf('::error title=%s::%s\n', exampleFiles{i}, strrep(e.message, newline, ' '));
    end
end

fprintf('\n=== Examples Summary ===\n');
fprintf('  Total:  %d\n', nExamples);
fprintf('  Passed: %d\n', passed);
fprintf('  Failed: %d\n', failed);

if failed > 0
    fprintf('\n  Failed examples:\n');
    for i = 1:length(failedNames)
        fprintf('    - %s\n', failedNames{i});
    end
    fprintf('\n');
    error('sid:examplesFailed', '%d of %d examples failed.', failed, nExamples);
else
    fprintf('\n  ALL EXAMPLES PASSED\n\n');
end
