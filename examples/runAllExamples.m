%% runAllExamples - Run all examples and report pass/fail summary
%
% Validates that every example script runs without error. Designed for
% headless CI execution (all figures are closed after each example).
%
% Usage:
%   run('examples/runAllExamples.m')

fprintf('=== sid-matlab Examples ===\n\n');

% Add paths
runner__thisDir = fileparts(mfilename('fullpath'));
runner__rootDir = fileparts(runner__thisDir);
addpath(runner__rootDir);
addpath(fullfile(runner__rootDir, 'internal'));

runner__exampleFiles = {
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

runner__nExamples = length(runner__exampleFiles);
runner__passed = 0;
runner__failed = 0;
runner__failedNames = {};

for runner__k = 1:runner__nExamples
    fprintf('Running %s...\n', runner__exampleFiles{runner__k});
    try
        run(fullfile(runner__thisDir, [runner__exampleFiles{runner__k} '.m']));
        close all;
        runner__passed = runner__passed + 1;
        fprintf('  %s: OK\n', runner__exampleFiles{runner__k});
    catch runner__e
        close all;
        runner__failed = runner__failed + 1;
        runner__failedNames{end+1} = runner__exampleFiles{runner__k}; %#ok<SAGROW>
        fprintf('  *** %s: FAILED ***\n', runner__exampleFiles{runner__k});
        fprintf('      Error: %s\n', runner__e.message);
        % Emit GitHub Actions annotation
        fprintf('::error title=%s::%s\n', ...
            runner__exampleFiles{runner__k}, ...
            strrep(runner__e.message, newline, ' '));
    end
end

fprintf('\n=== Examples Summary ===\n');
fprintf('  Total:  %d\n', runner__nExamples);
fprintf('  Passed: %d\n', runner__passed);
fprintf('  Failed: %d\n', runner__failed);

if runner__failed > 0
    fprintf('\n  Failed examples:\n');
    for runner__k = 1:length(runner__failedNames)
        fprintf('    - %s\n', runner__failedNames{runner__k});
    end
    fprintf('\n');
    error('sid:examplesFailed', '%d of %d examples failed.', ...
        runner__failed, runner__nExamples);
else
    fprintf('\n  ALL EXAMPLES PASSED\n\n');
end
