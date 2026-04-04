%% runAllTests - Master test runner for the sid-matlab package
%
% Runs all test scripts and reports pass/fail summary.
%
% Usage:
%   cd tests
%   runAllTests
%
% Or from the project root:
%   run('tests/runAllTests.m')

fprintf('=== sid-matlab Test Suite ===\n\n');

% Add paths
runner__thisDir = fileparts(mfilename('fullpath'));
runner__rootDir = fileparts(runner__thisDir);
addpath(runner__rootDir);
addpath(fullfile(runner__rootDir, 'internal'));

runner__testFiles = {
    'test_sidHannWin'
    'test_sidCov'
    'test_sidDFT'
    'test_sidWindowedDFT'
    'test_sidUncertainty'
    'test_sidValidate'
    'test_sidFreqBT'
    'test_sidFreqETFE'
    'test_sidFreqBTFDR'
    'test_sidPlotting'
    'test_validation'
    'test_crossMethod'
    'test_sidSpectrogram'
    'test_sidFreqMap'
    'test_sidSpectrogramPlot'
    'test_sidMapPlot'
    'test_sidLTVdisc'
    'test_sidLTVdiscTune'
    'test_sidLTVdiscVarLen'
    'test_sidLTVdiscUncertainty'
    'test_sidLTVdiscFrozen'
    'test_compareSpa'
    'test_compareEtfe'
    'test_compareSpafdr'
    'test_compareWelch'
    'test_multiTrajectory'
    'test_compareMultiTraj'
    'test_sidDetrend'
    'test_sidResidual'
    'test_sidCompare'
    'test_sidModelOrder'
    'test_sidLTIfreqIO'
    'test_sidLTVStateEst'
    'test_sidLTVdiscIO'
};

runner__nTests = length(runner__testFiles);
runner__passed = 0;
runner__failed = 0;
runner__failedNames = {};

for runner__k = 1:runner__nTests
    try
        run(fullfile(runner__thisDir, [runner__testFiles{runner__k} '.m']));
        runner__passed = runner__passed + 1;
    catch runner__e
        runner__failed = runner__failed + 1;
        runner__failedNames{end+1} = runner__testFiles{runner__k}; %#ok<SAGROW>
        fprintf('  *** %s: FAILED ***\n', runner__testFiles{runner__k});
        fprintf('      Error: %s\n', runner__e.message);
        % Emit GitHub Actions annotation so the error appears in CI check-run output
        fprintf('::error title=%s::%s\n', ...
            runner__testFiles{runner__k}, ...
            strrep(runner__e.message, newline, ' '));
    end
end

fprintf('\n=== Test Summary ===\n');
fprintf('  Total:  %d\n', runner__nTests);
fprintf('  Passed: %d\n', runner__passed);
fprintf('  Failed: %d\n', runner__failed);

if runner__failed > 0
    fprintf('\n  Failed tests:\n');
    for runner__k = 1:length(runner__failedNames)
        fprintf('    - %s\n', runner__failedNames{runner__k});
    end
    fprintf('\n');
    error('sid:testsFailed', '%d of %d test suites failed.', ...
        runner__failed, runner__nTests);
else
    fprintf('\n  ALL TESTS PASSED\n\n');
end
