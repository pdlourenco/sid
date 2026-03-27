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
thisDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(thisDir);
addpath(rootDir);
addpath(fullfile(rootDir, 'internal'));

testFiles = {
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
    'test_compareSpa'
    'test_compareEtfe'
    'test_compareSpafdr'
};

nTests = length(testFiles);
passed = 0;
failed = 0;
failedNames = {};

for i = 1:nTests
    try
        run(fullfile(thisDir, [testFiles{i} '.m']));
        passed = passed + 1;
    catch e
        failed = failed + 1;
        failedNames{end+1} = testFiles{i}; %#ok<SAGROW>
        fprintf('  *** %s: FAILED ***\n', testFiles{i});
        fprintf('      Error: %s\n', e.message);
    end
end

fprintf('\n=== Test Summary ===\n');
fprintf('  Total:  %d\n', nTests);
fprintf('  Passed: %d\n', passed);
fprintf('  Failed: %d\n', failed);

if failed > 0
    fprintf('\n  Failed tests:\n');
    for i = 1:length(failedNames)
        fprintf('    - %s\n', failedNames{i});
    end
    fprintf('\n');
    error('sid:testsFailed', '%d of %d test suites failed.', failed, nTests);
else
    fprintf('\n  ALL TESTS PASSED\n\n');
end
