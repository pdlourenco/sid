%% runAllTests - Master test runner for the sid package
%
% Discovers and runs all test_*.m scripts in this directory, then reports
% a pass/fail summary. New tests are picked up automatically — no manifest
% to edit.
%
% Usage:
%   cd tests
%   runAllTests
%
% Or from the project root:
%   run('matlab/tests/runAllTests.m')

fprintf('=== sid Test Suite ===\n\n');

% Add paths — resolve this script's directory robustly (mfilename may
% return only the base name when invoked via run() in some environments).
runner__thisDir = fileparts(mfilename('fullpath'));
if isempty(runner__thisDir)
    runner__thisDir = fileparts(which(mfilename));
end
runner__matlabDir = fileparts(runner__thisDir);
runner__sidDir = fullfile(runner__matlabDir, 'sid');
addpath(runner__sidDir);
% MATLAB ignores addpath on directories named 'private'. Copy to a
% temporary non-private-named directory for test-only access. Uses only
% cross-platform MATLAB builtins (copyfile, rmdir) — no OS calls.
runner__privateDir = fullfile(runner__sidDir, 'private');
runner__shimDir = fullfile(runner__thisDir, 'private_test_shim');
if exist(runner__shimDir, 'dir')
    rmdir(runner__shimDir, 's');
end
mkdir(runner__shimDir);
copyfile(fullfile(runner__privateDir, '*.m'), runner__shimDir);
addpath(runner__shimDir);

% Auto-discover test files matching test_*.m
runner__listing = dir(fullfile(runner__thisDir, 'test_*.m'));
if isempty(runner__listing)
    error('sid:noTests', 'No test_*.m files found in %s', runner__thisDir);
end
runner__testFiles = sort({runner__listing.name});
% Strip .m extension (2 characters)
for runner__k = 1:length(runner__testFiles)
    runner__testFiles{runner__k} = runner__testFiles{runner__k}(1:end-2);
end

runner__nFiles = length(runner__testFiles);
runner__filesPassed = 0;
runner__filesFailed = 0;
runner__failedNames = {};
runner__totalCases = 0;
runner__passedCases = 0;

try
    for runner__k = 1:runner__nFiles
        % Each test file sets runner__nPassed via the template convention.
        % Reset before each file so the workspace variable is clean.
        runner__nPassed = 0;
        try
            run(fullfile(runner__thisDir, ...
                [runner__testFiles{runner__k} '.m']));
            runner__filesPassed = runner__filesPassed + 1;
            runner__passedCases = runner__passedCases + runner__nPassed;
            runner__totalCases  = runner__totalCases  + runner__nPassed;
        catch runner__e
            runner__filesFailed = runner__filesFailed + 1;
            runner__failedNames{end+1} = runner__testFiles{runner__k}; %#ok<SAGROW>
            fprintf('  *** %s: FAILED ***\n', runner__testFiles{runner__k});
            fprintf('      Error: %s\n', runner__e.message);
            % Credit any cases that passed before the failure, plus 1 for
            % the case that failed.
            runner__passedCases = runner__passedCases + runner__nPassed;
            runner__totalCases  = runner__totalCases  + runner__nPassed + 1;
            % Emit GitHub Actions annotation so the error appears in CI check-run output
            fprintf('::error title=%s::%s\n', ...
                runner__testFiles{runner__k}, ...
                strrep(runner__e.message, newline, ' '));
        end
    end
catch runner__fatalErr
    % Ensure shim cleanup even on unexpected errors
    rmpath(runner__shimDir);
    rmdir(runner__shimDir, 's');
    rethrow(runner__fatalErr);
end

% Clean up temporary test shim directory
rmpath(runner__shimDir);
rmdir(runner__shimDir, 's');

fprintf('\n=== Test Summary ===\n');
fprintf('  Files:  %d passed, %d failed (%d total)\n', ...
    runner__filesPassed, runner__filesFailed, runner__nFiles);
fprintf('  Cases:  %d passed, %d failed (%d total)\n', ...
    runner__passedCases, runner__totalCases - runner__passedCases, ...
    runner__totalCases);

if runner__filesFailed > 0
    fprintf('\n  Failed tests:\n');
    for runner__k = 1:length(runner__failedNames)
        fprintf('    - %s\n', runner__failedNames{runner__k});
    end
    fprintf('\n');
    error('sid:testsFailed', '%d of %d test files failed (%d/%d cases passed).', ...
        runner__filesFailed, runner__nFiles, ...
        runner__passedCases, runner__totalCases);
else
    fprintf('\n  ALL TESTS PASSED\n\n');
end
