%% runAllTests - Master test runner for the sid-matlab package
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

fprintf('=== sid-matlab Test Suite ===\n\n');

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

runner__nTests = length(runner__testFiles);
runner__passed = 0;
runner__failed = 0;
runner__failedNames = {};

try
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
