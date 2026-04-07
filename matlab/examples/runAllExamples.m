%% runAllExamples - Run all examples and report pass/fail summary
%
% Discovers and runs all example*.m scripts in this directory, then reports
% a pass/fail summary. Designed for headless CI execution (all figures are
% closed after each example). New examples are picked up automatically —
% no manifest to edit.
%
% Usage:
%   run('matlab/examples/runAllExamples.m')

fprintf('=== sid-matlab Examples ===\n\n');

% Add paths — resolve this script's directory robustly (mfilename may
% return only the base name when invoked via run() in some environments).
runner__thisDir = fileparts(mfilename('fullpath'));
if isempty(runner__thisDir)
    runner__thisDir = fileparts(which(mfilename));
end
runner__matlabDir = fileparts(runner__thisDir);
runner__sidDir = fullfile(runner__matlabDir, 'sid');
addpath(runner__sidDir);

% Auto-discover example files matching example*.m
runner__listing = dir(fullfile(runner__thisDir, 'example*.m'));
if isempty(runner__listing)
    error('sid:noExamples', 'No example*.m files found in %s', runner__thisDir);
end
runner__exampleFiles = sort({runner__listing.name});
% Strip .m extension (2 characters)
for runner__k = 1:length(runner__exampleFiles)
    runner__exampleFiles{runner__k} = runner__exampleFiles{runner__k}(1:end-2);
end

runner__nExamples = length(runner__exampleFiles);
runner__passed = 0;
runner__failed = 0;
runner__failedNames = {};
runner__completedSections = 0;
runner__totalSections = 0;

for runner__k = 1:runner__nExamples
    fprintf('Running %s...\n', runner__exampleFiles{runner__k});
    runner__nCompleted = 0;
    try
        run(fullfile(runner__thisDir, [runner__exampleFiles{runner__k} '.m']));
        close all;
        runner__passed = runner__passed + 1;
        runner__completedSections = runner__completedSections + runner__nCompleted;
        runner__totalSections = runner__totalSections + runner__nCompleted;
        fprintf('  %s: OK (%d sections)\n', runner__exampleFiles{runner__k}, ...
            runner__nCompleted);
    catch runner__e
        close all;
        runner__failed = runner__failed + 1;
        runner__failedNames{end+1} = runner__exampleFiles{runner__k}; %#ok<SAGROW>
        runner__completedSections = runner__completedSections + runner__nCompleted;
        runner__totalSections = runner__totalSections + runner__nCompleted + 1;
        fprintf('  *** %s: FAILED (after %d sections) ***\n', ...
            runner__exampleFiles{runner__k}, runner__nCompleted);
        fprintf('      Error: %s\n', runner__e.message);
        % Emit GitHub Actions annotation
        fprintf('::error title=%s::%s\n', ...
            runner__exampleFiles{runner__k}, ...
            strrep(runner__e.message, newline, ' '));
    end
end

fprintf('\n=== Examples Summary ===\n');
fprintf('  Files:    %d passed, %d failed (%d total)\n', ...
    runner__passed, runner__failed, runner__nExamples);
fprintf('  Sections: %d completed, %d failed (%d total)\n', ...
    runner__completedSections, ...
    runner__totalSections - runner__completedSections, ...
    runner__totalSections);

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
