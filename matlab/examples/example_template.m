%% example_template - Template for sid-matlab example scripts
%
% Copy this file to create new examples. Follow the structure below
% so the example runner (runAllExamples.m) can track section progress.
%
% Convention:
%   - Initialize runner__nCompleted = 0 at the top
%   - Increment runner__nCompleted after each section completes
%   - Print 'Section N completed: description' for every section
%   - End with 'exampleName: N/N sections completed' summary line
%   - Each %% heading (after the file header) is one section
%
% NOTE: This file is named example_template.m (with underscore) so the
% runner's dir('example*.m') pattern picks it up. Rename when copying.

runner__nCompleted = 0;

%% Section 1: Generate example data
% Describe what this section does and why.
N = 256;
fs = 1;
t = (0:N-1)' / fs;
u = randn(N, 1);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Generate example data');

%% Section 2: Run identification
% Call sid functions here.
% result = sidFreqBT(u, y, fs);

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Run identification');

%% Section 3: Plot results
% Visualise the output.
% figure; plot(f, abs(G));

runner__nCompleted = runner__nCompleted + 1;
fprintf('  Section %d completed: %s.\n', ...
    runner__nCompleted, 'Plot results');

fprintf('example_template: %d/%d sections completed\n', ...
    runner__nCompleted, runner__nCompleted);
