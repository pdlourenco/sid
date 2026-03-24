function sidInstall()
%SIDINSTALL Add the sid package to the MATLAB/Octave path.
%
%   sidInstall
%
%   Adds the sid root folder and the internal subfolder to the path.
%   Run this once per session, or add the following to your startup.m:
%
%       run('/path/to/sid/sidInstall.m')
%

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    addpath(fullfile(rootDir, 'internal'));
    fprintf('sid: added to path.\n');
    fprintf('  %s\n', rootDir);
    fprintf('  %s\n', fullfile(rootDir, 'internal'));
end
