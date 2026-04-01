function sidInstall()
% SIDINSTALL Add the sid package to the MATLAB/Octave path.
%
%   sidInstall
%
%   Adds the sid root folder and the internal subfolder to the path.
%   Run this once per session, or add the following to your startup.m:
%
%       run('/path/to/sid/sidInstall.m')
%
%   EXAMPLES:
%     sidInstall
%     result = sidFreqBT(y, u);
%
%   See also: sidFreqBT, sidFreqETFE, sidLTVdisc
%
%   Changelog:
%   2026-03-24: First version by Pedro Lourenço.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid-matlab
%  -----------------------------------------------------------------------

    rootDir = fileparts(mfilename('fullpath'));
    addpath(rootDir);
    addpath(fullfile(rootDir, 'internal'));
    fprintf('sid: added to path.\n');
    fprintf('  %s\n', rootDir);
    fprintf('  %s\n', fullfile(rootDir, 'internal'));
end
