function dead = sidDeadInputChannels(u, epsTol)
% SIDDEADINPUTCHANNELS Flag individual (near-)constant input channels.
%
%   dead = sidDeadInputChannels(u)
%   dead = sidDeadInputChannels(u, epsTol)
%
%   A single constant channel among otherwise active channels does NOT make
%   the whole input degenerate (sidInputExcitationDegenerate stays false
%   because a healthy channel dominates), and it may not push cond(Phi_u)
%   past 1/eps at every frequency -- yet that channel's transfer-function
%   column is unidentifiable. Callers use this to warn without NaNing the
%   healthy channels. MATLAB twin of the Python dead_input_channels helper.
%
%   INPUTS:
%     u      - Input signal, (N x 1), (N x nu), or (N x nu x L).
%     epsTol - Relative tolerance, as in sidInputExcitationDegenerate.
%              Default 1e-10.
%
%   OUTPUTS:
%     dead - (1 x nu) logical row; true where the channel is (near-)constant.
%
%   EXAMPLES:
%     dead = sidDeadInputChannels([randn(50,1), ones(50,1)]);  % -> [false true]
%
%   SPECIFICATION:
%     SPEC.md §10.3 — Input Excitation Requirements
%
%   See also: sidInputExcitationDegenerate, sidRegularizeResponse
%
%   Changelog:
%   2026-07-23: First version by Pedro Lourenço.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenço, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid
%  -----------------------------------------------------------------------

    if nargin < 2 || isempty(epsTol)
        epsTol = 1e-10;
    end
    if isempty(u)
        dead = false(1, 0);
        return;
    end

    u = sidCollapseChannels(u);
    varAc = var(u, 1, 1);      % per-channel variance about the mean (÷N)
    meanSq = mean(u.^2, 1);    % per-channel mean square
    scale = max(meanSq);
    dead = (varAc <= epsTol * scale) | (meanSq <= realmin);
end
