function tf = sidInputExcitationDegenerate(u, epsTol)
% SIDINPUTEXCITATIONDEGENERATE Whole-signal check for un-exciting input.
%
%   tf = sidInputExcitationDegenerate(u)
%   tf = sidInputExcitationDegenerate(u, epsTol)
%
%   Returns true when the input carries no usable excitation: a constant
%   input (zero variance about its own mean) or an identically-zero input
%   cannot identify any dynamics on the (0, pi] grid. This is an ABSOLUTE
%   check on the sample variance relative to the mean square, so it fires
%   even though a constant input has a nonzero Phi_u under the biased-
%   covariance convention (which the relative per-frequency floor of §10.2
%   can never detect). It is the MATLAB twin of the Python
%   input_excitation_degenerate helper, so both languages apply §10.3
%   identically.
%
%   INPUTS:
%     u      - Input signal, (N x 1), (N x nu), or (N x nu x L). [] (time
%              series mode) is never degenerate in this sense.
%     epsTol - Relative tolerance; max_ch(var) <= epsTol*max_ch(meanSq)
%              counts as constant. Default 1e-10.
%
%   OUTPUTS:
%     tf - Logical scalar; true when u is constant or identically zero.
%
%   EXAMPLES:
%     tf = sidInputExcitationDegenerate(ones(100, 1));   % -> true
%
%   SPECIFICATION:
%     SPEC.md §10.3 — Input Excitation Requirements
%
%   See also: sidRegularizeResponse, sidDeadInputChannels
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
        tf = false;
        return;
    end

    u = sidCollapseChannels(u);
    varAc = var(u, 1, 1);      % per-channel variance about the mean (÷N)
    meanSq = mean(u.^2, 1);    % per-channel mean square
    maxMs = max(meanSq);
    tf = (max(varAc) <= epsTol * maxMs) || (maxMs <= realmin);
end
