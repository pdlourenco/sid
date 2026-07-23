function uc = sidCollapseChannels(u)
% SIDCOLLAPSECHANNELS Reshape input data to (samples x channels).
%
%   uc = sidCollapseChannels(u)
%
%   Normalises an input array to a 2-D (samples x channels) matrix for the
%   whole-signal degenerate-excitation checks: a column vector becomes
%   (N x 1), an (N x nu) matrix is returned unchanged, and an (N x nu x L)
%   multi-trajectory array is pooled across trajectories into (N*L x nu)
%   so the per-channel statistics see every sample. Mirrors the reshape at
%   the top of the Python input_excitation_degenerate / dead_input_channels
%   helpers.
%
%   INPUTS:
%     u - Input signal, (N x 1), (N x nu), or (N x nu x L).
%
%   OUTPUTS:
%     uc - (M x nu) matrix with all trajectories stacked along the rows.
%
%   EXAMPLES:
%     uc = sidCollapseChannels(randn(100, 2, 5));  % -> (500 x 2)
%
%   SPECIFICATION:
%     SPEC.md §10.3 — Input Excitation Requirements
%
%   See also: sidInputExcitationDegenerate, sidDeadInputChannels
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

    if ndims(u) == 3
        [N, nu, L] = size(u);
        uc = reshape(permute(u, [1 3 2]), N * L, nu);
    elseif isvector(u)
        uc = u(:);
    else
        uc = u;
    end
end
