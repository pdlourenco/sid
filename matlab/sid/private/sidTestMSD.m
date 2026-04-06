function [Ad, Bd] = sidTestMSD(m, k_spring, c_damp, F, Ts)
% SIDTESTMSD 3-mass spring-damper chain ZOH discretization.
%
%   [Ad, Bd] = sidTestMSD(m, k_spring, c_damp, F, Ts)
%
%   Discretises a chain of 3 masses connected by springs and dampers
%   using exact zero-order hold (matrix exponential).
%
%   State vector: [x1; x2; x3; v1; v2; v3] (positions, velocities).
%   Chain topology: wall--k1--m1--k2--m2--k3--m3
%
%   INPUTS:
%     m        - (3 x 1) masses (kg).
%     k_spring - (3 x 1) spring constants (N/m).
%     c_damp   - (3 x 1) damping coefficients (N s/m).
%     F        - (3 x q) force input distribution matrix.
%     Ts       - Sample time (s).
%
%   OUTPUTS:
%     Ad - (6 x 6) discrete dynamics matrix.
%     Bd - (6 x q) discrete input matrix.
%
%   EXAMPLES:
%     m = [1; 1; 1]; k = [10; 10; 10]; c = [0.1; 0.1; 0.1];
%     F = eye(3); Ts = 0.01;
%     [Ad, Bd] = sidTestMSD(m, k, c, F, Ts);
%
%   See also: test_sidLTVStateEst, test_sidLTIfreqIO, test_sidLTVdiscIO
%
%   Changelog:
%   2026-04-01: First version by Pedro Lourenço.
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

    M_mat = diag(m);
    K_mat = [k_spring(1) + k_spring(2), -k_spring(2), 0
             -k_spring(2), k_spring(2) + k_spring(3), ...
             -k_spring(3)
             0, -k_spring(3), k_spring(3)];
    C_mat = [c_damp(1) + c_damp(2), -c_damp(2), 0
             -c_damp(2), c_damp(2) + c_damp(3), ...
             -c_damp(3)
             0, -c_damp(3), c_damp(3)];

    n_mass = length(m);
    Ac = [zeros(n_mass), eye(n_mass)
          -M_mat \ K_mat, -M_mat \ C_mat];
    Bc = [zeros(n_mass, size(F, 2))
          M_mat \ F];

    n_state = 2 * n_mass;
    Ad = expm(Ac * Ts);
    Bd = Ac \ (Ad - eye(n_state)) * Bc;
end
