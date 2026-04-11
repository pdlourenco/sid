function [Ad, Bd] = util_msd(m, k_spring, c_damp, F, Ts)
% UTIL_MSD  LTI n-mass spring-damper chain (exact ZOH discretization).
%
%   [Ad, Bd] = util_msd(m, k_spring, c_damp, F, Ts)
%
%   Builds the continuous-time state-space model for a chain of n masses
%   connected by springs and dampers
%
%       wall --k1,c1-- m1 --k2,c2-- m2 --k3,c3-- ... --kn,cn-- mn
%
%   and discretizes it using the matrix exponential (zero-order hold).
%
%   State vector: [x1; x2; ...; xn; v1; v2; ...; vn]
%                 (n positions followed by n velocities).
%
%   INPUTS:
%     m        - (n x 1) or (1 x n) masses (kg). n is inferred from length(m).
%     k_spring - (n x 1) spring constants (N/m). k_spring(1) is the
%                wall-to-m1 spring; k_spring(i), i >= 2, is the spring
%                between m_{i-1} and m_i.
%     c_damp   - (n x 1) damping coefficients (N s/m), same indexing as
%                k_spring.
%     F        - (n x q) force input distribution matrix.
%     Ts       - Sample time (s).
%
%   OUTPUTS:
%     Ad - (2n x 2n) discrete dynamics matrix.
%     Bd - (2n x q)  discrete input matrix.
%
%   This helper is the canonical spring-mass-damper plant simulator
%   for both the example suite (matlab/examples/example*.m) and the
%   MATLAB unit tests (matlab/tests/test_*.m), which pick it up via
%   addpath(matlab/examples) in runAllTests.m. See spec/EXAMPLES.md
%   section 2.1 for the binding contract.
%
%   EXAMPLE:
%     m = [1; 1]; k = [100; 80]; c = [0.5; 0.5];
%     F = [1; 0]; Ts = 0.01;
%     [Ad, Bd] = util_msd(m, k, c, F, Ts);
%
%   See also: util_msd_ltv, util_msd_nl
%
%   Changelog:
%   2026-04-10: First version by Pedro Lourenco.
%
%  -----------------------------------------------------------------------
%   Copyright (c) 2026 Pedro Lourenco, All rights reserved.
%   This code is released under the MIT License. See LICENSE file in the
%   project root for full license information.
%
%   This function is part of the Open Source System Identification
%   Toolbox (SID).
%   For full documentation and examples, visit
%   https://github.com/pdlourenco/sid
%  -----------------------------------------------------------------------

    m = m(:);
    k_spring = k_spring(:);
    c_damp = c_damp(:);

    n_mass = length(m);
    if length(k_spring) ~= n_mass
        error('util_msd:shape', 'k_spring must have length n = %d', n_mass);
    end
    if length(c_damp) ~= n_mass
        error('util_msd:shape', 'c_damp must have length n = %d', n_mass);
    end
    if size(F, 1) ~= n_mass
        error('util_msd:shape', 'F must have %d rows, got %d', ...
            n_mass, size(F, 1));
    end

    K_mat = zeros(n_mass, n_mass);
    C_mat = zeros(n_mass, n_mass);
    for i = 1:n_mass
        if i < n_mass
            K_mat(i, i)   = k_spring(i) + k_spring(i + 1);
            C_mat(i, i)   = c_damp(i)   + c_damp(i + 1);
            K_mat(i, i+1) = -k_spring(i + 1);
            K_mat(i+1, i) = -k_spring(i + 1);
            C_mat(i, i+1) = -c_damp(i + 1);
            C_mat(i+1, i) = -c_damp(i + 1);
        else
            K_mat(i, i) = k_spring(i);
            C_mat(i, i) = c_damp(i);
        end
    end

    M_mat = diag(m);
    Ac = [zeros(n_mass), eye(n_mass)
          -M_mat \ K_mat, -M_mat \ C_mat];
    Bc = [zeros(n_mass, size(F, 2))
          M_mat \ F];

    n_state = 2 * n_mass;
    Ad = expm(Ac * Ts);
    Bd = Ac \ (Ad - eye(n_state)) * Bc;
end
