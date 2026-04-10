function x_traj = util_msd_nl(m, k_lin, k_cubic, c_damp, F, Ts, u, x0, substeps)
% UTIL_MSD_NL  RK4 simulation of an n-mass chain with Duffing cubic stiffness.
%
%   x = util_msd_nl(m, k_lin, k_cubic, c_damp, F, Ts, u)
%   x = util_msd_nl(m, k_lin, k_cubic, c_damp, F, Ts, u, x0)
%   x = util_msd_nl(m, k_lin, k_cubic, c_damp, F, Ts, u, x0, substeps)
%
%   Simulates the nonlinear n-mass chain
%
%       M x'' + C x' + K_lin(x) + K_cube(x) = F u(t)
%
%   where spring i produces the restoring force
%
%       f_spring_i = k_lin(i) * delta_i + k_cubic(i) * delta_i^3
%
%   with delta_0 = x(1) (wall spring) and delta_i = x(i+1) - x(i) for
%   i >= 1. Damping is linear in the velocity differences. Inputs are held
%   constant over each sample interval (zero-order hold).
%
%   Integration uses classic fixed-step RK4 at step Ts/substeps.
%
%   INPUTS:
%     m        - (n x 1) masses (kg).
%     k_lin    - (n x 1) linear spring constants (N/m).
%     k_cubic  - (n x 1) cubic stiffness coefficients (N/m^3). Pass zeros
%                for a purely linear plant.
%     c_damp   - (n x 1) damping coefficients (N s/m).
%     F        - (n x q) force input distribution matrix.
%     Ts       - Sample time (s).
%     u        - (N x q) input signal. Row k is the ZOH value over
%                [k*Ts, (k+1)*Ts].
%     x0       - (optional, 2n x 1) initial state [positions; velocities].
%                Defaults to zeros.
%     substeps - (optional) RK4 sub-steps per sample Ts. Default 1.
%
%   OUTPUTS:
%     x_traj - ((N+1) x 2n) state trajectory. Row k is the state at
%              time (k-1)*Ts, so row 1 is x0 and row N+1 is the final
%              state.
%
%   EXAMPLE:
%     rng(0);
%     N = 500; u = randn(N, 1);
%     x = util_msd_nl([1], [100], [1000], [0.5], [1], 0.01, u, [], 4);
%
%   See also: util_msd, util_msd_ltv
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
%  -----------------------------------------------------------------------

    if nargin < 8 || isempty(x0)
        x0 = [];
    end
    if nargin < 9 || isempty(substeps)
        substeps = 1;
    end

    m = m(:);
    k_lin = k_lin(:);
    k_cubic = k_cubic(:);
    c_damp = c_damp(:);

    n_mass = length(m);
    if length(k_lin) ~= n_mass
        error('util_msd_nl:shape', 'k_lin must have length n = %d', n_mass);
    end
    if length(k_cubic) ~= n_mass
        error('util_msd_nl:shape', 'k_cubic must have length n = %d', n_mass);
    end
    if length(c_damp) ~= n_mass
        error('util_msd_nl:shape', 'c_damp must have length n = %d', n_mass);
    end
    if size(F, 1) ~= n_mass
        error('util_msd_nl:shape', 'F must have %d rows, got %d', ...
            n_mass, size(F, 1));
    end
    if isvector(u)
        u = u(:);
    end
    if size(u, 2) ~= size(F, 2)
        error('util_msd_nl:shape', ...
            'u has %d inputs but F has %d columns', size(u, 2), size(F, 2));
    end
    if substeps < 1
        error('util_msd_nl:substeps', 'substeps must be >= 1, got %d', substeps);
    end

    N_samples = size(u, 1);
    n_state = 2 * n_mass;

    if isempty(x0)
        x = zeros(n_state, 1);
    else
        x0 = x0(:);
        if length(x0) ~= n_state
            error('util_msd_nl:shape', ...
                'x0 must have length %d, got %d', n_state, length(x0));
        end
        x = x0;
    end

    inv_m = 1.0 ./ m;
    h = Ts / substeps;

    x_traj = zeros(N_samples + 1, n_state);
    x_traj(1, :) = x(:)';
    for k = 1:N_samples
        uk = u(k, :)';
        for s = 1:substeps
            k1 = rhs_(x,            uk);
            k2 = rhs_(x + 0.5*h*k1, uk);
            k3 = rhs_(x + 0.5*h*k2, uk);
            k4 = rhs_(x + h*k3,     uk);
            x = x + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4);
        end
        x_traj(k + 1, :) = x(:)';
    end

    function dxdt = rhs_(state, force_vec)
        pos = state(1:n_mass);
        vel = state(n_mass+1:end);
        net_force = zeros(n_mass, 1);

        % Spring 1: between wall and mass 1
        delta0 = pos(1);
        dvel0  = vel(1);
        f0 = k_lin(1) * delta0 + k_cubic(1) * delta0^3 + c_damp(1) * dvel0;
        net_force(1) = net_force(1) - f0;

        % Springs 2..n: between mass i-1 and mass i
        for i = 2:n_mass
            delta = pos(i) - pos(i - 1);
            dvel  = vel(i) - vel(i - 1);
            fi = k_lin(i) * delta + k_cubic(i) * delta^3 + c_damp(i) * dvel;
            net_force(i - 1) = net_force(i - 1) + fi;
            net_force(i)     = net_force(i) - fi;
        end

        % External force (ZOH over the sample interval)
        net_force = net_force + F * force_vec;

        acc = inv_m .* net_force;
        dxdt = [vel; acc];
    end

end
