function [Ad, Bd] = util_msd_ltv(m, k_spring, c_damp, F, Ts, N)
% UTIL_MSD_LTV  LTV n-mass spring-damper chain with per-step parameters.
%
%   [Ad, Bd] = util_msd_ltv(m, k_spring, c_damp, F, Ts)
%   [Ad, Bd] = util_msd_ltv(m, k_spring, c_damp, F, Ts, N)
%
%   Builds a time-varying discrete state-space model for the n-mass chain
%   topology of util_msd. Any of m, k_spring, c_damp may be a constant
%   vector (n x 1) or a per-step trajectory (n x N). F may be (n x q) or
%   (n x q x N). N is inferred from the first time-varying input; if every
%   input is constant, N must be passed explicitly.
%
%   INPUTS:
%     m        - (n x 1) or (n x N) masses (kg).
%     k_spring - (n x 1) or (n x N) spring constants (N/m).
%     c_damp   - (n x 1) or (n x N) damping coefficients (N s/m).
%     F        - (n x q) or (n x q x N) force input distribution matrix.
%     Ts       - Sample time (s).
%     N        - (optional) number of time steps; required only when
%                every other input is time-invariant.
%
%   OUTPUTS:
%     Ad - (2n x 2n x N) discrete dynamics sequence.
%     Bd - (2n x q  x N) discrete input sequence.
%
%   EXAMPLE:
%     N = 200;
%     m = [1; 1]; c = [0.5; 0.5];
%     k_tv = zeros(2, N);
%     k_tv(1, :) = linspace(100, 150, N);   % time-varying k1
%     k_tv(2, :) = 80;
%     F = [1; 0];
%     [Ad, Bd] = util_msd_ltv(m, k_tv, c, F, 0.01);
%
%   See also: util_msd, util_msd_nl
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

    if nargin < 6
        N = [];
    end

    % Infer n_mass from m. Accept row or column vector for LTI, or
    % (n, N_m) matrix for LTV.
    if isvector(m)
        m_col = m(:);
        n_mass = length(m_col);
        m_is_tv = false;
    elseif ismatrix(m) && size(m, 2) > 1
        n_mass = size(m, 1);
        m_col = m;
        m_is_tv = true;
    else
        error('util_msd_ltv:shape', ...
            'm must be a vector or (n, N) matrix');
    end

    % Helper: classify a parameter as LTI vector (length n) or LTV
    % (n x N) matrix. Row vectors of length n are reshaped to columns.
    function [y, is_tv] = classify_(x, name)
        if isvector(x) && length(x) == n_mass
            y = x(:);
            is_tv = false;
        elseif ismatrix(x) && size(x, 1) == n_mass && size(x, 2) > 1
            y = x;
            is_tv = true;
        else
            error('util_msd_ltv:shape', ...
                '%s must be length-%d vector or %d-row matrix, got %s', ...
                name, n_mass, n_mass, mat2str(size(x)));
        end
    end

    [k_spring, k_is_tv] = classify_(k_spring, 'k_spring');
    [c_damp,   c_is_tv] = classify_(c_damp,   'c_damp');

    % F: (n, q) for LTI or (n, q, N) for LTV.
    if size(F, 1) ~= n_mass
        error('util_msd_ltv:shape', ...
            'F must have %d rows, got %d', n_mass, size(F, 1));
    end
    F_is_tv = ndims(F) == 3;

    if m_is_tv
        m = m_col;  % (n, N)
    else
        m = m_col;  % (n, 1)
    end

    % Infer N from the first time-varying input.
    N_candidates = [];
    if m_is_tv
        N_candidates(end+1) = size(m, 2);
    end
    if k_is_tv
        N_candidates(end+1) = size(k_spring, 2);
    end
    if c_is_tv
        N_candidates(end+1) = size(c_damp, 2);
    end
    if F_is_tv
        N_candidates(end+1) = size(F, 3);
    end

    if ~isempty(N_candidates)
        if any(N_candidates ~= N_candidates(1))
            error('util_msd_ltv:N', ...
                'Time-varying inputs disagree on N: %s', mat2str(N_candidates));
        end
        N_inferred = N_candidates(1);
        if ~isempty(N) && N ~= N_inferred
            error('util_msd_ltv:N', ...
                'Explicit N=%d conflicts with inferred N=%d', N, N_inferred);
        end
        N = N_inferred;
    else
        if isempty(N)
            error('util_msd_ltv:N', ...
                ['All inputs are time-invariant; pass N explicitly to ', ...
                 'replicate the LTI result across time steps.']);
        end
    end

    q = size(F, 2);

    Ad = zeros(2 * n_mass, 2 * n_mass, N);
    Bd = zeros(2 * n_mass, q, N);

    % Fast path: nothing varies in time -> compute once, replicate.
    if ~m_is_tv && ~k_is_tv && ~c_is_tv && ~F_is_tv
        [Ad_const, Bd_const] = util_msd(m, k_spring, c_damp, F, Ts);
        for k = 1:N
            Ad(:, :, k) = Ad_const;
            Bd(:, :, k) = Bd_const;
        end
        return;
    end

    for k = 1:N
        if m_is_tv
            m_k = m(:, k);
        else
            m_k = m;
        end
        if k_is_tv
            k_k = k_spring(:, k);
        else
            k_k = k_spring;
        end
        if c_is_tv
            c_k = c_damp(:, k);
        else
            c_k = c_damp;
        end
        if F_is_tv
            F_k = F(:, :, k);
        else
            F_k = F;
        end
        [Ad(:, :, k), Bd(:, :, k)] = util_msd(m_k, k_k, c_k, F_k, Ts);
    end
end
