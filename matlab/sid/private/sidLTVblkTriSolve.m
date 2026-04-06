function [w, Lbd] = sidLTVblkTriSolve(S, U, Theta)
% SIDLTVBLKTRISOLVE Generic block tridiagonal forward-backward solver.
%
%   [w, Lbd] = sidLTVblkTriSolve(S, U, Theta)
%
%   Solves a symmetric block tridiagonal linear system using Gaussian
%   elimination (forward pass) followed by back-substitution (backward
%   pass). Supports non-uniform block sizes via cell arrays.
%
%   The system has the form:
%
%       [ S{1}    U{1}                          ] [ w{1}   ]   [ Theta{1}   ]
%       [ U{1}'   S{2}    U{2}                  ] [ w{2}   ]   [ Theta{2}   ]
%       [         U{2}'   S{3}   ...             ] [ ...    ] = [ ...        ]
%       [                  ...   ...   U{K-1}    ] [ ...    ]   [ ...        ]
%       [                        U{K-1}'  S{K}   ] [ w{K}   ]   [ Theta{K}   ]
%
%   where the sub-diagonal blocks are U{k}' (symmetric Hessian).
%
%   INPUTS:
%     S     - Cell array of length K. S{k} is a (m_k x m_k) symmetric
%             positive definite diagonal block.
%     U     - Cell array of length K-1. U{k} is a (m_k x m_{k+1})
%             super-diagonal block coupling step k to step k+1.
%             The sub-diagonal is U{k}' by symmetry.
%     Theta - Cell array of length K. Theta{k} is a (m_k x n_rhs)
%             right-hand side.
%
%   OUTPUTS:
%     w   - Cell array of length K. w{k} is a (m_k x n_rhs) solution.
%     Lbd - Cell array of length K. Lbd{k} is a (m_k x m_k) forward
%           Schur complement (stored for potential reuse in uncertainty).
%
%   EXAMPLES:
%     K = 5; m = 3;
%     S = arrayfun(@(k) eye(m) * (k+1), 1:K, 'UniformOutput', false);
%     U = arrayfun(@(k) 0.1*eye(m), 1:K-1, 'UniformOutput', false);
%     Theta = arrayfun(@(k) ones(m, 1), 1:K, 'UniformOutput', false);
%     [w, Lbd] = sidLTVblkTriSolve(S, U, Theta);
%
%   ALGORITHM:
%     Forward pass (Gaussian elimination):
%       Lbd{1} = S{1},  Y{1} = Lbd{1} \ Theta{1}
%       Lbd{k} = S{k} - U{k-1}' * (Lbd{k-1} \ U{k-1})
%       Y{k}   = Lbd{k} \ (Theta{k} - U{k-1}' * Y{k-1})
%     Backward pass (back-substitution):
%       w{K} = Y{K}
%       w{k} = Y{k} - Lbd{k} \ (U{k} * w{k+1})
%     Complexity: O(K * m^3) where m is the typical block size.
%
%   SPECIFICATION:
%     docs/cosmic_output.md — Appendices A and B
%
%   See also: sidLTVcosmicSolve, sidLTVStateEst
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

    K = length(S);

    Lbd = cell(K, 1);
    Y   = cell(K, 1);
    w   = cell(K, 1);

    % ---- Forward pass ----
    Lbd{1} = S{1};
    Y{1}   = Lbd{1} \ Theta{1};

    for k = 2:K
        rc = rcond(Lbd{k-1});
        if rc < eps
            warning('sid:singularLbd', ...
                ['Block tridiagonal forward pass: Lbd{%d} is near-singular ' ...
                 '(rcond=%.2e). Results may be unreliable.'], k-1, rc);
        end
        LbdInvU = Lbd{k-1} \ U{k-1};
        Lbd{k}  = S{k} - U{k-1}' * LbdInvU;
        Y{k}    = Lbd{k} \ (Theta{k} - U{k-1}' * Y{k-1});
    end

    % ---- Backward pass ----
    w{K} = Y{K};

    for k = K-1:-1:1
        w{k} = Y{k} - Lbd{k} \ (U{k} * w{k+1});
    end

end
