function [bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, X_val, U_val, varargin)
%SIDLTVDISCTUNE Validation-based lambda tuning for sidLTVdisc.
%
%   [bestResult, bestLambda, allLosses] = sidLTVdiscTune(X_train, U_train, X_val, U_val)
%   [bestResult, bestLambda, allLosses] = sidLTVdiscTune(..., 'LambdaGrid', grid)
%
%   Runs sidLTVdisc for each candidate lambda, propagates the identified
%   model on validation initial conditions, and selects the lambda that
%   minimizes trajectory prediction RMSE.
%
%   INPUTS:
%     X_train - Training state data, (N+1 x p x L_train)
%     U_train - Training input data, (N x q x L_train)
%     X_val   - Validation state data, (N+1 x p x L_val)
%     U_val   - Validation input data, (N x q x L_val)
%
%   NAME-VALUE OPTIONS:
%     'LambdaGrid'   - Vector of candidate lambda values.
%                      Default: logspace(-3, 15, 50).
%     'Precondition'  - Passed through to sidLTVdisc. Default: false.
%     'Algorithm'     - Passed through to sidLTVdisc. Default: 'cosmic'.
%
%   OUTPUTS:
%     bestResult  - sidLTVdisc result struct at optimal lambda.
%     bestLambda  - Optimal scalar lambda value.
%     allLosses   - (nGrid x 1) trajectory prediction RMSE at each lambda.
%
%   TRAJECTORY PREDICTION LOSS:
%     For each validation trajectory l, the model is propagated from
%     x_l(0) using the identified A(k) and B(k):
%
%       x_hat(k+1) = A(k) x_hat(k) + B(k) u_l(k)
%
%     The loss is the average RMSE across validation trajectories:
%
%       L(lambda) = (1/L_val) sum_l sqrt( (1/N) sum_k ||x_hat(k) - x(k)||^2 )
%
%   EXAMPLE:
%     grid = logspace(0, 10, 30);
%     [best, bestLam, losses] = sidLTVdiscTune(Xtr, Utr, Xval, Uval, ...
%                                              'LambdaGrid', grid);
%     semilogx(grid, losses); xlabel('\lambda'); ylabel('RMSE');
%
%   See also: sidLTVdisc

    % ---- Parse options ----
    lambdaGrid = logspace(-3, 15, 50);
    extraArgs = {};

    k = 1;
    while k <= length(varargin)
        if ischar(varargin{k})
            switch lower(varargin{k})
                case 'lambdagrid'
                    lambdaGrid = varargin{k+1};
                    k = k + 2;
                case {'precondition', 'algorithm'}
                    extraArgs = [extraArgs, varargin(k), varargin(k+1)]; %#ok<AGROW>
                    k = k + 2;
                otherwise
                    error('sid:unknownOption', 'Unknown option: %s', varargin{k});
            end
        else
            error('sid:badInput', 'Expected option name at position %d.', k);
        end
    end

    lambdaGrid = lambdaGrid(:)';
    nGrid = length(lambdaGrid);

    % Dimensions
    N = size(X_val, 1) - 1;
    p = size(X_val, 2);
    L_val = size(X_val, 3);

    % ---- Grid search ----
    allLosses = zeros(nGrid, 1);

    for j = 1:nGrid
        % Identify model on training data
        res = sidLTVdisc(X_train, U_train, 'Lambda', lambdaGrid(j), extraArgs{:});

        % Propagate on validation trajectories
        allLosses(j) = trajectoryRMSE(res.A, res.B, X_val, U_val, N, p, L_val);
    end

    % ---- Select best ----
    [~, bestIdx] = min(allLosses);
    bestLambda = lambdaGrid(bestIdx);

    % ---- Re-run at optimal lambda ----
    bestResult = sidLTVdisc(X_train, U_train, 'Lambda', bestLambda, extraArgs{:});
end


function rmse = trajectoryRMSE(A, B, X_val, U_val, N, p, L_val)
%TRAJECTORYRMSE Average trajectory prediction RMSE over validation set.

    totalRMSE = 0;

    for l = 1:L_val
        x_hat = zeros(N + 1, p);
        x_hat(1, :) = X_val(1, :, l);   % initial condition from data

        for k = 1:N
            x_hat(k+1, :) = (A(:, :, k) * x_hat(k, :)' + ...
                             B(:, :, k) * reshape(U_val(k, :, l), [], 1))';
        end

        % RMSE for this trajectory
        x_true = X_val(:, :, l);
        err = x_hat - x_true;
        totalRMSE = totalRMSE + sqrt(mean(sum(err.^2, 2)));
    end

    rmse = totalRMSE / L_val;
end
