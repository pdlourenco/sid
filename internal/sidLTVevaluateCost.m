function [cost, fidelity, reg] = sidLTVevaluateCost(A, B, D, Xl, lambda, N, p, q)
% SIDLTVEVALUATECOST Compute COSMIC cost function value.
%
%   cost      = fidelity + reg
%   fidelity  = (1/2) Sigma_k ||D(k) C(k) - X'(k)||^2_F
%   reg       = (1/2) Sigma_k lambda_k ||C(k) - C(k-1)||^2_F
%   Handles both 3D array D (uniform) and cell array D (variable-length).

    fidelity = 0;
    priorVec = zeros(N - 1, 1);
    useCell = iscell(D);

    for k = 1:N
        % Data fidelity: ||D(k)*C(k) - X'(k)||^2_F
        % C(k) = [A(k)'; B(k)'] so D(k)*C(k) = D(k)*[A'; B']
        Ck = [A(:, :, k)'; B(:, :, k)'];
        if useCell
            residual = D{k} * Ck - Xl{k};
        else
            residual = D(:, :, k) * Ck - Xl(:, :, k);
        end
        fidelity = fidelity + norm(residual, 'fro')^2;

        % Regularization: ||C(k) - C(k-1)||^2_F
        if k < N
            Ck1 = [A(:, :, k+1)'; B(:, :, k+1)'];
            priorVec(k) = norm(Ck - Ck1, 'fro')^2;
        end
    end

    fidelity = 0.5 * fidelity;
    reg      = 0.5 * lambda' * priorVec;
    cost     = fidelity + reg;
end
