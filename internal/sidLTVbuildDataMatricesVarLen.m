function [D, Xl] = sidLTVbuildDataMatricesVarLen(X, U, N, p, q, L, horizons)
% SIDLTVBUILDDATAMATRICESVARLEN Construct D(k) and Xl(k) for variable-length trajectories.
%
%   At each time step k, only trajectories with horizon > k contribute.
%   D{k}  has size (|L(k)| x p+q), Xl{k} has size (|L(k)| x p).
%   Each is normalized by 1/sqrt(|L(k)|) to keep the cost well-scaled.

    D  = cell(N, 1);
    Xl = cell(N, 1);

    for k = 0:N-1
        % Active trajectories at step k: those with horizon > k
        active = find(horizons > k);
        Lk = length(active);

        if Lk == 0
            % No trajectories active — empty matrices
            D{k+1}  = zeros(0, p + q);
            Xl{k+1} = zeros(0, p);
            continue;
        end

        sqrtLk = sqrt(Lk);
        Dk  = zeros(Lk, p + q);
        Xlk = zeros(Lk, p);

        for ii = 1:Lk
            l = active(ii);
            Dk(ii, :)  = [X{l}(k+1, :), U{l}(k+1, :)] / sqrtLk;
            Xlk(ii, :) = X{l}(k+2, :) / sqrtLk;
        end

        D{k+1}  = Dk;
        Xl{k+1} = Xlk;
    end
end
