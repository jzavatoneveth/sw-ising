function [ Jmat ] = BuildPeriodicFourConnectedInteractionMatrix(N)
%BuildPeriodicFourConnectedInteractionMatrix: Symmetric interaction matrix
%for the Ising model on a two-dimensional square lattice with
%nearest-neighbor interactions and periodic boundary conditions.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License 

tic;
n = N^2;
Jmat = zeros(N, N, N, N);

for n1 = 1:N
    for n2 = 1:N
        if n1 == 1
            Jmat(n1, n2, N, n2) = 1;
            Jmat(n1, n2, n1+1, n2) = 1;
        elseif n1 == N
            Jmat(n1, n2, n1-1, n2) = 1;
            Jmat(n1, n2, 1, n2) = 1;
        else
            Jmat(n1, n2, n1-1, n2) = 1;
            Jmat(n1, n2, n1+1, n2) = 1;
        end
        
        if n2 == 1
            Jmat(n1, n2, n1, N) = 1;
            Jmat(n1, n2, n1, n2+1) = 1;
        elseif n2 == N
            Jmat(n1, n2, n1, n2-1) = 1;
            Jmat(n1, n2, n1, 1) = 1;
        else
            Jmat(n1, n2, n1, n2-1) = 1;
            Jmat(n1, n2, n1, n2+1) = 1;
        end
        
    end
end
Jmat = reshape(Jmat, n, n);

% Make the connectivity matrix sparse for the sake of efficiency
Jmat = sparse(Jmat);

fprintf('Built four-connected interaction matrix in %f seconds.\n', toc);

end

