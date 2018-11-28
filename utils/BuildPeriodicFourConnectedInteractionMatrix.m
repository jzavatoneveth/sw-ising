function [ Jmat ] = BuildPeriodicFourConnectedInteractionMatrix(N)
%BuildPeriodicFourConnectedInteractionMatrix: Symmetric interaction matrix
%for the Ising model on a two-dimensional square lattice with
%nearest-neighbor interactions and periodic boundary conditions.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License 

tic;

% Get the total lattice size
sz = N^2;

% Compute a mesh of indices
[ri, ci] = meshgrid((1:N)');
ri = ri(:);
ci = ci(:);

% Compute the linear indices of each element
lia = ci + (ri - 1) * N;

% Define a function to handle boundary values
wf = @(x) x .* (x>0 & x < N+1) + N * double(x==0) + double(x==N+1);

% Compute the linear indices of the four-connected neighbors
a = wf(ci-1) + (wf(ri) - 1) * N;
b = wf(ci+1) + (wf(ri) - 1) * N;
c = ci + (wf(ri-1) - 1) * N;
d = ci + (wf(ri+1) - 1) * N;

% Form the sparse adjacency matrix
ii = [a;b;c;d];
jj = repmat(lia, 4, 1);
v = ones(length(jj), 1);
Jmat = sparse(ii, jj, v, sz, sz);

fprintf('Built four-connected interaction matrix in %f seconds.\n', toc);

end

