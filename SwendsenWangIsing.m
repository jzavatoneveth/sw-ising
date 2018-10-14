function [E, M, x, PRNGState, X] = SwendsenWangIsing( N, T, J, nIter, displayIter )
%SWENDSENWANGISING Ising model simulation using the Swendsen-Wang algorithm.
%   This implementation uses MATLAB's sparse matrix API for maximum
%   computational efficiency and flexibility when the interaction matrix is
%   sparse.
%
%   Syntaxes:
%       [ E, M, x, PRNGState ] = SwendsenWangIsing( N, T, J, nIter )
%       [ ... ]                = SwendsenWangIsing( ..., displayIter )
%       [ ..., X ]             = SwendsenWangIsing( ... )
%
%   Inputs:
%       'N' - An integer indicating the desired total dimensionality of the
%               problem. For instance, if a 100 x 100 two-dimensional
%               lattice is desired, N should be 10000.
%
%       'T' - The thermodynamic temperature at which the system is
%               simulated. T must be nonzero.
%
%       'J' - The interaction matrix of the quadratic form defined by the
%               Hamiltonian. J must be a symmetric, non-negative matrix of
%               size NxN.
%
%       'nIter' - The number of iterations for which the simulation will be
%               run.
%
%       'displayIter' - The interval, in iterations, at which a status
%               update is printed to the terminal. If displayIter < 0, no
%               updates are displayed. If this input is ommitted, it
%               defaults to 0.
%
%   Outputs:
%       'E' - An nIter by 1 vector of the energy at each iteration.
%
%       'M' - An nIter by 1 vector of the magnetization at each iteration.
%
%       'x' - An N by 1 vector of the final spin state.
%
%       'PRNGState' - The state of the MATLAB PRNG as initialized at the
%               start of the simulation.
%
%       'X' - An nIter by N matrix containing the spin states at every
%               iteration. If nargout = 4, only the current spin state is
%               stored, increasing the speed and memory efficiency of the
%               simulation.
%
%   References:
%       [1] Swendsen, Robert H., and Jian-Sheng Wang. "Nonuniversal
%           critical dynamics in Monte Carlo simulations." Physical Review
%           Letters 58, no. 2 (1987): 86.
%
%       [2] Wang, Jian-Sheng, and Robert H. Swendsen. "Cluster Monte Carlo
%           algorithms." Physica A: Statistical Mechanics and its
%           Applications 167, no. 3 (1990): 565-579.
%
%       [3] Gilbert, John R., Cleve Moler, and Robert Schreiber. "Sparse
%           matrices in MATLAB: Design and implementation." SIAM Journal on
%           Matrix Analysis and Applications 13, no. 1 (1992): 333-356.
%
%
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License


%% Validate inputs

% Call local input validation function
ValidateSWIsingInputs(N, T, J, nIter);

% Check whether displayIter has been specified
if nargin == 4
    displayIter = 0;
end

% Check whether the full spin state history has been requested
if nargout > 4
    saveSpins = true;
else
    saveSpins = false;
end

%% Initialize the PRNG and spin state

% Initialize the MATLAB PRNG
% Note that MATLAB uses a global PRNG state. This call seeds the PRNG using
% the current time, and sets the PRNG to the Mersenne twister generator.
PRNGState = rng('shuffle', 'twister');

% Initialize the spin state
x = (-1).^round(rand(N,1));

%% Format the interaction matrix

% Make sure the interaction matrix is sparse
J = sparse(J);

% As the interactions are symmetric, take the lower trangular portion of J
J = tril(J, -1);

% Evaluate the link probabilities
prob = spfun(@(x) 1 - exp(-2 * x / T), J);

%% Run the simulation

% Allocate arrays to store the energy, magnetization, and (if desired) spin
% state at each iteration
E = zeros(nIter, 1);
M = zeros(nIter, 1);
if saveSpins
    X = zeros(nIter, N);
end

% If desired, store the spin state
if saveSpins
    X(1,:) = x;
end

% Compute the magnetization and energy
M(1) = mean(x);
E(1) = -(x' * J * x)/N;

% Print an update
if displayIter > 0
    fprintf('temperature %f, iteration %05d, energy %+f\n', T, 1, E(1));
end

% Iterate
tic;
for iter = 2:nIter

    % Find which spins are aligned
    sameSpin = tril(x' .* J .* x, -1) > 0;

    % Form the adjacency matrix of the graph
    A = (sprand(prob) < prob) .* sameSpin;
    A = A + speye(N);

    % Form the graph
    G = graph(A, 'lower');

    % Find connected components in the graph
    C = conncomp(G, 'OutputForm', 'cell');

    % Flip each cluster with probability 1/2
    for ind = 1:length(C)
        if rand(1) < 0.5
            x(C{ind}) = -x(C{ind});
        end
    end

    % If desired, store the spin state
    if saveSpins
        X(iter,:) = x;
    end

    % Compute the magnetization and energy
    M(iter) = mean(x);
    E(iter) = -(x' * J * x)/N;

    % Print an update
    if displayIter > 0 && ~mod(iter, displayIter)
        fprintf('temperature %f, iteration %05d, energy %+f\n', T, iter, E(iter));
    end
end

% Print an update
if displayIter > 0
    fprintf('Completed simulation in %f seconds.\n', toc);
end

end

function ValidateSWIsingInputs(N, T, J, nIter)
% Local utility function to check inputs

% Ensure that the dimensionality is a scalar
if numel(N) > 1
    error('SwendsenWangIsing:Problem dimensionality must be a scalar.');
end

% Ensure that the dimensionality is at least one
if (rem(N, 1) > 0) || N < 1
    error('SwendsenWangIsing:Problem dimensionality must an integer greater than or equal to one.');
end

% Ensure that the temperature is a scalar
if numel(T) > 1
    error('SwendsenWangIsing:Temperature must be a scalar.');
end

% Ensure that the temperature is nonzero
if T==0
    error('SwendsenWangIsing:Temperature must be non-zero.');
end

% Ensure that the interaction matrix is of the correct size
if ~isequal(size(J), [N, N])
    error('SwendsenWangIsing:Interaction matrix must have size equal to the dimensionality of the problem.');
end

% Ensure that all interactions are non-negative
if any(J(:) < 0)
    error('SwendsenWangIsing:Interaction strengths must be non-negative.');
end

% Ensure that some interactions are nonzero
if ~any(J(:))
    error('SwendsenWangIsing:Some interaction strengths must be nonzero.')
end

% Ensure that the interaction matrix is symmetric
if ~issymmetric(J)
    error('SwendsenWangIsing:Interactions must be symmetric.');
end

% Ensure that the number of iterations is a scalar
if numel(nIter) > 1
    error('SwendsenWangIsing:Number of iterations must be a scalar.');
end

% Ensure that the number of iterations is at least one
if (rem(nIter, 1) > 0) || nIter < 1
    error('SwendsenWangIsing:Number of iterations must an integer greater than or equal to one.');
end

end
