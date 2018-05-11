%% RunSW
% A script demonstrating use of SwendsenWangIsing to simulate the Ising
% model on a two-dimensional grid with period boundary conditions and
% nearest-neighbor interactions.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License 

%% Set parameters

% Add the utilities folder to the MATLAB path
addpath('utils');

% Grid size (NxN)
N = 25;

% Temperature range
T = [0.1:0.01:3.5, 3.6:0.1:10];

% Number of iterations
nIter = 4000;

% Interval at which to display updates
displayIter = 100;

% Number of iterations to ignore in computing statistics
nBurnin = 200;

% Build the interaction matrix
[ J ] = BuildPeriodicFourConnectedInteractionMatrix(N);

%% Simulate temperatures in serial

% Allocate containers
E_iter = zeros(nIter, length(T));
M_iter = zeros(nIter, length(T));
spin = zeros(N, N, length(T));
acorr = zeros(floor(N/2+1), length(T));

% Start a timer
timeAll = tic;

% Iterate in parallel
for ind = 1:length(T)
    fprintf('Working on temperature %f.\n', T(ind));
    
    % Run the simulation
    [E_iter(:,ind), M_iter(:,ind), x] = SwendsenWangIsing( N^2, T(ind), J, nIter, displayIter );
    
    % Reshape the spin state for output
    spin(:,:,ind) = reshape(x, N, N);
    
    % Compute autocorrelation
    [ acorr(:,ind) ] = Compute2dPeriodicIsingAutocorrelation( spin(:,:,ind), N );
end

% Get the elapsed time
elapsed = toc(timeAll);

% Print a status update
fprintf('\n\nSimulated %d temperatures in %f seconds.\n\n', length(T), elapsed);

%% Plot the results

MakeDemoPlots;