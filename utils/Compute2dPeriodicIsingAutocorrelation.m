function [ a ] = Compute2dPeriodicIsingAutocorrelation( spin, N )
%COMPUTE2DPERIODICISINGAUTOCORRELATION: Simple correlation computation
%utility.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License

% We only want to consider floor(N/2+1) offsets due to the periodic boundary
% conditions
numOffset = floor(N/2+1);

%% Compute autocorrelation along rows of grid

% Subtract off the mean spin
spin0 = spin - mean(spin, 2);

% Compute the correlation at each offset
rowCorr = zeros(m,N);
for n=1:numOffset
    rowCorr(n,:) = sum(circshift(spin0, [0 n-1]).*spin0, 2);
end
rowCorr = sum(rowCorr, 2);

%% Compute autocorrlation along columns of grid

% Subtract off the mean spin
spin0 = spin - mean(spin, 1);

% Compute the correlation at each offset
colCorr = zeros(m,N);
for n=1:numOffset
    colCorr(n,:) = sum(circshift(spin0, [n-1 0]).*spin0, 1);
end
colCorr = sum(colCorr,2);

%% Average autocorrelations along axes and normalize by N^2
a = (rowCorr+colCorr)/(2*N^2);

end
