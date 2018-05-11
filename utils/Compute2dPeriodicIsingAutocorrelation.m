function [ a ] = Compute2dPeriodicIsingAutocorrelation( spin, N )
%COMPUTE2DPERIODICISINGAUTOCORRELATION: Simple correlation computation
%utility.
%
%   Copyright (c) 2018 Jacob Zavatone-Veth, MIT License 

% We only want to consider floor(N/2+1) offsets due to the period boundary
% conditions
m = floor(N/2+1);

% Compute autocorrelation along rows of grid
meanrow = mean(spin, 2);
spin0 = spin - meanrow;
corrmat = zeros(m,N);
for kk=1:m
    corrmat(kk,:) = sum(circshift(spin0, [0 kk-1]).*spin0, 2);
end
corrtot1 = sum(corrmat, 2);

% Compute autocorrlation along columns of grid
meancolumn = mean(spin, 1);
spin0 = spin - meancolumn;
corrmat = zeros(m,N);
for kk=1:m
    corrmat(kk,:) = sum(circshift(spin0, [kk-1 0]).*spin0, 1);
end
corrtot2 = sum(corrmat,2);

% Average autocorrelations along axes and normalize by N^2
a = (corrtot1+corrtot2)/(2*N^2);

end