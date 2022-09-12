function [A, n] = fitLogDistanceModel( rssi, dist )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
tbl = table(dist, rssi, 'VariableNames',{'distances','rssi'});
beta0 = [2, 0];
modelfun = @(b,dist)10.*b(1).*log10(dist) + b(2);
mdl = fitnlm(tbl,modelfun,beta0);
b = mdl.Coefficients{1:2,{'Estimate'}};
n = b(1);
A = b(2);

end

