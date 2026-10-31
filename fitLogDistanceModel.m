function [referenceDbm, n] = fitLogDistanceModel(rssiDbm, dist_m, d0_m)
%FITLOGDISTANCEMODEL Fit signed dBm = referenceDbm - 10*n*log10(d/d0).
% Base MATLAB least squares; a measured trend, not isolated LOS calibration.
% API change: pass signed dBm, not magnitudes of negative dBm readings.
if nargin < 3
    d0_m = 1;
end
validateattributes(d0_m, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(dist_m, {'numeric'}, {'vector','real','finite','positive'});
validateattributes(rssiDbm, {'numeric'}, {'vector','real','finite'});
dist_m = dist_m(:);
rssiDbm = rssiDbm(:);
assert(numel(dist_m) == numel(rssiDbm) && numel(unique(dist_m)) >= 2, ...
    'Supply equal-sized vectors with at least two distinct distances.');
b = [ones(size(dist_m)), -10*log10(dist_m/d0_m)] \ rssiDbm;
referenceDbm = b(1);
n = b(2);
end
