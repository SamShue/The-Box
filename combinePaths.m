function [receivedDbm, h] = combinePaths(lengths_m, gamma, frequency_Hz, referenceDbm, n, d0_m)
%COMBINEPATHS Coherently sum scalar narrowband path amplitudes.
% gamma holds COMPLEX FIELD coefficients; include 1 for the LOS path.
% referenceDbm is a single-path reference at d0_m. n=2 is free-space spreading.
validateattributes(lengths_m, {'numeric'}, {'vector','real','finite','positive'});
validateattributes(gamma, {'numeric'}, {'vector','finite'});
validateattributes(frequency_Hz, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(referenceDbm, {'numeric'}, {'scalar','real','finite'});
validateattributes(n, {'numeric'}, {'scalar','real','finite','positive'});
validateattributes(d0_m, {'numeric'}, {'scalar','real','finite','positive'});
assert(numel(lengths_m) == numel(gamma), 'Each path needs one coefficient.');
lengths_m = lengths_m(:).';
gamma = gamma(:).';
k = 2*pi*frequency_Hz/299792458;
% Remove common propagation phase. h is relative to the first path's phase;
% received power is invariant to this change of reference.
h = sum(gamma .* (d0_m./lengths_m).^(n/2) .* ...
    exp(-1i*k*(lengths_m-lengths_m(1))));
receivedDbm = referenceDbm+20*log10(abs(h)); % Exact cancellation gives -Inf.
end
