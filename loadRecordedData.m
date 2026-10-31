function data = loadRecordedData(filename, encoding)
%LOADRECORDEDDATA Load historical workbooks without merging the rotations.
% Default ASSUMPTION: hex bytes are positive magnitudes of negative dBm.
% Alternative: 'signed-byte' for two's-complement dBm bytes.
if nargin < 2
    encoding = 'negative-magnitude';
end
encoding = validatestring(encoding, {'negative-magnitude','signed-byte'});
T = readcell(filename);
header = lower(string(T{1,1}));
if contains(header, 'inches')
    scale = 0.0254;
elseif contains(header, 'ft')
    scale = 0.3048;
else
    error('Unrecognized distance unit in workbook header.');
end
isDistance = cellfun(@(x) isnumeric(x) && isscalar(x) && ...
    isfinite(x) && x > 0, T(:,1));
assert(any(isDistance), 'No positive numeric distances found.');
data.distance_m = cell2mat(T(isDistance,1))*scale;
data.samplesDbm = decode(T(isDistance,2:11), encoding);
data.rotationMeanDbm = [mean(data.samplesDbm(:,1:5),2), ...
    mean(data.samplesDbm(:,6:10),2)];
data.rotationPowerMeanDbm = 10*log10([ ...
    mean(10.^(data.samplesDbm(:,1:5)/10),2), ...
    mean(10.^(data.samplesDbm(:,6:10)/10),2)]);
% Mean logarithmic RSSI and mean received power are different statistics.
data.meanDbm = mean(data.samplesDbm,2);
data.powerMeanDbm = 10*log10(mean(10.^(data.samplesDbm/10),2));
referenceRow = find(strcmpi(strtrim(string(T(:,1))), '1 meter'));
assert(isscalar(referenceRow), 'Expected one separately recorded 1 meter row.');
data.referenceSamplesDbm = decode(T(referenceRow,2:11), encoding);
data.referenceDbm = mean(data.referenceSamplesDbm);
data.referenceDistance_m = 1;
data.encoding = encoding;
end

function dbm = decode(cells, encoding)
values = reshape(hex2dec(string(cells)), size(cells));
assert(all(isfinite(values(:)) & values(:) >= 0 & values(:) <= 255), ...
    'RSSI samples must be hexadecimal bytes.');
if strcmp(encoding, 'negative-magnitude')
    dbm = -values;
else
    values(values >= 128) = values(values >= 128)-256;
    dbm = values;
end
end
