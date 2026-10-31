function testBoxModel
%TESTBOXMODEL Physical invariants and regression checks; base MATLAB only.
f = 2.4e9;
reference = -33;
% Direct path obeys inverse-square power; equal coherent fields give +6 dB.
p = combinePaths(2, 1, f, reference, 2, 1);
assert(abs(p-(reference-20*log10(2))) < 1e-10);
p = combinePaths([1,1], [1,1], f, reference, 2, 1);
assert(abs(p-(reference+20*log10(2))) < 1e-10);
p = combinePaths([1,1], [1,-1], f, reference, 2, 1);
assert(isinf(p) && p < 0);
% Known relative propagation phase, amplitude imbalance, and complex phase.
lambda = 299792458/f;
[p,h] = combinePaths([1,1+lambda/2], [1,1], f, reference, 2, 1);
expected = 1-1/(1+lambda/2);
assert(abs(h-expected) < 1e-12 && abs(p-(reference+20*log10(expected))) < 1e-10);
[p,h] = combinePaths([1,1], [1,1i], f, reference, 2, 1);
assert(abs(h-(1+1i)) < 1e-12 && abs(p-(reference+10*log10(2))) < 1e-10);
p1 = combinePaths([1,2,3], [1,.2i,-.1], f, reference, 2, 1);
p2 = combinePaths([3,1,2], [-.1,1,.2i]*exp(1i*.7), f, reference, 2, 1);
assert(abs(p1-p2) < 1e-10); % Reordering and common phase leave power unchanged.
p = combinePaths([1,5,10], [1,0,0], f, reference, 2, 1);
assert(abs(p-reference) < 1e-10);

% Unequal distances to planes and unequal heights: midpoint projection fails.
room = [12,14,3]; tx = [2,3,.7]; rx = [8,9,2.1];
[lengths,points] = boxReflectionPaths(room,tx,rx);
assert(abs(lengths(1)-norm(tx-rx)) < 1e-12);
for surface = 1:6
    dimension = ceil(surface/2);
    wall = 0;
    if mod(surface,2) == 0
        wall = room(dimension);
    end
    q = points(surface,:);
    assert(abs(q(dimension)-wall) < 1e-12);
    assert(all(q >= 0 & q <= room));
    incoming = (q-tx)/norm(q-tx);
    outgoing = (rx-q)/norm(rx-q);
    normal = zeros(1,3); normal(dimension) = 1;
    assert(norm(outgoing-(incoming-2*dot(incoming,normal)*normal)) < 1e-12);
    assert(abs(lengths(surface+1)-(norm(tx-q)+norm(q-rx))) < 1e-12);
end
% Reciprocity and original symmetric midpoint construction.
[reverseLengths,reversePoints] = boxReflectionPaths(room,rx,tx);
assert(max(abs(reverseLengths-lengths)) < 1e-12);
assert(max(abs(reversePoints(:)-points(:))) < 1e-12);
rx = [8,3,.7];
[~,points] = boxReflectionPaths(room,tx,rx);
for surface = 1:6
    dimension = ceil(surface/2);
    midpoint = (tx+rx)/2;
    midpoint(dimension) = 0;
    if mod(surface,2) == 0
        midpoint(dimension) = room(dimension);
    end
    assert(norm(points(surface,:)-midpoint) < 1e-12);
end
expectFailure(@() boxReflectionPaths(room,[0,1,1],rx));
expectFailure(@() boxReflectionPaths(room,tx,tx));
expectFailure(@() combinePaths([1,2],1,f,reference,2,1));

% Recover known signed-dBm parameters, including non-unit reference distance.
d = [1;2;4;8];
[A,n] = fitLogDistanceModel(-40-17*log10(d/2),d,2);
assert(abs(A+40) < 1e-10 && abs(n-1.7) < 1e-10);
expectFailure(@() fitLogDistanceModel([-30,-40],[1,1]));
repo = fileparts(mfilename('fullpath'));
a = loadRecordedData(fullfile(repo,'dataCollection_9_10_22.xlsx'));
b = loadRecordedData(fullfile(repo,'dataCollection_10_9_22.xlsx'));
assert(isequal(size(a.samplesDbm),[25,10]) && isequal(size(b.samplesDbm),[37,10]));
assert(max(abs(a.distance_m([1,end])-[.9144;8.2296])) < 1e-12);
assert(max(abs(b.distance_m([1,end])-[3.9624;4.8768])) < 1e-12);
assert(a.referenceDbm == -33 && b.referenceDbm == -34);
assert(all(a.samplesDbm(:) < 0) && all(b.samplesDbm(:) < 0));
assert(abs(max(abs(diff(a.rotationMeanDbm,1,2)))-11) < 1e-10);
assert(all(a.powerMeanDbm >= a.meanDbm-1e-12));
assert(isequal(size(a.rotationPowerMeanDbm),[25,2]));
combined = 10*log10(mean(10.^(a.rotationPowerMeanDbm/10),2));
assert(max(abs(combined-a.powerMeanDbm)) < 1e-12);
assert(all(a.powerMeanDbm >= min(a.rotationPowerMeanDbm,[],2)-1e-12));
assert(all(a.powerMeanDbm <= max(a.rotationPowerMeanDbm,[],2)+1e-12));
[A,n] = fitLogDistanceModel(a.meanDbm,a.distance_m);
assert(abs(A+38.638603) < 1e-6 && abs(n-1.257285) < 1e-6);
% Exercise the alternative byte interpretation against known workbook bytes.
signed = loadRecordedData(fullfile(repo,'dataCollection_9_10_22.xlsx'),'signed-byte');
assert(isequal(signed.samplesDbm,-a.samplesDbm));
fprintf('All Box model checks passed.\n');
end

function expectFailure(action)
failed = false;
try
    action();
catch
    failed = true;
end
assert(failed, 'Expected invalid input to be rejected.');
end
