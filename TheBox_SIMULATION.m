% LOS plus six single-reflection rays. See BACKGROUND.md for assumptions.
clear;
close all;
repoDirectory = fileparts(mfilename('fullpath'));

% Historical September room coordinates (meters).
room_m = [12.7127, 14.0462, 3];
txPos_m = [2.2479, 7.8486, 1];
rxHeight_m = 1;
data = loadRecordedData(fullfile(repoDirectory, 'dataCollection_9_10_22.xlsx'), ...
    'negative-magnitude'); % Assumption: verify original radio encoding.

% Set actual channel center frequency when known; 2.400 GHz is assumed.
frequency_Hz = 2.4e9;
pathExponent = 2; % Per-ray free-space spreading, not a fitted room trend.
% Illustrative coefficients, NOT measured material properties.
% Each has 3% power reflectance and a pi-radian field phase reversal.
% Order: x=0, x=L, y=0, y=W, z=0, z=H. Edit magnitude AND phase per surface.
powerReflectance = 0.03*ones(1,6);
reflectionPhase_rad = pi*ones(1,6);
gamma = [1, sqrt(powerReflectance).*exp(1i*reflectionPhase_rad)];
% The separate 1 m reading can contain multipath: this scale is provisional.
referenceDbm = data.referenceDbm;
d0_m = data.referenceDistance_m;

% Fitted measured trend is a comparison, not the ray model calibration.
[trendReferenceDbm, trendExponent] = fitLogDistanceModel(data.meanDbm, data.distance_m);
distance_m = linspace(min(data.distance_m), max(data.distance_m), ...
    ceil((max(data.distance_m)-min(data.distance_m))/0.005)+1);
predictedDbm = zeros(size(distance_m));
for index = 1:numel(distance_m)
    rxPos_m = [txPos_m(1)+distance_m(index), txPos_m(2), rxHeight_m];
    lengths_m = boxReflectionPaths(room_m, txPos_m, rxPos_m);
    predictedDbm(index) = combinePaths(lengths_m, gamma, frequency_Hz, ...
        referenceDbm, pathExponent, d0_m);
end
losDbm = referenceDbm-10*pathExponent*log10(distance_m/d0_m);
trendDbm = trendReferenceDbm-10*trendExponent*log10(distance_m/d0_m);

figure('Name','The Box: predicted and measured received power');
plot(distance_m, predictedDbm, 'DisplayName','LOS + six reflections');
hold on;
plot(distance_m, losDbm, '--', 'DisplayName','Direct path only');
plot(distance_m, trendDbm, ':', 'DisplayName','Fitted measured trend (comparison only)');
scatter(data.distance_m, data.rotationMeanDbm(:,1), 25, 'DisplayName','Rotation 1 mean');
scatter(data.distance_m, data.rotationMeanDbm(:,2), 25, 'DisplayName','Rotation 2 mean');
xlabel('Tx-Rx distance (m)');
ylabel('Received power (dBm; assumed RSSI decoding)');
title('Illustrative surface coefficients and provisional 1 m calibration');
grid on;
legend('Location','best');

% Show one receiver position without slowing the calculation with animation.
rxPos_m = [txPos_m(1)+median(data.distance_m), txPos_m(2), rxHeight_m];
[~, points_m] = boxReflectionPaths(room_m, txPos_m, rxPos_m);
figure('Name','Single-bounce image-source geometry');
plotBoxModel(room_m(1), room_m(2), room_m(3), txPos_m, rxPos_m);
plot3([txPos_m(1),rxPos_m(1)], [txPos_m(2),rxPos_m(2)], ...
    [txPos_m(3),rxPos_m(3)], 'k--');
for surface = 1:6
    vertices = [txPos_m; points_m(surface,:); rxPos_m];
    plot3(vertices(:,1), vertices(:,2), vertices(:,3));
    text(points_m(surface,1),points_m(surface,2),points_m(surface,3), ...
        sprintf(' %d',surface));
end
fprintf('Measured trend: reference = %.2f dBm, n = %.3f (not used for rays).\n', ...
    trendReferenceDbm, trendExponent);
fprintf('Separate 1 m reference: %.2f dBm; assumed carrier: %.6f GHz.\n', ...
    referenceDbm, frequency_Hz/1e9);
