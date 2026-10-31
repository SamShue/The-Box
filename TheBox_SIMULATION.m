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

distance_m = linspace(min(data.distance_m), max(data.distance_m), ...
    ceil((max(data.distance_m)-min(data.distance_m))/0.005)+1);
predictedDbm = zeros(size(distance_m));
for index = 1:numel(distance_m)
    rxPos_m = [txPos_m(1)+distance_m(index), txPos_m(2), rxHeight_m];
    lengths_m = boxReflectionPaths(room_m, txPos_m, rxPos_m);
    predictedDbm(index) = combinePaths(lengths_m, gamma, frequency_Hz, ...
        referenceDbm, pathExponent, d0_m);
end

% Evaluate at EXACT measurement locations, not interpolated plotting samples.
measuredRayDbm = zeros(size(data.distance_m));
for index = 1:numel(data.distance_m)
    rxPos_m = [txPos_m(1)+data.distance_m(index), txPos_m(2), rxHeight_m];
    lengths_m = boxReflectionPaths(room_m, txPos_m, rxPos_m);
    measuredRayDbm(index) = combinePaths(lengths_m, gamma, frequency_Hz, ...
        referenceDbm, pathExponent, d0_m);
end
% Average physical power across all ten packets (five per orientation).
% Orientation bounds use the same linear-power averaging within each group.
measuredMeanDbm = data.powerMeanDbm;
orientationMagnitude = -data.rotationPowerMeanDbm;
orientationLower = min(orientationMagnitude,[],2);
orientationUpper = max(orientationMagnitude,[],2);
errorDb = measuredRayDbm-measuredMeanDbm;
absoluteErrorDb = abs(errorDb);
rmseDb = sqrt(mean(errorDb.^2));
maeDb = mean(absoluteErrorDb);
biasDb = mean(errorDb);
errorSummary = table(rmseDb,maeDb,biasDb, ...
    'VariableNames',{'RMSE_dB','MAE_dB','Bias_dB'});
disp(errorSummary);

% Positive RSSI magnitudes for display; all calculations use signed dBm.
figure('Name','The Box: orientation-averaged measured versus predicted RSSI', ...
    'Position',[80 80 1200 650]);
rssiAxes = axes;
hold on;
fill([data.distance_m; flipud(data.distance_m)], ...
    [orientationLower; flipud(orientationUpper)], [0.7 0.7 0.7], ...
    'FaceAlpha',0.3,'EdgeColor','none', ...
    'DisplayName','Orientation spread (not a confidence interval)');
plot(distance_m,-predictedDbm,'-','Color',[0 0.447 0.741], ...
    'LineWidth',1.5,'DisplayName','Predicted RSSI (LOS + six reflections)');
plot(data.distance_m,-measuredMeanDbm,'-d','Color',[0.2 0.2 0.2], ...
    'LineWidth',1.5,'MarkerSize',5, ...
    'DisplayName','Measured RSSI (mean linear power, both orientations)');
xlim([min(data.distance_m),max(data.distance_m)]);
xticks(data.distance_m);
xticklabels(compose('%.3f',data.distance_m));
xtickangle(60);
xlabel('Tx-Rx distance (m)');
ylabel('RSSI (-dBm)');
title('Measured versus predicted RSSI');
subtitle({'Lower is stronger; shading spans the two orientation power means.', ...
    'Illustrative reflection coefficients and provisional 1 m calibration.'});
grid on;
legend('Location','best');

figure('Name','The Box: absolute error against orientation-averaged power', ...
    'Position',[100 100 1200 550]);
errorAxes = axes;
bar(data.distance_m,absoluteErrorDb,'FaceColor',[0 0.447 0.741]);
ylim([0 max(1,ceil(max(absoluteErrorDb(isfinite(absoluteErrorDb)))))]);
xticks(data.distance_m);
xticklabels(compose('%.3f',data.distance_m));
xtickangle(60);
xlabel('Tx-Rx distance (m)');
ylabel('Absolute prediction error (dB)');
title('Absolute error against measured mean linear power');
subtitle(sprintf('RMSE %.2f dB | MAE %.2f dB | Bias %+.2f dB (predicted - measured)', ...
    rmseDb,maeDb,biasDb));
grid on;

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
fprintf('Separate 1 m reference: %.2f dBm; assumed carrier: %.6f GHz.\n', ...
    referenceDbm, frequency_Hz/1e9);
