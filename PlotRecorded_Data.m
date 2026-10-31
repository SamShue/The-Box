% Dense October sweep: mean linear power with observed orientation spread.
repoDirectory = fileparts(mfilename('fullpath'));
data = loadRecordedData(fullfile(repoDirectory, 'dataCollection_10_9_22.xlsx'), ...
    'negative-magnitude');
figure('Name','Recorded RSSI averaged over orientations');
orientationMagnitude = -data.rotationPowerMeanDbm;
lower = min(orientationMagnitude,[],2);
upper = max(orientationMagnitude,[],2);
fill([data.distance_m;flipud(data.distance_m)], [lower;flipud(upper)], ...
    [0.7 0.7 0.7],'FaceAlpha',0.3,'EdgeColor','none', ...
    'DisplayName','Orientation spread (not a confidence interval)');
hold on;
plot(data.distance_m,-data.powerMeanDbm,'-d','Color',[0.2 0.2 0.2], ...
    'LineWidth',1.5,'DisplayName','Mean linear power, both orientations');
xlabel('Tx-Rx distance (m)');
ylabel('RSSI (-dBm)');
title('Measured RSSI averaged over the two sampled orientations');
subtitle('Lower is stronger; shading spans the two orientation power means.');
grid on;
legend('Location','best');
