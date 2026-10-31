% Dense October sweep. Keep antenna rotations separate.
repoDirectory = fileparts(mfilename('fullpath'));
data = loadRecordedData(fullfile(repoDirectory, 'dataCollection_10_9_22.xlsx'), ...
    'negative-magnitude');
figure('Name','Recorded RSSI by rotation');
plot(data.distance_m, data.rotationMeanDbm(:,1), 'o-', 'DisplayName','Rotation 1 mean');
hold on;
plot(data.distance_m, data.rotationMeanDbm(:,2), 's-', 'DisplayName','Rotation 2 mean');
xlabel('Tx-Rx distance (m)');
ylabel('Received power (dBm; assumed RSSI decoding)');
grid on;
legend('Location','best');
