% Collected data
T = readcell('dataCollection_10_9_22.xlsx');
RealDist_m = cell2mat(T(2:38,1)).*0.0254;
RealRSSI_nDBm = mean(hex2dec(string(T(2:38,2:11))), 2);

% Plot measured data points
scatter(RealDist_m, RealRSSI_nDBm); hold on;
plot(RealDist_m, RealRSSI_nDBm);