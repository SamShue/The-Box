clc;
clear all;
close all;

% Anchor Node @: 
% 25' 9" from wall depression near door (7.8486 m)
% 7' 4.5" from wall on the left while entering room (2.2479 m)
% 
% Room Measurements: 
% 46' 1" from wall depression near door to opposite wall (14.0462 m)
% 41' 8.5" from right to left wall while entering room (12.7127 m)

% Room dimensions
roomLength_m = 12.7127;
roomWidth_m = 14.0462;
roomHeight_m = 3;

% Collected data
% RealRSSI_nDBm = [32, 33, 38, 44, 44, 39, 40, 41]';
% RealDist_m =  [1, 2, 3, 4, 5, 6, 7, 8]';
T = readcell('dataCollection_9_10_22.xlsx');
RealDist_m = cell2mat(T(2:26,1)).*0.3048;
RealRSSI_nDBm = mean(hex2dec(string(T(2:26,2:11))), 2);

% Tx Rx positions
txPos_m = [2.2479,7.8486, 1];
rxPos_m = [3.2479, 7.8486, 1];

% Log-Distance Path Loss Model Parameters
% A = 32;     % 1 meter reference RSSI value
% n = 1.14;    % Path-loss exponent
[A, n] = fitLogDistanceModel(RealRSSI_nDBm, RealDist_m);
% A = hex2dec('21'); % recorded separately

% Transmission settings
Frequency_Hz = 2400000000; % 2.4 GHz
Wavelength_m = physconst('LightSpeed')/Frequency_Hz;

plotBoxModel(roomLength_m, roomWidth_m, roomHeight_m, txPos_m, rxPos_m);

% Get all points associated with cube
ii = 1;
for x = 0:roomLength_m:roomLength_m
    for y = 0:roomWidth_m:roomWidth_m
        for z = 0:roomHeight_m:roomHeight_m
            p(ii,:) = [x, y ,z];
            ii = ii + 1;
        end
    end
end

% Get points that define each wall of the cube
% Each wall is defined by the collection of points that contain the same
% value along x, y, or z, so long as the cube is right angle aligned with
% the coordinate plane it exists within
w1 = p(p(:,1) == 0,:);
w2 = p(p(:,1) == roomLength_m,:);
w3 = p(p(:,2) == 0,:);
w4 = p(p(:,2) == roomWidth_m,:);
w5 = p(p(:,3) == 0,:);
w6 = p(p(:,3) == roomHeight_m,:);

% put all wall cordinates into cell array
w = {w1, w2, w3, w4, w5, w6};
% w = {w3, w4, w5, w6};
w = {w4}

distIncrements = 1:0.01:8;
positionIncrements = txPos_m(1) + distIncrements;
losRssi = A + 10*n*log10(distIncrements);
losPhase = real(exp(-j*((2*pi*Frequency_Hz)/physconst('LightSpeed'))*distIncrements));

for kk = 1:length(positionIncrements)
    rxPos_m(1) = positionIncrements(kk);

    % Get RSSI value from left wall reflection
    % Get the vectors that define the plane of the first wall
    % Point Projection onto Plane Reference:
    % https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d
    
    % % Define plane by normal vector
    % a = w1(2,:) - w1(1,:);
    % b = w1(3,:) - w1(1,:);
    % n = cross(a,b);
    % n = n/norm(n);
    % % Project point onto plane
    % v = txPos_m - w1(1,:);
    % dist = dot(v,n);
    % projected_point = txPos_m - dist.*n;
    % % Plot projected point
    % scatter3(projected_point(1), projected_point(2), projected_point(3));

    clf;
    hold on;
    
    % Loop thru each wall projecting points
    for ii = 1:length(w)
        wallPoints = cell2mat(w(ii));

        a = wallPoints(2,:) - wallPoints(1,:);
        b = wallPoints(3,:) - wallPoints(1,:);
        normalVec = cross(a,b);
        normalVec = normalVec/norm(normalVec);
        % Project tx point onto plane
        midpoint = (txPos_m + rxPos_m)/2;
        v = midpoint - wallPoints(1,:);
        dist = dot(v,normalVec);
        
        projectedPoint = midpoint - dist.*normalVec;

        
        % save multipath travel distance
        dist = norm(txPos_m - projectedPoint) + norm(projectedPoint - rxPos_m);
%         % divide and get remainder and scale by 2pi
%         phase(ii) = rem(dist(ii)/Wavelength_m, 1)*2*pi;
%         % use phase to get signal interference
%         interference(ii) = cos(phase(ii));
        
        phase = real(exp(-j*((2*pi*Frequency_Hz)/physconst('LightSpeed'))*dist));
        phaseDifference = losPhase(kk) - phase;

        % use phase difference to get signal interference
        interference = phaseDifference/2.0;

        % reflectance coefficient
        reflectance = 0.03;

        % get RSSI
        rssi(ii) = (A + 10*n*log10(dist)).*interference.*reflectance;

        % Plot line connecting projection points
        hold on;
        xpoints = [txPos_m(1), projectedPoint(1), rxPos_m(1)];
        ypoints = [txPos_m(2), projectedPoint(2), rxPos_m(2)];
        zpoints = [txPos_m(3), projectedPoint(3), rxPos_m(3)];
        plot3(xpoints, ypoints, zpoints)
        % Plot projected point
%         scatter3(projectedPoint(1), projectedPoint(2), projectedPoint(3));
%         scatter3(midpoint(1), midpoint(2), midpoint(3));
        text(projectedPoint(1),projectedPoint(2), projectedPoint(3)+roomHeight_m*0.1, int2str(ii));
    end
    plotBoxModel(roomLength_m, roomWidth_m, roomHeight_m, txPos_m, rxPos_m);
    pause(0.1);
    
    rssiMultipath(kk) = sum(rssi);
end

figure()
ActualRssi = losRssi + rssiMultipath;
plot(distIncrements, ActualRssi)
hold on;
% Plot measured data points
scatter(RealDist_m, RealRSSI_nDBm);
% Plot calibrated model
distances = 1:0.1:8;
plot(distances, (A + 10*n*log10(distances)))


% figure()
% hold on;
% plot(distIncrements, losRssi)
% plot(distIncrements, rssiMultipath)