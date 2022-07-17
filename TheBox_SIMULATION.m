clc;
clear all;
close all;

% Anchor Node @: 
% 25' 9" from wall depression near door
% 7' 4.5" from wall on the left while entering room
% 
% Room Measurements: 
% 46' 1" from wall depression near door to opposite wall
% 41' 8.5" from right to left wall while entering room

% Room dimensions
roomLength_m = 10;
roomWidth_m = 20;
roomHeight_m = 3;

% Tx Rx positions
txPos_m = [5,5,1];
rxPos_m = [5,15,1];

% Transmission settings
Frequency_Hz = 2400000; % 2.4 GHz
Wavelength_m = physconst('LightSpeed')/Frequency_Hz;

% Log-Distance Path Loss Model Parameters
A = 45;     % 1 meter reference RSSI value
n = 2.1;    % Path-loss exponent
d = norm(txPos_m - rxPos_m);    % distance between

% Render environment for fun
hold on;
plotcube([roomLength_m roomWidth_m roomHeight_m],[0  0  0],.0,[0 0 1]);
scatter3(txPos_m(1),txPos_m(2),txPos_m(3));
text(txPos_m(1),txPos_m(2),txPos_m(3)+roomHeight_m*0.1,'tx');
scatter3(rxPos_m(1),rxPos_m(2),rxPos_m(3));
text(rxPos_m(1),rxPos_m(2),rxPos_m(3)+roomHeight_m*0.1,'rx');

% Get RSSI value from line of sight measuremnt
RSSI_los = A + 10*n*log(d);

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

w = {w1, w2, w3, w4, w5, w6};

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

% Loop thru each wall projecting points
for ii = 1:length(w)
    points = cell2mat(w(ii));
    a = points(2,:) - points(1,:);
    b = points(3,:) - points(1,:);
    n = cross(a,b);
    n = n/norm(n);
    % Project tx point onto plane
    v = txPos_m - points(1,:);
    dist1 = dot(v,n);
    % Project rx point onto plane
    v = rxPos_m - points(1,:);
    dist2 = dot(v,n);
    % Save the closer point
    if(dist1 < dist2)
        projectedPoint = txPos_m - dist1.*n;
    else
        projectedPoint = rxPos_m - dist2.*n;
    end
    % Plot line connecting projection points
    xpoints = [txPos_m(1), projectedPoint(1), rxPos_m(1)];
    ypoints = [txPos_m(2), projectedPoint(2), rxPos_m(2)];
    zpoints = [txPos_m(3), projectedPoint(3), rxPos_m(3)];
    plot3(xpoints, ypoints, zpoints)
    % Plot projected point
    scatter3(projectedPoint(1), projectedPoint(2), projectedPoint(3));
end
RSSI_left_wall = A + 10*n*log(d);
