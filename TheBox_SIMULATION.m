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

% Render environment for fun
hold on;
plotcube([roomLength_m roomWidth_m roomHeight_m],[0  0  0],.0,[0 0 1]);
scatter3(txPos_m(1),txPos_m(2),txPos_m(3));
text(txPos_m(1),txPos_m(2),txPos_m(3)+roomHeight_m*0.1,'tx');
scatter3(rxPos_m(1),rxPos_m(2),rxPos_m(3));
text(rxPos_m(1),rxPos_m(2),rxPos_m(3)+roomHeight_m*0.1,'rx');

% Get RSSI value from line of sight measuremnt

