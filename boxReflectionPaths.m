function [lengths_m, points_m] = boxReflectionPaths(room_m, tx_m, rx_m)
%BOXREFLECTIONPATHS LOS plus six single-bounce paths in an axis-aligned room.
% Order: LOS, x=0, x=L, y=0, y=W, z=0, z=H.
% points_m has one reflection point per surface (six rows).
validateattributes(room_m, {'numeric'}, {'real','finite','positive','numel',3});
validateattributes(tx_m, {'numeric'}, {'real','finite','numel',3});
validateattributes(rx_m, {'numeric'}, {'real','finite','numel',3});
room_m = reshape(room_m,1,3);
tx_m = reshape(tx_m,1,3);
rx_m = reshape(rx_m,1,3);
assert(all(tx_m > 0 & tx_m < room_m) && all(rx_m > 0 & rx_m < room_m), ...
    'Tx and Rx must be strictly inside the room.');
assert(norm(tx_m-rx_m) > 0, 'Tx and Rx must be distinct.');
lengths_m = zeros(1,7);
lengths_m(1) = norm(tx_m-rx_m);
points_m = zeros(6,3);
for surface = 1:6
    axisIndex = ceil(surface/2);
    wall = 0;
    if mod(surface,2) == 0
        wall = room_m(axisIndex);
    end
    imageTx = tx_m;
    imageTx(axisIndex) = 2*wall-tx_m(axisIndex);
    direction = rx_m-imageTx;
    fraction = (wall-imageTx(axisIndex))/direction(axisIndex);
    point = imageTx+fraction*direction;
    tolerance = 1e-10*max(room_m);
    assert(fraction > 0 && fraction < 1 && ...
        all(point >= -tolerance & point <= room_m+tolerance), ...
        'Reflection point lies outside the finite surface.');
    points_m(surface,:) = point;
    lengths_m(surface+1) = norm(direction);
end
end
