function plotBoxModel(l, w, h, tx, rx)
%PLOTBOX Summary of this function goes here
%   Detailed explanation goes here

    % Render environment for fun
    hold on;
    plotcube([l w h],[0  0  0],.0,[0 0 1]);
    scatter3(tx(1),tx(2),tx(3));
    text(tx(1),tx(2),tx(3)+h*0.1,'tx');
    scatter3(rx(1),rx(2),rx(3));
    text(rx(1),rx(2),rx(3)+h*0.1,'rx');

end

