function phandle = plot_AB()

% plots the border of alberta
% to use: figure; plot_AB()

hold on

filename = './data/Canada_Borders/PROVINCE.SHP';
provshape = shaperead(filename,'UseGeoCoords',true);
border = plot(provshape(1).Lon, provshape(1).Lat,'k'); %AB
border.LineWidth = 2;
xlim([-120.01 -109.99])
ylim([48.99 60.01])