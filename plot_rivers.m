function phandle = plot_rivers()

% plots rivers in canada, from HYDAT database
% to use: figure; plot_rivers()

hold on
filename = './data/Rivers_Lakes/RiversAndLakes_7.5m.shp';
rivermap = shaperead(filename);
for kk = 1:length(rivermap)
    plot(rivermap(kk).X,rivermap(kk).Y,'k')
end