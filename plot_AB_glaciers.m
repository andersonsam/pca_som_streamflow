function phandle = plot_AB_glaciers(size,colour)

% Syntax: phandle = plot_AB_glaciers(size,colour)
% size == size of circle, double
% colour == either rgb triplet, or colour like 'k'
% plots the mean lat/lon of each glacier in Alberta; faster than plot_AB_rgi()

load('AB_rgi.mat')

hold on
phandle = scatter(glacierLon,glacierLat,size,colour,'filled');