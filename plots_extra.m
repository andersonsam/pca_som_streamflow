%% pca results

figure
for kk = 1:4
    subplot(2,4,2*kk-1)
    plot(eigvecsF(:,kk),'k')
    title(sprintf('Flow Mode #%d \n Frac Var = %0.1f%%',kk,100*fracVarF(kk)),'fontsize',14)
    ylim([min(min(eigvecsF(:,1:4))), max(max(eigvecsF(:,1:4)))])
    
    subplot(2,4,2*kk)
    plot(PCsF(:,kk))
    title(sprintf('PCs Mode #%d',kk),'fontsize',14)
end
    
figure
for kk = 1:2
    subplot(1,2,kk)
    plot(eigvecsF(:,kk),'k','LineWidth',2)
    title(sprintf('Flow Mode #%d \n Percentage Variance = %0.1f%%',kk,100*fracVarF(kk)),'fontsize',14)
    xlim([1,31])
    ylim([min(min(eigvecsF(:,1:2))), max(max(eigvecsF(:,1:2)))])
    xlabel('Day of August')
    ylabel('Mode Strength')
    set(gca,'FontSize',30)
end

%% SOM clusters in space

figure('Position',[0,0,600,600])
plot_rivers()
feval(plot_prov)
stationplot = scatter(stationLon(dataInds),stationLat(dataInds),150,'filled','LineWidth',2,'MarkerEdgeColor','k');
set(stationplot,'CData',cStations);

cmap = cbmus;
colormap(cmap);
cbstationvals = colorbar; 
cbstationvals.Location = 'eastoutside';
cbstationvals.TickLabels = [];
cbstationvals.Ticks = [];
cbstationvals.Ticks = linspace(0,1,length(cbmus));
cbstationvals.TickLabels = 1:length(cbmus);

titlestr = sprintf('Clustering of August PC Trajectory');
title(titlestr)
xlabel('Longitude')
ylabel('Latitude')
set(gca,'FontSize',24)

%% plot % glaciers vs PC1

figure
s = scatter(meanPC1,stdPC1,100,'filled','LineWidth',2,'MarkerEdgeColor','k');
s.CData = cStations;
xlabel('Mean PC1')
ylabel('Standard Deviation PC1')
title('Mean PC1 vs Standard Deviation PC1')
set(gca,'FontSize',14)

%% Visualize MLR model coefficients

figure %make boxplots of a values for each predictand
for kk = 1:length(predictandsAll(1,:)) %for each predictand
    subplot(1,2,kk)
    bpdata = horzcat(models1.a{:,kk})';
    boxplot(bpdata(:,2:end))
    title(predictandsStr{kk})
end  

figure %make kde of a values for each predictand
for kk = 1:length(predictandsAll(1,:)) %for each predictand
    subplot(1,2,kk), hold on
    kdedata = horzcat(models1.a{:,kk});
    for jj = 2:length(kdedata(:,kk)) %for each predictor, make kde; don't include a(1) values (ie: intercepts)
        [kde,XI] = ksdensity(kdedata(jj,:));
        p = plot(XI,kde);
        p.LineWidth = 3;
    end
    title(predictandsStr{kk})
end

%% modelled vs measured mean and std PC1

figure
for kk = 1:length(predictandsAll(1,:)) %for each predictand
    subplot(1,2,kk), hold on
    scatter(predictandsAll(:,kk),predictands_MLR(:,kk),100,cStations,'filled','LineWidth',2,'MarkerEdgeColor','k')
    refline(1,0)
    xlabel('Measured')
    ylabel('Modelled')
    title(predictandsStr{kk})
end

figure
for kk = 1:length(predictandsAll(1,:)) %for each predictand
    subplot(1,2,kk), hold on
    scatter(predictandsAll(stationsWithoutGlaciers,kk),predictands_MLR(stationsWithoutGlaciers,kk),150,cStations(stationsWithoutGlaciers,:),'filled','LineWidth',2,'MarkerEdgeColor','k')
    refline(1,0)
    xlabel('Measured')
    ylabel('Modelled')
    title(predictandsStr{kk})
    set(gca,'fontsize',36)
    if kk==1
        xlim([-3,3])
        ylim([-3,3])
    end
end
    
%% plot change in meanQ and stdQ

figure, hold on
scatter(station_Q_mean_glaciers,station_Q_mean_MLR,100,cStations,'filled')
scatter(intake_Q_mean_glaciers,intake_Q_mean_MLR,100,intake_dstdQ,'filled')
scatter(dam_Q_mean_glaciers,dam_Q_mean_MLR,100,dam_dstdQ,'filled','Marker','^')
refline(1,0)
colorbar
xlabel('Glaciers')
ylabel('MLR')

figure, hold on
scatter(station_dmeanQ,station_dstdQ,100,cStations,'filled')
scatter(intake_dmeanQ,intake_dstdQ,75,'ks')
scatter(dam_dmeanQ,dam_dstdQ,'k^')
xlabel('$\Delta Q_{mean}$')
ylabel('$\Delta\sigma_{Q}$')
set(gca,'FontSize',18)