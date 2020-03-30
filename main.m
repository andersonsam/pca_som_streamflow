%
% Sam Anderson
% sanderson@eoas.ubc.ca
%
% This code reproduces results from the paper "Identification of local water resource vulnerability to rapid deglaciation"
%
% Note: Requires somtoolbox (http://www.cis.hut.fi/somtoolbox/) for
% self-organizing maps functions

clear
clc

data = load_data(); %loads streamflow data, climate data, stream gauge data

%unpack data, which has two layers: outer layer is a 'group' (i.e. stream gauge variables), inner layer is the variables of that group (i.e. station elevation)
fields_outer = fieldnames(data);
for kk = 1:length(fields_outer)
    
    field_outer = fields_outer{kk};
    fields_inner = fieldnames(data.(field_outer));
    
    for jj = 1:length(fields_inner)
        
        field_inner = fields_inner{jj};
        varname = genvarname(field_inner);
        eval([varname '= data.(field_outer).(field_inner);']);
        
    end
    
end

%% smooth streamflow data

smooth_window = 30;  
for kk = 1:length(flow_years_norm(:,1)) %for each station/year, smooth
    flow_years_norm_smooth(kk,:) = movmean(flow_years_norm(kk,:), smooth_window);
end

%% Do PCA on August flow data

summerInds = 213:243; %day numbers of August
pcaInputF = flow_years_norm_smooth(:,summerInds); %input only august to PCA (F for 'flow')
[eigvecsF, PCsF, eigvalsF] = pca(pcaInputF); %do PCA on normalized, smoothed August streamflow
fracVarF = eigvalsF/sum(eigvalsF); %calculate the fraction of variance explained by each PCA mode
  
%% Convert PCs into complex coordinates to cluster using SOM

complex_coord = PCsF(:,1) + 1i*PCsF(:,2);
complex_coord = reshape(complex_coord,[],length(stationLat))';

%% Cluster using a SOM

somInput = complex_coord;
[en,nx_som,ny_som] = get_numberOfPatterns(somInput); %determine the number of nodes/size of map to use
[sM, sT, sMap, bmus] = do_SOM(somInput, nx_som, ny_som); %make SOM with this map size

%determine RGB colours to use for...
cbmus = som_colorcode(sM,'rgb2'); %... each cluster
for kk = 1:length(somInput(:,1))
    cStations(kk,:) = cbmus(bmus(kk),:); %... and each station
end
    
%% Name PC1/2, meanPC1/2, and stdPC1/2 for ease of use

PC1 = real(complex_coord);
PC2 = imag(complex_coord);
meanPC1 = mean(real(complex_coord),2);
meanPC2 = mean(imag(complex_coord),2);
stdPC1 = std(real(complex_coord),0,2);
stdPC2 = std(imag(complex_coord),0,2);

%% Calculate min seasonal flow / max seasonal flow

stationMinFlow = nanmin(all_flowseason,[],2);
stationMaxFlow = nanmax(all_flowseason,[],2);
stationPercentMin = stationMinFlow./stationMaxFlow;
stationPercentMinLog = log10(stationPercentMin);
stationPercentMinLog(isinf(stationPercentMinLog)) = min(stationPercentMinLog(~isinf(stationPercentMinLog)));

%% prep data for regression

%interpolate intake/dam fields 
intakePercentMinLog = griddata(stationLon,stationLat,stationPercentMinLog,intakeLon,intakeLat);
damPercentMinLog = griddata(stationLon,stationLat,stationPercentMinLog,damLon,damLat);

lonERA = double(lonERA);
latERA = double(latERA);
for kk = 1:length(lonERA) %for each lat
    for jj = 1:length(latERA) %for each lon

        %inverse distance weighting, power 2
        dlat = latERA(jj) - stationLat;
        dlon = lonERA(kk) - stationLon;
        d2 = dlat.^2 + dlon.^2;
        w = 1./sqrt(d2);
        w = w.^2;        
        provincePercentMinLog(kk,jj) = sum(w.*stationPercentMinLog) / sum(w);
    
    end
end    

stationElevation = double(stationElevation);
intakeElevation = double(intakeElevation);
damElevation = double(damElevation);

%prep data for regression
data = [stationTmean, stationPmean, stationPJJAmean, stationTJJAmean, stationEJJAmean, stationElevation, log10(stationDrainageArea), stationStreamOrder, stationPercentMinLog];
intakeData = [intakeTmean, intakePmean, intakePJJAmean, intakeTJJAmean, intakeEJJAmean, intakeElevation, log10(intakeDrainageArea), intakeStreamOrder, intakePercentMinLog];
damData = [damTmean, damPmean, damPJJAmean, damTJJAmean, damEJJAmean, damElevation, log10(damDrainageArea), damStreamOrder, damPercentMinLog];

%normalize all variables
data_norm = zeros(size(data));
damData_norm = zeros(size(damData));
intakeData_norm = zeros(size(intakeData));
for kk = 1:length(data(1,:))
    data_norm(:,kk) = (data(:,kk) - mean(data(:,kk)))/std(data(:,kk));  
    intakeData_norm(:,kk) = (intakeData(:,kk) - nanmean(data(:,kk)))/nanstd(data(:,kk));
    damData_norm(:,kk) = (damData(:,kk) - nanmean(data(:,kk)))/nanstd(data(:,kk));
end


%% Use stepwise regression to see which columns in data are statistically significant in predicting meanPC1 and stdPC1

inds = stationsWithoutGlaciers;
predictandsAll = [meanPC1,stdPC1];
predictandsStr = {'Mean PC1','Std PC1'};

for kk = 1:length(predictandsAll(1,:)) %for each predictand
    
    predictand = predictandsAll(stationsWithoutGlaciers,kk);
    predictors = data_norm(stationsWithoutGlaciers,:);

    [a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(predictors,predictand,'penter',0.01,'premove',0.01); %stepwise to calculate which predictors to keep
    predictorsKeep{kk,1} = find(INMODEL==1)';
    
end

%% remove predictors which have high VIF (VIF > vif_threshold)

vif_threshold = 5;
predictorsRemoved = cell(1,2); %keep track of predictors removed in each model

disp('VIF test')
for model = 1:length(predictandsAll(1,:)) %for each model (meanPC1, stdPC1)
    
    disp(sprintf('\n Model #%d',model))
    
    keep = predictorsKeep{model};
    
    VIF = vif(data_norm(:,keep))
    [maxVIF, indMaxVIF] = max(VIF);
    
    jj = 1;
    while maxVIF > vif_threshold %if vif too high, remove predictor
        disp(sprintf('\t \t Removed predictor %d, VIF = %0.2f',predictorsKeep{model}(indMaxVIF),maxVIF))
        predictorsRemoved{model}(jj) = predictorsKeep{model}(indMaxVIF);
        predictorsKeep{model}(indMaxVIF) = [];
        keep = predictorsKeep{model};
        VIF = vif(data_norm(:,keep));
        [maxVIF, indMaxVIF] = max(VIF);
        jj = jj+1;
    end
    
end

%% bootstrap predictors + loop through all combinations of possible predictors to see which are most often used 

clear models r_best_ThisNumPredictors ind_best_ThisNumPredictors predictors_best_ThisNumPredictors

numIterations = 100;
fracCalib = 1;
nCalib = floor(fracCalib*length(stationLon(stationsWithoutGlaciers)));

for kk = 1:length(predictandsAll(1,:)) %for each predictand
    
    clear r_best_ThisNumPredictors predictors_best_ThisNumPredictors
    
    disp(['Predictand ',num2str(kk),'/',num2str(length(predictandsAll(1,:)))])
    
    predictands = predictandsAll(:,kk);
    
    for nn = 1:numIterations %for each bootstrapping iteration
        
        disp(['    Iteration ',num2str(nn),'/',num2str(numIterations)])
        
        indsCalib = randi([1,length(stationsWithoutGlaciers)],[nCalib,1]);
        indsValid = 1:length(stationsWithoutGlaciers);
        
        predictandsCalib = predictands(stationsWithoutGlaciers(indsCalib),:);
        predictandsValid = predictands(stationsWithoutGlaciers(indsValid),:);
        
        for ii = 1:length(data_norm(1,:)) - length(predictorsRemoved{kk}) %for each number of predictors
            
            if length(predictorsRemoved{kk})>0 %if some predictor has been removed, don't use it here
                predictorInds = nchoosek(1:length(data_norm(1,:))-length(predictorsRemoved{kk}),ii);
                dummy = 1:length(data_norm(1,:));
                dummy(predictorsRemoved{kk}) = [];
                if ii==1
                    predictors = dummy(predictorInds)';
                else
                    predictors = dummy(predictorInds);
                end
            else %if no predictors removed, use all of them here
                predictors = nchoosek(1:length(data_norm(1,:)),ii);
            end
            clear r_test
            clear p_test
            
            for jj = 1:length(predictors(:,1)) %for each combination of this number of predictors
                
                predictorsCalib = data_norm(stationsWithoutGlaciers(indsCalib),predictors(jj,:));
                predictorsValid = data_norm(stationsWithoutGlaciers(indsValid),predictors(jj,:));
                
                XCalib = [ones(size(predictorsCalib(:,1))), predictorsCalib];
                YCalib = predictandsCalib;
                
                XValid = [ones(size(predictorsValid(:,1))), predictorsValid];
                YValid = predictandsValid;
                
                a_test = XCalib\YCalib;
                
                Y_reg = XValid*a_test;
                
                [r0_test, p0_test] = corrcoef(Y_reg,YValid);
                r_test(jj,1) = r0_test(1,2);
                p_test(jj,1) = p0_test(1,2);
                
            end
            
            [r_best_ThisNumPredictors(ii,1), ind_best_ThisNumPredictors(ii,1)] = max(r_test);
            predictors_best_ThisNumPredictors{ii,1} = predictors(ind_best_ThisNumPredictors(ii,1),:);
            
        end
        
        [r_best_ThisIter, ind_best_ThisIter] = max(r_best_ThisNumPredictors);
        predictors_best_ThisIter = predictors_best_ThisNumPredictors{ind_best_ThisIter};
        
        X_best_ThisIter = [ones(size(stationsWithoutGlaciers(indsCalib)')), data_norm(stationsWithoutGlaciers(indsCalib),predictors_best_ThisIter)];
        a_best_ThisIter = X_best_ThisIter\YCalib;
        
        models.a{nn,kk} = a_best_ThisIter(2:end);
        models.predictors{nn,kk} = predictors_best_ThisIter;
        models.r(nn,kk) = r_best_ThisIter;
        
    end
    
end

%% Visualize predictors used

for kk = 1:length(predictandsAll(1,:)) %for each predictand
    predictorsUsed{kk,1} = horzcat(models.predictors{:,kk});
end

figure
for kk = 1:length(predictandsAll(1,:)) %make histogram of predictors used for each predictand
    subplot(1,2,kk), hold on
    h=histogram(predictorsUsed{kk,1},'BinMethod','integers'); %bin count = frequency of use of each predictor in bootstrapping + all combinations approach
    for ii = 1:length(predictorsKeep{kk}) %colour red = predictors kept by stepwise regression
        jj = predictorsKeep{kk}(ii);
        x = h.BinEdges(jj:jj+1); 
        x = [x fliplr(x)]; 
        y = [0 0 repmat(h.Values(jj), 1, 2)];
        h_patch = patch(x, y, 'r');
    end
    xlabel('Predictor')
    ylabel('Count')
    if kk==1
        leg = legend('Bootstrapping','Bootstrapping and Stepwise');
    end
    set(gca,'FontSize',16)
end

figure %make boxplot of correlation values for each predictand
boxplot(models.r, 'Labels',{'Mean PC1','Std PC1'});
ylabel('Correlation')
set(gca,'FontSize',18)

%% remove any predictors that were not: chosen by stepwise AND used most often in calib/valid

%do this manually
%e.g. if the first predictor of the second predictand should be removed: predictorsKeep{2}(1) = [];

%% now, bootstrapping to determine regression parameters

clear models r_best_ThisNumPredictors ind_best_ThisNumPredictors predictors_best_ThisNumPredictors
numIterations = 1000;

for kk = 1:length(predictandsAll(1,:)) %for each predictand
    
    disp(['Predictand ',num2str(kk),'/',num2str(length(predictandsAll(1,:)))])
    
    predictands = predictandsAll(:,kk);
    predictors = predictorsKeep{kk};
    
    for nn = 1:numIterations %for each bootstrapping iteration
        
        disp(['    Iteration ',num2str(nn),'/',num2str(numIterations)])
        
        indsCalib = randi([1,length(stationsWithoutGlaciers)],[nCalib,1]);
        indsValid = 1:length(stationsWithoutGlaciers);
        
        predictandsCalib = predictands(stationsWithoutGlaciers(indsCalib),:);
        predictandsValid = predictands(stationsWithoutGlaciers(indsValid),:);
        
        clear r_test p_test
        
        predictorsCalib = data_norm(stationsWithoutGlaciers(indsCalib),predictors);
        predictorsValid = data_norm(stationsWithoutGlaciers(indsValid),predictors);
        
        XCalib = [ones(size(predictorsCalib(:,1))), predictorsCalib];
        YCalib = predictandsCalib;
        
        XValid = [ones(size(predictorsValid(:,1))), predictorsValid];
        YValid = predictandsValid;
        
        a_ThisIter = XCalib\YCalib;
        
        Y_reg = XValid*a_ThisIter;
        
        [r0_test, p0_test] = corrcoef(Y_reg,YValid);
        
        r_ThisIter(jj,1) = r0_test(1,2);
        p_ThisIter(jj,1) = p0_test(1,2);
        
        models.a{nn,kk} = a_ThisIter(1:end); 
        models.predictors{nn,kk} = predictors;
        models.r(nn,kk) = r_ThisIter(jj,1);
        
    end
    
end

%% use MLR to predict mean PC1 and std PC1 for stream gauges, intakes, and dams

for kk = 1:length(predictandsAll(1,:)) %for each predictand, find average 'a' values (coefficients)
    aAllIters = horzcat(models.a{:,kk});
    models.aMean{:,kk} = mean(aAllIters,2);
end

for kk = 1:length(predictandsAll(1,:)) %for meanPC1 and stdPC1, compute modelled meanPC1 and stdPC1 values for stream gauges, communities, and dams
    
   %stream gauges
   X = [ones(size(data_norm(:,1))), data_norm(:,predictorsKeep{kk})];
   predictands_MLR(:,kk) = X*models.aMean{:,kk};
   
   %intakes
   X = [ones(size(intakeData_norm(:,1))), intakeData_norm(:,predictorsKeep{kk})];
   intakePredictands_MLR(:,kk) = X*models.aMean{:,kk};
   
   %dams
   X = [ones(size(damData_norm(:,1))), damData_norm(:,predictorsKeep{kk})];
   damPredictands_MLR(:,kk) = X*models.aMean{:,kk};
      
end

meanPC1_MLR = predictands_MLR(:,1);
stdPC1_MLR = predictands_MLR(:,2);

%% make exponential fit to PC1 vs %glaciers data

%fit to equation: PC1 = m*log(PG) + b
%equivalent to: log(PG) = (1/m)*(PC1 - b)
%equivalent to: log(PG) = (1/m)*PC1 - b/m
%equivalent to: PG = exp((1/m)*PC1 - b/m)
%equivalent to: PG = alpha*exp(beta*PC1)

X = [ones(size(stationsWithGlaciers')), meanPC1(stationsWithGlaciers)];
Y = log(stationPercentGlaciers(stationsWithGlaciers));

a = regress(Y,X);
a = X\Y;
alpha = exp(a(1));
beta = a(2);

PC1_regress = linspace(-10,10,100); %for plotting a regression line
percentGlaciers_regress = alpha*exp(beta*PC1_regress); %for plotting a regression line

%% estimate mean PC1 values from percentage glaciation for intake/dam locations

intakePercentGlaciers(isnan(intakePercentGlaciers)) = 0;
intakePC1_regress = NaN*zeros(size(intakePercentGlaciers));
intakePC1_regress(intakePercentGlaciers>0) = log(intakePercentGlaciers(intakePercentGlaciers>0)/alpha)/beta;

damPC1_regress = NaN*zeros(size(damPercentGlaciers));
damPC1_regress(damPercentGlaciers>0) = log(damPercentGlaciers(damPercentGlaciers>0)/alpha)/beta;

%% for readability, extract individual variables that we care about as (intake/dam)(meanPC1/stdPC1)_(MLR/glaciers)

%modelled mean PC1 and std PC1 for intakes + dams using MLR (case: no glaciers)
intakeMeanPC1_MLR = intakePredictands_MLR(:,1);
intakeStdPC1_MLR = intakePredictands_MLR(:,2);

damMeanPC1_MLR = intakePredictands_MLR(:,1);
damStdPC1_MLR = intakePredictands_MLR(:,2);

%modelled mean PC1 and std PC1 for intakes + dams using glacier fit (case: sufficiently glaciated)
PGthreshold = 0.01; 

intakeMeanPC1_glaciers = intakeMeanPC1_MLR;
intakeMeanPC1_glaciers(intakePercentGlaciers>PGthreshold) = intakePC1_regress(intakePercentGlaciers>PGthreshold);
intakeStdPC1_glaciers = intakeStdPC1_MLR;
intakeStdPC1_glaciers(intakePercentGlaciers>PGthreshold) = mean(stdPC1(stationPercentGlaciers>PGthreshold));

damMeanPC1_glaciers = intakeMeanPC1_MLR;
damMeanPC1_glaciers(damPercentGlaciers>PGthreshold) = damPC1_regress(damPercentGlaciers>PGthreshold);
damStdPC1_glaciers = damStdPC1_MLR;
damStdPC1_glaciers(damPercentGlaciers>PGthreshold) = mean(stdPC1(stationPercentGlaciers>PGthreshold));

%% reconstruct discharge for stations

N = 1e4; %number of iterations (times to pull PC1 values from a distribution with meanPC1 and stdPC1 for reconstructing streamflow)

%initialize variables
Q_kde = zeros(length(stationLon),100);
Q_x = zeros(size(Q_kde));
Q_kde_MLR = zeros(size(Q_kde));
Q_x_MLR = zeros(size(Q_kde));

Q_glaciers = zeros(length(stationLon),N*31);
Q_MLR = zeros(size(Q_glaciers));

for kk = 1:length(stationLon) %for each station, make Q kde for observed/modelled
    
    disp(['Discharge for Station ',num2str(kk),'/',num2str(length(stationLon))])

    PC1 = normrnd(meanPC1(kk),stdPC1(kk),N,1);
    PC1_MLR = normrnd(meanPC1_MLR(kk),stdPC1_MLR(kk),N,1);

    Q0 = PC1*eigvecsF(:,1)' + mean(pcaInputF,1);
    Q_glaciers(kk,:) = reshape(Q0,1,[]);
    [Q_kde(kk,:),Q_x(kk,:)] = ksdensity(Q_glaciers(kk,:));

    Q0_MLR = PC1_MLR*eigvecsF(:,1)' + mean(pcaInputF,1);
    Q_MLR(kk,:) = reshape(Q0_MLR,1,[]);
    [Q_kde_MLR(kk,:),Q_x_MLR(kk,:)] = ksdensity(Q_MLR(kk,:));
    
end

%% reconstruct discharge for communities

%initialize variables
intake_Q_glaciers = zeros(length(cityLon),N*31);
intake_Q_MLR = zeros(size(intake_Q_glaciers));

for kk = 1:length(intakeLon) %for each station, make Q kde for observed/modelled
    
    disp(['Discharge for Community ',num2str(kk),'/',num2str(length(intakeLon))])
    
    if intakeStdPC1_glaciers(kk)>0 && ~isnan(intakePredictands_MLR(kk,1)) %if there are glaciers in headwaters, reconstruct streamflow
    
        disp('    CALCULATED')
        %before: _glaciers
        intakePC1_glaciers = normrnd(intakeMeanPC1_glaciers(kk),intakeStdPC1_glaciers(kk),N,1);
        
        %after: _MLR
        intakePC1_MLR = normrnd(intakeMeanPC1_MLR(kk),intakeStdPC1_MLR(kk),N,1);
        
        Q0_glaciers = intakePC1_glaciers*eigvecsF(:,1)' + mean(pcaInputF,1);
        intake_Q_glaciers(kk,:) = reshape(Q0_glaciers,1,[]);
        [intake_Q_kde_glaciers(kk,:),intake_Q_x_glaciers(kk,:)] = ksdensity(intake_Q_glaciers(kk,:));
        
        Q0_MLR = intakePC1_MLR*eigvecsF(:,1)' + mean(pcaInputF,1);
        intake_Q_MLR(kk,:) = reshape(Q0_MLR,1,[]);
        [intake_Q_kde_MLR(kk,:),intake_Q_x_MLR(kk,:)] = ksdensity(intake_Q_MLR(kk,:));
        
    else %if there are not glaciers in headwaters, don't reconstruct streamflow
        
        intake_Q_kde_glaciers(kk,:) = NaN*ones(1,100);
        intake_Q_x_glaciers(kk,:) = NaN*ones(1,100);
        intake_Q_glaciers(kk,:) = NaN;
        intake_Q_kde_MLR(kk,:) = NaN*ones(1,100);
        intake_Q_x_MLR(kk,:) = NaN*ones(1,100);
        intake_Q_MLR(kk,:) = NaN;
       
    end
end

 %% reconstruct discharge for dams

 %initialize variables
dam_Q_glaciers = zeros(length(damLon),N*31);
dam_Q_MLR = zeros(size(dam_Q_glaciers));

for kk = 1:length(damLon) %for each station, make Q kde for observed/modelled
    
    disp(['Discharge for Dam ',num2str(kk),'/',num2str(length(damLon))])
    
    if damStdPC1_glaciers(kk)>0 && ~isnan(damPredictands_MLR(kk,1)) %if there are glaciers in headwaters, reconstruct streamflow
    
        %before: _glaciers
        damPC1_glaciers = normrnd(damMeanPC1_glaciers(kk),damStdPC1_glaciers(kk),N,1);
        
        %after: _MLR
        damPC1_MLR = normrnd(damMeanPC1_MLR(kk),damStdPC1_MLR(kk),N,1);
        
        Q0_glaciers = damPC1_glaciers*eigvecsF(:,1)' + mean(pcaInputF,1);
        dam_Q_glaciers(kk,:) = reshape(Q0_glaciers,1,[]);
        [dam_Q_kde_glaciers(kk,:),dam_Q_x_glaciers(kk,:)] = ksdensity(dam_Q_glaciers(kk,:));
        
        Q0_MLR = damPC1_MLR*eigvecsF(:,1)' + mean(pcaInputF,1);
        dam_Q_MLR(kk,:) = reshape(Q0_MLR,1,[]);
        [dam_Q_kde_MLR(kk,:),dam_Q_x_MLR(kk,:)] = ksdensity(dam_Q_MLR(kk,:));
        
    else %if there are not glaciers in headwaters, do not reconstruct streamflow
        
        dam_Q_kde_glaciers(kk,:) = NaN*ones(1,100);
        dam_Q_x_glaciers(kk,:) = NaN*ones(1,100);
        dam_Q_glaciers(kk,:) = NaN;
        dam_Q_kde_MLR(kk,:) = NaN*ones(1,100);
        dam_Q_x_MLR(kk,:) = NaN*ones(1,100);
        dam_Q_MLR(kk,:) = NaN;
       
    end
end

%% interpolate reconstructed discharge onto same x-axis (for comparison)

%stations
Nx = 1000;
Q_x_interp = linspace(min(min([Q_x,Q_x_MLR])),max(max([Q_x,Q_x_MLR])),Nx);
Q_kde_interp = zeros(length(stationLon),Nx);
Q_kde_MLR_interp = zeros(size(Q_kde_interp));
for kk = 1:length(stationLon) %for each station, interpolate onto same x-axis
    Q_kde_interp(kk,:) = interp1(Q_x(kk,:),Q_kde(kk,:),Q_x_interp);
    if ~isnan(sum(Q_kde_MLR(kk,:)))
        Q_kde_MLR_interp(kk,:) = interp1(Q_x_MLR(kk,:),Q_kde_MLR(kk,:),Q_x_interp);
    end
end
Q_kde_interp(isnan(Q_kde_interp))=0;
Q_kde_MLR_interp(isnan(Q_kde_MLR_interp))=0;

%communities
intake_Q_x_interp = linspace(nanmin(nanmin([intake_Q_x_glaciers,intake_Q_x_MLR])),nanmax(nanmax([intake_Q_x_glaciers,intake_Q_x_MLR])),Nx);
for kk = 1:length(intakeLon) %for each community, interpolate onto same x-axis
    if ~isnan(sum(intake_Q_kde_glaciers(kk,:)))
        intake_Q_kde_glaciers_interp(kk,:) = interp1(intake_Q_x_glaciers(kk,:),intake_Q_kde_glaciers(kk,:),intake_Q_x_interp);
    end
    if ~isnan(sum(intake_Q_kde_MLR(kk,:)))
        intake_Q_kde_MLR_interp(kk,:) = interp1(intake_Q_x_MLR(kk,:),intake_Q_kde_MLR(kk,:),intake_Q_x_interp);
    end
end
intake_Q_kde_glaciers_interp(isnan(intake_Q_kde_glaciers_interp))=0;
intake_Q_kde_MLR_interp(isnan(intake_Q_kde_MLR_interp))=0;

%dams
dam_Q_x_interp = linspace(nanmin(nanmin([dam_Q_x_glaciers,dam_Q_x_MLR])),nanmax(nanmax([dam_Q_x_glaciers,dam_Q_x_MLR])),Nx);
for kk = 1:length(damLon) %for each community, interpolate onto same x-axis
    if ~isnan(sum(dam_Q_kde_glaciers(kk,:)))
        dam_Q_kde_glaciers_interp(kk,:) = interp1(dam_Q_x_glaciers(kk,:),dam_Q_kde_glaciers(kk,:),dam_Q_x_interp);
    end
    if ~isnan(sum(dam_Q_kde_MLR(kk,:)))
        dam_Q_kde_MLR_interp(kk,:) = interp1(dam_Q_x_MLR(kk,:),dam_Q_kde_MLR(kk,:),dam_Q_x_interp);
    end
end
dam_Q_kde_glaciers_interp(isnan(dam_Q_kde_glaciers_interp))=0;
dam_Q_kde_MLR_interp(isnan(dam_Q_kde_MLR_interp))=0;


%% calculate mean and std Q

station_Q_mean_glaciers = mean(Q_glaciers,2);
station_Q_std_glaciers = std(Q_glaciers,[],2);

station_Q_mean_MLR = mean(Q_MLR,2);
station_Q_std_MLR = std(Q_MLR,[],2);

station_dmeanQ = station_Q_mean_MLR - station_Q_mean_glaciers;
station_dstdQ = station_Q_std_MLR - station_Q_std_glaciers;

intake_Q_mean_glaciers = mean(intake_Q_glaciers,2);
intake_Q_std_glaciers = std(intake_Q_glaciers,[],2);

intake_Q_mean_MLR = mean(intake_Q_MLR,2);
intake_Q_std_MLR = std(intake_Q_MLR,[],2);

intake_dmeanQ = intake_Q_mean_MLR - intake_Q_mean_glaciers;
intake_dstdQ = intake_Q_std_MLR - intake_Q_std_glaciers;

dam_Q_mean_glaciers = mean(dam_Q_glaciers,2);
dam_Q_std_glaciers = std(dam_Q_glaciers,[],2);

dam_Q_mean_MLR = mean(dam_Q_MLR,2);
dam_Q_std_MLR = std(dam_Q_MLR,[],2);

dam_dmeanQ = dam_Q_mean_MLR - dam_Q_mean_glaciers;
dam_dstdQ = dam_Q_std_MLR - dam_Q_std_glaciers;

%% compute RMSE

for kk = 1:length(stationLon) %for each station
    stationRMSE(kk,1) = sqrt(mean((Q_kde_interp(kk,:) - Q_kde_MLR_interp(kk,:)).^2));
    [r0,p0] = corrcoef(Q_kde_interp(kk,:)',Q_kde_MLR_interp(kk,:)');
    stationQcorr(kk,1) = r0(1,2);
    stationQp(kk,1) = p0(1,2);
end

meanRMSE_noglacier = mean(stationRMSE(stationsWithoutGlaciers));
stdRMSE_noglacier = std(stationRMSE(stationsWithoutGlaciers));
RMSE_glacier = stationRMSE(stationsWithGlaciers);

%%
%%
%%
%%
%%
%% Make figures

figFS = 20; %fontsize to use

minLat = 48.99; %for bounding box around AB
maxLat = 60.01;
minLon = -120.01;
maxLon = -109.98;
glacierRGB = [124,106,238]/256; %colour of glaciers

indLL_city = find(strcmp(cityName,'Lake Louise')); %indices of key locations
indHinton_city = find(strcmp(cityName,'Hinton'));
indBighorn_dam = find(strcmp(damName,'Bighorn Dam'));

%% Figure 1: visualize rivers, communities, and connections

% find point in river nearest to via
for kk = 1:length(cityLon) %for each community    
    if isnan(viaLat(kk)) %if there is no via lat/lon, then the water is sourced from the community lat/lon
            viaLat(kk) = cityLat(kk);
            viaLon(kk) = cityLon(kk);
    end  
end

riverInd = zeros(size(cityLon));
intakeInd = zeros(size(cityLon));
dist_river2via_min = zeros(size(cityLon));
for kk = 1:length(cityLon) %for each community
    
    for jj = 1:length(riverName) %find the river/headwater index that matches this location     
        if strcmp(riverName{jj},source{kk})
            riverInd(kk) = jj; %index in riverName of the river that is the source for community kk
        end
    end
    
    if riverInd(kk)~=0
        d2lat = (flowRoute{riverInd(kk)}(:,1) - viaLat(kk)).^2;
        d2lon =  (flowRoute{riverInd(kk)}(:,2) - viaLon(kk)).^2;
        dist = sqrt(d2lat + d2lon);
    end
    
    [dist_river2via_min(kk,1), intakeInd(kk,1)] = min(dist);
    
end

intakeLon = NaN*zeros(size(riverInd));
intakeLat = NaN*zeros(size(intakeLon));
for kk = 1:length(cityLon) %for each community
    if riverInd(kk)>0
        intakeLon(kk) = flowRoute{riverInd(kk)}(intakeInd(kk),2);
        intakeLat(kk) = flowRoute{riverInd(kk)}(intakeInd(kk),1);
    end
end

inds = find(riverInd>0);
colours = parula(length(riverName));
glacierRGB = [124,106,238]/256;

figure('Position',[0,0,800,1200])

plot_rivers()
plot_AB() %plot border
plot_AB_glaciers(10,glacierRGB)

%for each glacier-fed river, plot river
for kk = 1:length(riverName)
    inds = strcmp(source,riverName{kk});
    s1 = scatter(flowRoute{kk}(:,2),flowRoute{kk}(:,1),20,'filled');
    set(s1,'CData',colours(kk,:))
end

%plot non-glacier-fed-river-sourced communitites
a = 0.2; %transparency of non-glacier-fed locations
s = 100; %size
for kk = 1:5 %for each source type
    inds = find(sourceType == kk);
    if kk == 1 %if groundwater
        scatter(cityLon(inds), cityLat(inds), s, 'filled', 'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',a, 'MarkerEdgeAlpha',0.5,'Marker','d')
    elseif kk == 2 %if river
        scatter(cityLon(inds), cityLat(inds), s, 'filled', 'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',a, 'MarkerEdgeAlpha',0.5,'Marker','o')
    elseif kk == 3 %if lake
        scatter(cityLon(inds), cityLat(inds), s, 'filled', 'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',a, 'MarkerEdgeAlpha',0.5,'Marker','s')
    elseif kk == 4 %if missing
        %do nothing
    elseif kk == 5 %if combination
        scatter(cityLon(inds), cityLat(inds), s, 'filled', 'LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerFaceAlpha',a, 'MarkerEdgeAlpha',0.5,'Marker','s')
    end
end

%for each glacier-fed-river-sourced community, plot line connecting community with intake location
for kk = 1:length(cityLon)
    if riverInd(kk)>0 %if sourcing from a glacier-fed river, plot line
        l = line([cityLon(kk), intakeLon(kk)],[cityLat(kk),intakeLat(kk)]);
        l.Color = 'k';
    end
end

%for each glacier-fed river, plot communities
for kk = 1:length(riverName)
    inds = strcmp(source,riverName{kk});
    s3 = scatter(flowRoute{kk}(intakeInd(find(riverInd==kk)),2),flowRoute{kk}(intakeInd(find(riverInd==kk)),1),10,'k','filled');
    s2 = scatter(cityLon(inds),cityLat(inds),s,'filled','LineWidth',1,'MarkerEdgeColor','k');
    set(s2,'CData',colours(kk,:))
end

%plot Calgary on top (coloured as glacier-fed, since it receives from both Bow (glacier-fed) and Elbow (non-glacier-fed) rivers
calgaryInd = 41;
s0 = scatter(cityLon(calgaryInd),cityLat(calgaryInd),s,'filled','LineWidth',1,'MarkerEdgeColor','k');
set(s0,'CData',colours(2,:))

%plot Edmonton on top
edmontonInd = 86;
s0 = scatter(cityLon(edmontonInd),cityLat(edmontonInd),s,'filled','LineWidth',1,'MarkerEdgeColor','k');
set(s0,'CData',colours(4,:))

%plot dams
scatter(damLon,damLat,s,'k^','filled')

%% Figure 2: PCA and SOM results

%figure 2, panel a: eigenvectors

figure('Position',[0,0,500,1200])

subplot(2,1,1), hold on %first, plot eigenvectors
for kk = 1:2
    p(kk) = plot(eigvecsF(:,kk),'k','LineWidth',2);
end
p(2).LineStyle = '--';
leg = legend('Mode 1','Mode 2');
leg.Location = 'best';
leg.String = leg.String(1:2);

minX = 1;
maxX = 31;
minY = min(min(eigvecsF(:,1:2)));
maxY = max(max(eigvecsF(:,1:2)));
fracX = -0.18;
fracY = 1.1;
textX = minX + (maxX - minX)*fracX;
textY  = minY + (maxY - minY)*fracY;
text(textX,textY,'a)','fontsize',figFS)

xlim([minX,maxX])
ylim([minY, maxY])
xlabel('Day of August')
ylabel('Normalized Streamflow')
title('August Streamflow Modes')
set(gca,'FontSize',figFS)

%figure 2, panel b: SOM clusters

C = som_colorcode(sM,'rgb2');

subplot(2,1,2), hold on

plot_rivers()
plot_AB_glaciers(10,glacierRGB);
plot_AB()
stationplot = scatter(stationLon,stationLat,150,'filled','LineWidth',2,'MarkerEdgeColor','k');
set(stationplot,'CData',cStations);
cmap = cbmus;
colormap(gca,cmap);
cbstationvals = colorbar;
cbstationvals.Location = 'eastoutside';
cbstationvals.TickLabels = [];
cbstationvals.Ticks = [];
cbstationvals.Ticks = linspace(1/12,11/12,length(cbmus));
cbstationvals.TickLabels = 1:length(cbmus);
xlabel('Longitude')
ylabel('Latitude')
ylabel(cbstationvals,'Cluster','fontsize',figFS)
set(gca,'FontSize',figFS)
titlestr = sprintf('Streamflow Clusters');
title(titlestr,'FontSize',figFS)
minX = minLon;
maxX = maxLon;
minY = minLat;
maxY = maxLat;
fracX = -0.225;
textX = minX + (maxX - minX)*fracX;
textY  = minY + (maxY - minY)*fracY;
text(textX,textY,'b)','fontsize',figFS)
xlim([minLon,maxLon])
ylim([minLat,maxLat])
%%
%figure 2, panel c: SOM patterns

figure('Position',[0,0,1000,1200]), hold on
h = som_cplane(sM,C);
h.FaceAlpha = 0.7;
titlestr = sprintf('SOM With Clusters in PC-Space');
title(titlestr)
set(gca,'FontSize',figFS)

xdata = real(sM.codebook);
ydata = imag(sM.codebook);

xmids = h.XData(1,:);
ymids = (h.YData(2,:) + h.YData(3,:))/2;

w = 1; %width of hexagon to take up
h = 1; %height of hexagon to take up

xscale = w/(max(max(xdata)) - min(min(xdata)));
yscale = h/(max(max(ydata)) - min(min(ydata)));

axW = 0.14 * 1.25;
axH = 0.14 * 1.15;
axOffset = 0.128;
axLeftEdge = 0.25;
axPos = [axLeftEdge,0.685,axW,axH;
    axLeftEdge + axOffset, 0.445, axW, axH;
    axLeftEdge,0.205,axW,axH;
    axLeftEdge + 2*axOffset,0.685,axW,axH;
    axLeftEdge + 3*axOffset, 0.445, axW, axH;
    axLeftEdge + 2*axOffset,0.205,axW,axH];

for kk = 1:en %for each cluster, scatter data
    
    otherRows = find([1:en]~=kk);
    
    xmid = xmids(kk);
    ymid = ymids(kk);
    
    xdiff = mean(mean(xdata))*xscale - xmid+0.1;
    ydiff = mean(mean(ydata))*yscale - ymid;
    
    xnew = xdata*xscale;
    
    g = 1;
    grey = [g,g,g];
    circleSize = 80;
    axes('Position',axPos(kk,:)), hold on
    if kk==1 %if first cluster, label c) for whole figure
        minX = min(min(xdata))-1;
        maxX = max(max(xdata))+1;
        minY = min(min(ydata))-1;
        maxY = max(max(ydata))+1;
        fracX = -0.5;
        fracY = 1.92;
        textX = minX + (maxX - minX)*fracX;
        textY  = minY + (maxY - minY)*fracY;
        text(textX,textY,'c)','fontsize',figFS)
    end
    box on
    s = scatter(reshape(xdata(otherRows,:),[],1),reshape(ydata(otherRows,:),[],1),'Marker','^','MarkerFaceColor',grey,'MarkerEdgeColor','k','SizeData',circleSize); %plot all points
    s = scatter(xdata(kk,:),ydata(kk,:),'MarkerFaceColor',cbmus(kk,:),'MarkerEdgeColor','k','SizeData',circleSize); %plot all points
    xlim([min(min(xdata))-1,max(max(xdata))+1])
    ylim([min(min(ydata))-1,max(max(ydata))+1])
    xlabel('PC1')
    ylabel('PC2')
    %titlestr = sprintf('$Cluster \pound$%d',kk);
    titlestr = sprintf('$Cluster %#$%d',kk);
    %#5.0f
    titlestr = titlestr(1:9);
    title(sprintf('Cluster #%1.0d',kk))
    set(gca,'FontSize',14)
    set(gca,'XTick',[-4:4:8])
    set(gca,'XTickLabel',{-4:4:8})
    set(gca,'YTick',[-4:4:4])
    set(gca,'YTickLabel',{-4:4:4})
    set(gca,'XColor',[0,0,0])
    set(gca,'YColor',[0,0,0])
    
    grid on
    
end

%% Figure 3: MLR + OLR results

figure('Position',[0,0,1000,1500])

%figure 3, panel a: stdPC1 vs meanPC1, observed

subplot(2,2,1)
scatter(meanPC1,stdPC1,100,cStations,'filled','LineWidth',2,'MarkerEdgeColor','k')

minX = min([meanPC1;meanPC1_MLR])-0.2;
maxX = max([meanPC1;meanPC1_MLR])+0.2;
minY = min([stdPC1;stdPC1_MLR])-0.2;
maxY = max([stdPC1;stdPC1_MLR])+0.2;
fracX = -0.15;
fracY = 1.08;
textX = minX + (maxX - minX)*fracX;
textY  = minY + (maxY - minY)*fracY;
text(textX,textY,'a)','fontsize',figFS)

xlim([minX,maxX])
ylim([minY,maxY])
xlabel('PC_1')
ylabel('\sigma_{PC_1}')
title('Measured: With Glaciers')
set(gca,'fontsize',figFS)
grid on

%figure 3, panel b: stdPC1 vs meanPC1, modelled

subplot(2,2,3)
scatter(meanPC1_MLR,stdPC1_MLR,100,cStations,'filled','LineWidth',2,'MarkerEdgeColor','k')
fracX = -0.15;
fracY = 1.08;
textX = minX + (maxX - minX)*fracX;
textY  = minY + (maxY - minY)*fracY;
text(textX,textY,'b)','fontsize',figFS)
xlim([minX,maxX])
ylim([minY,maxY])
xlabel('PC_1')
ylabel('\sigma_{PC_1}')
title('Modelled: Without Glaciers')
set(gca,'fontsize',figFS)
grid on

%figure 3, panel c: meanPC1 vs % glaciation

subplot(2,2,2), hold on
err = stdPC1;
errorbar(stationPercentGlaciers,meanPC1,err,'k','LineStyle','none')
s = scatter(stationPercentGlaciers,meanPC1,100,'filled','LineWidth',2,'MarkerEdgeColor','k');
s.CData = cStations;
plot(percentGlaciers_regress,PC1_regress,'k--','LineWidth',3)
minX = 0.0001;
maxX = 1;
minY = -4;
maxY = 11;
fracX = 0.0185;
fracY = 1.08;
textX = log10(minX) + (log10(maxX) - log10(minX))*log10(fracX);
textY  = minY + (maxY - minY)*fracY;
text(exp(textX),textY,'c)','fontsize',figFS)
ylim([minY,maxY])
xlabel('% Glaciation')
ylabel('PC_1')
title('Glacier Influence')
set(gca,'FontSize',figFS)
xlim([minX,maxX])
ylim([-3,11])
set(gca,'XScale','log')
cmap = cbmus;
colormap(gca,cmap);
cbstationvals = colorbar;
cbstationvals.Location = 'eastoutside';
cbstationvals.TickLabels = [];
cbstationvals.Ticks = [];
cbstationvals.Ticks = linspace(1/12,11/12,length(cbmus));
cbstationvals.TickLabels = 1:length(cbmus);
grid on

%figure 3, panel d: projected streamflow

textX = -2;
textY = 1.8;
textYd= -0.15;
textFirstTwoDiff = 0.0;
Ylim = 2.2;

subplot(2,2,4),hold on

RGB = brewermap([],'RdYlBu');
RGBdam = RGB(10,:);
RGBhinton = RGB(53,:);
RGBll = RGB(64,:);

lw = 5;
plot(dam_Q_x_interp,dam_Q_kde_glaciers_interp(indBighorn_dam,:),'color',RGBdam,'LineWidth',lw)
plot(intake_Q_x_interp,intake_Q_kde_glaciers_interp(indLL_city,:),'color',RGBll,'LineWidth',lw)
plot(intake_Q_x_interp,intake_Q_kde_glaciers_interp(indHinton_city,:),'color',RGBhinton,'LineWidth',lw)

plot(dam_Q_x_interp,dam_Q_kde_MLR_interp(indBighorn_dam,:),'color',RGBdam,'LineStyle','--','LineWidth',lw)
plot(intake_Q_x_interp,intake_Q_kde_MLR_interp(indLL_city,:),'color',RGBll,'LineStyle','--','LineWidth',lw)
plot(intake_Q_x_interp,intake_Q_kde_MLR_interp(indHinton_city,:),'color',RGBhinton,'LineStyle','--','LineWidth',lw)

xlabel('Normalized Discharge')
ylabel('Probability')
title('Historical and Projected Streamflow')

minX = -1.5;
maxX = 2.2;
minY = 0;
maxY = Ylim;
fracX = -0.15;
fracY = 1.08;
textX1 = minX + (maxX - minX)*fracX;
textY1 = minY + (maxY - minY)*fracY;
text(textX1,textY1,'d)','fontsize',figFS)

xlim([minX,maxX])
ylim([minY,maxY])

leg = legend('Bighorn Dam','Lake Louise','Hinton');
leg.Location = 'best';

set(gca,'fontsize',figFS)


%% percentage variance explained by modes

figure, hold on
plot(1:length(fracVarF),fracVarF*100,'k','linewidth',2)
scatter(1:length(fracVarF),fracVarF*100,50,'k','filled')
xlabel('Mode Number')
ylabel('% Variance Explained')
xlim([1,31])
title('Variance Explained by PCA Modes')
set(gca,'fontsize',figFS)
grid on

%% PC1-PC2-time figures (regular, coloured by SOM, and clusters)

yearRange = [1987,2010];

figure('Position',[0,0,1600,400])

s1 = subplot(1,3,1); hold on
for kk = 1:length(stationLon)
    plot3(real(complex_coord(kk,:)),imag(complex_coord(kk,:)),yearRange(1):yearRange(2),'linewidth',2)
end
xlabel('PC_1')
ylabel('PC_2')
zlabel('Year')
title('August Flow in PC-Space')
set(gca,'fontsize',figFS)

minX = -8;
maxX = 12;
minY = -5;
maxY = 8;
minZ = yearRange(1);
maxZ = yearRange(2);
fracX = -0.15;
fracY = 1.1;
fracZ = 1.35;
textX = minX + (maxX - minX)*fracX;
textY = minY + (maxY - minY)*fracY;
textZ = minZ + (maxZ - minZ)*fracZ;
text(textX,textY,textZ,'a)','fontsize',figFS)

xlim([minX,maxX])
ylim([minY,maxY])
zlim([minZ,maxZ])
grid on
v = [-5 -8 4];
view(v)

s2 = subplot(1,3,2); hold on
for kk = 1:length(stationLon)
    plot3(real(complex_coord(kk,:)),imag(complex_coord(kk,:)),yearRange(1):yearRange(2),'color',cStations(kk,:),'linewidth',2)
end
xlabel('PC_1')
ylabel('PC_2')
zlabel('Year')
title('August Flow in PC-Space')
set(gca,'fontsize',figFS)
text(textX,textY,textZ,'b)','fontsize',figFS)
xlim([minX,maxX])
ylim([minY,maxY])
zlim([minZ,maxZ])
grid on
v = [-5 -8 4];
view(v)

s3 = subplot(1,3,3); hold on
grid on
for kk = 1:en
    plot3(real(sM.codebook(kk,:)),imag(sM.codebook(kk,:)),1987:2010,'color',cbmus(kk,:),'linewidth',2)
end
xlabel('PC_1')
ylabel('PC_2')
zlabel('Year')
xlim([minX,maxX])
ylim([minY,maxY])
zlim([minZ,maxZ])
text(textX,textY,textZ,'c)','fontsize',figFS)
title('SOM Patterns')
set(gca,'fontsize',figFS)
originalPosition = get(s3,'Position');
cmap = cbmus;
colormap(gca,cmap);
cbstationvals = colorbar;
cbstationvals.Location = 'eastoutside';
cbstationvals.TickLabels = [];
cbstationvals.Ticks = [];
cbstationvals.Ticks = linspace(1/12,11/12,length(cbmus));
cbstationvals.TickLabels = 1:length(cbmus);
ylabel(cbstationvals,'Cluster')
set(s3,'Position',originalPosition)
v = [-5 -8 4];
view(v)

%% SF: show AB map w/ clustered gauges if there are fewer/more clusters

% initilizing SOM chose size of the map for SOM
data = complex_coord;
ny_som = 3;
nx_som = 1;
en = nx_som*ny_som;
[sM, sT, sMap, bmus] = do_SOM(data, nx_som, ny_som);

cbmus = som_colorcode(sM,'rgb2'); %rgb colours of clusters
for kk = 1:length(data(:,1))
    cStations(kk,:) = cbmus(bmus(kk),:); %rgb colours of stations
end

figure('Position',[0,0,1600,600])

subplot(1,2,1)
plot_rivers()
plot_AB_glaciers(10,glacierRGB);
plot_AB()
stationplot = scatter(stationLon,stationLat,150,'filled','LineWidth',2,'MarkerEdgeColor','k');
set(stationplot,'CData',cStations);

cmap = cbmus;
colormap(gca,cmap);
cbstationvals = colorbar;
cbstationvals.Location = 'eastoutside';
cbstationvals.TickLabels = [];
cbstationvals.Ticks = [];
cbstationvals.Ticks = linspace(1/6,5/6,length(cbmus));
cbstationvals.TickLabels = 1:length(cbmus);

titlestr = sprintf('3x1 SOM Clusters');
title(titlestr)
xlabel('Longitude')
ylabel('Latitude')
set(gca,'FontSize',24)

xlim([-120.01,-109.98])
ylim([48.99,60.01])

minX = -120;
maxX = -110;
minY = 49;
maxY = 60;
fracX = -0.2;
fracY = 1.08;
textX = minX + (maxX - minX)*fracX;
textY = minY + (maxY - minY)*fracY;
text(textX,textY,'a)','fontsize',figFS+4)

ny_som = 3;
nx_som = 3;
en = nx_som*ny_som;
[sM, sT, sMap, bmus] = do_SOM(data, nx_som, ny_som);

cbmus = som_colorcode(sM,'rgb2'); %rgb colours of clusters
for kk = 1:length(data(:,1))
    cStations(kk,:) = cbmus(bmus(kk),:); %rgb colours of stations (according to their cluster)
end


subplot(1,2,2)
plot_rivers()
plot_AB_glaciers(10,glacierRGB);
plot_AB()
stationplot = scatter(stationLon,stationLat,150,'filled','LineWidth',2,'MarkerEdgeColor','k');
set(stationplot,'CData',cStations);

cmap = cbmus;
colormap(gca,cmap);
cbstationvals = colorbar;
cbstationvals.Location = 'eastoutside';
cbstationvals.TickLabels = [];
cbstationvals.Ticks = [];
cbstationvals.Ticks = linspace(1/20,19/20,length(cbmus));
cbstationvals.TickLabels = 1:length(cbmus);

titlestr = sprintf('3x3 SOM Clusters');
title(titlestr)
xlabel('Longitude')
ylabel('Latitude')
set(gca,'FontSize',24)

xlim([-120.01,-109.98])
ylim([48.99,60.01])

text(textX,textY,'b)','fontsize',figFS+4)

%% SF: Show predictors used in MLR

%% first, interpolate DEM

load('./data/AB_topo_lowres.mat')

%%

filename = './data/Canada_Borders/PROVINCE.SHP';
provshape = shaperead(filename,'UseGeoCoords',true);
abLat = provshape(1).Lat;
abLon = provshape(1).Lon;

minX = -120;
maxX = -110;
minY = 49;
maxY = 60;
fracX = -0.2;
fracY = 1.2;
textX = minX + (maxX - minX)*fracX;
textY = minY + (maxY - minY)*fracY;

for kk = 1:length(lonERA)
    for jj  = 1:length(latERA)
        inProv(kk,jj) = inpolygon(lonERA(kk),latERA(jj),abLon,abLat);
    end
end
inProv = double(inProv);
inProv(inProv~=1) = NaN;

for kk = 1:length(demLonInterp)
    for jj  = 1:length(demLatInterp)
        inProvDEM(kk,jj) = inpolygon(demLonInterp(kk),demLatInterp(jj),abLon,abLat);
    end
end
inProvDEM = double(inProvDEM);
inProvDEM(inProvDEM~=1) = NaN;

minLat = 48.99;
maxLat = 60.01;
minLon = -120.01;
maxLon = -109.98;

%%

figure('Position',[0,0,800,800])

subplot(3,2,1)
levels = [floor(min(min((TJJAmean.*inProv)' - 273.15))):0.7:ceil(max(max((TJJAmean.*inProv)' - 273.15)))];
contourf(lonERA,latERA,(TJJAmean.*inProv)' - 273.15,levels);
plot_AB;
set(gca,'YDir','normal');
colormap(gca,redblue);
ylabel('Latitude')
title(sprintf('\\_     \n T_{JAS}'))
title('T_{JJA}')
text(textX,textY,'a)','fontsize',figFS)
set(gca,'FontSize',figFS)
cb = colorbar;
ylabelstr = ['Temperature [',char(176),'C]'];
ylabel(cb,ylabelstr)
set(cb,'fontsize',figFS-4)
xlim([minLon,maxLon])
ylim([minLat,maxLat])
colourStr = '*RdBu';
colormap(gca,brewermap([],colourStr))

subplot(3,2,2)
contourf(lonERA,latERA,(PJJAmean.*inProv)');
plot_AB;
set(gca,'YDir','normal');
colormap(gca,bluered);
title(sprintf('\\_     \n P_{JAS}'))
title('P_{JJA}')
text(textX,textY,'b)','fontsize',figFS)
set(gca,'FontSize',figFS)
cb = colorbar;
ylabel(cb,'Precipitation [mm]')
set(cb,'fontsize',figFS-4)
xlim([minLon,maxLon])
ylim([minLat,maxLat])
colourStr = 'RdBu';
colormap(gca,brewermap([],colourStr))

subplot(3,2,3)
levels = [floor(nanmin(nanmin((Tmean.*inProv)' - 273.15))):ceil(nanmax(nanmax((Tmean.*inProv)' - 273.15)))];
contourf(lonERA,latERA,(Tmean.*inProv)' - 273.15,levels);
plot_AB;
set(gca,'YDir','normal');
colormap(gca,redblue);
ylabel('Latitude')
title(sprintf('\\_     \n T_{year}'))
title('T_{year}')
text(textX,textY,'c)','fontsize',figFS)
set(gca,'FontSize',figFS)
cb = colorbar;
ylabelstr = ['Temperature [',char(176),'C]'];
ylabel(cb,ylabelstr)
set(cb,'fontsize',figFS-4)
xlim([minLon,maxLon])
ylim([minLat,maxLat])
colourStr = '*RdBu';
colormap(gca,brewermap([],colourStr))

subplot(3,2,4)
contourf(lonERA,latERA,(provincePercentMinLog.*inProv)');
plot_AB;
set(gca,'YDir','normal');
colormap(gca,bluered);
title('log[Q_{min}/Q_{max}]')
text(textX,textY,'d)','fontsize',figFS)
set(gca,'FontSize',figFS)
cb = colorbar;
ylabelstr = '$log\bigg[\frac{Q_{min}}{Q_{max}}\bigg]$';
ylabelstr = 'log[Q_{min}/Q_{max}]';
ylabel(cb,ylabelstr)
set(cb,'fontsize',figFS-4)
xlim([minLon,maxLon])
ylim([minLat,maxLat])
colourStr = 'RdBu';
colormap(gca,brewermap([],colourStr))

subplot(3,2,5)
levels = linspace(nanmin(nanmin(EJJAmean.*inProv(2:end,:))),nanmax(nanmax(EJJAmean.*inProv(2:end,:))),8);
contourf(lonERA(2:end),latERA,(EJJAmean.*inProv(2:end,:))',levels);
plot_AB;
set(gca,'YDir','normal');
colormap(gca,bluered);
title(sprintf('\\_     \n E_{JJA}'))
title('E_{JJA}')
text(textX,textY,'e)','fontsize',figFS)
set(gca,'FontSize',figFS)
cb = colorbar;
xlabel('Longitude')
ylabel('Latitude')
ylabel(cb,'Evaporation [m.w.e.]')
set(cb,'fontsize',figFS-4)
xlim([minLon,maxLon])
ylim([minLat,maxLat])
colourStr = 'RdBu';
colormap(gca,brewermap([],colourStr))

subplot(3,2,6)
data = demInterp.*inProvDEM';
[nr,nc] = size(data);
surf(demLonInterp,demLatInterp,data,'edgecolor','none');
view(2)
plot_AB;
set(gca,'YDir','normal');
demcmap(data)
title(sprintf('\\_     \n E_{JJA}'))
title('Elevation')
text(textX,textY,'f)','fontsize',figFS)
set(gca,'FontSize',figFS)
cb = colorbar;
xlabel('Longitude')
ylabel(cb,'Elevation [m]')
set(cb,'fontsize',figFS-4)
xlim([minLon,maxLon])
ylim([minLat,maxLat])

%% SF: KDEs of MLR coefficients

labels = {'$\overline{T}_{year}$','$\overline{P}_{JJA}$','$\overline{T}_{JJA}$','$\overline{E}_{JJA}$','Elevation','$log\bigg[\frac{Q_{min}}{Q_{max}}\bigg]$'};
inds1 = [1,2,3,4,6];
inds2 = [2,4,5,6];
figure('Position',[0,0,1200,400]) %make kde of a values for each predictand
for kk = 1:2 %for each predictand
    
    subplot(1,2,kk), hold on
    
    kdedata = horzcat(models.a{:,kk});
    for jj = 2:length(kdedata(:,kk)) %for each predictor, make kde; don't include a(1) values (ie: intercepts)
        [kde,XI] = ksdensity(kdedata(jj,:));
        p = plot(XI,kde);
        p.LineWidth = 4;
    end
    xlabel('Regression Coefficient')
    ylabel('Probability')
    if kk == 1
        title(sprintf('\\_  \n PC_1'))
        title('PC_1')
        title('$\overline{PC_1}$','Interpreter','latex')
        xlim([-1,2])
        ylim([0,7])
        minX = -1;
        maxX = 1.1;
        minY = 0;
        maxY = 7;
        fracX = -0.15;
        fracY = 1.15;
        textX = minX + (maxX - minX)*fracX;
        textY = minY + (maxY - minY)*fracY;
        text(textX,textY,'a)','fontsize',figFS)
        leg = legend(labels(inds1),'Interpreter','latex','location','best');
    elseif kk==2
        title('\sigma_{PC_1}')
        xlim([-1,1.5])
        ylim([0,8])
        minX = -1;
        maxX = 1.1;
        minY = 0;
        maxY = 8;
        fracX = -0.15;
        fracY = 1.15;
        textX = minX + (maxX - minX)*fracX;
        textY = minY + (maxY - minY)*fracY;
        text(textX,textY,'b)','fontsize',figFS)
        leg = legend(labels(inds2),'Interpreter','latex','location','best');
    end
    set(gca,'fontsize',figFS)
end

%% show minflow/maxflow seasonal hydrograph

RGB = brewermap([],'RdYlBu');
RGB1 = RGB(10,:);
RGB2 = RGB(53,:);
RGBll = RGB(64,:);

lw = 3;

ind1 = 134;
ind2 = 19;

sf1 = movmean(all_flowseason(ind1,:),30);
sf2 = movmean(all_flowseason(ind2,:),30);

min1 = min(sf1);
min2 = min(sf2);
max1 = max(sf1);
max2 = max(sf2);

time = linspace(1,365);

dash_min1 = min1*ones(size(time));
dash_min2 = min2*ones(size(time));
dash_max1 = max1*ones(size(time));
dash_max2 = max2*ones(size(time));

pmin1 = [5,min1-8];
pmax1 = [5,max1+8];
pmin2 = [5,min2+8];
pmax2 = [5,max2-8];

figure('Position',[0,0,1000,500])

subplot(1,2,1), hold on
plot(sf1,'color',RGB1,'linewidth',lw)
plot(time,dash_min1,'k--','linewidth',lw-1)
plot(time,dash_max1,'k--','linewidth',lw-1)
text(pmin1(1),pmin1(2),'Q_{min, Lesser Slave}','fontsize',figFS)
text(pmax1(1),pmax1(2),'Q_{max, Lesser Slave}','fontsize',figFS)
xlabel('Day of Year')
ylabel('Q [m^3/s]')
title(stationName(ind1))
xlim([1,365])
ylim([0,160])
yticks([0,50,100,150])
set(gca,'fontsize',figFS)

subplot(1,2,2), hold on
plot(sf2,'color',RGBll,'linewidth',lw)
plot(time,dash_min2,'k--','linewidth',lw-1)
plot(time,dash_max2,'k--','linewidth',lw-1)
text(pmin2(1),pmin2(2),'Q_{min, Bow}','fontsize',figFS)
text(pmax2(1),pmax2(2),'Q_{max, Bow}','fontsize',figFS)
xlabel('Day of Year')
ylabel('Q [m^3/s]')
title(stationName(ind2))
xlim([1,365])
ylim([0,160])
yticks([0,50,100,150])
set(gca,'fontsize',figFS)

%% show normalization steps: raw, normalized, smoothed, seasonal

ind = 103;
lw = 2;
RGB = brewermap([],'RdYlBu');
RGB1 = RGB(10,:);
RGB2 = RGB(53,:);
RGB3 = RGB(64,:);

textx = -95;
yfactor = 0.2;

figure('Position',[0,0,800,600])

subplot(2,2,1), hold on %raw
plot(flow_years(24*(ind-1)+1,:),'color',RGB3,'linewidth',lw)
xlabel('Day of Year')
ylabel('Q [m^3/s]')
title('Raw')
ymin = 0;
ymax = 500;
texty = ymax + yfactor*(ymax-ymin);
xlim([0,365])
ylim([ymin,ymax])
text(textx,texty,'a)','fontsize',figFS)
set(gca,'fontsize',figFS)

subplot(2,2,2), hold on %normalized
plot(flow_years_norm(24*(ind-1)+1,:),'color',RGB3,'linewidth',lw)
xlabel('Day of Year')
ylabel('Q')
title('Normalized')
ymin = -1;
ymax = 4;
texty = ymax + yfactor*(ymax-ymin);
xlim([0,365])
ylim([ymin,ymax])
text(textx,texty,'b)','fontsize',figFS)
set(gca,'fontsize',figFS)

subplot(2,2,3), hold on %smoothed
plot(flow_years_norm_smooth(24*(ind-1)+1,:),'color',RGB3,'linewidth',lw)
xlabel('Day of Year')
ylabel('Q')
title('Normalized and Smoothed')
ymin = -1;
ymax = 3;
texty = ymax + yfactor*(ymax-ymin);
xlim([0,365])
ylim([ymin,ymax])
text(textx,texty,'c)','fontsize',figFS)
set(gca,'fontsize',figFS)

subplot(2,2,4), hold on %seasonal
plot(all_flowseason(ind,:),'color',RGB3,'linewidth',lw)
xlabel('Day of Year')
ylabel('Q [m^3/s]')
title('Seasonal')
ymin = 0;
ymax = 300;
texty = ymax + yfactor*(ymax-ymin);
xlim([0,365])
ylim([ymin,ymax])
text(textx,texty,'d)','fontsize',figFS)
set(gca,'fontsize',figFS)

%% physical interpretation of PC1 (and PC2?)

%calc slope of august flow

for kk = 1:length(pcaInputF(:,1))
    x1 = [1:31]';
    X = [ones(size(x1)),x1];
    y = pcaInputF(kk,:)';
    a = X\y;
    slope(kk) = a(2);
end

%%
% n = 100;
% gridx1 = linspace(min(PCsF(:,1)),max(PCsF(:,1)),n);
% x = mean(flow_years_norm_smooth(:,summerInds),2);
% gridx2 = linspace(min(x),max(x),n);
% [x1,x2] = meshgrid(gridx1,gridx2);
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% xdata = [PCsF(:,1),x];
% 
% kde1 = ksdensity(xdata,xi);
% k1 = reshape(kde1,[n,n]);
% 
% for kk = 1:length(PCsF(:,1)) %find kde value
%     dPC = abs(PCsF(kk,1) - gridx1);
%     dF = abs(x(kk) - gridx2);
%     [~,indPC] = min(dPC);
%     [~,indF] = min(dF);
%     z1(kk,1) = k1(indF,indPC);
% end
% 
% n = 100;
% gridx1 = linspace(min(PCsF(:,2)),max(PCsF(:,2)),n);
% x = slope';
% gridx2 = linspace(min(x),max(x),n);
% [x1,x2] = meshgrid(gridx1,gridx2);
% x1 = x1(:);
% x2 = x2(:);
% xi = [x1 x2];
% xdata = [PCsF(:,2),x];
% 
% kde2 = ksdensity(xdata,xi);
% k2 = reshape(kde2,[n,n]);
% 
% for kk = 1:length(PCsF(:,2)) %find kde value
%     dPC = abs(PCsF(kk,2) - gridx1);
%     dF = abs(x(kk) - gridx2);
%     [~,indPC] = min(dPC);
%     [~,indF] = min(dF);
%     z2(kk,1) = k2(indF,indPC);
% end

figure('Position',[0,0,1000,400])

subplot(1,2,1)
r = corr(PCsF(:,1),mean(flow_years_norm_smooth(:,summerInds),2))
s = scatter(mean(flow_years_norm_smooth(:,summerInds),2),PCsF(:,1));%,'filled');
text(-1.5,8,sprintf('Correlation > 0.99'),'fontsize',figFS)
ylabel('PC_1')
xlabel('Mean Normalized August Flow')
title('Interpretation of Mode #1')
set(gca,'fontsize',figFS)

subplot(1,2,2)
r = corr(PCsF(:,2),slope')
s = scatter(slope,PCsF(:,2));
text(-0.17,-5,sprintf('Correlation = %0.2f',r),'fontsize',figFS)
ylabel('PC_2')
xlabel('Slope of Normalized August Flow')
title('Interpretation of Mode #2')
set(gca,'fontsize',figFS)

%% large SOM figure

data = complex_coord;

sD = som_data_struct(data);
sM = som_make(sD); %this is the SOM
ny_som = sM.topol.msize(1);
nx_som = sM.topol.msize(2);
en=ny_som*nx_som;

Bmus = som_bmus(sM,sD); %best matching units of each node

% plot colored SOM (with actual hexagonal nodes)
%  and distance among neigbouring nodes
U = som_umat(sM);
Um = U(1:2:size(U,1),1:2:size(U,2));
C = som_colorcode(sM,'rgb2');

myhits = som_hits(sM,sD); %number of nodes in each BMU
[inbNaN dummy]=find(myhits == 0);

C0=C;
C0(inbNaN,1)=1;
C0(inbNaN,2)=1;
C0(inbNaN,3)=1;

%do PCA on the SOM nodes
[eigvec, pcs, eigvals] = pca(sM.codebook);
pc_topo = sqrt(pcs(:,1).^2 + pcs(:,2).^2);

pc_topo_mat = reshape(pc_topo,ny_som,nx_som);

if sum(sum(imag(pc_topo_mat)))==0 %if real pc matrix
    localmax_mat = imregionalmax(pc_topo_mat); %1s where local maxima exist, 0s elsewhere
else
    localmax_mat = imregionalmax(abs(pc_topo_mat));
end

localmax_lin = reshape(localmax_mat,length(pc_topo(:,1)),length(pc_topo(1,:))); %convert to a line for plotting

globalmin = min(pc_topo_mat(:)); %global minimum value in pc topography space
globalmin_mat = zeros(length(pc_topo_mat(:,1)),length(pc_topo_mat(1,:)));
globalmin_mat(find(pc_topo_mat==globalmin)) = 1;
globalmin_lin = reshape(globalmin_mat,length(pc_topo(:,1)),length(pc_topo(1,:))); %convert to a line for plotting

maxmin_mat = globalmin_mat + localmax_mat;
maxmin_lin = globalmin_lin + localmax_lin;

maxmin_patterns = sM.codebook(find(maxmin_lin==1),:); %patterns at the global min/local max locations

figure('Renderer','Painters') %visualize locations of maxima/minima
som_cplane(sM,maxmin_lin/2);
title('Locations of Local Maxima and Global Minima')
set(gca,'xticklabel',[])
colorbar

figure('Position',[0,0,1200,1000])
subplot(1,2,1)
som_cplane(sM,C0);
title('Large SOM')
set(gca,'xticklabel',[])
textx = -0.3;
texty = -0.6;
text(textx,texty,'a)','fontsize',figFS)
set(gca,'FontSize',figFS)

subplot(1,2,2)
som_cplane(sM,abs(pc_topo));
title('PC1^2 + PC2^2 of Large SOM Patterns')
set(gca,'xticklabel',[])
cb = colorbar;
ylabel(cb,'PC_1^2 + PC_2^2')
text(textx,texty-0.1,'b)','fontsize',figFS)
set(gca,'FontSize',figFS)

