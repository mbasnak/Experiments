%Group analysis code for the gain change experiment

clear all; close all; clc;

%% Load data

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data';

folderContents = dir(path);

for content = 1:length(folderContents)
   if (contains(folderContents(content).name,'60D05'))
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\gain_change_data.mat']);
   end
end

%% Clean and combine data

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%Make the moving variable a numerical variable for concatenation
for fly = 1:length(data)
    for window = 1:6
        data(fly).modelTable{1,window}.Moving = double(data(fly).modelTable{1,window}.Moving);
    end
end

for fly = 1:length(data)
   if data(fly).type_of_fly == 1
       type_1_data{fly} = data(fly).modelTable;
   else
       type_2_data{fly} = data(fly).modelTable;
   end
   all_flies_data{fly} = data(fly).modelTable;
end

%Remove empty cells
type_1_data = type_1_data(~cellfun('isempty',type_1_data));
type_2_data = type_2_data(~cellfun('isempty',type_2_data));
all_flies_data = all_flies_data(~cellfun('isempty',all_flies_data));

%Combine the tables into one
for window = 1:6
    flyNumber = [];
    allType1Data{window} = array2table(zeros(0,10),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','TotalMovement','Moving','YawSpeed','BumpMagnitude','BumpWidth','Rsq','HeadingOffset','BarOffset'});
    for fly = 1:length(type_1_data)
        flyNumber = [flyNumber,repelem(fly,length(type_1_data{1,fly}{1,window}.BumpMagnitude))];
        allType1Data{window} = [allType1Data{window};type_1_data{1,fly}{1,window}];
    end
    allType1Data{window} = addvars(allType1Data{window},flyNumber');
    allType1Data{window}.Properties.VariableNames{'Var11'} = 'Fly';
end

%For type 2 flies
for window = 1:6
    flyNumber = [];
    allType2Data{window} = array2table(zeros(0,10),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','TotalMovement','Moving','YawSpeed','BumpMagnitude','BumpWidth','Rsq','HeadingOffset','BarOffset'});
    for fly = 1:length(type_2_data)
        flyNumber = [flyNumber,repelem(fly,length(type_2_data{1,fly}{1,window}.BumpMagnitude))];
        allType2Data{window} = [allType2Data{window};type_2_data{1,fly}{1,window}];
    end
    allType2Data{window} = addvars(allType2Data{window},flyNumber');
    allType2Data{window}.Properties.VariableNames{'Var11'} = 'Fly';
end

%Combine all data
for window = 1:6
    flyNumber = [];
    allData{window} = array2table(zeros(0,10),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','TotalMovement','Moving','YawSpeed','BumpMagnitude','BumpWidth','Rsq','HeadingOffset','BarOffset'});
    for fly = 1:length(data)
        flyNumber = [flyNumber,repelem(fly,length(all_flies_data{1,fly}{1,window}.BumpMagnitude))];
        allData{window} = [allData{window};all_flies_data{1,fly}{1,window}];
    end
    allData{window} = addvars(allData{window},flyNumber');
    allData{window}.Properties.VariableNames{'Var11'} = 'Fly';
end

%% Compute model for bump parameters as a function of bar offset variability for all flies

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BM{window} = fitlme(allData{window},'BumpMagnitude~BarOffsetVariability+TotalMovement+(1|Fly)');
    mdl_BW{window} = fitlme(allData{window},'BumpWidth~BarOffsetVariability+TotalMovement+(1|Fly)');    
    %Model Rsquared
    Rsquared_BM(window) = mdl_BM{window}.Rsquared.Adjusted;
    Rsquared_BW(window) = mdl_BW{window}.Rsquared.Adjusted;    
end

figure,
subplot(1,2,1)
plot(Rsquared_BM,'-o')
title('Bump magnitude model fit with total mvt');
ylabel('Rsquared');
xlabel('window #');

subplot(1,2,2)
plot(Rsquared_BW,'-o')
title('Bump width model fit with total mvt');
ylabel('Rsquared');
xlabel('window #');

%Get window for max Rsquared model
[m I] = max(Rsquared_BM);
%mdl_BM{I}
I = I;

[m Iw] = max(Rsquared_BW);
Iw = Iw;


%% Plot bump parameters vs bar offset variability across flies for all window sizes

clear meanBin
clear meanBinw

%define bin limits
allBarOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(all_flies_data)
        allBarOffsetVariability = [allBarOffsetVariability;all_flies_data{1,fly}{1,window}.BarOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(all_flies_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(all_flies_data{1,fly}{1,window}.BumpWidth)];
    end
end
%Define bins
nbins = 10;
max_bin = prctile(allBarOffsetVariability,95,'all');
min_bin = prctile(allBarOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for window = 1:6
    
    figure,
    
    for fly = 1:length(all_flies_data)
        
        subplot(1,2,1)
        
        %Getting binned means
        allBumpMag = zscore(all_flies_data{1,fly}{1,window}.BumpMagnitude);
        bar_offset_variability = all_flies_data{1,fly}{1,window}.BarOffsetVariability;
        total_movement = all_flies_data{1,fly}{1,window}.TotalMovement;
        adj_rs = all_flies_data{1,fly}{1,window}.Rsq;
        mvt_thresh = 20;
        nbins = 10;
        
        for bin = 1:length(Bins)-1
            meanBin(fly,bin) = nanmean(allBumpMag((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (adj_rs > 0.5) & (bar_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel('Zscored bump magnitude'); xlabel('Bar offset variability');
        ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
        xlim([mvtAxes(1) mvtAxes(end)]);
        
        subplot(1,2,2)
        
        allHalfWidth = zscore(all_flies_data{1,fly}{1,window}.BumpWidth);
        for bin = 1:length(Bins)-1
            meanBinw(fly,bin) = nanmean(allHalfWidth((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (adj_rs > 0.5) & (bar_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel('Zscored bump half width'); xlabel('Bar offset variability');
        ylim([min(min(meanBin))-0.5 max(max(meanBinw))+0.5]);
        xlim([mvtAxes(1) mvtAxes(end)]);
    end
    
    plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
    subplot(1,2,1)
    plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)
    
    suptitle(['Window = ',num2str(window)]);
 
end


%% Plot bump parameters vs heading offset variability across flies for all window sizes

clear meanBin
clear meanBinw

%define bin limits
allHeadingOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(all_flies_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;all_flies_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(all_flies_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(all_flies_data{1,fly}{1,window}.BumpWidth)];
    end
end
%Define bins
nbins = 10;
max_bin = prctile(allHeadingOffsetVariability,95,'all');
min_bin = prctile(allHeadingOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for window = 1:6
    
    figure,
    
    for fly = 1:length(all_flies_data)
        
        subplot(1,2,1)
        
        %Getting binned means
        allBumpMag = zscore(all_flies_data{1,fly}{1,window}.BumpMagnitude);
        heading_offset_variability = all_flies_data{1,fly}{1,window}.HeadingOffsetVariability;
        total_movement = all_flies_data{1,fly}{1,window}.TotalMovement;
        adj_rs = all_flies_data{1,fly}{1,window}.Rsq;
        
        mvt_thresh = 20;
        nbins = 10;
        
        for bin = 1:length(Bins)-1
            meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (adj_rs > 0.5) & (heading_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel('Zscored bump magnitude'); xlabel('heading offset variability');
        ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
        xlim([mvtAxes(1) mvtAxes(end)]);
        
        subplot(1,2,2)
        
        allHalfWidth = zscore(all_flies_data{1,fly}{1,window}.BumpWidth);
        for bin = 1:length(Bins)-1
            meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (adj_rs > 0.5) & (heading_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel('Zscored bump half width'); xlabel('heading offset variability');
        ylim([min(min(meanBin))-0.5 max(max(meanBinw))+0.5]);
        xlim([mvtAxes(1) mvtAxes(end)]);
    end
    
    plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
    subplot(1,2,1)
    plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)
    
    suptitle(['Window = ',num2str(window)]);
 
end


%% Cluster data to see if it can be naturally divided into type I and type II flies

%I will use for each fly the following criteria:
%offset_var
%heading_var
%total_mvt
%mean bump mag
%mean bump width at half max

%create array
for fly = 1:length(data)
    
    type_of_fly(fly) = data(fly).type_of_fly;
    gof = data(fly).modelTableNG{1,1}.Rsq(1:1836);
    
    heading_var = data(fly).modelTableNG{1,1}.HeadingVariability(1:1836);
    clusterDataHeading(fly) = mean(heading_var(gof>0.5));
    
    BM = data(fly).modelTableNG{1,1}.BumpMagnitude(1:1836);
    clusterDataBM(fly) = mean(BM);
    
    BW = data(fly).modelTableNG{1,1}.BumpWidth(1:1836);
    clusterDataBW(fly) = mean(BW);
    
    offset = deg2rad(data(fly).modelTableNG{1,1}.Offset(1:1836));
    [~,clusterDataOffsetVar(fly)] = circ_std(offset);
    
end

%  Look at mean bump parameters in type I flies and type II flies in the first bout

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterDataBM,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
%ylim([0 1.8]);

subplot(1,3,2)
boxplot(clusterDataBW,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterDataBW,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterDataOffsetVar,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsPredictors_nothresh.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_strategy_nothresh.png')

%% Repeat, thresholding by goodness of fit

for fly = 1:length(data)
    
    type_of_fly(fly) = data(fly).type_of_fly;
    gof = data(fly).modelTableNG{1,1}.Rsq(1:1836);
    
    heading_var = data(fly).modelTableNG{1,1}.HeadingVariability(1:1836);
    clusterDataHeading(fly) = mean(heading_var(gof>0.5));
    
    BM = data(fly).modelTableNG{1,1}.BumpMagnitude(1:1836);
    clusterDataBM_thresh(fly) = mean(BM(gof>0.5));
    
    BW = data(fly).modelTableNG{1,1}.BumpWidth(1:1836);
    clusterDataBW_thresh(fly) = mean(BW(gof>0.5));
    
    offset = deg2rad(data(fly).modelTableNG{1,1}.Offset(1:1836));
    [~,clusterDataOffsetVar_thresh(fly)] = circ_std(offset(gof>0.5));
    
end

%  Look at mean bump parameters in type I flies and type II flies in the first bout

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM_thresh,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterDataBM_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
%ylim([0 1.8]);

subplot(1,3,2)
boxplot(clusterDataBW_thresh,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterDataBW_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar_thresh,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterDataOffsetVar_thresh,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsPredictors_thresh.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_strategy_thresh.png')

%% Define if the fly is learning or not based on the ratio of heading offset variabilities
    
heading_offset_var_ratio = [];

for window = 1:6
    
    for fly = 1:length(data)
        
        quarter_period = floor(length(data(fly).modelTable{1,window}.HeadingOffsetVariability)/4);
        heading_offset_var_ratio(fly,window) = mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(1:quarter_period))/mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(end-quarter_period:end));        
        
    end
    
end

mean_heading_offset_var_ratio = mean(heading_offset_var_ratio,2);

learning = ones(1,length(data));
learning(mean_heading_offset_var_ratio < 1.15) = 2;


%  Look at mean bump parameters in type I flies and type II flies in the first bout

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM,learning,'color','k')
hold on
scatter(learning,clusterDataBM,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});

subplot(1,3,2)
boxplot(clusterDataBW,learning,'color','k')
hold on
scatter(learning,clusterDataBW,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar,learning,'color','k')
hold on
scatter(learning,clusterDataOffsetVar,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);
xticks([1 2]);
xticklabels({'learners','non-learners'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsLearningPredictors_nothresh.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_learning_nothresh.png')

%% Repeat thresholding by gof


figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM_thresh,learning,'color','k')
hold on
scatter(learning,clusterDataBM_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});

subplot(1,3,2)
boxplot(clusterDataBW_thresh,learning,'color','k')
hold on
scatter(learning,clusterDataBW_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar_thresh,learning,'color','k')
hold on
scatter(learning,clusterDataOffsetVar_thresh,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);
xticks([1 2]);
xticklabels({'learners','non-learners'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsLearningPredictors_thresh.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_learning_thresh.png')

%% Defining the learners manually by what I see in the raw data, wihtout thresholding


learning_manual = 2*ones(1,length(data));
learning_manual(1,[1,2,3,6,8,13,15]) = 1;

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM,learning_manual,'color','k')
hold on
scatter(learning_manual,clusterDataBM,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});

subplot(1,3,2)
boxplot(clusterDataBW,learning_manual,'color','k')
hold on
scatter(learning_manual,clusterDataBW,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar,learning_manual,'color','k')
hold on
scatter(learning_manual,clusterDataOffsetVar,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);
xticks([1 2]);
xticklabels({'learners','non-learners'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsLearningPredictors_nothresh_manual.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_learning_nothresh_manual.png')


%% Repeat thresholding

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM_thresh,learning_manual,'color','k')
hold on
scatter(learning_manual,clusterDataBM_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});

subplot(1,3,2)
boxplot(clusterDataBW_thresh,learning_manual,'color','k')
hold on
scatter(learning_manual,clusterDataBW_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar_thresh,learning_manual,'color','k')
hold on
scatter(learning_manual,clusterDataOffsetVar_thresh,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);
xticks([1 2]);
xticklabels({'learners','non-learners'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsLearningPredictors_thresh_manual.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_learning_thresh_manual.png')

%% Classify the flies using the ratio of circular standard deviations of the heading offset (rather than the mean of heading offset var)

for fly = 1:length(data)
    
    
    %get offset, gof and total movement
    heading_offset_initial_quarter = deg2rad(data(fly).modelTable{1,1}.HeadingOffset(1:quarter_period));
    heading_offset_final_quarter = deg2rad(data(fly).modelTable{1,1}.HeadingOffset(end-quarter_period:end));
    rsq_initial_quarter = data(fly).modelTable{1,1}.Rsq(1:quarter_period);
    rsq_final_quarter = data(fly).modelTable{1,1}.Rsq(end-quarter_period:end);
    mvt_initial_quarter = data(fly).modelTable{1,1}.TotalMovement(1:quarter_period);
    mvt_final_quarter = data(fly).modelTable{1,1}.TotalMovement(end-quarter_period:end);    
    
    %get offset var
    [~,heading_offset_var_initial_quarter] = circ_std(heading_offset_initial_quarter);
    [~,heading_offset_var_initial_quarter_thresh] = circ_std(heading_offset_initial_quarter(rsq_initial_quarter > 0.5 & mvt_initial_quarter > 25));
    [~,heading_offset_var_final_quarter] = circ_std(heading_offset_final_quarter);
    [~,heading_offset_var_final_quarter_thresh] = circ_std(heading_offset_final_quarter(rsq_final_quarter > 0.5 & mvt_final_quarter > 25));
    
    %get offset var ratios
    heading_offset_var_ratio2(fly) = heading_offset_var_initial_quarter/heading_offset_var_final_quarter;
    heading_offset_var_ratio2_thresh(fly) = heading_offset_var_initial_quarter_thresh/heading_offset_var_final_quarter_thresh;
       
end

learning = ones(1,length(data));
learning(heading_offset_var_ratio2_thresh < 1.15) = 2;


%  Look at mean bump parameters in type I flies and type II flies in the first bout

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM,learning,'color','k')
hold on
scatter(learning,clusterDataBM,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});

subplot(1,3,2)
boxplot(clusterDataBW,learning,'color','k')
hold on
scatter(learning,clusterDataBW,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar,learning,'color','k')
hold on
scatter(learning,clusterDataOffsetVar,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);
xticks([1 2]);
xticklabels({'learners','non-learners'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsLearningPredictors_nothresh_bettermethod.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_learning_nothresh_bettermethod.png')

%% Repeat thresholding by gof

%  Look at mean bump parameters in type I flies and type II flies in the first bout

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterDataBM_thresh,learning,'color','k')
hold on
scatter(learning,clusterDataBM_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});

subplot(1,3,2)
boxplot(clusterDataBW_thresh,learning,'color','k')
hold on
scatter(learning,clusterDataBW_thresh,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
xticks([1 2]);
xticklabels({'learners','non-learners'});
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterDataOffsetVar_thresh,learning,'color','k')
hold on
scatter(learning,clusterDataOffsetVar_thresh,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);
xticks([1 2]);
xticklabels({'learners','non-learners'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsLearningPredictors_thresh_bettermethod.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_learning_thresh_bettermethod.png')



%% Plot parameters in the first bout vs ratio of heading offset variability

for window = 1:6
    
    ratio_offset_var = [];
    
    for fly = 1:length(data)
        
        half_period = floor(length(data(fly).modelTable{1,window}.HeadingOffsetVariability)/2);
        
        ratio_offset_var = [ratio_offset_var,mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(1:half_period))/mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(half_period:end))];        
                
        pre_bump_mag = data(fly).modelTableNG{1,window}.BumpMagnitude(1:1836);
        pre_bump_width = data(fly).modelTableNG{1,window}.BumpWidth(1:1836);
        pre_bump_offset = deg2rad(data(fly).modelTableNG{1,window}.Offset(1:1836));
        
        %threhsolds
        pre_total_mvt = data(fly).modelTableNG{1,window}.TotalMovement(1:1836);
        gof = data(fly).modelTableNG{1,window}.Rsq(1:1836);
        
        pre_mean_BM(fly) = mean(pre_bump_mag(gof >0.5));
        pre_mean_BW(fly) = mean(pre_bump_width(gof >0.5)); 
        [~,mean_offset_var(fly)] = circ_std(pre_bump_offset(gof>0.5));
        
    end
    
    figure('Position',[100 100 1400 600]),
    subplot(1,3,1)
    plot(ratio_offset_var,pre_mean_BM,'ko')
    xlabel('Mean initial offset variability / mean final offset variability');
    ylabel('Bump magnitude in the preceding block')
    ylim([0.5 3]);
    
    subplot(1,3,2)
    plot(ratio_offset_var,pre_mean_BW,'ko')
    xlabel('Mean initial offset variability / mean final offset variability');
    ylabel('Bump width in the preceding block')
    ylim([1.5 2.5]);
    
    subplot(1,3,3)
    plot(ratio_offset_var,mean_offset_var,'ko')
    xlabel('Mean initial offset variability / mean final offset variability');
    ylabel('Offset variability in the preceding block')
    
    suptitle(['Window = ',num2str(window)]);
    
    %saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\corrBumpParOffsetVar_',num2str(window),'.png'])
    
end


%% Repeat using the first and last quarter of the data instead of half to compute the variability ratio

for window = 1:6
    
    ratio_offset_var = [];
    
    for fly = 1:length(data)
        
        quarter_period = floor(length(data(fly).modelTable{1,window}.HeadingOffsetVariability)/4);
        ratio_offset_var = [ratio_offset_var,mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(1:quarter_period))/mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(end-quarter_period:end))];        
        
        pre_bump_mag = data(fly).modelTableNG{1,window}.BumpMagnitude(1:1836);
        pre_bump_width = data(fly).modelTableNG{1,window}.BumpWidth(1:1836);
        pre_bump_offset = deg2rad(data(fly).modelTableNG{1,window}.Offset(1:1836));
        
        %threhsolds
        pre_total_mvt = data(fly).modelTableNG{1,window}.TotalMovement(1:1836);
        gof = data(fly).modelTableNG{1,window}.Rsq(1:1836);
        
        pre_mean_BM(fly) = mean(pre_bump_mag(gof >0.5));
        pre_mean_BW(fly) = mean(pre_bump_width(gof >0.5)); 
        [~,mean_offset_var(fly)] = circ_std(pre_bump_offset(gof>0.5));
        
    end
    
    figure('Position',[100 100 1400 600]),
    subplot(1,3,1)
    plot(ratio_offset_var,pre_mean_BM,'ko')
    xlabel('Mean initial offset variability / mean final offset variability');
    ylabel('Bump magnitude in the preceding block')
    ylim([0.5 3]);
    
    subplot(1,3,2)
    plot(ratio_offset_var,pre_mean_BW,'ko')
    xlabel('Mean initial offset variability / mean final offset variability');
    ylabel('Bump width in the preceding block')
    ylim([1.5 2.5]);
    
    subplot(1,3,3)
    plot(ratio_offset_var,mean_offset_var,'ko')
    xlabel('Mean initial offset variability / mean final offset variability');
    ylabel('Offset variability in the preceding block')
    
    suptitle(['Window = ',num2str(window)]);
    
    %saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\corrBumpParOffsetVar_',num2str(window),'.png'])
    
end

%% Relate to fly strategy

for window = 1:6
    
    ratio_offset_var = [];
    
    for fly = 1:length(data)
        quarter_period = floor(length(data(fly).modelTable{1,window}.HeadingOffsetVariability)/2);
        ratio_offset_var = [ratio_offset_var,mean(data(fly).modelTable{1,window}.BarOffsetVariability(1:end))/mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(1:end))];        
        pre_bump_mag = data(fly).modelTableNG{1,window}.BumpMagnitude(1:1836);
        pre_bump_width = data(fly).modelTableNG{1,window}.BumpWidth(1:1836);
        pre_total_mvt = data(fly).modelTableNG{1,window}.TotalMovement(1:1836);
        pre_mean_BM(fly) = mean(pre_bump_mag(pre_total_mvt>25));
        pre_mean_BW(fly) = mean(pre_bump_width(pre_total_mvt>25)); 
        [~,mean_offset_var(fly)] = circ_std(deg2rad(data(fly).modelTableNG{1,window}.Offset(1:1836)));
    end
    figure('Position',[100 100 1000 800]),
    subplot(1,3,1)
    plot(ratio_offset_var,pre_mean_BM,'ko')
    xlabel('Mean bar offset variability / mean heading offset variability');
    ylabel('Bump magnitude in the preceding block')
    ylim([0.5 3]);
    
    subplot(1,3,2)
    plot(ratio_offset_var,pre_mean_BW,'ko')
    xlabel('Mean bar offset variability / mean heading offset variability');
    ylabel('Bump width in the preceding block')
    ylim([1.5 2.5]);
    
    subplot(1,3,3)
    plot(ratio_offset_var,mean_offset_var,'ko')
    xlabel('Mean bar offset variability / meanheading offset variability');
    ylabel('Offset variability in the preceding block')
    
    suptitle(['Window = ',num2str(window)]);
    
    
end