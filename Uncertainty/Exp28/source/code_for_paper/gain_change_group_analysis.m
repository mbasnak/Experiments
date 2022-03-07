%Group analysis code for the gain change experiment

clear all; close all;

%% Load data

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data';

folderContents = dir(path);

for content = 1:length(folderContents)
   if (contains(folderContents(content).name,'60D05') & ~contains(folderContents(content).name,'20210205'))
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
end

%Remove empty cells
type_1_data = type_1_data(~cellfun('isempty',type_1_data));
type_2_data = type_2_data(~cellfun('isempty',type_2_data));

%Combine the tables into one
for window = 1:6
    flyNumber = [];
    allType1Data{window} = array2table(zeros(0,7),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','TotalMovement','Moving','YawSpeed','BumpMagnitude','BumpWidth'});
    for fly = 1:length(type_1_data)
        flyNumber = [flyNumber,repelem(fly,length(type_1_data{1,fly}{1,window}.BumpMagnitude))];
        allType1Data{window} = [allType1Data{window};type_1_data{1,fly}{1,window}];
    end
    allType1Data{window} = addvars(allType1Data{window},flyNumber');
    allType1Data{window}.Properties.VariableNames{'Var8'} = 'Fly';
end

%For type 2 flies
for window = 1:6
    flyNumber = [];
    allType2Data{window} = array2table(zeros(0,7),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','TotalMovement','Moving','YawSpeed','BumpMagnitude','BumpWidth'});
    for fly = 1:length(type_2_data)
        flyNumber = [flyNumber,repelem(fly,length(type_2_data{1,fly}{1,window}.BumpMagnitude))];
        allType2Data{window} = [allType2Data{window};type_2_data{1,fly}{1,window}];
    end
    allType2Data{window} = addvars(allType2Data{window},flyNumber');
    allType2Data{window}.Properties.VariableNames{'Var8'} = 'Fly';
end

%% Compute model for bump magnitude for all the type 1 flies

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BM_type1{window} = fitlme(allType1Data{window},'BumpMagnitude~HeadingOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BM_type1(window) = mdl_BM_type1{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BM_type1,'-o')
title('Bump magnitude model fit with total mvt');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
[m I1] = max(Rsquared_BM_type1);
mdl_BM_type1{I1}

%% Repeat with yaw speed instead of total movement

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BM_type1{window} = fitlme(allType1Data{window},'BumpMagnitude~HeadingOffsetVariability+YawSpeed+(1|Fly)');
    %Model Rsquared
    Rsquared_BM_type1(window) = mdl_BM_type1{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BM_type1,'-o')
title('Bump magnitude model fit with yaw speed');
ylabel('Rsquared');
xlabel('window #');

% %show coefficients for max Rsquared model
% [m I1] = max(Rsquared_BM_type1);
% mdl_BM_type1{I1}

%% Repeat for bump width at half max

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BW_type1{window} = fitlme(allType1Data{window},'BumpWidth~HeadingOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BW_type1(window) = mdl_BW_type1{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BW_type1,'-o')
title('Bump width model fit with total mvt');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
[m Iw1] = max(Rsquared_BW_type1);
mdl_BW_type1{Iw1}


%% Repeat using yaw speed instead of total movement 

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BW_type1{window} = fitlme(allType1Data{window},'BumpWidth~HeadingOffsetVariability+YawSpeed+(1|Fly)');
    %Model Rsquared
    Rsquared_BW_type1(window) = mdl_BW_type1{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BW_type1,'-o')
title('Bump width model fit with yaw speed');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
% [m Iw1] = max(Rsquared_BW_type1);
% mdl_BW_type1{Iw1}

%% BM for type 2 flies

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BM_type2{window} = fitlme(allType2Data{window},'BumpMagnitude~HeadingOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BM_type2(window) = mdl_BM_type2{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BM_type2,'-o')
title('Bump magnitude');
ylabel('Rsquared');
xlabel('window #');
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\modelFit_type1.png')

%show coefficients for max Rsquared model
[m I2] = max(Rsquared_BM_type2);
mdl_BM_type2{I2}

%% Repeat for bump width at half max

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BW_type2{window} = fitlme(allType2Data{window},'BumpWidth~HeadingOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BW_type2(window) = mdl_BW_type2{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BW_type2,'-o')
title('Bump width');
ylabel('Rsquared');
xlabel('window #');
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\modelFit_type2.png')

%show coefficients for max Rsquared model
[m Iw2] = max(Rsquared_BW_type2);
mdl_BW_type2{Iw2}

%% BM for type 2 flies with bar offset variability

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BM_type2{window} = fitlme(allType2Data{window},'BumpMagnitude~BarOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BM_type2(window) = mdl_BM_type2{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BM_type2,'-o')
title('Bump magnitude');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
[m I2] = max(Rsquared_BM_type2);
mdl_BM_type2{I2}

%% Repeat for bump width at half max

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BW_type2{window} = fitlme(allType2Data{window},'BumpWidth~BarOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BW_type2(window) = mdl_BW_type2{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BW_type2,'-o')
title('Bump width');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
[m Iw2] = max(Rsquared_BW_type2);
mdl_BW_type2{Iw2}

%% Plot bump magnitude and half width as a function of offset variability per fly and on average

figure,

%define bin limits
allHeadingOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_1_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;type_1_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;type_1_data{1,fly}{1,window}.BumpMagnitude];
        allBumpWidth = [allBumpWidth;type_1_data{1,fly}{1,window}.BumpWidth];    
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

for fly = 1:length(type_1_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = type_1_data{1,fly}{1,I1}.BumpMagnitude;
    heading_offset_variability = type_1_data{1,fly}{1,I1}.HeadingOffsetVariability;
    total_movement = type_1_data{1,fly}{1,I1}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    ylim([0 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump magnitude');
    
    subplot(1,2,2)
    
    allHalfWidth = type_1_data{1,fly}{1,Iw1}.BumpWidth;
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Heading offset variability');
    ylim([0 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump half width');
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)
    
%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParamVsHeadingVariability.png')

%% Repeat for type 2 flies with bar offset variability

clear meanBin
clear meanBinw
    
figure,

%define bin limits
allBarOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_2_data)
        allBarOffsetVariability = [allBarOffsetVariability;type_2_data{1,fly}{1,window}.BarOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;type_2_data{1,fly}{1,window}.BumpMagnitude];
        allBumpWidth = [allBumpWidth;type_2_data{1,fly}{1,window}.BumpWidth];    
    end
end
%Define bins
max_bin = prctile(allBarOffsetVariability,95,'all');
min_bin = prctile(allBarOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_2_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = type_2_data{1,fly}{1,I2}.BumpMagnitude;
    bar_offset_variability = type_2_data{1,fly}{1,I2}.BarOffsetVariability;
    total_movement = type_2_data{1,fly}{1,I2}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;   
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
    ylim([0 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump magnitude');
    
    subplot(1,2,2)
    
    allHalfWidth = type_2_data{1,fly}{1,Iw2}.BumpWidth;
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Bar offset variability');
    ylim([0 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump half width');
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParamVsBarVariability.png')

%% Repeat z-scoring the bump pars data for type 1 flies

figure,

%define bin limits
allHeadingOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_1_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;type_1_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(type_1_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(type_1_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allHeadingOffsetVariability,95,'all');
min_bin = prctile(allHeadingOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_1_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(type_1_data{1,fly}{1,I1}.BumpMagnitude);
    heading_offset_variability = type_1_data{1,fly}{1,I1}.HeadingOffsetVariability;
    total_movement = type_1_data{1,fly}{1,I1}.TotalMovement;
    mvt_thresh = 20;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Zscored bump magnitude'}); xlabel('Heading offset variability');
    ylim([-1.5 1.5]);
    xlim([0.3 mvtAxes(end)]);
    
    subplot(1,2,2)
    
    allHalfWidth = zscore(type_1_data{1,fly}{1,Iw1}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Zscored bump width'}); xlabel('Heading offset variability');
    ylim([-1.5 1.5]);
    xlim([0.3 mvtAxes(end)]);
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsHeadingVariability.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\zBumpParamVsHeadingVariability_type1.png')

%% Plot for every window size

clear meanBin
clear meanBinw

%define bin limits
allHeadingOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_1_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;type_1_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(type_1_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(type_1_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allHeadingOffsetVariability,95,'all');
min_bin = prctile(allHeadingOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for window = 1:6
    
    figure,
    
    for fly = 1:length(type_1_data)
        
        subplot(1,2,1)
        
        %Getting binned means
        allBumpMag = zscore(type_1_data{1,fly}{1,window}.BumpMagnitude);
        heading_offset_variability = type_1_data{1,fly}{1,window}.HeadingOffsetVariability;
        total_movement = type_1_data{1,fly}{1,window}.TotalMovement;
        mvt_thresh = 20;
        nbins = 10;
        
        for bin = 1:length(Bins)-1
            meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel({'Zscored bump magnitude'}); xlabel('Heading offset variability');
        ylim([-1.5 1.5]);
        xlim([0.3 mvtAxes(end)]);
        
        subplot(1,2,2)
        
        allHalfWidth = zscore(type_1_data{1,fly}{1,Iw1}.BumpWidth);
        for bin = 1:length(Bins)-1
            meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel({'Zscored bump width'}); xlabel('Heading offset variability');
        ylim([-1.5 1.5]);
        xlim([0.3 mvtAxes(end)]);
    end
    
    plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
    subplot(1,2,1)
    plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)
    
    suptitle(['Window = ',num2str(window)]);
    
    saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsHeadingVariability_',num2str(window),'.png'])
    saveas(gcf,['C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\zBumpParamVsHeadingVariability_type1_',num2str(window),'.png'])
    
end

%% Repeat for type 2 flies

clear meanBin
clear meanBinw

figure,

%define bin limits
allBarOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_2_data)
        allBarOffsetVariability = [allBarOffsetVariability;type_2_data{1,fly}{1,window}.BarOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(type_2_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(type_2_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allBarOffsetVariability,95,'all');
min_bin = prctile(allBarOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_2_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(type_2_data{1,fly}{1,I2}.BumpMagnitude);
    bar_offset_variability = type_2_data{1,fly}{1,I2}.BarOffsetVariability;
    total_movement = type_2_data{1,fly}{1,I2}.TotalMovement;
    mvt_thresh = 20;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Zscored bump magnitude'); xlabel('Bar offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    
    subplot(1,2,2)
    
    allHalfWidth = zscore(type_2_data{1,fly}{1,Iw2}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
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

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsBarVariability.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\zBumpParamVsHeadingVariability_type2.png')

%% Repeat for every window size

clear meanBin
clear meanBinw

%define bin limits
allBarOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_2_data)
        allBarOffsetVariability = [allBarOffsetVariability;type_2_data{1,fly}{1,window}.BarOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(type_2_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(type_2_data{1,fly}{1,window}.BumpWidth)];
    end
end
%Define bins
max_bin = prctile(allBarOffsetVariability,95,'all');
min_bin = prctile(allBarOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for window = 1:6
    
    figure,
    
    for fly = 1:length(type_2_data)
        
        subplot(1,2,1)
        
        %Getting binned means
        allBumpMag = zscore(type_2_data{1,fly}{1,window}.BumpMagnitude);
        bar_offset_variability = type_2_data{1,fly}{1,window}.BarOffsetVariability;
        total_movement = type_2_data{1,fly}{1,window}.TotalMovement;
        mvt_thresh = 20;
        nbins = 10;
        
        for bin = 1:length(Bins)-1
            meanBin(fly,bin) = nanmean(allBumpMag((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
        end
        
        plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
        hold on
        ylabel('Zscored bump magnitude'); xlabel('Bar offset variability');
        ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
        xlim([mvtAxes(1) mvtAxes(end)]);
        
        subplot(1,2,2)
        
        allHalfWidth = zscore(type_2_data{1,fly}{1,Iw2}.BumpWidth);
        for bin = 1:length(Bins)-1
            meanBinw(fly,bin) = nanmean(allHalfWidth((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
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

    %save
    saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsBarVariability_',num2str(window),'.png'])
    saveas(gcf,['C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\zBumpParamVsHeadingVariability_type2_',num2str(window),'.png'])
    
end

%% Repeat combining heading offset variability in type 1 flies and bar offset variability in type 2 flies

clear meanBin
clear meanBinw

figure,

for fly = 1:length(type_1_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(type_1_data{1,fly}{1,I1}.BumpMagnitude);
    bar_offset_variability = type_1_data{1,fly}{1,I1}.BarOffsetVariability;
    total_movement = type_1_data{1,fly}{1,I1}.TotalMovement;
    mvt_thresh = 20;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    
    subplot(1,2,2)
    allHalfWidth = zscore(type_1_data{1,fly}{1,Iw1}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on

end

for fly = 1:length(type_2_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(type_2_data{1,fly}{1,I2}.BumpMagnitude);
    bar_offset_variability = type_2_data{1,fly}{1,I2}.BarOffsetVariability;
    total_movement = type_2_data{1,fly}{1,I2}.TotalMovement;
    mvt_thresh = 20;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly+length(type_1_data),bin) = nanmean(allBumpMag((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly+length(type_1_data),:),'color',[.5 .5 .5])
    hold on
    ylabel('Zscored bump magnitude'); xlabel('Offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    
    subplot(1,2,2)
    
    allHalfWidth = zscore(type_2_data{1,fly}{1,Iw2}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly+length(type_1_data),bin) = nanmean(allHalfWidth((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly+length(type_1_data),:),'color',[.5 .5 .5])
    hold on
    ylabel('Zscored bump half width'); xlabel('Offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsVariability.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\zBumpParamVsHeadingVariability_all_flies.png')

%% Repeat for every window?

%% Look at heading offset variability in type II flies

clear meanBin
clear meanBinw

figure,

%define bin limits
allHeadingOffsetVariability2 = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_2_data)
        allHeadingOffsetVariability2 = [allHeadingOffsetVariability2;type_2_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(type_2_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(type_2_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allHeadingOffsetVariability2,95,'all');
min_bin = prctile(allHeadingOffsetVariability2,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_2_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(type_2_data{1,fly}{1,I2}.BumpMagnitude);
    heading_offset_variability2 = type_2_data{1,fly}{1,I2}.HeadingOffsetVariability;
    total_movement = type_2_data{1,fly}{1,I2}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability2 > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability2 < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump magnitude');
    
    subplot(1,2,2)
    
    allHalfWidth = zscore(type_2_data{1,fly}{1,Iw2}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability2 > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability2 < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Heading offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump half width');
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

%save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsHeadingVariability2.png')

%% Now look at the correlation between bump parameters and bar offset variability for type I flies

clear meanBin
clear meanBinw

figure,

%define bin limits
allBarOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(type_2_data)
        allBarOffsetVariability1 = [allBarOffsetVariability;type_1_data{1,fly}{1,window}.BarOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(type_1_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(type_1_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allBarOffsetVariability1,95,'all');
min_bin = prctile(allBarOffsetVariability1,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_1_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(type_1_data{1,fly}{1,I1}.BumpMagnitude);
    bar_offset_variability1 = type_1_data{1,fly}{1,I1}.BarOffsetVariability;
    total_movement = type_1_data{1,fly}{1,I1}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((bar_offset_variability1 > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability1 < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump magnitude');
    
    subplot(1,2,2)
    
    allHalfWidth = zscore(type_1_data{1,fly}{1,Iw1}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((bar_offset_variability1 > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability1 < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Bar offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump half width');
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

%save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsBarVariability1.png')

%% Look at correlation between bar and heading offset in type I and type II flies in the best fit window

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
for fly = 1:length(type_1_data)
   bar_offset_variability_I = type_1_data{1,fly}{1,I1}.BarOffsetVariability; 
end
plot(heading_offset_variability,bar_offset_variability_I,'o')
corr_type_I = corrcoef(heading_offset_variability,bar_offset_variability_I);
hold on
text(0.8,1.2,['Corr = ',num2str(round(corr_type_I(1,2),2))],'fontweight','bold','fontsize',12);
title('Type I flies');
ylabel('Bar offset variability');
xlabel('Heading offset variability');

subplot(1,2,2)
for fly = 1:length(type_2_data)
   heading_offset_variability_II = type_2_data{1,fly}{1,I2}.HeadingOffsetVariability; 
end
plot(bar_offset_variability,heading_offset_variability_II,'ro')
corr_type_II = corrcoef(heading_offset_variability_II,bar_offset_variability);
hold on
text(0.8,1.2,['Corr = ',num2str(round(corr_type_II(1,2),2))],'fontweight','bold','fontsize',12);
title('Type II flies');
xlabel('Heading offset variability');

%Save figure
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\BarVsHeadingOffset.png')


%% Get relationship between bump parameters and offset during normal gain bouts

for fly = 1:length(data)
       NG_data{fly} = data(fly).modelTableNG;
end

%Remove empty cells
NG_data = NG_data(~cellfun('isempty',NG_data));

%Combine the tables into one
for window = 1:6
    flyNumber = [];
    allNGData{window} = array2table(zeros(0,5),'VariableNames', {'HeadingOffsetVariability','TotalMovement','BumpMagnitude','BumpWidth','HeadingVariability'});
    for fly = 1:length(NG_data)
        flyNumber = [flyNumber,repelem(fly,length(NG_data{1,fly}{1,window}.BumpMagnitude))];
        allNGData{window} = [allNGData{window};NG_data{1,fly}{1,window}];
    end
    allNGData{window} = addvars(allNGData{window},flyNumber');
    allNGData{window}.Properties.VariableNames{'Var6'} = 'Fly';
end

%% Compute model for bump magnitude for the normal gain period

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BM{window} = fitlme(allNGData{window},'BumpMagnitude~HeadingOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BM(window) = mdl_BM{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BM,'-o')
title('Bump magnitude');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
[m I] = max(Rsquared_BM);
mdl_BM{I}

%% Repeat for bump width at half max

%Fit mixed linear model using fly number as a random variable
for window = 1:6
    mdl_BW{window} = fitlme(allNGData{window},'BumpWidth~HeadingOffsetVariability+TotalMovement+(1|Fly)');
    %Model Rsquared
    Rsquared_BW(window) = mdl_BW{window}.Rsquared.Adjusted;
end

figure,
plot(Rsquared_BW,'-o')
title('Bump width');
ylabel('Rsquared');
xlabel('window #');

%show coefficients for max Rsquared model
[m Iw] = max(Rsquared_BW);
mdl_BW{Iw}

%% Plot bump magnitude and half width as a function of offset variability per fly and on average

figure,

%define bin limits
allHeadingOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(NG_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;NG_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;NG_data{1,fly}{1,window}.BumpMagnitude];
        allBumpWidth = [allBumpWidth;NG_data{1,fly}{1,window}.BumpWidth];    
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

for fly = 1:length(NG_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = NG_data{1,fly}{1,I}.BumpMagnitude;
    heading_offset_variability = NG_data{1,fly}{1,I}.HeadingOffsetVariability;
    total_movement = NG_data{1,fly}{1,I}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    ylim([0 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1)+binWidth mvtAxes(end)]);
    title('Bump magnitude');
    
    subplot(1,2,2)
    
    allHalfWidth = NG_data{1,fly}{1,Iw}.BumpWidth;
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Heading offset variability');
    ylim([0 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1)+binWidth mvtAxes(end)]);
    title('Bump half width');
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)
    
%Save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParamVsHeadingVariabilityNG.png')

%% Repeat zscoring data

figure,

%define bin limits
allHeadingOffsetVariability = [];
allBumpMagnitude = [];
allBumpWidth = [];
for window = 1:6
    for fly = 1:length(NG_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;NG_data{1,fly}{1,window}.HeadingOffsetVariability];
        allBumpMagnitude = [allBumpMagnitude;zscore(NG_data{1,fly}{1,window}.BumpMagnitude)];
        allBumpWidth = [allBumpWidth;zscore(NG_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allHeadingOffsetVariability,95,'all');
min_bin = prctile(allHeadingOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(NG_data)
    
    subplot(1,2,1)
    
    %Getting binned means
    allBumpMag = zscore(NG_data{1,fly}{1,I}.BumpMagnitude);
    heading_offset_variability = NG_data{1,fly}{1,I}.HeadingOffsetVariability;
    total_movement = NG_data{1,fly}{1,I}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allBumpMag((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1)+binWidth mvtAxes(end)]);
    title('Zscored bump magnitude');
    
    subplot(1,2,2)
    
    allHalfWidth = zscore(NG_data{1,fly}{1,Iw}.BumpWidth);
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = nanmean(allHalfWidth((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Heading offset variability');
    ylim([min(min(meanBinw))-0.5 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1)+binWidth mvtAxes(end)]);
    title('Zscored bump half width');
end

plot(mvtAxes,nanmean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsHeadingVariabilityNG.png')

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
    clusterDataHeading{fly,1} = mean(data(fly).modelTableNG{1,I}.HeadingVariability(1:1836));
    clusterDataOffset{fly,1} = mean(data(fly).modelTableNG{1,I}.HeadingOffsetVariability(1:1836));
    clusterDataMvt{fly,1} = mean(data(fly).modelTableNG{1,I}.TotalMovement(1:1836));
    clusterDataBM{fly,1} = mean(data(fly).modelTableNG{1,I}.BumpMagnitude(1:1836));
    clusterDataBW{fly,1} = mean(data(fly).modelTableNG{1,I}.BumpWidth(1:1836));
    clusterDataAll{fly,1} = [mean(data(fly).modelTableNG{1,I}.HeadingOffsetVariability(1:1836)),mean(data(fly).modelTableNG{1,I}.HeadingVariability(1:1836)),mean(data(fly).modelTableNG{1,I}.TotalMovement(1:1836)),mean(data(fly).modelTableNG{1,I}.BumpMagnitude(1:1836)),mean(data(fly).modelTableNG{1,I}.BumpWidth(1:1836))];
    clusterDataMvtVars{fly,1} = [mean(data(fly).modelTableNG{1,I}.HeadingVariability(1:1836)),mean(data(fly).modelTableNG{1,I}.TotalMovement(1:1836))];
    clusterDataBumpPar{fly,1} = [mean(data(fly).modelTableNG{1,I}.BumpMagnitude(1:1836)),mean(data(fly).modelTableNG{1,I}.BumpWidth(1:1836))];    
end

clusterdataHeading = cell2mat(clusterDataHeading);
%run kmeans clustering with the data, defining two clusters
idx = kmeans(clusterdataHeading,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataOffset = cell2mat(clusterDataOffset);
idx = kmeans(clusterdataOffset,2)
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataMvt = cell2mat(clusterDataMvt);
idx = kmeans(clusterdataMvt,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataBM = cell2mat(clusterDataBM);
idx = kmeans(clusterdataBM,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataBW = cell2mat(clusterDataBW);
idx = kmeans(clusterdataBW,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataAll = cell2mat(clusterDataAll);
idx = kmeans(clusterdataAll,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataMvtVars = cell2mat(clusterDataMvtVars);
idx = kmeans(clusterdataMvtVars,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


clusterdataBumpPar = cell2mat(clusterDataBumpPar);
idx = kmeans(clusterdataBumpPar,2);
opposite_idx(idx==1) = 2;
opposite_idx(idx==2) = 1;
%compare the result from the kmeans with the type of fly
differences1 = type_of_fly - idx';
differences2 = type_of_fly - opposite_idx;
%percentage correct
percent_correct1 = (length(idx)-sum(abs(differences1)))*100/length(idx);
percent_correct2 = (length(idx)-sum(abs(differences2)))*100/length(idx);
percent_correct = max(percent_correct1,percent_correct2)


%%  Look at mean bump parameters in type I flies and type II flies in the first bout

figure('Position',[100 100 1200 800]),
subplot(1,3,1)
boxplot(clusterdataBM,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterdataBM,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
%ylim([0 1.8]);

subplot(1,3,2)
boxplot(clusterdataBW,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterdataBW,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 4]);

subplot(1,3,3)
boxplot(clusterdataOffset,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterdataOffset,[],[.5 .5 .5],'filled')
ylabel('Offset variability');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 2]);

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsPredictors.png')
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\InvertedGain-Experiment\predictors_of_strategy.png')

%% Plot bump parameters in the first bout vs ratio of offset variabilities

for window = 1:6
    %focus on last 1/4 of the block
    ratio_offset_var_final = [];
    
    for fly = 1:length(data)
        start_of_period = floor(length(data(fly).modelTable{1,window}.BarOffsetVariability)/4);
        ratio_offset_var_final = [ratio_offset_var_final,mean(data(fly).modelTable{1,window}.BarOffsetVariability(end-start_of_period:end))/mean(data(fly).modelTable{1,window}.HeadingOffsetVariability(end-start_of_period:end))];        
    end
    figure('Position',[100 100 1000 800]),
    subplot(1,2,1)
    plot(ratio_offset_var_final(type_of_fly == 1),clusterdataBM(type_of_fly == 1),'ko')
    hold on
    plot(ratio_offset_var_final(type_of_fly == 2),clusterdataBM(type_of_fly == 2),'ko')
    xlabel('Mean bar offset variability / mean heading offset variability');
    ylabel('Bump magnitude in the preceding block')
    ylim([0.5 3]);
    
    subplot(1,2,2)
    plot(ratio_offset_var_final(type_of_fly == 1),clusterdataBW(type_of_fly == 1),'ko')
    hold on
    plot(ratio_offset_var_final(type_of_fly == 2),clusterdataBW(type_of_fly == 2),'ko')
    xlabel('Mean bar offset variability / mean heading offset variability');
    ylabel('Bump width in the preceding block')
    ylim([1.5 2.5]);
    
    suptitle(['Window = ',num2str(window)]);
    
    saveas(gcf,['Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\corrBumpParOffsetVar_',num2str(window),'.png'])
    
end


%% Using the fist half of the inverted gain bout to compute the bump parameters
%Relate bump parameters in first half of inverted gain bout with learning

ratio_offset_var_final = [];
bump_mag_pre = [];
bump_width_pre = [];

for fly = 1:length(data)
    start_of_period = floor(length(data(fly).modelTable{1,1}.BarOffsetVariability)/4);
    end_of_period = floor(length(data(fly).modelTable{1,1}.BarOffsetVariability)/2);
    ratio_offset_var_final = [ratio_offset_var_final,mean(data(fly).modelTable{1,1}.BarOffsetVariability(end-start_of_period:end))/mean(data(fly).modelTable{1,1}.HeadingOffsetVariability(end-start_of_period:end))];
    bump_mag_pre = [bump_mag_pre,mean(data(fly).modelTable{1,1}.BumpMagnitude(1:end_of_period))];
    bump_width_pre = [bump_width_pre,mean(data(fly).modelTable{1,1}.BumpWidth(1:end_of_period))];    
end
figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(ratio_offset_var_final,bump_mag_pre,'ko')
xlabel('Mean bar offset variability / mean heading offset variability');
ylabel('Bump magnitude in the first half')
ylim([0.5 3]);

subplot(1,2,2)
plot(ratio_offset_var_final,bump_width_pre,'ko')
xlabel('Mean bar offset variability / mean heading offset variability');
ylabel('Bump width in the first half')
ylim([1.5 2.5]);

%% Divide the inverted gain portion into quartiles and look at bump parameter evolution

quartile_width = floor(length(data(1).modelTable{1,1}.BumpMagnitude)/4);
quartile_changes = [1:quartile_width:length(data(1).modelTable{1,1}.BumpMagnitude)];
warning('off','all')

for quartile = 1:4
    for fly = 1:14
        BM(fly,:) = data(fly).modelTable{1,1}.BumpMagnitude(quartile_changes(quartile):quartile_changes(quartile+1));
        mvt_q = data(fly).modelTable{1,1}.TotalMovement(quartile_changes(quartile):quartile_changes(quartile+1));
        BMq(fly,quartile) = nanmean(BM(fly,mvt_q > mvt_thresh));
        
        BW(fly,:) = data(fly).modelTable{1,1}.BumpWidth(quartile_changes(quartile):quartile_changes(quartile+1));
        BWq(fly,quartile) = nanmean(BW(fly,mvt_q > mvt_thresh));
    end
end

%% Plot parameter evolution for type 1 flies

figure,
subplot(1,2,1)
plot(BMq(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BMq(type_of_fly'==1,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Bump Magnitude');
xlim([0 5]);
ylim([0 1.5]);

subplot(1,2,2)
plot(BWq(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BWq(type_of_fly'==1,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Bump Width');
xlim([0 5]);
ylim([0 3.5]);
            
%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\BumpParQuartilesType1.png')

%% Plot parameter evolution for type 2 flies

figure,
subplot(1,2,1)
plot(BMq(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BMq(type_of_fly'==2,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Bump Magnitude');
xlim([0 5]);
%ylim([0 1.5]);

subplot(1,2,2)
plot(BWq(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BWq(type_of_fly'==2,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Bump Width');
xlim([0 5]);
ylim([0 3.5]);

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\BumpParQuartilesType2.png')

%% Zscore the data

for quartile = 1:4
    for fly = 1:14
        zBM = zscore(data(fly).modelTable{1,1}.BumpMagnitude);
        BMz(fly,:) = zBM(quartile_changes(quartile):quartile_changes(quartile+1));
        mvt_q = data(fly).modelTable{1,1}.TotalMovement(quartile_changes(quartile):quartile_changes(quartile+1));
        BMqz(fly,quartile) = nanmean(BMz(fly,mvt_q > mvt_thresh));
        
        zBW = zscore(data(fly).modelTable{1,1}.BumpWidth);
        BWz(fly,:) = zBW(quartile_changes(quartile):quartile_changes(quartile+1));
        BWqz(fly,quartile) = nanmean(BWz(fly,mvt_q > mvt_thresh));
    end
end

%% Plot parameter evolution for type 1 flies

figure,
subplot(1,2,1)
plot(BMqz(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BMqz(type_of_fly'==1,:)),'-ko','lineWidth',2,'MarkerFaceColor','k')
xticks([1:4])
xlabel('Quartile');
ylabel('Zscored bump magnitude')
xlim([0 5]);

subplot(1,2,2)
plot(BWqz(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BWqz(type_of_fly'==1,:)),'-ko','lineWidth',2,'MarkerFaceColor','k')
xticks([1:4])
xlabel('Quartile');
ylabel('Zscored bump width');
xlim([0 5]);
%ylim([0 3.5]);

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParQuartilesType1.png')
% 
% %% Repeat breaking into different block sizes
% 
% quintile_width = floor(length(data(1).modelTable{1,1}.BumpMagnitude)/5);
% quintile_changes = [1:quintile_width:length(data(1).modelTable{1,1}.BumpMagnitude)];
% 
% for quintile = 1:3
%     for fly = 1:14
%         zBM2 = zscore(data(fly).modelTable{1,1}.BumpMagnitude);
%         BMz2(fly,:) = zBM2(quintile_changes(quintile):quintile_changes(quintile+1));
%         mvt_q2 = data(fly).modelTable{1,1}.TotalMovement(quintile_changes(quintile):quintile_changes(quintile+1));
%         BMqz2(fly,quintile) = nanmean(BMz2(fly,mvt_q2 > mvt_thresh));
%         
%         zBW2 = zscore(data(fly).modelTable{1,1}.BumpWidth);
%         BWz2(fly,:) = zBW2(quintile_changes(quintile):quintile_changes(quintile+1));
%         BWqz2(fly,quintile) = nanmean(BWz2(fly,mvt_q2 > mvt_thresh));
%     end
% end
% 
% figure,
% subplot(1,2,1)
% plot(BMqz2(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
% hold on
% plot(mean(BMqz2(type_of_fly'==1,:)),'-ko','lineWidth',2,'MarkerFaceColor','k')
% xticks([1:4])
% xlabel('Quartile');
% ylabel('Zscored bump magnitude')
% xlim([0 6]);
% 
% subplot(1,2,2)
% plot(BWqz2(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
% hold on
% plot(mean(BWqz2(type_of_fly'==1,:)),'-ko','lineWidth',2,'MarkerFaceColor','k')
% xticks([1:4])
% xlabel('Quartile');
% ylabel('Zscored bump width');
% xlim([0 6]);
% 
% 
% 
% thirds_width = floor(length(data(1).modelTable{1,1}.BumpMagnitude)/3);
% third_changes = [1:2448:length(data(1).modelTable{1,1}.BumpMagnitude)];
% %third_changes = [third_changes,7348];
% 
% for third = 1:3
%     for fly = 1:14
%         zBM3 = zscore(data(fly).modelTable{1,1}.BumpMagnitude);
%         BMz3(fly,:) = zBM3(third_changes(third):third_changes(third+1));
%         mvt_q3 = data(fly).modelTable{1,1}.TotalMovement(third_changes(third):third_changes(third+1));
%         BMqz3(fly,third) = nanmean(BMz3(fly,mvt_q3 > mvt_thresh));
% %         
%         zBW3 = zscore(data(fly).modelTable{1,1}.BumpWidth);
%         BWz3(fly,:) = zBW3(third_changes(third):third_changes(third+1));
%         BWqz3(fly,third) = nanmean(BWz3(fly,mvt_q3 > mvt_thresh));
%     end
% end
% 
% figure,
% subplot(1,2,1)
% plot(BMqz3(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
% hold on
% plot(mean(BMqz3(type_of_fly'==1,:)),'-ko','lineWidth',2,'MarkerFaceColor','k')
% xticks([1:3])
% xlabel('Third');
% ylabel('Zscored bump magnitude')
% xlim([0 4]);
% 
% subplot(1,2,2)
% plot(BWqz3(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
% hold on
% plot(mean(BWqz3(type_of_fly'==1,:)),'-ko','lineWidth',2,'MarkerFaceColor','k')
% xticks([1:3])
% xlabel('Third');
% ylabel('Zscored bump width');
% xlim([0 4]);


%% Plot parameter evolution for type 2 flies

figure,
subplot(1,2,1)
plot(BMqz(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BMqz(type_of_fly'==2,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Zscored bump Magnitude');
xlim([0 5]);
%ylim([0 1.5]);

subplot(1,2,2)
plot(BWqz(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BWqz(type_of_fly'==2,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Zscored bump Width');
xlim([0 5]);
%ylim([0 3.5]);

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParQuartilesType2.png')


%% Divide the inverted gain portion into quartiles and look at bump parameter evolution

quartile_width = floor(length(data(1).modelTable{1,1}.BumpMagnitude)/4);
quartile_changes = [1:quartile_width:length(data(1).modelTable{1,1}.BumpMagnitude)];
warning('off','all')

for quartile = 1:4
    for fly = 1:14
        BM(fly,:) = data(fly).modelTable{1,1}.BumpMagnitude(quartile_changes(quartile):quartile_changes(quartile+1));
        mvt_q = data(fly).modelTable{1,1}.TotalMovement(quartile_changes(quartile):quartile_changes(quartile+1));
        BMq(fly,quartile) = nanmean(BM(fly,mvt_q > mvt_thresh));
        
        BW(fly,:) = data(fly).modelTable{1,1}.BumpWidth(quartile_changes(quartile):quartile_changes(quartile+1));
        BWq(fly,quartile) = nanmean(BW(fly,mvt_q > mvt_thresh));
    end
end

%% Plot parameter evolution for type 1 flies

%Get bump parameters for the two normal gain bouts
for fly = 1:14
    BM_NG_1(fly,:) = data(fly).modelTableNG{1,1}.BumpMagnitude(1:1837);
    BM_NG_2(fly,:) = data(fly).modelTableNG{1,1}.BumpMagnitude(1838:end);
    
    BW_NG_1(fly,:) = data(fly).modelTableNG{1,1}.BumpWidth(1:1837);
    BW_NG_2(fly,:) = data(fly).modelTableNG{1,1}.BumpWidth(1838:end);
end

figure('Position',[100 100 1000 600]),
%Bump magnitude
subplot(1,13,1)
boxplot(mean(BM_NG_1(type_of_fly'==1,:)'),'-o','color',[.5 .5 .5])
ylim([0 1.5]);

subplot(1,13,[2:5])
plot(BMq(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BMq(type_of_fly'==1,:)),'-ko','lineWidth',2)
xticks([1:4])
xlabel('Quartile');
title('Bump Magnitude');
xlim([0 5]);
ylim([0 1.5]);
set(gca,'yticklabel',{[]});

subplot(1,13,6)
boxplot(mean(BM_NG_2(type_of_fly'==1,:)'),'-o','color',[.5 .5 .5])
ylim([0 1.5]);
set(gca,'yticklabel',{[]});

%Bump width
subplot(1,13,8)
boxplot(mean(BW_NG_1(type_of_fly'==1,:)'),'-o','color',[.5 .5 .5])
ylim([0 4]);

subplot(1,13,[9:12])
plot(BWq(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BWq(type_of_fly'==1,:)),'-ko','lineWidth',2)
xticks([1:4])
xlabel('Quartile');
title('Bump Width');
xlim([0 5]);
ylim([0 4]);
set(gca,'yticklabel',{[]});

subplot(1,13,13)
boxplot(mean(BW_NG_2(type_of_fly'==1,:)'),'-o','color',[.5 .5 .5])
ylim([0 4]);
set(gca,'yticklabel',{[]});


%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\BumpParEvolutionType1.png')

%% Repeat for type 2 flies

figure('Position',[100 100 1000 600]),
%Bump magnitude
subplot(1,13,1)
boxplot(mean(BM_NG_1(type_of_fly'==2,:)'),'-o','color',[.5 .5 .5])
ylim([0 1.5]);

subplot(1,13,[2:5])
plot(BMq(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BMq(type_of_fly'==2,:)),'-ko','lineWidth',2)
xticks([1:4])
xlabel('Quartile');
title('Bump Magnitude');
xlim([0 5]);
ylim([0 1.5]);
set(gca,'yticklabel',{[]});

subplot(1,13,6)
boxplot(mean(BM_NG_2(type_of_fly'==2,:)'),'-o','color',[.5 .5 .5])
ylim([0 1.5]);
set(gca,'yticklabel',{[]});

%Bump width
subplot(1,13,8)
boxplot(mean(BW_NG_1(type_of_fly'==2,:)'),'-o','color',[.5 .5 .5])
ylim([0 4]);

subplot(1,13,[9:12])
plot(BWq(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(BWq(type_of_fly'==2,:)),'-ko','lineWidth',2)
xticks([1:4])
xlabel('Quartile');
title('Bump Width');
xlim([0 5]);
ylim([0 4]);
set(gca,'yticklabel',{[]});

subplot(1,13,13)
boxplot(mean(BW_NG_2(type_of_fly'==2,:)'),'-o','color',[.5 .5 .5])
ylim([0 4]);
set(gca,'yticklabel',{[]});


%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\BumpParEvolutionType2.png')

%% Relationship between mvt data and heading offset variability data for type 1 flies

figure,

%define bin limits
allHeadingOffsetVariability = [];
allTotalMovement = [];

for window = 1:6
    for fly = 1:length(type_1_data)
        allHeadingOffsetVariability = [allHeadingOffsetVariability;type_1_data{1,fly}{1,window}.HeadingOffsetVariability];
        allTotalMovement = [allTotalMovement;zscore(type_1_data{1,fly}{1,window}.TotalMovement)];
    end
end
%Define bins
max_bin = prctile(allHeadingOffsetVariability,95,'all');
min_bin = prctile(allHeadingOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_1_data)
    
    %Getting binned means
    allTotalMvt = zscore(type_1_data{1,fly}{1,I1}.TotalMovement);
    heading_offset_variability = type_1_data{1,fly}{1,I1}.HeadingOffsetVariability;
    total_movement = type_1_data{1,fly}{1,I1}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allTotalMvt((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    
    subplot(1,2,1)
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('zscored total movement'); xlabel('Heading offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    
    subplot(1,2,2)
    allYawS = zscore(type_1_data{1,fly}{1,I1}.YawSpeed);
    for bin = 1:length(Bins)-1
        meanBinYaw(fly,bin) = nanmean(allYawS((heading_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (heading_offset_variability < Bins(bin+1))));
    end
    plot(mvtAxes,meanBinYaw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('zscored yaw speed'); xlabel('Heading offset variability');
    ylim([min(min(meanBinYaw))-0.5 max(max(meanBinYaw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    

end
plot(mvtAxes,nanmean(meanBinYaw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\MvtVsHeadingOffset.png')

%% Repeat for type 2 flies

clear meanBin
clear meanBinw
clear meanBinYaw

figure,

%define bin limits
allBarOffsetVariability = [];
allTotalMovement = [];
for window = 1:6
    for fly = 1:length(type_2_data)
        allBarOffsetVariability = [allBarOffsetVariability;type_2_data{1,fly}{1,window}.BarOffsetVariability];
        allTotalMovement = [allTotalMovement;zscore(type_2_data{1,fly}{1,window}.TotalMovement)];
        allBumpWidth = [allBumpWidth;zscore(type_2_data{1,fly}{1,window}.BumpWidth)];    
    end
end
%Define bins
max_bin = prctile(allBarOffsetVariability,95,'all');
min_bin = prctile(allBarOffsetVariability,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

for fly = 1:length(type_2_data)
       
    %Getting binned means
    allTotalMvt = zscore(type_2_data{1,fly}{1,I2}.TotalMovement);
    bar_offset_variability = type_2_data{1,fly}{1,I2}.BarOffsetVariability;
    total_movement = type_2_data{1,fly}{1,I2}.TotalMovement;
    mvt_thresh = 50;
    nbins = 10;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = nanmean(allTotalMvt((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    
    subplot(1,2,1)
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Zscoed total movement'); xlabel('Bar offset variability');
    ylim([min(min(meanBin))-0.5 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);

    subplot(1,2,2)
    allYawS = zscore(type_2_data{1,fly}{1,I2}.YawSpeed);
    for bin = 1:length(Bins)-1
        meanBinYaw(fly,bin) = nanmean(allYawS((bar_offset_variability > Bins(bin)) & (total_movement > mvt_thresh) & (bar_offset_variability < Bins(bin+1))));
    end
    plot(mvtAxes,meanBinYaw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('zscored yaw speed'); xlabel('Bar offset variability');
    ylim([min(min(meanBinYaw))-0.5 max(max(meanBinYaw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    

end
plot(mvtAxes,nanmean(meanBinYaw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,nanmean(meanBin),'k','linewidth',2)

%save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\MvtVsBarOffset.png')


%% Zscore the movement data

for quartile = 1:4
    for fly = 1:14
        zTM = zscore(data(fly).modelTable{1,1}.TotalMovement);
        TMz(fly,:) = zTM(quartile_changes(quartile):quartile_changes(quartile+1));
        mvt_q = data(fly).modelTable{1,1}.TotalMovement(quartile_changes(quartile):quartile_changes(quartile+1));
        TMqz(fly,quartile) = nanmean(TMz(fly,mvt_q > mvt_thresh));
        
        zYS = zscore(data(fly).modelTable{1,1}.YawSpeed);
        YSz(fly,:) = zYS(quartile_changes(quartile):quartile_changes(quartile+1));
        YSqz(fly,quartile) = nanmean(YSz(fly,mvt_q > mvt_thresh));
        
    end
end

%% Plot total movement evolution for all flies

figure,
subplot(1,2,1)
plot(TMqz(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(TMqz(type_of_fly'==1,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Type 1 flies');
xlim([0 5]);
%ylim([0 1.5]);
ylabel('Zscored total movement');

subplot(1,2,2)
plot(TMqz(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(TMqz(type_of_fly'==2,:)),'-ko','lineWidth',2)
%set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Type 2 flies');
xlim([0 5]);

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zTotalMvtQuartiles.png')


%% Plot yaw speed evolution for all flies

figure,
subplot(1,2,1)
plot(YSqz(type_of_fly'==1,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(YSqz(type_of_fly'==1,:)),'-ko','lineWidth',2)
xticks([1:4])
xlabel('Quartile');
title('Type 1 flies');
xlim([0 5]);
ylabel('Zscored yaw speed');

subplot(1,2,2)
plot(YSqz(type_of_fly'==2,:)','-o','color',[.5 .5 .5])
hold on
plot(mean(YSqz(type_of_fly'==2,:)),'-ko','lineWidth',2)
xticks([1:4])
xlabel('Quartile');
title('Type 2 flies');
xlim([0 5]);

%Save
%saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zYawSpeedQuartiles.png')


%% Plot distribution of angular speed under normal and inverted gain conditions

%Combine movement variables
allYawSpeedNG = [];
allYawSpeedIG = [];
for fly = 1:length(data)
   allYawSpeedNG = [allYawSpeedNG,data(fly).modelTableNG{1,1}.TotalMovement];
   allYawSpeedIG = [allYawSpeedIG,data(fly).modelTable{1,1}.TotalMovement];       
end


figure,
subplot(2,1,1)
histogram(allYawSpeedNG);
title('Normal gain');
xlim([0 400]);

subplot(2,1,2)
histogram(allYawSpeedIG);
title('Inverted gain');
xlim([0 350]);
xlabel('Angular speed (deg/s)');