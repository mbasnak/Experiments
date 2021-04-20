%group analysis code for the gain change experiment

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');

folderContents = dir(path);

for content = 1:length(folderContents)
   if (contains(folderContents(content).name,'60D05') & ~contains(folderContents(content).name,'20210205'))
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\gain_change_data.mat']);
   end
end


%% Clean and combine data

%remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

for fly = 1:length(data)
   if data(fly).type_of_fly == 1
       %type_1_data{fly} = data(fly).modelTable;
       %create data set restricted to the inverted gain bout
       type_1_data{fly} = data(fly).modelTable(1837:9183,:);
   else
       %type_2_data{fly} = data(fly).modelTable;
       %create data set restricted to the inverted gain bout
       type_2_data{fly} = data(fly).modelTable(1837:9183,:);
   end
end

%remove empty cells
type_1_data = type_1_data(~cellfun('isempty',type_1_data));
type_2_data = type_2_data(~cellfun('isempty',type_2_data));

%combine the tables into one
allType1Data = array2table(zeros(0,12),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','SmoothBarOffsetVariability','SmoothHeadingOffsetVariability','TotalMovement','ZscoredTotalMovement','Time','BumpMagnitude','BumpWidth','ZBumpMagnitude','ZBumpWidth','HeadingVariability'});
flyNumber= [];
for fly = 1:length(type_1_data)
    flyNumber = [flyNumber,repelem(fly,length(type_1_data{1,fly}.BumpMagnitude))];
    allType1Data = [allType1Data;type_1_data{1,fly}]; 
end
allType1Data = addvars(allType1Data,flyNumber');
allType1Data.Properties.VariableNames{'Var13'} = 'Fly';

%for type 2 flies
allType2Data = array2table(zeros(0,12),'VariableNames', {'BarOffsetVariability','HeadingOffsetVariability','SmoothBarOffsetVariability','SmoothHeadingOffsetVariability','TotalMovement','ZscoredTotalMovement','Time','BumpMagnitude','BumpWidth','ZBumpMagnitude','ZBumpWidth','HeadingVariability'});
flyNumber= [];
for fly = 1:length(type_2_data)
    flyNumber = [flyNumber,repelem(fly,length(type_2_data{1,fly}.BumpMagnitude))];
    allType2Data = [allType2Data;type_2_data{1,fly}]; 
end
allType2Data = addvars(allType2Data,flyNumber');
allType2Data.Properties.VariableNames{'Var13'} = 'Fly';



%% Compute model for bump magnitude for all the type 1 flies

%fit mixed linear model using fly number as a random variable
mdl_BM_type1 = fitlme(allType1Data,'BumpMagnitude~HeadingOffsetVariability+ZscoredTotalMovement+(1|Fly)')
%mdl_zBM_type1 = fitlme(allType1Data,'ZBumpMagnitude~HeadingOffsetVariability+ZscoredTotalMovement+Time+(1|Fly)')


%% Repeat for bump width at half max

mdl_BW_type1 = fitlme(allType1Data,'BumpWidth~HeadingOffsetVariability+ZscoredTotalMovement+(1|Fly)')
%mdl_zBW_type1 = fitlme(allType1Data,'ZBumpWidth~HeadingOffsetVariability+ZscoredTotalMovement+Time+(1|Fly)')

%% BM for type 2 flies

mdl_BM_type2 = fitlme(allType2Data,'BumpMagnitude~BarOffsetVariability+ZscoredTotalMovement+(1|Fly)')
%mdl_zBM_type2 = fitlme(allType2Data,'ZBumpMagnitude~BarOffsetVariability+ZscoredTotalMovement+Time+(1|Fly)')

%% Repeat for bump width at half max

mdl_BW_type2 = fitlme(allType2Data,'BumpWidth~BarOffsetVariability+ZscoredTotalMovement+(1|Fly)')
%mdl_zBW_type2 = fitlme(allType2Data,'ZBumpWidth~BarOffsetVariability+ZscoredTotalMovement+Time+(1|Fly)')

%% Plot bump magnitude and half width as a function of offset variability per fly and on average

%bin the data
figure,
for fly = 1:length(type_1_data)
    subplot(1,2,1)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allBumpMag = type_1_data{1,fly}.BumpMagnitude;
    heading_offset_variability = type_1_data{1,fly}.HeadingOffsetVariability;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = median(allBumpMag((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    ylim([0 1.5]);
    title('Bump magnitude');

subplot(1,2,2)
%for fly = 1:length(type_1_data)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allHalfWidth = type_1_data{1,fly}.BumpWidth;
    
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = median(allHalfWidth((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Heading offset variability');
    ylim([2 3.5]);
    title('Bump half width');
end
plot(mvtAxes,mean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,mean(meanBin),'k','linewidth',2)

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParamVsHeadingVariability.png')

%% Repeat for type 2 flies with bar offset variability

%bin the data
figure,
for fly = 1:length(type_2_data)
    subplot(1,2,1)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allBumpMag = type_2_data{1,fly}.BumpMagnitude;
    bar_offset_variability = type_2_data{1,fly}.BarOffsetVariability;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = median(allBumpMag((bar_offset_variability > Bins(bin)) & (bar_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
    ylim([0 1.5]);
    title('Bump magnitude');

    subplot(1,2,2)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allHalfWidth = type_2_data{1,fly}.BumpWidth;
    
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = median(allHalfWidth((bar_offset_variability > Bins(bin)) & (bar_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Bar offset variability');
    ylim([2 3.5]);
    title('Bump half width');
end
plot(mvtAxes,mean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,mean(meanBin),'k','linewidth',2)

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParamVsBarVariability.png')

%% Repeat z-scoring the bump mag data for type 1 flies

figure,
for fly = 1:length(type_1_data)
    subplot(1,2,1)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allBumpMag = zscore(type_1_data{1,fly}.BumpMagnitude);
    heading_offset_variability = type_1_data{1,fly}.HeadingOffsetVariability;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = median(allBumpMag((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Zscored bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    title('Zscored bump magnitude');

    subplot(1,2,2)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allHalfWidth = zscore(type_1_data{1,fly}.BumpWidth);
    
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = median(allHalfWidth((heading_offset_variability > Bins(bin)) & (heading_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Zscored bump half width (EB wedges)'); xlabel('Heading offset variability');
    %ylim([2 3.5]);
    title('Zscored bump half width');
end
plot(mvtAxes,mean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,mean(meanBin),'k','linewidth',2)

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsHeadingVariability.png')


%% Repeat for type 2 flies

figure,
for fly = 1:length(type_2_data)
    subplot(1,2,1)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allBumpMag = zscore(type_2_data{1,fly}.BumpMagnitude);
    bar_offset_variability = type_2_data{1,fly}.BarOffsetVariability;
    
    for bin = 1:length(Bins)-1
        meanBin(fly,bin) = median(allBumpMag((bar_offset_variability > Bins(bin)) & (bar_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBin(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel({'Zscored bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
    %ylim([1 3]);
    title('Bump magnitude');

    subplot(1,2,2)
    nbins = 10;
    max_bin = 1;
    min_bin = 0.2;
    binWidth = max_bin/nbins;
    Bins = [0.2:binWidth:max_bin];
    
    %getting binned means
    allHalfWidth = zscore(type_2_data{1,fly}.BumpWidth);
    
    for bin = 1:length(Bins)-1
        meanBinw(fly,bin) = median(allHalfWidth((bar_offset_variability > Bins(bin)) & (bar_offset_variability < Bins(bin+1))));
    end
    
    %create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    mvtAxes(end) = mvtAxes(end-1)+binWidth; 
    plot(mvtAxes,meanBinw(fly,:),'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width (EB wedges)'); xlabel('Bar offset variability');
    %ylim([2 3.5]);
    title('Bump half width');
end
plot(mvtAxes,mean(meanBinw),'k','linewidth',2)
subplot(1,2,1)
plot(mvtAxes,mean(meanBin),'k','linewidth',2)

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zBumpParamVsBarVariability.png')

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
    clusterDataHeading{fly,1} = mean(data(fly).modelTable.HeadingVariability(1:1836));
    clusterDataOffset{fly,1} = mean(data(fly).modelTable.BarOffsetVariability(1:1836));
    clusterDataMvt{fly,1} = mean(data(fly).modelTable.TotalMovement(1:1836));
    clusterDataBM{fly,1} = mean(data(fly).modelTable.BumpMagnitude(1:1836));
    clusterDataBW{fly,1} = mean(data(fly).modelTable.BumpWidth(1:1836));
    clusterDataAll{fly,1} = [mean(data(fly).modelTable.BarOffsetVariability(1:1836)),mean(data(fly).modelTable.HeadingVariability(1:1836)),mean(data(fly).modelTable.TotalMovement(1:1836)),mean(data(fly).modelTable.BumpMagnitude(1:1836)),mean(data(fly).modelTable.BumpWidth(1:1836))];
    clusterDataMvtVars{fly,1} = [mean(data(fly).modelTable.HeadingVariability(1:1836)),mean(data(fly).modelTable.TotalMovement(1:1836))];
    clusterDataBumpPar{fly,1} = [mean(data(fly).modelTable.BumpMagnitude(1:1836)),mean(data(fly).modelTable.BumpWidth(1:1836))];    
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


%%  Look at mean bump width in type I flies and type II flies in the first bout

figure,
subplot(1,2,1)
boxplot(clusterdataBM,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterdataBM,[],[.5 .5 .5],'filled')
ylabel('Bump magnitude');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 1.8]);

subplot(1,2,2)
boxplot(clusterdataBW,type_of_fly,'color','k')
hold on
scatter(type_of_fly,clusterdataBW,[],[.5 .5 .5],'filled')
ylabel('Bump width');
set(findobj(gca,'type','line'),'linew',2)
xlabel('Type of fly');
ylim([0 4]);

%save
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\bumpParAsPredictors.png')