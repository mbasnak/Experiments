%Code for the group analysis of the virtual hallway

clear all; close all;

%prompt the user to open the folder for analysis
cd 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp16\data\experimental flies'
dirName = uigetdir();
cd(dirName)
folderNames = dir;

%load the .mat data files inside the corresponding genotype folder as cells
%of an array
for i = 1:length(folderNames)
    if contains(folderNames(i).name,'2019')
        Data{i} = load(strcat(folderNames(i).folder,'\',folderNames(i).name,'\dataFromAnalysis\data.mat'));
    else
        Data{i} = [];
    end
end

Data = Data(~cellfun('isempty',Data));


%extract the forward velocities
for i = 1:length(Data)
    for j = 1:length(Data{1,i}.smoothed)
        velocities{i,j} = Data{1,i}.smoothed{1,j}.xVel;
        angvelocities{i,j} = Data{1,i}.smoothed{1,j}.angularVel;
    end
    velAroundExpansion{i} = Data{1,i}.velAroundExpansionPoint;
    velAroundPulse{i} = Data{1,i}.velAroundPulse;   
end



%% Plot the fwd vel and ang speed vs the panels y dimension for every fly as a subplot in the same fig
%Not normalized data

% bin per distance step
for i = 1:size(velocities,1)
    for j = 1:size(Data{1,i}.smoothed,2)
        for k = 1:13
            stepMeans{i,j,k} = nanmean(velocities{i,j}(Data{1,i}.step{j,k}));
            stepMeansAng{i,j,k} = nanmean(angvelocities{i,j}(Data{1,i}.step{j,k}));
        end        
    end
end

stepMeansAng = cellfun(@abs,stepMeansAng,'UniformOutput',0);

% get the median and error for each animal
for i = 1:size(velocities,1)
    optoMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
    optoError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
    
    optoAngMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeansAng(i,Data{1,i}.WorkingOptoPulse,:))));
    probeAngMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeansAng(i,Data{1,i}.WorkingProbeTrials,:))));
    optoAngError(i,:) = squeeze(nanstd(cell2mat(stepMeansAng(i,Data{1,i}.WorkingOptoPulse,:))));
    probeAngError(i,:) = squeeze(nanstd(cell2mat(stepMeansAng(i,Data{1,i}.WorkingProbeTrials,:))));
end


%Plot one figure per animal
ydimension = [1:13];

figure('Position',[100 100 1400 800]),
for i = 1:size(velocities,1)
    subplot(2,ceil(size(velocities,1)/2),i)
    boundedline(ydimension,optoMeans(i,:),optoError(i,:)./sqrt(size(optoError,1)),'-ro','alpha');
    hold on
    boundedline(ydimension,probeMeans(i,:),probeError(i,:)./sqrt(size(optoError,1)),'-ko','alpha');
    line([12 12],[0  max(max(probeMeans(i,:)))+1],'Color' ,'r','lineWidth',2);
    ylim([0 max(max(probeMeans(i,:)))+1]); xlim([1 13]);
    title(['Fly ',num2str(i)])
end
suptitle('Median fwd vel vs distance all flies');
[ax,h1]=suplabel('Distance (y panel dimension)','x');
[ax,h2]=suplabel('Forward velocity (mm/s)','y');

saveas(gcf,[dirName,'\globalPlots\AllFliesFwdVelVsPanels.png'])


figure('Position',[100 100 1400 800]),
for i = 1:size(velocities,1)
    subplot(2,ceil(size(velocities,1)/2),i)
    boundedline(ydimension,optoAngMeans(i,:),optoAngError(i,:)./sqrt(size(optoAngError,1)),'-ro','alpha');
    hold on
    boundedline(ydimension,probeAngMeans(i,:),probeAngError(i,:)./sqrt(size(optoAngError,1)),'-ko','alpha');
    line([12 12],[0  max(max(probeAngMeans(i,:)))+10],'Color' ,'r','lineWidth',2);
    ylim([0  max(max(probeAngMeans(i,:)))+10]); xlim([1 13]);
    title(['Fly ',num2str(i)])
end
suptitle('Median ang speed vs distance all flies');
[ax,h1]=suplabel('Distance (y panel dimension)','x');
[ax,h2]=suplabel('Angular speed (deg/s)','y');

saveas(gcf,[dirName,'\globalPlots\AllFliesAngSpeedVsPanels.png'])

%% Normalize the data in different ways

%1) z-scoring
zscoreVel = cellfun(@zscore,velocities,'UniformOutput',0);
zscoreAngVel = cellfun(@zscore,angvelocities,'UniformOutput',0);

% bin per distance step
for i = 1:size(zscoreVel,1)
    for j = 1:size(Data{1,i}.smoothed,2)
        for k = 1:13
            stepMeans{i,j,k} = nanmean(zscoreVel{i,j}(Data{1,i}.step{j,k}));
            stepMeansAng{i,j,k} = nanmean(zscoreAngVel{i,j}(Data{1,i}.step{j,k}));
        end        
    end
end

stepMeansAng = cellfun(@abs,stepMeansAng,'UniformOutput',0);

% get the median and error for each animal
for i = 1:size(zscoreVel,1)
    optoMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
    optoError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
    
    optoAngMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeansAng(i,Data{1,i}.WorkingOptoPulse,:))));
    probeAngMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeansAng(i,Data{1,i}.WorkingProbeTrials,:))));
    optoAngError(i,:) = squeeze(nanstd(cell2mat(stepMeansAng(i,Data{1,i}.WorkingOptoPulse,:))));
    probeAngError(i,:) = squeeze(nanstd(cell2mat(stepMeansAng(i,Data{1,i}.WorkingProbeTrials,:))));
end


%Plot one figure per animal
ydimension = [1:13];

figure,
subplot(1,2,1)
for i = 1:size(optoMeans,1)
    h1 = plot(ydimension,optoMeans(i,:),'-o');
    hold on
    set(h1, 'markerfacecolor', get(h1, 'color'));
    line([12 12],[min(min([probeMeans;optoMeans]))-0.5 max(max([probeMeans;optoMeans]))+0.5],'Color' ,'r','lineWidth',2);

end
title('Opto trials'); xlabel('y dimension');
ylabel('z-scored forward velocity (mm/s)');
ylim([min(min([probeMeans;optoMeans]))-0.5 max(max([probeMeans;optoMeans]))+0.5]);

subplot(1,2,2)
for i = 1:size(probeMeans,1)
    h2 = plot(ydimension,probeMeans(i,:),'-o');
    hold on
    set(h2, 'markerfacecolor', get(h2, 'color'));
    line([12 12],[min(min([probeMeans;optoMeans]))-0.5 max(max([probeMeans;optoMeans]))+0.5],'Color' ,'r','LineStyle', ':','lineWidth',2);
end
title('Probe trials');
xlabel('y dimension');
ylim([min(min([probeMeans;optoMeans]))-0.5 max(max([probeMeans;optoMeans]))+0.5]);

saveas(gcf,[dirName,'\globalPlots\zscoredVelvsDistance.png'])

%% Plotting z-scored angular speed

figure,
subplot(1,2,1)
for i = 1:size(optoMeans,1)
    h1 = plot(ydimension,optoAngMeans(i,:),'-o');
    hold on
    set(h1, 'markerfacecolor', get(h1, 'color')); 
    line([12 12],[min(min([probeAngMeans;optoAngMeans]))-0.5 max(max([probeAngMeans;optoAngMeans]))+0.5],'Color' ,'r','lineWidth',2);

end
title('Opto trials'); xlabel('y dimension');
ylabel('z-scored angular speed (deg/s)');
ylim([min(min([probeAngMeans;optoAngMeans]))-0.5 max(max([probeAngMeans;optoAngMeans]))+0.5]);

subplot(1,2,2)
for i = 1:size(probeAngMeans,1)
    h2 = plot(ydimension,probeAngMeans(i,:),'-o');
    hold on
    set(h2, 'markerfacecolor', get(h2, 'color'));
    line([12 12],[min(min([probeAngMeans;optoAngMeans]))-0.5 max(max([probeAngMeans;optoAngMeans]))+0.5],'Color' ,'r','LineStyle', ':','lineWidth',2);
end
title('Probe trials');
xlabel('y dimension');
ylim([min(min([probeAngMeans;optoAngMeans]))-0.5 max(max([probeAngMeans;optoAngMeans]))+0.5]);

saveas(gcf,[dirName,'\globalPlots\zscoredAngSpeedvsDistance.png'])

%% Max-min normalization

%(value-min)/(max-min)
maxVels = cellfun(@max,velocities,'UniformOutput',0);
minVels = cellfun(@min, velocities,'UniformOutput',0);

for i = 1:size(velocities,1)
    for j = 1:size(velocities,2)
        normVel{i,j} = (velocities{i,j}-minVels{i,j})/(maxVels{i,j}-minVels{i,j});
    end
end

% bin per distance step
for i = 1:size(zscoreVel,1)
    for j = 1:size(Data{1,i}.smoothed,2)
        for k = 1:13
            stepMeans{i,j,k} = nanmean(normVel{i,j}(Data{1,i}.step{j,k}));
        end        
    end
end

% get the median and error for each animal
for i = 1:size(zscoreVel,1)
    optoMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
    optoError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
end


%Plot one figure per animal
ydimension = [1:13];

figure,
subplot(1,2,1)
for i = 1:size(optoMeans,1)
    h1 = plot(ydimension,optoMeans(i,:),'-o');
    hold on
    set(h1, 'markerfacecolor', get(h1, 'color'));
    line([12 12],[0 1],'Color' ,'r','lineWidth',2);

end
title('Opto trials'); xlabel('y dimension');
ylabel('Max-min normalized forward velocity (mm/s)');
ylim([0 1]);

subplot(1,2,2)
for i = 1:size(probeMeans,1)
    h2 = plot(ydimension,probeMeans(i,:),'-o');
    hold on
    set(h2, 'markerfacecolor', get(h2, 'color'));
    line([12 12],[0 1],'Color' ,'r','lineWidth',2,'lineStyle',':');
end
title('Probe trials');
xlabel('y dimension');
ylim([0 1]);

saveas(gcf,[dirName,'\globalPlots\MaxMinVelvsDistance.png'])

%% Profile plots around the pulse

ystim = 12;
  
figure('Position',[100 100 1400 800])

for i = 1:length(Data)
    subplot(1,length(Data),i)
    plot([ystim-1:ystim+1],cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingOptoPulse,ystim-1:ystim+1))),'color',[0.5 0.5 0.5])
    hold on
    errorbar([ystim-1,ystim,ystim+1],nanmean(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingOptoPulse,ystim-1:ystim+1)))),nanstd(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingOptoPulse,ystim-1:ystim+1))))/sqrt(length(Data{1,i}.WorkingOptoPulse)),'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
    xlim([ystim-2 ystim+2]);
    xticks([ystim-1 ystim ystim+1])
    xticklabels({num2str(ystim-1),num2str(ystim),num2str(ystim+1)})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)
    
    if i == 1
        ylabel('Min-max forward velocity', 'FontSize', 16);
    elseif i == ceil(length(Data)/2)
        xlabel('Distance (y panel dimension)', 'FontSize', 16);
    end   
    line([12 12],[0 1],'Color' ,'r','lineWidth',2);
end
%suptitle('Opto trials');
saveas(gcf,[dirName,'\globalPlots\ProfilesOptoAroundPulse.png'])


figure('Position',[100 100 1400 800]),
for i = 1:length(Data)
    subplot(1,length(Data),i)
    plot([ystim-1:ystim+1],cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingProbeTrials,ystim-1:ystim+1))),'color',[0.5 0.5 0.5])
    hold on
    errorbar([ystim-1,ystim,ystim+1],nanmean(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingProbeTrials,ystim-1:ystim+1)))),nanstd(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingProbeTrials,ystim-1:ystim+1))))/sqrt(length(Data{1,i}.WorkingProbeTrials)),'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
    xlim([ystim-2 ystim+2]);
    xticks([ystim-1 ystim ystim+1])
    xticklabels({num2str(ystim-1),num2str(ystim),num2str(ystim+1)})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)
    if i == 1
        ylabel('Min-max forward velocity', 'FontSize', 16);
    elseif i == ceil(length(Data)/2)
        xlabel('Distance (y panel dimension)', 'FontSize', 16);
    end      
    line([12 12],[0 1],'Color' ,'r','lineWidth',2,'lineStyle',':');   
end
    
%suptitle('Probe trials');
saveas(gcf,[dirName,'\globalplots\ProfilesProbeAroundPulse.png'])

%% For angular speed

maxSpeed = cellfun(@max,angvelocities,'UniformOutput',0);
minSpeed = cellfun(@min, angvelocities,'UniformOutput',0);

for i = 1:size(velocities,1)
    for j = 1:size(velocities,2)
        normSpeed{i,j} = (angvelocities{i,j}-minSpeed{i,j})/(maxSpeed{i,j}-minSpeed{i,j});
    end
end

% bin per distance step
for i = 1:size(zscoreVel,1)
    for j = 1:size(Data{1,i}.smoothed,2)
        for k = 1:13
            stepMeans{i,j,k} = nanmean(normSpeed{i,j}(Data{1,i}.step{j,k}));
        end        
    end
end

% get the median and error for each animal
for i = 1:size(zscoreVel,1)
    optoMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeMeans(i,:) = squeeze(nanmedian(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
    optoError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingOptoPulse,:))));
    probeError(i,:) = squeeze(nanstd(cell2mat(stepMeans(i,Data{1,i}.WorkingProbeTrials,:))));
end


%Plot one figure per animal
ydimension = [1:13];

figure,
subplot(1,2,1)
for i = 1:size(optoMeans,1)
    h1 = plot(ydimension,optoMeans(i,:),'-o');
    hold on
    set(h1, 'markerfacecolor', get(h1, 'color'));
    line([12 12],[0 1],'Color' ,'r','LineWidth', 2);

end
title('Opto trials'); xlabel('y dimension');
ylabel('Max-min normalized angular speed (deg/s)');
ylim([0 1]);

subplot(1,2,2)
for i = 1:size(probeMeans,1)
    h1 = plot(ydimension,probeMeans(i,:),'-o');
    hold on
    set(h1, 'markerfacecolor', get(h1, 'color'));
    line([12 12],[0 1],'Color' ,'r', 'lineWidth', 2, 'lineStyle', ':');
end
title('Probe trials');
xlabel('y dimension');
ylim([0 1]);

saveas(gcf,[dirName,'\globalPlots\MaxMinAngSpeedvsDistance.png'])


%% Profile plots for normalized angular speed

figure('Position',[100 100 1400 800])

for i = 1:length(Data)
    subplot(1,length(Data),i)
    plot([ystim-1:ystim+1],cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingOptoPulse,ystim-1:ystim+1))),'color',[0.5 0.5 0.5])
    hold on
    errorbar([ystim-1,ystim,ystim+1],nanmean(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingOptoPulse,ystim-1:ystim+1)))),nanstd(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingOptoPulse,ystim-1:ystim+1))))/sqrt(length(Data{1,i}.WorkingOptoPulse)),'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
    xlim([ystim-2 ystim+2]);
    xticks([ystim-1 ystim ystim+1])
    xticklabels({num2str(ystim-1),num2str(ystim),num2str(ystim+1)})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)
    
    if i == 1
        ylabel('Min-max angular speed', 'FontSize', 16);
    elseif i == ceil(length(Data)/2)
        xlabel('Distance (y panel dimension)', 'FontSize', 16);
    end
    line([12 12],[0 1],'Color' ,'r','lineWidth',2);
end
%suptitle('Opto trials');
saveas(gcf,[dirName,'\globalPlots\ProfilesOptoAroundPulseAng.png'])


figure('Position',[100 100 1400 800]),  
for i = 1:length(Data)
    subplot(1,length(Data),i)
    plot([ystim-1:ystim+1],cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingProbeTrials,ystim-1:ystim+1))),'color',[0.5 0.5 0.5])
    hold on
    errorbar([ystim-1,ystim,ystim+1],nanmean(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingProbeTrials,ystim-1:ystim+1)))),nanstd(cell2mat(squeeze(stepMeans(i,Data{1,i}.WorkingProbeTrials,ystim-1:ystim+1))))/sqrt(length(Data{1,i}.WorkingProbeTrials)),'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
    xlim([ystim-2 ystim+2]);
    xticks([ystim-1 ystim ystim+1])
    xticklabels({num2str(ystim-1),num2str(ystim),num2str(ystim+1)})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)
    if i == 1
        ylabel('Min-max angular speed', 'FontSize', 16);
    elseif i == ceil(length(Data)/2)
        xlabel('Distance (y panel dimension)', 'FontSize', 16);
    end      
    line([12 12],[0 1],'Color' ,'r','lineWidth',2,'lineStyle',':');
end
    
%suptitle('Probe trials');
saveas(gcf,[dirName,'\globalPlots\ProfilesProbeAroundPulseAng.png'])

%% Temporal analysis

%(value-min)/(max-min)
maxVelsAP = cellfun(@max,velAroundExpansion,'UniformOutput',0);
minVelsAP = cellfun(@min, velAroundExpansion,'UniformOutput',0);

for i = 1:size(velAroundExpansion,2)
    normVelP{i} = (velAroundExpansion{1,i}-minVelsAP{1,i})./(maxVelsAP{1,i}-minVelsAP{1,i});
end

% get the median and error for each animal
for i = 1:size(zscoreVel,1)
    APMeans(i,:) = nanmedian(normVelP{1,i},2);
    APError(i,:) = nanstd(normVelP{1,i},[],2);   
end

time = linspace(-0.5, 1, 39);

figure,
subplot(1,2,2)
for i = 1:size(APMeans,1)
    h1 = plot(time, APMeans(i,:),'LineWidth',1.5);
    hold on
    line([0 0],[0 1],'Color' ,'r','lineWidth',2,'LineStyle',':');

end
title('Probe trials'); xlabel('Time (s)');
ylabel('Max-min forward velocity (mm/s)');
xlim([-0.5 1]);


maxVelsAPO = cellfun(@max,velAroundPulse,'UniformOutput',0);
minVelsAPO = cellfun(@min, velAroundPulse,'UniformOutput',0);

for i = 1:size(velAroundPulse,2)
    normVelO{i} = (velAroundPulse{1,i}-minVelsAPO{1,i})./(maxVelsAPO{1,i}-minVelsAPO{1,i});
end

% get the median and error for each animal
for i = 1:size(zscoreVel,1)
    APMeansO(i,:) = nanmedian(normVelO{1,i},2);
    APErrorO(i,:) = nanstd(normVelO{1,i},[],2);   
end

time = linspace(-0.5, 1, 39);

subplot(1,2,1)
for i = 1:size(APMeansO,1)
    h1 = plot(time, APMeansO(i,:),'LineWidth',1.5);
    hold on
    line([0 0],[0 1],'Color' ,'r','lineWidth',2);

end
title('Opto trials'); xlabel('Time (s)');
ylabel('Max-min forward velocity (mm/s)');

saveas(gcf,[dirName,'\globalPlots\MaxMinFwdVelinTimeOpto.png'])

%% Distance estimate analysis with local min and local max


for i = 1:length(folderNames)
    if contains(folderNames(i).name,'2019')
        distanceEstimate{i} = load(strcat(folderNames(i).folder,'\',folderNames(i).name,'\dataFromAnalysis\distanceEstimates.mat'));
        angDistanceEstimate{i} = load(strcat(folderNames(i).folder,'\',folderNames(i).name,'\dataFromAnalysis\AngDistanceEstimates.mat'));        
    else
        distanceEstimate{i} = [];
        angDistanceEstimate{i} = [];
    end
end

distanceEstimate = distanceEstimate(~cellfun('isempty',distanceEstimate));
angDistanceEstimate = angDistanceEstimate(~cellfun('isempty',angDistanceEstimate));

allDistanceEstimates = [];
allAngDistanceEstimates = [];

for i = 1:length(distanceEstimate)
    allDistanceEstimates = [allDistanceEstimates,distanceEstimate{1,i}.estimate];
    allAngDistanceEstimates = [allAngDistanceEstimates,angDistanceEstimate{1,i}.angEstimate];
end


figure,
boxplot(allDistanceEstimates,'orientation','horizontal')
hold on
scatter(allDistanceEstimates,ones(1,length(allDistanceEstimates)).*(1+(rand(1,length(allDistanceEstimates))-0.5)/10),'k','filled')
xlim([1 12]);
xlabel('y dimension');
xline(12,'--r','lineWidth',2)
set(gca,'ytick',[])
title('Estimated reward distance (local minimum)');

saveas(gcf,[dirName,'\globalPlots\LocalMinEstimate.png'])


figure,
boxplot(allAngDistanceEstimates,'orientation','horizontal')
hold on
scatter(allAngDistanceEstimates,ones(1,length(allAngDistanceEstimates)).*(1+(rand(1,length(allAngDistanceEstimates))-0.5)/10),'k','filled')
xlim([1 12]);
xlabel('y dimension');
xline(12,'--r','lineWidth',2)
set(gca,'ytick',[])
title('Estimated reward distance (local maximum)');

saveas(gcf,[dirName,'\globalPlots\LocalMaxEstimate.png'])
