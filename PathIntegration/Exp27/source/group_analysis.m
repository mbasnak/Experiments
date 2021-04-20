%group analysis

close all; clear all;
%get the directory
pathName = uigetdir;

%import all the data
cd (pathName)
folderNames = dir;

for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'2021')==1
        data{folder} = load([folderNames(folder).folder,'\',folderNames(folder).name,'\ball\distanceData.mat']);
    end
end

%remove empty cells
data = data(~cellfun('isempty',data));


%% Make directory to save plots if it doesn't already exist
  
%Move to the analysis folder
cd(pathName)
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'globalPlots') == 0)
   mkdir(pathName,'globalPlots'); 
end
%List the contents of the 'plots' folder
cd([pathName,'\globalPlots\'])

%% get distance

%Convert the data to distance in relevant units
ydimension = [1:48]; %number of stimulus dimensions we're using  (check if this is 48 or 46!)
%convert the ydimensions to distance using the xgain
distance = 9*pi*(ydimension/92)/0.5;
%get reward location
rewardDimension = round(50*48/100);

%% combine data

allOptoMeans = [];
allProbeMeans = [];
for fly = 1:length(data)
    allOptoMeans = [allOptoMeans;data{1,fly}.optoMeans];
    allProbeMeans = [allProbeMeans;data{1,fly}.probeMeans];
end

figure,
subplot(1,2,1)
plot(distance,allOptoMeans','r');
hold on
plot(distance,mean(allOptoMeans),'r','linewidth',3)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Opto trials');
ylim([0 10]);
ylabel('Forward velocity (mm/s)');
xlabel('Distance (mm)');

subplot(1,2,2)
plot(distance,allProbeMeans','k')
hold on
plot(distance,mean(allProbeMeans),'k','linewidth',3)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');
title('Probe trials');
ylim([0 10]);

saveas(gcf,'meanVelPooledData.png');

%% zscoring the data

zscored_opto_means = zscore(allOptoMeans,[],2);
zscored_probe_means = zscore(allProbeMeans,[],2);

figure,
subplot(1,2,1)
plot(distance,zscored_opto_means','color',[0.5,0.5,0.5])
hold on
plot(distance,mean(zscored_opto_means),'r','linewidth',3)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('zscored forward velocity');
xlabel('Distance (mm)');
title('Opto trials');
ylim([-4,2]);

subplot(1,2,2)
plot(distance,zscored_probe_means','color',[0.5,0.5,0.5])
hold on
plot(distance,mean(zscored_probe_means),'k','linewidth',3)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');
title('Probe trials');
ylim([-4,2]);

saveas(gcf,'meanVelZscoredPooledData.png')

%% heatmap of the results

figure,
subplot(1,2,1)
imagesc(zscored_opto_means)
title('Opto trials');
colormap(gray)

subplot(1,2,2)
imagesc(zscored_probe_means)
title('Probe trials');
colormap(gray)

saveas(gcf,'heatmapVelZscoredPooledData.png')
