%Code to analyze the relationship between the fly's velocity and the bump
%magnitude across flies

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
%Get folder names
folderContents = dir(path);

%get what type of block it is
prompt = 'Is this for offset stabilizer bouts? Write 1 if yes ';
answer = input(prompt);

%Load the summary data of the folder that correspond to experimental flies
for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       if answer == 1
        data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\vel_bm_data.mat']);
       else
        data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\vel_bm_data_ov.mat']);           
       end
   end
end

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%% Plot relationship between forward velocity and bump magnitude

figure('Position',[200 200 800 600]),

subplot(1,2,1)
for fly = 1:size(data,2)
    nbins = 20;
    yaw_speed = abs(data(fly).yaw_speed);
    maxBinYS = min([max(yaw_speed),nanmean(yaw_speed)+3*nanstd(yaw_speed)]);
    binWidthYS = maxBinYS/nbins;
    YSBins = [0:binWidthYS:maxBinYS];
    
    %getting binned means
    for bin = 1:length(YSBins)-1
        meanBin(bin) = nanmean(data(fly).bump_magnitude((yaw_speed(1:length(data(fly).bump_magnitude)) > YSBins(bin)) & (yaw_speed(1:length(data(fly).bump_magnitude)) < YSBins(bin+1))));
    end
    
    %create axes for plot
    YSAxes = YSBins - binWidthYS;
    YSAxes = YSAxes(2:end);
    YSAxes(end) = YSAxes(end-1)+binWidthYS;
    
    %Plot
    plot(YSAxes,meanBin,'-o','color',[.5 .5 .5])
    ylabel('Mean bump magnitude'); xlabel('Yaw speed (deg/s)');
    ylim([0 1.2]);
    
    hold on
end

subplot(1,2,2)
for fly = 1:size(data,2)
    for_vel = abs(data(fly).for_vel);
    maxBinFV = min([max(for_vel),nanmean(for_vel)+3*nanstd(for_vel)]);
    binWidthFV = maxBinFV/nbins;
    FVBins = [0:binWidthFV:maxBinFV];
    
    %getting binned means
    for bin = 1:length(FVBins)-1
        meanBin(bin) = nanmean(data(fly).bump_magnitude((for_vel(1:length(data(fly).bump_magnitude)) > FVBins(bin)) & (for_vel(1:length(data(fly).bump_magnitude)) < FVBins(bin+1))));
    end
    
    %create axes for plot
    FVAxes = FVBins - binWidthFV;
    FVAxes = FVAxes(2:end);
    FVAxes(end) = FVAxes(end-1)+binWidthFV;
    
    %Plot
    plot(FVAxes,meanBin,'-o','color',[.5 .5 .5])
    ylabel('Mean bump magnitude'); xlabel('Forward velocity (mm/s)');
    ylim([0 1.2]);
    
    hold on
end

if answer == 1
    saveas(gcf,[path,'\globalPlots\vel_vs_BM.png']);
else
    saveas(gcf,[path,'\globalPlots\vel_vs_BM_ov.png']);
end

%% Combine all the data without zscoring

allYawSpeed = [];
allForVel = [];
allBM = [];

for fly = 1:length(data)
    allYawSpeed = [allYawSpeed,data(fly).yaw_speed];
    allForVel = [allForVel,data(fly).for_vel];
    allBM = [allBM,data(fly).bump_magnitude];
end


nbins = 20;
maxBinYS = min([max(allYawSpeed), nanmean(allYawSpeed)+3*nanstd(allYawSpeed)]);
maxBinFV = min([max(allForVel), nanmean(allForVel)+3*nanstd(allForVel)]);
binWidthYS = maxBinYS/nbins;
YSBins = [0:binWidthYS:maxBinYS];
binWidthFV = maxBinFV/nbins;
FVBins = [0:binWidthFV:maxBinFV];

%getting binned means
for ys_bin = 1:length(YSBins)-1
    for fv_bin = 1:length(FVBins)-1
        doubleBin(ys_bin,fv_bin) = mean(allBM((allYawSpeed(1:length(allBM)) > YSBins(ys_bin)) & (allYawSpeed(1:length(allBM)) < YSBins(ys_bin+1)) & (allForVel(1:length(allBM)) > FVBins(fv_bin)) & (allForVel(1:length(allBM)) < FVBins(fv_bin+1))));
    end
end

%flip the data such that high forward velocity values are at the top
binned_data = flip(doubleBin);

figure,
imagesc(binned_data)
xticks([1:4:20])
set(gca, 'XTickLabel', round(FVBins(1:4:20)))
xlabel('Forward velocity (mm/s)','fontsize',12,'fontweight','bold');
ylabel('Yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
yticks([1:4:20])
set(gca, 'YTickLabel', round(YSBins(20:-4:1)))
c = colorbar;
ylabel(c, 'Bump magnitude')

if answer == 1
    saveas(gcf,[path,'\globalPlots\vel_BM_heatmap.png']);
else
    saveas(gcf,[path,'\globalPlots\vel_BM_heatmap_ov.png']);
end

%% Zscoring the data

all_data = struct2cell(data);
%zscored_data = cellfun(@zscore,all_data,'UniformOutput',0);
zscored_data = cellfun(@zscore_omitting_nans,all_data,'UniformOutput',0);
fields = {'yaw_speed','for_vel','bump_magnitude'};
zscored_data = cell2struct(zscored_data,fields,1);

% %Plot the z-scored data distribution
% figure('Position',[100 100 800 1200]),
% for fly = 1:14
%     subplot(15,3,fly*3-2)
%     histogram(zscored_data(fly).yaw_speed);
%     
%     subplot(15,3,fly*3-1)
%     histogram(zscored_data(fly).for_vel);
%     
%     subplot(15,3,fly*3)
%     histogram(zscored_data(fly).bump_magnitude);
% end

allYawSpeed = [];
allForVel = [];
allBM = [];

for fly = 1:length(data)
    allYawSpeed = [allYawSpeed,zscored_data(fly).yaw_speed];
    allForVel = [allForVel,zscored_data(fly).for_vel];
    allBM = [allBM,zscored_data(fly).bump_magnitude];
end


nbins = 20;
maxBinYS = min([max(allYawSpeed), nanmean(allYawSpeed)+3*nanstd(allYawSpeed)]);
maxBinFV = min([max(allForVel), nanmean(allForVel)+3*nanstd(allForVel)]);
minBinYS = min(allYawSpeed);
minBinFV = min(allForVel);
binWidthYS = (maxBinYS-minBinYS)/nbins;
YSBins = [minBinYS:binWidthYS:maxBinYS];
binWidthFV = (maxBinFV-minBinFV)/nbins;
FVBins = [minBinFV:binWidthFV:maxBinFV];

%getting binned means
for ys_bin = 1:length(YSBins)-1
    for fv_bin = 1:length(FVBins)-1
        doubleBin(ys_bin,fv_bin) = mean(allBM((allYawSpeed(1:length(allBM)) > YSBins(ys_bin)) & (allYawSpeed(1:length(allBM)) < YSBins(ys_bin+1)) & (allForVel(1:length(allBM)) > FVBins(fv_bin)) & (allForVel(1:length(allBM)) < FVBins(fv_bin+1))));
    end
end

%flip the data such that high forward velocity values are at the top
binned_data = flip(doubleBin);

figure,
imagesc(binned_data)
xticks([1:4:20])
set(gca, 'XTickLabel', round(FVBins(1:4:20)))
xlabel('Zscored forward velocity (mm/s)','fontsize',12,'fontweight','bold');
ylabel('Zscored yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
yticks([1:4:20])
set(gca, 'YTickLabel', round(YSBins(20:-4:1)))
c = colorbar;
ylabel(c, 'Zscored bump magnitude')

if answer == 1
    saveas(gcf,[path,'\globalPlots\zscored_vel_BM_heatmap.png']);
else
    saveas(gcf,[path,'\globalPlots\zscored_vel_BM_heatmap_ov.png']);
end

%%
close all; clear all
