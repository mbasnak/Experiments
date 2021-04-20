%code to analyze more carefully the relationship between total movement and
%bump parameters

%so far we have fit linear models, but what if this relationship has a
%different form?

%I will pull the data from all the flies for this experiment to explore
%this relationship better

clear all; close all;

%list all of the folders in the directory
folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');

%initialize empty vectors
total_mvt = [];
bump_mag = [];
yaw_speed = [];
allBumpMag = [];
side_vel = [];
fwd_vel = [];
half_max_width_ind = [];
half_width = [];

%for the fly folders
for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')==1
        %load the file with the sessions info
        load([folderNames(folder).folder,'\',folderNames(folder).name,'\sessions_info.mat']);
        %list every file
        fileNames{folder} = dir(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis'));
        for file = 1:length(fileNames{1,folder})
            %if this is an analysis file matching an offset stabilizing or
            %gain change session, load it and use the data
            if (contains(fileNames{1,folder}(file).name,'analysis')==1 & contains(fileNames{1,folder}(file).name,num2str(sessions_info.offset_stabilizer))==1)
                load(fullfile(fileNames{1,folder}(file).folder,fileNames{1,folder}(file).name))
                
                %load the mvt data into one big vector
                total_mvt_ds = fillmissing(data.total_mvt_ds,'linear');
                zscored_total_mvt = zscore(total_mvt_ds);
                total_mvt = [total_mvt,zscored_total_mvt];
                
                %load the yaw mvt data into one big vector
                yaw_speed_ds = abs(fillmissing(data.vel_yaw_ds,'linear'));
                zscored_yaw_speed_ds = zscore(yaw_speed_ds);
                yaw_speed = [yaw_speed,zscored_yaw_speed_ds];
                
                %load the fwd mvt data into one big vector
                fwd_vel_ds = fillmissing(data.vel_for_deg_ds,'linear');
                zscored_fwd_vel_ds = zscore(fwd_vel_ds);
                fwd_vel = [fwd_vel,zscored_fwd_vel_ds'];
                
                %load the side mvt data into one big vector
                side_vel_ds = abs(fillmissing(data.vel_side_deg_ds,'linear'));
                zscored_side_vel_ds = zscore(side_vel_ds);
                side_vel = [side_vel,zscored_side_vel_ds'];
                
        
                %load the bump mag data into one big vector
                bump_mag_ind = max(data.mean_dff_EB)-min(data.mean_dff_EB);
                zscored_bump_mag = zscore(bump_mag_ind);
                bump_mag = [bump_mag,zscored_bump_mag];
                
                %add analysis for bump width at half max
                for timepoint = 1:size(data.mean_dff_EB,2)
                    
                    %linearly interpolate to have 1000 datapoints instead of 8
                    interp_ex_data = interp1([1:8],data.mean_dff_EB(:,timepoint),[1:7/1000:8]);
                    %Find the half max point
                    half_max = (max(data.mean_dff_EB(:,timepoint))-min(data.mean_dff_EB(:,timepoint)))/2 + min(data.mean_dff_EB(:,timepoint));
                    [ex_bump_mag_interp I_interp] = max(interp_ex_data);
                    %Find in each half the index closest to the half max
                    diff_data = abs(interp_ex_data-half_max);
                    [sortedVals,indexes] = sort(diff_data);
                    %remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
                    %index(1)
                    diff_indexes = abs(indexes-indexes(1));
                    indexes(diff_indexes<175 & diff_indexes>0)=NaN;
                    indexes = indexes(~isnan(indexes));
                    two_indexes = [indexes(1), indexes(2)];
                    I1 = min(two_indexes);
                    I2 = max(two_indexes);
                    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
                        half_max_w = I1+1000-I2;
                    else
                        half_max_w = I2-I1;
                    end
                    %convert to EB tiles
                    half_max_width_ind(timepoint) = half_max_w*8/1001;
                    
                end
                zscored_half_width = zscore(half_max_width_ind);
                half_width = [half_width,zscored_half_width];
                
            end
        end
    end

end

%% Plot the distributions of the zscored mvt data

figure('Position',[100 100 1400 500]),
subplot(2,5,[1:4])
plot(fwd_vel)
ylabel('Zscored fwd vel')

subplot(2,5,5)
histogram(fwd_vel)

subplot(2,5,[6:9])
plot(yaw_speed)
ylabel('Zscored yaw speed')

subplot(2,5,10)
histogram(yaw_speed)

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\zscored_mvt.png');
% 
% %% Plot the relationship
% 
% figure('Position',[100 100 1600 800]),
% subplot(1,4,1)
% binscatter(fwd_vel,bump_mag)
% correlation1 = corrcoef(fwd_vel,bump_mag);
% hold on
% text(3,4,['Corr = ',num2str(round(correlation1(2,1),2))]);
% %xlim([0,5]);
% %ylim([0,5]);
% title('Fwd vel');
% xlabel('Zscored fwd vel');
% ylabel('Bump magnitude (max-min)');
% 
% 
% subplot(1,4,2)
% binscatter(side_vel,bump_mag)
% correlation2 = corrcoef(side_vel,bump_mag);
% hold on
% text(3,4,['Corr = ',num2str(round(correlation2(2,1),2))]);
% %xlim([0,5]);
% %ylim([0,5]);
% title('Side speed');
% xlabel('Zscored side speed');
% ylabel('Bump magnitude (max-min)');
% 
% subplot(1,4,3)
% binscatter(yaw_speed,bump_mag)
% correlation3 = corrcoef(yaw_speed,bump_mag);
% hold on
% text(3,4,['Corr = ',num2str(round(correlation3(2,1),2))]);
% %xlim([0,5]);
% %ylim([0,5]);
% title('Yaw speed');
% xlabel('Zscored yaw speed');
% ylabel('Bump magnitude (max-min)');
% 
% subplot(1,4,4)
% binscatter(total_mvt,bump_mag)
% correlation4 = corrcoef(total_mvt,bump_mag);
% hold on
% text(3,4,['Corr = ',num2str(round(correlation4(2,1),2))]);
% %xlim([0,5]);
% %ylim([0,5]);
% title('Total movement');
% xlabel('Zscored total movement');
% ylabel('Bump magnitude (max-min)');

%% Bin the data, and focus on z-scored movement below 2 

figure('Position',[200 200 1400 600]),

nbins = 20;
max_fwd_bin = 2;
min_fwd_bin = -2;
binWidth = max_fwd_bin/nbins;
fwdBins = [min_fwd_bin:binWidth:max_fwd_bin]; 

%getting binned means 
for bin = 1:length(fwdBins)-1
    meanBin(bin) = mean(bump_mag((fwd_vel > fwdBins(bin)) & (fwd_vel < fwdBins(bin+1))));
    stdBin(bin) = std(bump_mag((fwd_vel > fwdBins(bin)) & (fwd_vel < fwdBins(bin+1))));
    errBin(bin) = stdBin(bin)./sqrt(length(bump_mag((fwd_vel > fwdBins(bin)) & (fwd_vel < fwdBins(bin+1)))));
end

%create axes for plot
mvtAxes = fwdBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,1)
boundedline(mvtAxes,meanBin,errBin)
ylabel('Mean z-scored bump magnitude'); xlabel('Zscored fwd vel');
ylim([-1 1]);
title('Forward velocity');


max_side_bin = 2;
min_side_bin = -2;
binWidth = max_side_bin/nbins;
sideBins = [min_side_bin:binWidth:max_side_bin]; 

%getting binned means 
for bin = 1:length(sideBins)-1
    meanSideBin(bin) = mean(bump_mag((side_vel > sideBins(bin)) & (side_vel < sideBins(bin+1))));
    stdSideBin(bin) = std(bump_mag((side_vel > sideBins(bin)) & (side_vel < sideBins(bin+1))));
    errSideBin(bin) = stdSideBin(bin)./sqrt(length(bump_mag((side_vel > sideBins(bin)) & (side_vel < sideBins(bin+1)))));
end

%create axes for plot
mvtAxes = sideBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,2)
boundedline(mvtAxes,meanSideBin,errSideBin,'m')
ylabel('Mean z-scored bump magnitude'); xlabel('Zscored side speed');
ylim([-1 1]);
title('Side speed');


max_yaw_bin = 2;
min_yaw_bin = -2;
binWidth = max_yaw_bin/nbins;
yawBins = [min_yaw_bin:binWidth:max_yaw_bin]; 

%getting binned means 
for bin = 1:length(sideBins)-1
    meanYawBin(bin) = mean(bump_mag((yaw_speed > yawBins(bin)) & (yaw_speed < yawBins(bin+1))));
    stdYawBin(bin) = std(bump_mag((yaw_speed > yawBins(bin)) & (yaw_speed < yawBins(bin+1))));
    errYawBin(bin) = stdYawBin(bin)./sqrt(length(bump_mag((yaw_speed > yawBins(bin)) & (yaw_speed < yawBins(bin+1)))));
end

%create axes for plot
mvtAxes = sideBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,3)
boundedline(mvtAxes,meanYawBin,errYawBin,'r')
ylabel('Mean z-scored bump magnitude'); xlabel('Zscored yaw speed');
ylim([-1 1]);
title('Yaw speed');



min_mvt_bin = -2;
max_mvt_bin = 2;
binWidth = max_mvt_bin/nbins;
mvtBins = [min_mvt_bin:binWidth:max_mvt_bin]; 

%getting binned means 
for bin = 1:length(mvtBins)-1
    meanMvtBin(bin) = mean(bump_mag((total_mvt > mvtBins(bin)) & (total_mvt < mvtBins(bin+1))));
    stdMvtBin(bin) = std(bump_mag((total_mvt > mvtBins(bin)) & (total_mvt < mvtBins(bin+1))));
    errMvtBin(bin) = stdMvtBin(bin)./sqrt(length(bump_mag((total_mvt > mvtBins(bin)) & (total_mvt < mvtBins(bin+1)))));
end

%create axes for plot
mvtAxes = mvtBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,4)
boundedline(mvtAxes,meanMvtBin,errMvtBin,'k')
ylabel('Mean z-scored bump magnitude'); xlabel('Zscored total movement');
ylim([-1 1]);
title('Total movement');

%save the data
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\Mvt_vs_bm.png');


%% repeat for bump width

figure('Position',[200 200 1400 600]),

nbins = 20;
max_fwd_bin = 2;
min_fwd_bin = -2;
binWidth = max_fwd_bin/nbins;
fwdBins = [min_fwd_bin:binWidth:max_fwd_bin]; 

%getting binned means 
for bin = 1:length(fwdBins)-1
    meanBin(bin) = mean(half_width((fwd_vel > fwdBins(bin)) & (fwd_vel < fwdBins(bin+1))));
    stdBin(bin) = std(half_width((fwd_vel > fwdBins(bin)) & (fwd_vel < fwdBins(bin+1))));
    errBin(bin) = stdBin(bin)./sqrt(length(half_width((fwd_vel > fwdBins(bin)) & (fwd_vel < fwdBins(bin+1)))));
end

%create axes for plot
mvtAxes = fwdBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,1)
boundedline(mvtAxes,meanBin,errBin)
ylabel('Mean z-scored half width'); xlabel('Zscored fwd vel');
ylim([-1 1]);
title('Forward velocity');


max_side_bin = 2;
min_side_bin = -2;
binWidth = max_side_bin/nbins;
sideBins = [min_side_bin:binWidth:max_side_bin]; 

%getting binned means 
for bin = 1:length(sideBins)-1
    meanSideBin(bin) = mean(half_width((side_vel > sideBins(bin)) & (side_vel < sideBins(bin+1))));
    stdSideBin(bin) = std(half_width((side_vel > sideBins(bin)) & (side_vel < sideBins(bin+1))));
    errSideBin(bin) = stdSideBin(bin)./sqrt(length(half_width((side_vel > sideBins(bin)) & (side_vel < sideBins(bin+1)))));
end

%create axes for plot
mvtAxes = sideBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,2)
boundedline(mvtAxes,meanSideBin,errSideBin,'m')
ylabel('Mean z-scored half width'); xlabel('Zscored side speed');
ylim([-1 1]);
title('Side speed');


max_yaw_bin = 2;
min_yaw_bin = -2;
binWidth = max_yaw_bin/nbins;
yawBins = [min_yaw_bin:binWidth:max_yaw_bin]; 

%getting binned means 
for bin = 1:length(sideBins)-1
    meanYawBin(bin) = mean(half_width((yaw_speed > yawBins(bin)) & (yaw_speed < yawBins(bin+1))));
    stdYawBin(bin) = std(half_width((yaw_speed > yawBins(bin)) & (yaw_speed < yawBins(bin+1))));
    errYawBin(bin) = stdYawBin(bin)./sqrt(length(half_width((yaw_speed > yawBins(bin)) & (yaw_speed < yawBins(bin+1)))));
end

%create axes for plot
mvtAxes = sideBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,3)
boundedline(mvtAxes,meanYawBin,errYawBin,'r')
ylabel('Mean z-scored half width'); xlabel('Zscored yaw speed');
ylim([-1 1]);
title('Yaw speed');



min_mvt_bin = -2;
max_mvt_bin = 2;
binWidth = max_mvt_bin/nbins;
mvtBins = [min_mvt_bin:binWidth:max_mvt_bin]; 

%getting binned means 
for bin = 1:length(mvtBins)-1
    meanMvtBin(bin) = mean(half_width((total_mvt > mvtBins(bin)) & (total_mvt < mvtBins(bin+1))));
    stdMvtBin(bin) = std(half_width((total_mvt > mvtBins(bin)) & (total_mvt < mvtBins(bin+1))));
    errMvtBin(bin) = stdMvtBin(bin)./sqrt(length(half_width((total_mvt > mvtBins(bin)) & (total_mvt < mvtBins(bin+1)))));
end

%create axes for plot
mvtAxes = mvtBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,4)
boundedline(mvtAxes,meanMvtBin,errMvtBin,'k')
ylabel('Mean z-scored bump width'); xlabel('Zscored total movement');
ylim([-1 1]);
title('Total movement');

%save the data
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\Mvt_vs_hw.png');


%% heatmap of fwd vel and ang speed

%I want to see if there is a systematic relationship between moments of
%high fwd vel and high ang speed

max_yaw_bin2 = 2;
min_yaw_bin2 = -1;
binWidth2 = max_yaw_bin2/nbins;
yawBins2 = [min_yaw_bin2:binWidth2:max_yaw_bin2]; 

%getting binned means 
for bin = 1:length(yawBins2)-1
    meanYawBin2(bin) = mean(fwd_vel((yaw_speed > yawBins2(bin)) & (yaw_speed < yawBins2(bin+1))));
    stdYawBin2(bin) = std(fwd_vel((yaw_speed > yawBins2(bin)) & (yaw_speed < yawBins2(bin+1))));
    errYawBin2(bin) = stdYawBin2(bin)./sqrt(length(fwd_vel((yaw_speed > yawBins2(bin)) & (yaw_speed < yawBins2(bin+1)))));
end

%create axes for plot
mvtAxes = yawBins2 - binWidth2;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth2;

%Plot
figure,
boundedline(mvtAxes,meanYawBin2,errYawBin2)
ylabel('Mean z-scored fwd vel'); xlabel('Zscored yaw speed');
ylim([-2 2]);

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\fwd_vs_yaw.png');
