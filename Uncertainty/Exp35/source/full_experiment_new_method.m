%Full experiment analysis
%Code to analyze the full experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Set colormap

folderNames = dir(path(1:53));
flyNames = struct();
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'60D05') & ~contains(folderNames(folder).name,'txt'))
        flyNames(folder).name = folderNames(folder).name;
    end
end
%Remove empty rows
flyNames = flyNames(~cellfun(@isempty,{flyNames.name}));

%Assign fly number
for fly = 1:length(flyNames)
    if strcmp(flyNames(fly).name,path(54:end))
        fly_ID = fly;
    end
end

colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
    0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880];

fly_color = colors_for_plots(fly_ID,:);

%% Analyze initial closed-loop panels

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.initial_cl_bar),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
pre_panels_heading = deg2rad(heading);
pre_panels_heading_thresh = pre_panels_heading(continuous_data.adj_rs>=0.5);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
pre_panels_offset = deg2rad(offset);
pre_panels_offset_above_thresh = pre_panels_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude')

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width')

suptitle('Initial trial with just panels');

saveas(gcf,[path,'\analysis\plots\initial_panels.png']);

%get mean bump parameters
meanBM_pre_panels = mean(continuous_data.bump_magnitude);
meanBW_pre_panels = mean(continuous_data.bump_width);
meanBM_thresh_pre_panels = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5));
meanBW_thresh_pre_panels = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5));

%get mean vel
mean_total_mvt_pre_panels = nanmean(continuous_data.total_mvt_ds);
mean_total_mvt_thresh_pre_panels = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5));

%store all bump param and mvt values
allBumpMag = continuous_data.bump_magnitude;
allBumpWidth = continuous_data.bump_width;
allTotalMvt = continuous_data.total_mvt_ds;
blockType = repelem(1,1,length(continuous_data.bump_magnitude));

%% Analyze inital open-loop panels

if isfield(sessions,'initial_ol_bar')
    
    for sid = 1:length(sessions.initial_ol_bar)
        
        load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.initial_ol_bar(sid)),'_tid_0.mat'])
        
        figure('Position',[100 100 1000 800]),
        subplot(4,5,[1:4])
        dff = continuous_data.dff_matrix';
        imagesc(flip(dff))
        colormap(flipud(gray))
        title('EPG activity');
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[6:9])
        bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
        %Remove wrapped lines to plot
        [x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
        plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
        hold on
        stim = wrapTo180(continuous_data.panel_angle);
        [x_out_stim,stim_to_plot] = removeWrappedLines(continuous_data.time,stim);
        heading = wrapTo180(-continuous_data.heading_deg);
        [x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
        plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
        plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
        legend('Bump estimate','Stim position','Fly position','Location','best');
        title('Bump and stim position');
        ylim([-180 180]);
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[11:14])
        offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle))));
        [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
        plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
        title('Stimulus offset');
        set(gca,'xticklabel',{[]})
        ylim([-180 180]);
        
        subplot(4,5,15)
        polarhistogram(deg2rad(offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        subplot(4,5,[16:19])
        mvt_offset = wrapTo180(continuous_data.offset);
        [x_out_mvt_offset,mvt_offset_to_plot] = removeWrappedLines(continuous_data.time,mvt_offset);
        plot(x_out_mvt_offset,mvt_offset_to_plot,'LineWidth',1.5,'color','k')
        title('Movement offset');
        xlabel('Time (sec)');
        ylim([-180 180]);
        
        subplot(4,5,20)
        polarhistogram(deg2rad(mvt_offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        suptitle('Initial OL panels');
        
        saveas(gcf,[path,'\analysis\plots\initial_OL_panels_trial',num2str(sid),'.png']);
        
    end
    
end

%%  Analyze initial closed-loop wind

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.initial_cl_wind),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
pre_wind_heading = deg2rad(heading);
pre_wind_heading_thresh = pre_wind_heading(continuous_data.adj_rs>=0.5);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
pre_wind_offset = deg2rad(offset);
pre_wind_offset_above_thresh = pre_wind_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude')

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width')

suptitle('Initial trial with just wind');

saveas(gcf,[path,'\analysis\plots\initial_wind.png']);


%get mean bump parameters
meanBM_pre_wind = mean(continuous_data.bump_magnitude);
meanBW_pre_wind = mean(continuous_data.bump_width);
meanBM_thresh_pre_wind = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5));
meanBW_thresh_pre_wind = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5));
%get mean vel
mean_total_mvt_pre_wind = nanmean(continuous_data.total_mvt_ds);
mean_total_mvt_thresh_pre_wind = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5));

%store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(2,1,length(continuous_data.bump_magnitude))];

%% Analyze inital open-loop wind

if isfield(sessions,'initial_old_wind')
    for sid = 1:length(sessions.initial_ol_wind)
        
        load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.initial_ol_wind(sid)),'_tid_0.mat'])
        
        figure('Position',[100 100 1000 800]),
        subplot(4,5,[1:4])
        dff = continuous_data.dff_matrix';
        imagesc(flip(dff))
        colormap(flipud(gray))
        title('EPG activity');
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[6:9])
        bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
        %Remove wrapped lines to plot
        [x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
        plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
        hold on
        stim = wrapTo180(rad2deg(continuous_data.motor_pos'));
        [x_out_stim,stim_to_plot] = removeWrappedLines(continuous_data.time,stim);
        heading = wrapTo180(-continuous_data.heading_deg);
        [x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
        plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
        plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
        legend('Bump estimate','Stim position','Fly position','Location','best');
        title('Bump and stim position');
        ylim([-180 180]);
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[11:14])
        offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',continuous_data.motor_pos')));
        [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
        plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
        title('Stimulus offset');
        set(gca,'xticklabel',{[]})
        ylim([-180 180]);
        
        subplot(4,5,15)
        polarhistogram(deg2rad(offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        subplot(4,5,[16:19])
        mvt_offset = wrapTo180(continuous_data.offset);
        [x_out_mvt_offset,mvt_offset_to_plot] = removeWrappedLines(continuous_data.time,mvt_offset);
        plot(x_out_mvt_offset,mvt_offset_to_plot,'LineWidth',1.5,'color','k')
        title('Movement offset');
        xlabel('Time (sec)');
        ylim([-180 180]);
        
        subplot(4,5,20)
        polarhistogram(deg2rad(mvt_offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        suptitle('Initial OL wind');
        
        saveas(gcf,[path,'\analysis\plots\initial_OL_wind_trial',num2str(sid),'.png']);
        
    end
end
%% Analyze the cue combination trial

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.cue_combination),'_tid_0.mat'])

%Cut the data for the fly for which fictrac failed
if contains(path, '20210922')
    frames = find(floor(continuous_data.time) == 600);
    first_frame = frames(1);
    continuous_data.dff_matrix = continuous_data.dff_matrix(1:first_frame,:);
    continuous_data.bump_pos = continuous_data.bump_pos(1:first_frame);
    continuous_data.heading_deg = continuous_data.heading_deg(1:first_frame);
    continuous_data.time = continuous_data.time(1:first_frame);
    continuous_data.heading = continuous_data.heading(1:first_frame);
    continuous_data.adj_rs = continuous_data.adj_rs(1:first_frame);
    continuous_data.bump_magnitude = continuous_data.bump_magnitude(1:first_frame);
    continuous_data.bump_width = continuous_data.bump_width(1:first_frame);
    continuous_data.total_mvt_ds = continuous_data.total_mvt_ds(1:first_frame);
    continuous_data.vel_yaw_ds = continuous_data.vel_yaw_ds(1:first_frame);
    continuous_data.vel_for_ds = continuous_data.vel_for_ds(1:first_frame);
end

figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]});

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
combined_heading = deg2rad(heading);
combined_heading_thresh = combined_heading(continuous_data.adj_rs>=0.5);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','Best');
set(gca,'xticklabel',{[]});
ylim([-180 180]);
xlim([0 x_out_heading(end)]);

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
combined_offset = deg2rad(offset);
combined_offset_above_thresh = combined_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset');
ylim([-180 180]);
xlim([0 x_out_offset(end)]);

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude');
xlim([0 continuous_data.time(end)]);

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width');
xlim([0 continuous_data.time(end)]);


suptitle('Trial with both cues');

saveas(gcf,[path,'\analysis\plots\cue_combination.png']);

%get mean bump parameters
meanBM_combined = mean(continuous_data.bump_magnitude);
meanBW_combined = mean(continuous_data.bump_width);
meanBM_thresh_combined = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5));
meanBW_thresh_combined = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5));
%get mean vel
mean_total_mvt_combined = nanmean(continuous_data.total_mvt_ds);
mean_total_mvt_thresh_combined = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5));


%store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(3,1,length(continuous_data.bump_magnitude))];

%% Relationship between bump parameters and velocity for the cue combination trial

nbins = 20;
yaw_speed = abs(continuous_data.vel_yaw_ds);
maxBinYS = min([max(yaw_speed),nanmean(yaw_speed)+3*nanstd(yaw_speed)]);
binWidthYS = maxBinYS/nbins;
YSBins = [0:binWidthYS:maxBinYS];
for_vel = abs(continuous_data.vel_for_ds');
maxBinFV = min([max(for_vel),nanmean(for_vel)+3*nanstd(for_vel)]);
binWidthFV = maxBinFV/nbins;
FVBins = [0:binWidthFV:maxBinFV];

%getting binned means
for ys_bin = 1:length(YSBins)-1
    for fv_bin = 1:length(FVBins)-1
        doubleBinBM(ys_bin,fv_bin) = nanmean(continuous_data.bump_magnitude((yaw_speed(1:length(continuous_data.bump_magnitude)) > YSBins(ys_bin)) & (yaw_speed(1:length(continuous_data.bump_magnitude)) < YSBins(ys_bin+1)) & (for_vel(1:length(continuous_data.bump_magnitude)) > FVBins(fv_bin)) & (for_vel(1:length(continuous_data.bump_magnitude)) < FVBins(fv_bin+1))));
        doubleBinBW(ys_bin,fv_bin) = nanmean(continuous_data.bump_width((yaw_speed(1:length(continuous_data.bump_width)) > YSBins(ys_bin)) & (yaw_speed(1:length(continuous_data.bump_width)) < YSBins(ys_bin+1)) & (for_vel(1:length(continuous_data.bump_width)) > FVBins(fv_bin)) & (for_vel(1:length(continuous_data.bump_width)) < FVBins(fv_bin+1))));        
    end
end

%flip the data such that high forward velocity values are at the top
binned_data_BM = flip(doubleBinBM);
binned_data_BW = flip(doubleBinBW);

figure('Position',[100 100 1400 500]),
subplot(1,2,1)
imagesc(binned_data_BM)
xticks([1:4:20])
set(gca, 'XTickLabel', round(FVBins(1:4:20)))
xlabel('Forward velocity (mm/s)','fontsize',12,'fontweight','bold');
ylabel('Yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
yticks([1:4:20])
set(gca, 'YTickLabel', round(YSBins(20:-4:1)))
c = colorbar;
ylabel(c, 'Bump magnitude')

subplot(1,2,2)
imagesc(binned_data_BW)
xticks([1:4:20])
set(gca, 'XTickLabel', round(FVBins(1:4:20)))
xlabel('Forward velocity (mm/s)','fontsize',12,'fontweight','bold');
ylabel('Yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
yticks([1:4:20])
set(gca, 'YTickLabel', round(YSBins(20:-4:1)))
c = colorbar;
ylabel(c, 'Bump with')

saveas(gcf,[path,'\analysis\plots\vel_vs_bm_heatmap.png']);

%% Analyze final closed-loop panels

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.final_cl_bar),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
post_panels_heading = deg2rad(heading);
post_panels_heading_thresh = post_panels_heading(continuous_data.adj_rs>=0.5);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
post_panels_offset = deg2rad(offset);
post_panels_offset_above_thresh = post_panels_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude')

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width')

suptitle('Final trial with just panels');

saveas(gcf,[path,'\analysis\plots\final_panels.png']);

%get mean bump parameters
meanBM_post_panels = mean(continuous_data.bump_magnitude);
meanBW_post_panels = mean(continuous_data.bump_width);
meanBM_thresh_post_panels = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5));
meanBW_thresh_post_panels = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5));
%get fly vel
mean_total_mvt_post_panels = nanmean(continuous_data.total_mvt_ds);
mean_total_mvt_thresh_post_panels = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5));

%store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(4,1,length(continuous_data.bump_magnitude))];

%% Analyze final open-loop panels

if isfield(sessions,'final_ol_bar')
    for sid = 1:length(sessions.final_ol_bar)
        
        load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.final_ol_bar(sid)),'_tid_0.mat'])
        
        figure('Position',[100 100 1000 800]),
        subplot(4,5,[1:4])
        dff = continuous_data.dff_matrix';
        imagesc(flip(dff))
        colormap(flipud(gray))
        title('EPG activity');
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[6:9])
        bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
        %Remove wrapped lines to plot
        [x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
        plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
        hold on
        stim = wrapTo180(continuous_data.panel_angle);
        [x_out_stim,stim_to_plot] = removeWrappedLines(continuous_data.time,stim);
        heading = wrapTo180(-continuous_data.heading_deg);
        [x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
        plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
        plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
        legend('Bump estimate','Stim position','Fly position','Location','best');
        title('Bump and stim position');
        ylim([-180 180]);
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[11:14])
        offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle))));
        [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
        plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
        title('Stimulus offset');
        set(gca,'xticklabel',{[]})
        ylim([-180 180]);
        
        subplot(4,5,15)
        polarhistogram(deg2rad(offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        subplot(4,5,[16:19])
        mvt_offset = wrapTo180(continuous_data.offset);
        [x_out_mvt_offset,mvt_offset_to_plot] = removeWrappedLines(continuous_data.time,mvt_offset);
        plot(x_out_mvt_offset,mvt_offset_to_plot,'LineWidth',1.5,'color','k')
        title('Movement offset');
        xlabel('Time (sec)');
        ylim([-180 180]);
        
        subplot(4,5,20)
        polarhistogram(deg2rad(mvt_offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        suptitle('Final OL panels');
        
        saveas(gcf,[path,'\analysis\plots\final_OL_panels_trial',num2str(sid),'.png']);
        
    end
end

%%  Analyze final closed-loop wind

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.final_cl_wind),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
post_wind_heading = deg2rad(heading);
post_wind_heading_thresh = post_wind_heading(continuous_data.adj_rs>=0.5);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
post_wind_offset = deg2rad(offset);
post_wind_offset_above_thresh = post_wind_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude')

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width')

suptitle('Final trial with just wind');

saveas(gcf,[path,'\analysis\plots\final_wind.png']);

%get mean bump parameters
meanBM_post_wind = mean(continuous_data.bump_magnitude);
meanBW_post_wind = mean(continuous_data.bump_width);
meanBM_thresh_post_wind = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5));
meanBW_thresh_post_wind = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5));
%get fly vel
mean_total_mvt_post_wind = nanmean(continuous_data.total_mvt_ds);
mean_total_mvt_thresh_post_wind = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5));

%store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(5,1,length(continuous_data.bump_magnitude))];

%% Analyze final open-loop wind

if isfield(sessions,'final_ol_wind')
    for sid = 1:length(sessions.final_ol_wind)
        
        load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.final_ol_wind(sid)),'_tid_0.mat'])
        
        figure('Position',[100 100 1000 800]),
        subplot(4,5,[1:4])
        dff = continuous_data.dff_matrix';
        imagesc(flip(dff))
        colormap(flipud(gray))
        title('EPG activity');
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[6:9])
        bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
        %Remove wrapped lines to plot
        [x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
        plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
        hold on
        stim = wrapTo180(rad2deg(continuous_data.motor_pos'));
        [x_out_stim,stim_to_plot] = removeWrappedLines(continuous_data.time,stim);
        heading = wrapTo180(-continuous_data.heading_deg);
        [x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
        plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
        plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
        legend('Bump estimate','Stim position','Fly position','Location','best');
        title('Bump and stim position');
        ylim([-180 180]);
        set(gca,'xticklabel',{[]})
        
        subplot(4,5,[11:14])
        offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',continuous_data.motor_pos')));
        [x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
        plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
        title('Stimulus offset');
        set(gca,'xticklabel',{[]})
        ylim([-180 180]);
        
        subplot(4,5,15)
        polarhistogram(deg2rad(offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        subplot(4,5,[16:19])
        mvt_offset = wrapTo180(continuous_data.offset);
        [x_out_mvt_offset,mvt_offset_to_plot] = removeWrappedLines(continuous_data.time,mvt_offset);
        plot(x_out_mvt_offset,mvt_offset_to_plot,'LineWidth',1.5,'color','k')
        title('Movement offset');
        xlabel('Time (sec)');
        ylim([-180 180]);
        
        subplot(4,5,20)
        polarhistogram(deg2rad(mvt_offset),'FaceColor','k')
        set(gca,'ThetaZeroLocation','top');
        Ax = gca; % current axes
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
        
        suptitle('Final OL wind');
        
        saveas(gcf,[path,'\analysis\plots\final_OL_wind_trial',num2str(sid),'.png']);
        
    end
end

%%
close all
%% Offset evolution


if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_offset,'FaceColor',fly_color)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_offset,'FaceColor',fly_color)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset,'FaceColor',fly_color)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_offset,'FaceColor',fly_color)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_offset,'FaceColor',fly_color)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_offset,'FaceColor',fly_color)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_offset,'FaceColor',fly_color)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset,'FaceColor',fly_color)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_offset,'FaceColor',fly_color)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_offset,'FaceColor',fly_color)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution');
end

saveas(gcf,[path,'\analysis\plots\offset_evolution.png']);


%% Thresholded offset evolution

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_offset_above_thresh,'FaceColor',fly_color)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_offset_above_thresh,'FaceColor',fly_color)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset_above_thresh,'FaceColor',fly_color)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_offset_above_thresh,'FaceColor',fly_color)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_offset_above_thresh,'FaceColor',fly_color)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution (above threshold)');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_offset_above_thresh,'FaceColor',fly_color)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_offset_above_thresh,'FaceColor',fly_color)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset_above_thresh,'FaceColor',fly_color)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_offset_above_thresh,'FaceColor',fly_color)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_offset_above_thresh,'FaceColor',fly_color)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution (above threshold)');
end

saveas(gcf,[path,'\analysis\plots\offset_evolution_thresh.png']);

%% Plot circ_std of offset per block

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    offset_var = [circ_std(pre_wind_offset_above_thresh),circ_std(pre_panels_offset_above_thresh),circ_std(combined_offset_above_thresh),circ_std(post_wind_offset_above_thresh),circ_std(post_panels_offset_above_thresh)];
    
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    offset_var = [circ_std(pre_panels_offset_above_thresh),circ_std(pre_wind_offset_above_thresh),circ_std(combined_offset_above_thresh),circ_std(post_panels_offset_above_thresh),circ_std(post_wind_offset_above_thresh)];
     
end

figure,
plot(offset_var,'-ko');
ylabel('circ_std offset');
xlabel('Block #');
xlim([0 6]);

saveas(gcf,[path,'\analysis\plots\circ_std_offset_thresh.png']);


%% Correct offset variation for the duration of the trial

%To do this, I will randomly sample timepoints in the longest trial to
%match the length of the shortest ones
min_frames = min([length(pre_wind_offset_above_thresh),length(pre_panels_offset_above_thresh),length(combined_offset_above_thresh),length(post_wind_offset_above_thresh),length(post_panels_offset_above_thresh)]);

pre_wind_offset_r = randsample(pre_wind_offset_above_thresh,min_frames);
pre_panels_offset_r = randsample(pre_panels_offset_above_thresh,min_frames);
combined_offset_r = randsample(combined_offset_above_thresh,min_frames);
post_wind_offset_r = randsample(post_wind_offset_above_thresh,min_frames);
post_panels_offset_r = randsample(post_panels_offset_above_thresh,min_frames);

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    offset_var_r = [circ_std(pre_wind_offset_r),circ_std(pre_panels_offset_r),circ_std(combined_offset_r),circ_std(post_wind_offset_r),circ_std(post_panels_offset_r)];
    
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    offset_var_r = [circ_std(pre_panels_offset_r),circ_std(pre_wind_offset_r),circ_std(combined_offset_r),circ_std(post_panels_offset_r),circ_std(post_wind_offset_r)];
     
end

figure,
plot(offset_var_r,'-ko');
ylabel('circ_std offset');
xlabel('Block #');
xlim([0 6]);

saveas(gcf,[path,'\analysis\plots\circ_std_offset_thresh_equal_length.png']);

%% Plot circ_mean of offset per block

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    offset_mean = [circ_mean(pre_wind_offset_above_thresh),circ_mean(pre_panels_offset_above_thresh),circ_mean(combined_offset_above_thresh),circ_mean(post_wind_offset_above_thresh),circ_mean(post_panels_offset_above_thresh)];
    
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    offset_mean = [circ_mean(pre_panels_offset_above_thresh),circ_mean(pre_wind_offset_above_thresh),circ_mean(combined_offset_above_thresh),circ_mean(post_panels_offset_above_thresh),circ_mean(post_wind_offset_above_thresh)];
     
end

figure,
subplot(2,1,1)
plot(rad2deg(offset_mean),'-ko');
ylabel('circ_mean offset');
xlabel('Block #');
xlim([0 6]);
ylim([-180 180]);

subplot(2,1,2)
plot([nan,diff(rad2deg(offset_mean))],'-ko');
ylabel('diff(circ_mean offset)');
xlabel('Block #');
xlim([0 6]);
ylim([-180 180]);
yline(0,'r','linestyle','--')

saveas(gcf,[path,'\analysis\plots\circ_mean_offset_thresh.png']);


%% Heading evolution


if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_heading,'FaceColor',fly_color)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_heading,'FaceColor',fly_color)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading,'FaceColor',fly_color)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_heading,'FaceColor',fly_color)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_heading,'FaceColor',fly_color)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_heading,'FaceColor',fly_color)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_heading,'FaceColor',fly_color)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading,'FaceColor',fly_color)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_heading,'FaceColor',fly_color)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_heading,'FaceColor',fly_color)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');
end

saveas(gcf,[path,'\analysis\plots\heading_evolution.png']);

%% Thresholded heading evolution

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_heading_thresh,'FaceColor',fly_color)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_heading_thresh,'FaceColor',fly_color)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading_thresh,'FaceColor',fly_color)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_heading_thresh,'FaceColor',fly_color)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_heading_thresh,'FaceColor',fly_color)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_heading_thresh,'FaceColor',fly_color)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_heading_thresh,'FaceColor',fly_color)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading_thresh,'FaceColor',fly_color)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_heading_thresh,'FaceColor',fly_color)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_heading_thresh,'FaceColor',fly_color)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');
end

saveas(gcf,[path,'\analysis\plots\heading_evolution_thresh.png']);

%% Plot circ_mean of heading per block

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    heading_mean = [circ_mean(pre_wind_heading_thresh),circ_mean(pre_panels_heading_thresh),circ_mean(combined_heading_thresh),circ_mean(post_wind_heading_thresh),circ_mean(post_panels_heading_thresh)];
    
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    heading_mean = [circ_mean(pre_panels_heading_thresh),circ_mean(pre_wind_heading_thresh),circ_mean(combined_heading_thresh),circ_mean(post_panels_heading_thresh),circ_mean(post_wind_heading_thresh)];
     
end

figure,
subplot(2,1,1)
plot(rad2deg(heading_mean),'-ko');
ylabel('circ_mean heading');
xlabel('Block #');
xlim([0 6]);
ylim([-180 180]);

subplot(2,1,2)
plot([nan,diff(rad2deg(heading_mean))],'-ko');
ylabel('diff(circ_mean heading)');
xlabel('Block #');
xlim([0 6]);
ylim([-180 180]);
yline(0,'r','linestyle','--')

saveas(gcf,[path,'\analysis\plots\circ_mean_heading.png']);


%%  Bump parameter evolution using threshold

if sessions.initial_cl_wind < sessions.initial_cl_bar
    allBM_thresh = [meanBM_thresh_pre_wind;meanBM_thresh_pre_panels;meanBM_thresh_combined;meanBM_thresh_post_wind;meanBM_thresh_post_panels];
    allBW_thresh = [meanBW_thresh_pre_wind;meanBW_thresh_pre_panels;meanBW_thresh_combined;meanBW_thresh_post_wind;meanBW_thresh_post_panels];
    all_total_mvt_thresh = [mean_total_mvt_pre_wind;mean_total_mvt_pre_panels;mean_total_mvt_combined;mean_total_mvt_post_wind;mean_total_mvt_post_panels];
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    allBM_thresh = [meanBM_thresh_pre_panels;meanBM_thresh_pre_wind;meanBM_thresh_combined;meanBM_thresh_post_panels;meanBM_thresh_post_wind];  
    allBW_thresh = [meanBW_thresh_pre_panels;meanBW_thresh_pre_wind;meanBW_thresh_combined;meanBW_thresh_post_panels;meanBW_thresh_post_wind]; 
    all_total_mvt_thresh = [mean_total_mvt_pre_panels;mean_total_mvt_pre_wind;mean_total_mvt_combined;mean_total_mvt_post_panels;mean_total_mvt_post_wind]; 
end

figure('Position',[100 100 1000 600]),
subplot(1,2,1)
yyaxis left
plot(allBM_thresh,'-o')
ylim([0 3]);
ylabel('Bump magnitude');
yyaxis right
plot(all_total_mvt_thresh,'-o')
xlim([0 6]);
ylim([0 300]);
ylabel('Total movement (deg/s)');
xlabel('Block #');

subplot(1,2,2)
yyaxis left
plot(allBW_thresh,'-o')
ylim([0 3.5]);
ylabel('Bump width');
yyaxis right
plot(all_total_mvt_thresh,'-o')
xlim([0 6]);
ylim([0 300]);

%% Save variables

save([path,'\analysis\data.mat'],'offset_var','offset_var_r','offset_mean','heading_mean','allBM_thresh','allBW_thresh','all_total_mvt_thresh')

%% Clear

close all; clear all