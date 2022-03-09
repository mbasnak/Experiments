%Full experiment analysis
%Code to analyze the full experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Set colormap

folderNames = dir(path(1:69));
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
    %if strcmp(flyNames(fly).name,path(54:end))
    if strcmp(flyNames(fly).name,path(71:end)) %70 for low reliability/71 for high reliability (code this)
        fly_ID = fly;
    end
end

%Set colors for individual fly plots
colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
    0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880;...
    0 0.4470 0.7410; 0.75, 0.1, 0.75; 0.75, 0.75, 0;...
    0 0.75 0.75; 0.75 0 0.75; 1 0.2 0.2; 0.2 1 0.2; 0.2 0.2 1;...
    0.25 0.25 0.25; 0.25 0.45 0.15; 0.70 0.15 0.50];

%Assign color for this fly's plots
fly_color = colors_for_plots(fly_ID,:);

%% Make directory to save plots

%Move to the analysis folder
contents = dir([path,'\analysis']);

%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'continuous_plots') == 0)
   mkdir(path,'\analysis\continuous_plots'); 
end

%% Analyze initial closed-loop panels

%% Full experiment plot

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
pre_panels_heading_thresh = pre_panels_heading(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
pre_panels_offset_above_thresh = deg2rad(offset(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));
[~, offset_var_pre_panels_offset_above_thresh] = circ_std(pre_panels_offset_above_thresh);
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

saveas(gcf,[path,'\analysis\continuous_plots\initial_panels_full_experiment.png']);

%% Bump parameters

%Get mean bump parameters
meanBM_thresh_pre_panels = nanmean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
meanBW_thresh_pre_panels = nanmean(continuous_data.bump_width(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
%Get mean vel
mean_total_mvt_thresh_pre_panels = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));

%Store all bump param and mvt values
allBumpMag = continuous_data.bump_magnitude;
allBumpWidth = continuous_data.bump_width;
allTotalMvt = continuous_data.total_mvt_ds;
blockType = repelem(1,1,length(continuous_data.bump_magnitude));

%% Probability of being stopped

p_stopped_pre_panels = sum(continuous_data.total_mvt_ds <= 25)/length(continuous_data.total_mvt_ds);

sec_to_frames = length(continuous_data.time)/continuous_data.time(end);
%analyze probability of being stopped in shorter timescales
%Divide into ~1 min bouts and plot
bout_boundaries = 1:floor(60*sec_to_frames):length(continuous_data.bump_magnitude);
for bout = 1:length(bout_boundaries)-1
    mvt_bout = continuous_data.total_mvt_ds(bout_boundaries(bout):bout_boundaries(bout+1));
    p_stop_pre_panels(bout) = sum(mvt_bout <= 25)/length(mvt_bout);
end
figure,
plot(p_stop_pre_panels,'-o')
ylim([0 1]); ylabel('P(stopped)');
xlabel('# 60 sec bout');

saveas(gcf,[path,'\analysis\continuous_plots\prob_stopping_pre_panels.png']);

%%  Analyze initial closed-loop wind

%% Full experiment plot

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
pre_wind_heading_thresh = pre_wind_heading(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
pre_wind_offset_above_thresh = deg2rad(offset(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
[~, offset_var_pre_wind] = circ_std(pre_wind_offset_above_thresh);
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

saveas(gcf,[path,'\analysis\continuous_plots\initial_wind_full_experiment.png']);

%% Bump parameters

%Get mean bump parameters
meanBM_thresh_pre_wind = nanmean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
meanBW_thresh_pre_wind = nanmean(continuous_data.bump_width(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
%Get mean vel
mean_total_mvt_thresh_pre_wind = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));

%Store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(2,1,length(continuous_data.bump_magnitude))];


%% Probability of being stopped

p_stopped_pre_wind = sum(continuous_data.total_mvt_ds <= 25)/length(continuous_data.total_mvt_ds);

sec_to_frames = length(continuous_data.time)/continuous_data.time(end);
%analyze probability of being stopped in shorter timescales
%Divide into ~1 min bouts and plot
bout_boundaries = 1:floor(60*sec_to_frames):length(continuous_data.bump_magnitude);
for bout = 1:length(bout_boundaries)-1
    mvt_bout = continuous_data.total_mvt_ds(bout_boundaries(bout):bout_boundaries(bout+1));
    p_stop_pre_wind(bout) = sum(mvt_bout <= 25)/length(mvt_bout);
end
figure,
plot(p_stop_pre_wind,'-o')
ylim([0 1]); ylabel('P(stopped)');
xlabel('# 60 sec bout');

saveas(gcf,[path,'\analysis\continuous_plots\prob_stopping_pre_wind.png']);

%% Analyze the cue combination trial

%% Full experiment analysis

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
combined_heading_thresh = combined_heading(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25);
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
combined_offset_above_thresh = deg2rad(offset(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
[~, offset_var_combined] = circ_std(combined_offset_above_thresh);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset');
ylim([-180 180]);
xlim([0 x_out_offset(end-1)]);

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

saveas(gcf,[path,'\analysis\continuous_plots\cue_combination_full_experiment.png']);

%% Bump parameters

%Get mean bump parameters
meanBM_thresh_combined = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
meanBW_thresh_combined = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));

%Get mean vel
mean_total_mvt_thresh_combined = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));

%Store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(3,1,length(continuous_data.bump_magnitude))];


%% Probability of being stopped

p_stopped_cue_combination = sum(continuous_data.total_mvt_ds <= 25)/length(continuous_data.total_mvt_ds);

sec_to_frames = length(continuous_data.time)/continuous_data.time(end);

%analyze probability of being stopped in shorter timescales
%Divide into ~1 min bouts and plot
bout_boundaries = 1:floor(60*sec_to_frames):length(continuous_data.bump_magnitude);
for bout = 1:length(bout_boundaries)-1
    mvt_bout = continuous_data.total_mvt_ds(bout_boundaries(bout):bout_boundaries(bout+1));
    p_stop_cue_combination(bout) = sum(mvt_bout <= 25)/length(mvt_bout);
end
figure,
plot(p_stop_cue_combination,'-o')
ylim([0 1]); ylabel('P(stopped)');
xlabel('# 60 sec bout');

saveas(gcf,[path,'\analysis\continuous_plots\prob_stopping_cue_combination.png']);

%% Analyze final closed-loop panels

%% Full experiment analysis

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
post_panels_heading_thresh = post_panels_heading(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
post_panels_offset_above_thresh = deg2rad(offset(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds>25));
[~, offset_var_post_panels] = circ_std(post_panels_offset_above_thresh);
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

saveas(gcf,[path,'\analysis\continuous_plots\final_panels_full_experiment.png']);

%% Bump parameters

%Get mean bump parameters
meanBM_thresh_post_panels = nanmean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));
meanBW_thresh_post_panels = nanmean(continuous_data.bump_width(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));

%Get fly vel
mean_total_mvt_thresh_post_panels = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));

%Store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(4,1,length(continuous_data.bump_magnitude))];


%% Probability of being stopped

p_stopped_post_panels = sum(continuous_data.total_mvt_ds <= 25)/length(continuous_data.total_mvt_ds);

sec_to_frames = length(continuous_data.time)/continuous_data.time(end);

%analyze probability of being stopped in shorter timescales
%Divide into ~1 min bouts and plot
bout_boundaries = 1:floor(60*sec_to_frames):length(continuous_data.bump_magnitude);
for bout = 1:length(bout_boundaries)-1
    mvt_bout = continuous_data.total_mvt_ds(bout_boundaries(bout):bout_boundaries(bout+1));
    p_stop_post_panels(bout) = sum(mvt_bout <= 25)/length(mvt_bout);
end
figure,
plot(p_stop_post_panels,'-o')
ylim([0 1]); ylabel('P(stopped)');
xlabel('# 60 sec bout');

saveas(gcf,[path,'\analysis\continuous_plots\prob_stopping_post_panels.png']);

%%  Analyze final closed-loop wind

%% Full experiment analysis

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
post_wind_heading_thresh = post_wind_heading(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
post_wind_offset_above_thresh = deg2rad(offset(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));
[~, offset_var_post_wind] = circ_std(post_wind_offset_above_thresh);
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

saveas(gcf,[path,'\analysis\continuous_plots\final_wind_full_experiment.png']);

%% Bump parameters

%Get mean bump parameters
meanBM_thresh_post_wind = nanmean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));
meanBW_thresh_post_wind = nanmean(continuous_data.bump_width(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));

%Get fly vel
mean_total_mvt_thresh_post_wind = nanmean(continuous_data.total_mvt_ds(continuous_data.adj_rs>=0.5 & continuous_data.total_mvt_ds > 25));

%Store all bump param and mvt values
allBumpMag = [allBumpMag,continuous_data.bump_magnitude];
allBumpWidth = [allBumpWidth,continuous_data.bump_width];
allTotalMvt = [allTotalMvt,continuous_data.total_mvt_ds];
blockType = [blockType,repelem(5,1,length(continuous_data.bump_magnitude))];


%% Probability of being stopped

p_stopped_post_wind = sum(continuous_data.total_mvt_ds <= 25)/length(continuous_data.total_mvt_ds);

sec_to_frames = length(continuous_data.time)/continuous_data.time(end);

%analyze probability of being stopped in shorter timescales
%Divide into ~1 min bouts and plot
bout_boundaries = 1:floor(60*sec_to_frames):length(continuous_data.bump_magnitude);
for bout = 1:length(bout_boundaries)-1
    mvt_bout = continuous_data.total_mvt_ds(bout_boundaries(bout):bout_boundaries(bout+1));
    p_stop_post_wind(bout) = sum(mvt_bout <= 25)/length(mvt_bout);
end
figure,
plot(p_stop_post_wind,'-o')
ylim([0 1]); ylabel('P(stopped)');
xlabel('# 60 sec bout');

saveas(gcf,[path,'\analysis\continuous_plots\prob_stopping_post_wind.png']);

%% Thresholded offset evolution

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_offset_above_thresh,15,'FaceColor',fly_color)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_offset_above_thresh,15,'FaceColor',fly_color)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset_above_thresh,15,'FaceColor',fly_color)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_offset_above_thresh,15,'FaceColor',fly_color)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_offset_above_thresh,15,'FaceColor',fly_color)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution (above threshold)');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_offset_above_thresh,15,'FaceColor',fly_color)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_offset_above_thresh,15,'FaceColor',fly_color)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset_above_thresh,15,'FaceColor',fly_color)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_offset_above_thresh,15,'FaceColor',fly_color)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_offset_above_thresh,15,'FaceColor',fly_color)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution (above threshold)');
end

saveas(gcf,[path,'\analysis\continuous_plots\offset_evolution_thresh.png']);


%% Offset variability

%1) Determine the num frames of shorter session
min_frames = min([length(pre_wind_offset_above_thresh),length(pre_panels_offset_above_thresh),length(combined_offset_above_thresh),length(post_wind_offset_above_thresh),length(post_panels_offset_above_thresh)]);

%2) Randomly sample timepoints in the longest sessions to
%match the length of the shorter one
pre_wind_offset_r = randsample(pre_wind_offset_above_thresh,min_frames);
pre_panels_offset_r = randsample(pre_panels_offset_above_thresh,min_frames);
combined_offset_r = randsample(combined_offset_above_thresh,min_frames);
post_wind_offset_r = randsample(post_wind_offset_above_thresh,min_frames);
post_panels_offset_r = randsample(post_panels_offset_above_thresh,min_frames);

%3) Compute offset variability
[~,offset_var_pre_wind_offset_r] = circ_std(pre_wind_offset_r);
[~,offset_var_pre_panels_offset_r] = circ_std(pre_panels_offset_r);
[~,offset_var_combined_offset_r] = circ_std(combined_offset_r);
[~,offset_var_post_wind_offset_r] = circ_std(post_wind_offset_r);
[~,offset_var_post_panels_offset_r] = circ_std(post_panels_offset_r);

%Combine data
if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    offset_var_r = [offset_var_pre_wind_offset_r,offset_var_pre_panels_offset_r,offset_var_combined_offset_r,offset_var_post_wind_offset_r,offset_var_post_panels_offset_r];
    
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    offset_var_r = [offset_var_pre_panels_offset_r,offset_var_pre_wind_offset_r,offset_var_combined_offset_r,offset_var_post_panels_offset_r,offset_var_post_wind_offset_r];
     
end

%Plot variability per block
figure,
plot(offset_var_r,'-ko');
ylabel('circ_std offset');
xlabel('Block #');
xlim([0 6]); ylim([0 1.5]);

saveas(gcf,[path,'\analysis\continuous_plots\circ_std_offset_thresh_equal_length.png']);

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
for block = 1:length(offset_mean)-1
    offset_diff(block) = abs(rad2deg(circ_dist(offset_mean(block+1),offset_mean(block))));
end
plot([nan,offset_diff],'-ko');
ylabel('diff(circ_mean offset)');
xlabel('Block #');
xlim([0 6]);
ylim([0 180]);
yline(0,'r','linestyle','--')

saveas(gcf,[path,'\analysis\continuous_plots\circ_mean_offset_thresh.png']);


%% Overlay offset mean over offset evolution plot

if sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(1),offset_mean(1)],[0,rl(2)],'k','linewidth',2)
    title('Initial visual offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,2)
    polarhistogram(pre_wind_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(2),offset_mean(2)],[0,rl(2)],'k','linewidth',2)
    title('Initial wind offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,3)
    polarhistogram(combined_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(3),offset_mean(3)],[0,rl(2)],'k','linewidth',2)
    title('Cue combination offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,4)
    polarhistogram(post_panels_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(4),offset_mean(4)],[0,rl(2)],'k','linewidth',2)
    title('Final visual offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,5)
    polarhistogram(post_wind_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(5),offset_mean(5)],[0,rl(2)],'k','linewidth',2)
    title('Final wind offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    
else
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(1),offset_mean(1)],[0,rl(2)],'k','linewidth',2)
    title('Initial wind offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,2)
    polarhistogram(pre_panels_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(2),offset_mean(2)],[0,rl(2)],'k','linewidth',2)
    title('Initial visual offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,3)
    polarhistogram(combined_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(3),offset_mean(3)],[0,rl(2)],'k','linewidth',2)
    title('Cue combination offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,4)
    polarhistogram(post_wind_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(4),offset_mean(4)],[0,rl(2)],'k','linewidth',2)
    title('Final wind offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
    
    subplot(1,5,5)
    polarhistogram(post_panels_offset_above_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([offset_mean(5),offset_mean(5)],[0,rl(2)],'k','linewidth',2)
    title('Final visual offset','fontsize',12);
    set(gca,'ThetaZeroLocation','top');
    Ax = gca;
    Ax.RTickLabel = [];
        
end

saveas(gcf,[path,'\analysis\continuous_plots\offset_evo_with_mean.png']);


%% Thresholded heading evolution

if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_heading_thresh,15,'FaceColor',fly_color)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_heading_thresh,15,'FaceColor',fly_color)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading_thresh,15,'FaceColor',fly_color)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_heading_thresh,15,'FaceColor',fly_color)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_heading_thresh,15,'FaceColor',fly_color)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_heading_thresh,15,'FaceColor',fly_color)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_heading_thresh,15,'FaceColor',fly_color)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading_thresh,15,'FaceColor',fly_color)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_heading_thresh,15,'FaceColor',fly_color)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_heading_thresh,15,'FaceColor',fly_color)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');
end

saveas(gcf,[path,'\analysis\continuous_plots\heading_evolution_thresh.png']);

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
for block = 1:length(offset_mean)-1
    heading_diff(block) = abs(rad2deg(circ_dist(heading_mean(block+1),heading_mean(block))));
end
plot([nan,heading_diff],'-ko');
ylabel('diff(circ_mean heading)');
xlabel('Block #');
xlim([0 6]);
ylim([0 180]);
yline(0,'r','linestyle','--')

saveas(gcf,[path,'\analysis\continuous_plots\circ_mean_heading.png']);


%% Overlay heading mean over offset evolution plot

if sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(1),heading_mean(1)],[0,rl(2)],'k','linewidth',2)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(2),heading_mean(2)],[0,rl(2)],'k','linewidth',2)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(3),heading_mean(3)],[0,rl(2)],'k','linewidth',2)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(4),heading_mean(4)],[0,rl(2)],'k','linewidth',2)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(5),heading_mean(5)],[0,rl(2)],'k','linewidth',2)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('heading evolution (above threshold)');
    
else
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(1),heading_mean(1)],[0,rl(2)],'k','linewidth',2)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(2),heading_mean(2)],[0,rl(2)],'k','linewidth',2)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(3),heading_mean(3)],[0,rl(2)],'k','linewidth',2)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(4),heading_mean(4)],[0,rl(2)],'k','linewidth',2)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_heading_thresh,15,'FaceColor',fly_color)
    hold on
    rl = rlim;
    polarplot([heading_mean(5),heading_mean(5)],[0,rl(2)],'k','linewidth',2)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution (above threshold)');
    
end

saveas(gcf,[path,'\analysis\continuous_plots\heading_evo_with_mean.png']);

%%  Bump parameter evolution using threshold

if sessions.initial_cl_wind < sessions.initial_cl_bar
    allBM_thresh = [meanBM_thresh_pre_wind;meanBM_thresh_pre_panels;meanBM_thresh_combined;meanBM_thresh_post_wind;meanBM_thresh_post_panels];
    allBW_thresh = [meanBW_thresh_pre_wind;meanBW_thresh_pre_panels;meanBW_thresh_combined;meanBW_thresh_post_wind;meanBW_thresh_post_panels];
    all_total_mvt_thresh = [mean_total_mvt_thresh_pre_wind;mean_total_mvt_thresh_pre_panels;mean_total_mvt_thresh_combined;mean_total_mvt_thresh_post_wind;mean_total_mvt_thresh_post_panels];
elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    allBM_thresh = [meanBM_thresh_pre_panels;meanBM_thresh_pre_wind;meanBM_thresh_combined;meanBM_thresh_post_panels;meanBM_thresh_post_wind];  
    allBW_thresh = [meanBW_thresh_pre_panels;meanBW_thresh_pre_wind;meanBW_thresh_combined;meanBW_thresh_post_panels;meanBW_thresh_post_wind]; 
    all_total_mvt_thresh = [mean_total_mvt_thresh_pre_panels;mean_total_mvt_thresh_pre_wind;mean_total_mvt_thresh_combined;mean_total_mvt_thresh_post_panels;mean_total_mvt_thresh_post_wind]; 
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

saveas(gcf,[path,'\analysis\continuous_plots\bump_par_evolution_thresh.png']);

%% Link plasticity to initial differences between bar and wind offset

initial_cue_diff = rad2deg(circ_dist(circ_mean(pre_wind_offset_above_thresh),circ_mean(pre_panels_offset_above_thresh)));
bar_offset_diff = rad2deg(circ_dist(circ_mean(pre_panels_offset_above_thresh),circ_mean(post_panels_offset_above_thresh)));
wind_offset_diff = rad2deg(circ_dist(circ_mean(pre_wind_offset_above_thresh),circ_mean(post_wind_offset_above_thresh)));

%% Probability of being stopped per block

if sessions.initial_cl_wind > sessions.initial_cl_bar
    p_stopped = [p_stopped_pre_panels,p_stopped_pre_wind,p_stopped_cue_combination,p_stopped_post_panels,p_stopped_post_wind];
else
    p_stopped = [p_stopped_pre_wind,p_stopped_pre_panels,p_stopped_cue_combination,p_stopped_post_wind,p_stopped_post_panels];
end

%% Determine and save configuration

if sessions.initial_cl_wind > sessions.initial_cl_bar
    configuration = 1;
else
    configuration = 2;
end

%% Combine certain variables in table

summary_data = table(allBumpMag',allBumpWidth',allTotalMvt',blockType','VariableNames',{'BumpMag','BumpWidth','TotalMvt','BlockType'});


%% Save variables

save([path,'\analysis\data.mat'],'summary_data','offset_var_r','offset_mean','offset_mean_final','heading_mean','allBM_thresh','allBW_thresh','all_total_mvt_thresh','allBM_thresh_final','allBW_thresh_final','all_total_mvt_thresh_final','initial_cue_diff','bar_offset_diff','wind_offset_diff','initial_cue_diff_last_part','bar_offset_diff_last_part','wind_offset_diff_last_part','configuration','mean_bout_BM_pre_panels','mean_bout_BW_pre_panels','mean_bout_BM_pre_wind','mean_bout_BW_pre_wind','mean_bout_BM_cue_combination','mean_bout_BW_cue_combination','mean_bout_BM_post_panels','mean_bout_BW_post_panels','mean_bout_BM_post_wind','mean_bout_BW_post_wind','p_stop_pre_panels','p_stop_pre_wind','p_stop_cue_combination','p_stop_post_panels','p_stop_post_wind','p_stopped')

%% Clear

close all; clear all