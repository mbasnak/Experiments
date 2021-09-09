%Code to analyze the cue conflict experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Analyze initial closed-loop panels

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.pre_panels),'_tid_0.mat'])

%Heatmap
figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
pre_panels_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
pre_panels_offset = deg2rad(offset);
%store the offset minus the points where the fit is low
pre_panels_offset_above_thresh = pre_panels_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

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


saveas(gcf,[path,'\analysis\plots\pre_panels.png']);

%%  Analyze initial closed-loop wind

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.pre_wind),'_tid_0.mat'])

%Initial closed-loop panels
figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
pre_wind_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
pre_wind_offset = deg2rad(offset);
pre_wind_offset_above_thresh = pre_wind_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

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

saveas(gcf,[path,'\analysis\plots\pre_wind.png']);

%% Analyze the cue combination trial

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.combined),'_tid_0.mat'])

%Initial closed-loop panels
figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
combined_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
combined_offset = deg2rad(offset);
combined_offset_above_thresh = combined_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

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

suptitle('Trial with both cues');

saveas(gcf,[path,'\analysis\plots\cue_combination.png']);

%% Analyze initial closed-loop panels

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.post_panels),'_tid_0.mat'])

%Initial closed-loop panels
figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
post_panels_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
post_panels_offset = deg2rad(offset);
post_panels_offset_above_thresh = post_panels_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

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

saveas(gcf,[path,'\analysis\plots\post_panels.png']);

%%  Analyze final closed-loop wind

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.post_wind),'_tid_0.mat'])

%Initial closed-loop panels
figure('Position',[100 100 1200 800]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
post_wind_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(5,1,3)
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%store offset for later
post_wind_offset = deg2rad(offset);
post_wind_offset_above_thresh = post_wind_offset(continuous_data.adj_rs>=0.5);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

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

saveas(gcf,[path,'\analysis\plots\post_wind.png']);


%% Offset evolution

figure('Position',[100 100 1400 400]),
subplot(1,5,1)
polarhistogram(pre_panels_offset)
title('Initial panels offset');

subplot(1,5,2)
polarhistogram(pre_wind_offset)
title('Initial wind offset');

subplot(1,5,3)
polarhistogram(combined_offset)
title('Cue combination offset');

subplot(1,5,4)
polarhistogram(post_panels_offset)
title('Final panels offset');

subplot(1,5,5)
polarhistogram(post_wind_offset)
title('Final wind offset');

suptitle('Offset evolution');

saveas(gcf,[path,'\analysis\plots\offset_evolution.png']);

%% Offset evolution with only the datapoints where the fit was above threshold

figure('Position',[100 100 1400 400]),
subplot(1,5,1)
polarhistogram(pre_panels_offset_above_thresh)
title('Initial panels offset');

subplot(1,5,2)
polarhistogram(pre_wind_offset_above_thresh)
title('Initial wind offset');

subplot(1,5,3)
polarhistogram(combined_offset_above_thresh)
title('Cue combination offset');

subplot(1,5,4)
polarhistogram(post_panels_offset_above_thresh)
title('Final panels offset');

subplot(1,5,5)
polarhistogram(post_wind_offset_above_thresh)
title('Final wind offset');

suptitle('Offset evolution, fit above threshold');

saveas(gcf,[path,'\analysis\plots\offset_evolution_above_thresh.png']);

%% Heading evolution

figure('Position',[100 100 1400 400]),
subplot(1,5,1)
polarhistogram(pre_panels_heading)
set(gca,'ThetaZeroLocation','top')
title('Initial panels heading');

subplot(1,5,2)
polarhistogram(pre_wind_heading)
set(gca,'ThetaZeroLocation','top')
title('Initial wind heading');

subplot(1,5,3)
polarhistogram(combined_heading)
set(gca,'ThetaZeroLocation','top')
title('Cue combination heading');

subplot(1,5,4)
polarhistogram(post_panels_heading)
set(gca,'ThetaZeroLocation','top')
title('Final panels heading');

subplot(1,5,5)
polarhistogram(post_wind_heading)
set(gca,'ThetaZeroLocation','top')
title('Final wind heading');

suptitle('Heading evolution');

saveas(gcf,[path,'\analysis\plots\heading_evolution.png']);

%%
close all; clear all;