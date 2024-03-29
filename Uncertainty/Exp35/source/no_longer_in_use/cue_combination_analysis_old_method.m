%Code to analyze the cue conflict experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Analyze initial closed-loop panels

load([path,'\analysis\analysis_sid_',num2str(sessions.pre_panels),'_tid_0.mat'])

figure,
subplot(3,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
pre_panels_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Initial trial with just panels');

saveas(gcf,[path,'\analysis\plots\pre_panels_old.png']);

%%  Analyze initial closed-loop wind

load([path,'\analysis\analysis_sid_',num2str(sessions.pre_wind),'_tid_0.mat'])

figure,
subplot(3,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
pre_wind_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Initial trial with just wind');

saveas(gcf,[path,'\analysis\plots\pre_wind_old.png']);

%% Analyze the cue combination trial

load([path,'\analysis\analysis_sid_',num2str(sessions.combined),'_tid_0.mat'])

figure,
subplot(3,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
combined_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Trial with both cues');

saveas(gcf,[path,'\analysis\plots\cue_combination_old.png']);

%% Analyze initial closed-loop panels

load([path,'\analysis\analysis_sid_',num2str(sessions.post_panels),'_tid_0.mat'])

figure,
subplot(3,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
post_panels_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Final trial with just panels');

saveas(gcf,[path,'\analysis\plots\post_panels_old.png']);

%%  Analyze final closed-loop wind

load([path,'\analysis\analysis_sid_',num2str(sessions.post_wind),'_tid_0.mat'])

figure,
subplot(3,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
post_wind_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Final trial with just wind');

saveas(gcf,[path,'\analysis\plots\post_wind_old.png']);


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

saveas(gcf,[path,'\analysis\plots\offset_evolution_old.png']);

%% Clear

close all; clear all