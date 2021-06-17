%Code to analyze the cue conflict experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp34\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Analyze initial closed-loop panels

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.pre_panels),'_tid_0.mat'])

%Initial closed-loop panels
figure,
subplot(3,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(continuous_data.offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Initial trial with just panels');

saveas(gcf,[path,'\analysis\plots\pre_panels.png']);

%%  Analyze initial closed-loop wind

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.pre_wind),'_tid_0.mat'])

%Initial closed-loop panels
figure,
subplot(3,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(continuous_data.offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Initial trial with just wind');

saveas(gcf,[path,'\analysis\plots\pre_wind.png']);

%% Analyze the cue combination trial

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.combined),'_tid_0.mat'])

%Initial closed-loop panels
figure,
subplot(3,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(continuous_data.offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Trial with both cues');

saveas(gcf,[path,'\analysis\plots\cue_combination.png']);

%% Analyze initial closed-loop panels

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.post_panels),'_tid_0.mat'])

%Initial closed-loop panels
figure,
subplot(3,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(continuous_data.offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Final trial with just panels');

saveas(gcf,[path,'\analysis\plots\post_panels.png']);

%%  Analyze final closed-loop wind

load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.post_wind),'_tid_0.mat'])

%Initial closed-loop panels
figure,
subplot(3,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-continuous_data.heading_deg);
[x_out_heading,heading_to_plot] = removeWrappedLines(continuous_data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
%xlim([1 length(continuous_data.bump_pos)]);
ylim([-180 180]);

subplot(3,1,3)
offset = wrapTo180(continuous_data.offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
title('Offset')
ylim([-180 180]);

suptitle('Final trial with just wind');

saveas(gcf,[path,'\analysis\plots\post_wind.png']);
