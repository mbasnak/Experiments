%Full experiment analysis
%Code to analyze the full experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Analyze initial closed-loop panels

load([path,'\analysis\analysis_sid_',num2str(sessions.initial_cl_bar),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(4,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
pre_panels_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
pre_panels_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.bump_magnitude,'LineWidth',1.5,'color','k')
title('Bump magnitude');
xlabel('Time (sec)');

suptitle('Initial trial with just panels');

saveas(gcf,[path,'\analysis\plots\initial_panels_old.png']);

%% Analyze inital open-loop panels

for sid = 1:length(sessions.initial_ol_bar)
    
    load([path,'\analysis\analysis_sid_',num2str(sessions.initial_ol_bar(sid)),'_tid_0.mat'])
    
    figure,
    subplot(3,1,1)
    dff = data.dff_matrix;
    imagesc(dff)
    colormap(flipud(gray))
    title('EPG activity');
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,2)
    bump_pos = wrapTo180(rad2deg(data.phase));
    %Remove wrapped lines to plot
    [x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
    plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
    hold on
    stim = wrapTo180(data.panel_angle);
    [x_out_stim,stim_to_plot] = removeWrappedLines(data.time,stim);
    plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
    legend('Bump estimate','Stim position','Location','best');
    title('Bump and stim position');
    ylim([-180 180]);
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,3)
    offset = wrapTo180(rad2deg(circ_dist(data.phase',deg2rad(data.panel_angle))));
    [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
    plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
    title('Offset');
    xlabel('Time (sec)');
    ylim([-180 180]);
    
    suptitle('Initial OL panels');
    
    saveas(gcf,[path,'\analysis\plots\initial_OL_panels_old_trial',num2str(sid),'.png']);
    
end

%%  Analyze initial closed-loop wind

load([path,'\analysis\analysis_sid_',num2str(sessions.initial_cl_wind),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(4,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
pre_wind_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
pre_wind_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.bump_magnitude,'LineWidth',1.5,'color','k')
title('Bump magnitude');
xlabel('Time (sec)');

suptitle('Initial trial with just wind');

saveas(gcf,[path,'\analysis\plots\initial_wind_old.png']);


%% Analyze inital open-loop wind

for sid = 1:length(sessions.initial_ol_wind)
    
    load([path,'\analysis\analysis_sid_',num2str(sessions.initial_ol_wind(sid)),'_tid_0.mat'])
    
    figure,
    subplot(3,1,1)
    dff = data.dff_matrix;
    imagesc(dff)
    colormap(flipud(gray))
    title('EPG activity');
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,2)
    bump_pos = wrapTo180(rad2deg(data.phase));
    %Remove wrapped lines to plot
    [x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
    plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
    hold on
    stim = wrapTo180(rad2deg(data.motor_pos));
    [x_out_stim,stim_to_plot] = removeWrappedLines(data.time,stim');
    plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
    legend('Bump estimate','Stim position','Location','Best');
    title('Bump and stim position');
    ylim([-180 180]);
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,3)
    offset = wrapTo180(rad2deg(circ_dist(data.phase',data.motor_pos')));
    [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
    plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
    title('Offset')
    ylim([-180 180]);
    
    suptitle('Initial OL wind');
    
    saveas(gcf,[path,'\analysis\plots\initial_OL_wind_old_trial',num2str(sid),'.png']);

end

%% Analyze the cue combination trial

load([path,'\analysis\analysis_sid_',num2str(sessions.cue_combination),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]});

subplot(4,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
combined_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','Best');
set(gca,'xticklabel',{[]});
ylim([-180 180]);

subplot(4,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
combined_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);

subplot(4,1,4)
plot(data.time,data.bump_magnitude,'LineWidth',1.5,'color','k')
title('Bump magnitude');
xlabel('Time (sec)');

suptitle('Trial with both cues');

saveas(gcf,[path,'\analysis\plots\cue_combination_old.png']);

%% Analyze final closed-loop panels

load([path,'\analysis\analysis_sid_',num2str(sessions.final_cl_bar),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(4,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
post_panels_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
post_panels_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.bump_magnitude,'LineWidth',1.5,'color','k')
title('Bump magnitude');
xlabel('Time (sec)');

suptitle('Final trial with just panels');

saveas(gcf,[path,'\analysis\plots\final_panels_old.png']);


%% Analyze final open-loop panels

for sid = 1:length(sessions.final_ol_bar)
    
    load([path,'\analysis\analysis_sid_',num2str(sessions.final_ol_bar(sid)),'_tid_0.mat'])
    
    figure,
    subplot(3,1,1)
    dff = data.dff_matrix;
    imagesc(dff)
    colormap(flipud(gray))
    title('EPG activity');
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,2)
    bump_pos = wrapTo180(rad2deg(data.phase));
    %Remove wrapped lines to plot
    [x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
    plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
    hold on
    stim = wrapTo180(data.panel_angle);
    [x_out_stim,stim_to_plot] = removeWrappedLines(data.time,stim);
    plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
    legend('Bump estimate','Stim position','Location','best');
    title('Bump and stim position');
    ylim([-180 180]);
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,3)
    offset = wrapTo180(rad2deg(circ_dist(data.phase',deg2rad(data.panel_angle))));
    [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
    plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
    title('Offset');
    xlabel('Time (sec)');
    ylim([-180 180]);
    
    suptitle('Initial trial with just panels');
    
    saveas(gcf,[path,'\analysis\plots\final_OL_panels_old_trial',num2str(sid),'.png']);
    
end


%%  Analyze final closed-loop wind

load([path,'\analysis\analysis_sid_',num2str(sessions.final_cl_wind),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff_matrix;
imagesc(dff)
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(4,1,2)
bump_pos = wrapTo180(rad2deg(data.phase));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
heading = wrapTo180(-data.heading_deg);
post_wind_heading = deg2rad(heading);
[x_out_heading,heading_to_plot] = removeWrappedLines(data.time,heading);
plot(x_out_heading,heading_to_plot,'LineWidth',1.5)
title('Bump and fly position');
legend('Bump estimate','Fly position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
offset = wrapTo180(rad2deg(circ_dist(data.phase',-data.heading)));
%store offset for later
post_wind_offset = deg2rad(offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.bump_magnitude,'LineWidth',1.5,'color','k')
title('Bump magnitude');
xlabel('Time (sec)');

suptitle('Final trial with just wind');

saveas(gcf,[path,'\analysis\plots\final_wind_old.png']);


%% Analyze final open-loop wind

for sid = 1:length(sessions.final_ol_wind)
    
    load([path,'\analysis\analysis_sid_',num2str(sessions.final_ol_wind(sid)),'_tid_0.mat'])
    
    figure,
    subplot(3,1,1)
    dff = data.dff_matrix;
    imagesc(dff)
    colormap(flipud(gray))
    title('EPG activity');
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,2)
    bump_pos = wrapTo180(rad2deg(data.phase));
    %Remove wrapped lines to plot
    [x_out_bump,bump_to_plot] = removeWrappedLines(data.time,bump_pos');
    plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
    hold on
    stim = wrapTo180(rad2deg(data.motor_pos));
    [x_out_stim,stim_to_plot] = removeWrappedLines(data.time,stim');
    plot(x_out_stim,stim_to_plot,'LineWidth',1.5)
    legend('Bump estimate','Stim position','Location','Best');
    title('Bump and stim position');
    ylim([-180 180]);
    set(gca,'xticklabel',{[]})
    
    subplot(3,1,3)
    offset = wrapTo180(rad2deg(circ_dist(data.phase',data.motor_pos')));
    [x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
    plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
    title('Offset')
    ylim([-180 180]);
    
    suptitle('Final trial with just wind');
    
    saveas(gcf,[path,'\analysis\plots\final_OL_wind_old_trial',num2str(sid),'.png']);
    
end

%% Offset evolution


if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_offset)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_offset)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_offset)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_offset)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_offset)
    title('Initial panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_offset)
    title('Initial wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_offset)
    title('Cue combination offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_offset)
    title('Final panels offset');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_offset)
    title('Final wind offset');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Offset evolution');
end

saveas(gcf,[path,'\analysis\plots\offset_evolution_old.png']);



%% Heading evolution


if sessions.initial_cl_wind < sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_wind_heading)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_panels_heading)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_wind_heading)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_panels_heading)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');

elseif sessions.initial_cl_wind > sessions.initial_cl_bar
    
    figure('Position',[100 100 1400 400]),
    
    subplot(1,5,1)
    polarhistogram(pre_panels_heading)
    title('Initial panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,2)
    polarhistogram(pre_wind_heading)
    title('Initial wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,3)
    polarhistogram(combined_heading)
    title('Cue combination heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,4)
    polarhistogram(post_panels_heading)
    title('Final panels heading');
    set(gca,'ThetaZeroLocation','top');
    
    subplot(1,5,5)
    polarhistogram(post_wind_heading)
    title('Final wind heading');
    set(gca,'ThetaZeroLocation','top');
    
    suptitle('Heading evolution');
end

saveas(gcf,[path,'\analysis\plots\heading_evolution_old.png']);

%% Clear

close all; clear all