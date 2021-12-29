%Full experiment analysis
%Code to analyze the full experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp38\data\third_version\');
load(fullfile(path,file))

warning('off','all');

%% Make directory to save plots

%Move to the analysis folder
cd(path)
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(path,'\plots'); 
end
%List the contents of the 'plots' folder
cd([path,'\plots\'])

%% Set colormap

folderNames = dir(path(1:67));
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
    if strcmp(flyNames(fly).name,path(68:end-10))
        fly_ID = fly;
    end
end

colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
    0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880;...
    0 0.4470 0.7410; 0.75, 0.1, 0.75; 0.75, 0.75, 0];

fly_color = colors_for_plots(fly_ID,:);


%% Determine when the stimuli are on

panels_on = continuous_data.fr_y_ds>7;
wind_on = continuous_data.wind_valve>2;

figure,
subplot(2,1,1)
plot(panels_on)
ylim([-1 2]);
title('Panels on');

subplot(2,1,2)
plot(wind_on)
ylim([-1 2]);
title('Wind on');


%% Determine stimulus configuration

if mode(panels_on(1:100)) == 1
    configuration = 1; %bar first
else
    configuration = 2; %wind first
end

%% Look for change in stimuli

%Find the frames where the stimuli change
panels_change = abs(diff(panels_on));
panels_change_frames = find(panels_change>0.5);
wind_change = abs(diff(wind_on));
wind_change_frames = find(wind_change>0.5);

%Conversion factors
sec_to_frames = length(continuous_data.dff_matrix)/continuous_data.time(end);
frames_to_sec = continuous_data.time(end)/length(continuous_data.dff_matrix);


%% Find the jump frames

if configuration == 1
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1200*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(1800*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2400*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(3000*length(continuous_data.dff_matrix)/continuous_data.time(end))];
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1500*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2100*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2700*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(3300*length(continuous_data.dff_matrix)/continuous_data.time(end))];   
else
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1200*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(1800*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2400*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(3000*length(continuous_data.dff_matrix)/continuous_data.time(end))];
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1500*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2100*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2700*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(3300*length(continuous_data.dff_matrix)/continuous_data.time(end))];   
end

%Correct bar jump frames
abs_diff_bar_signal = abs(diff(unwrap(deg2rad(continuous_data.panel_angle))));
for jump = 1:length(coded_bar_jump_frames)
    frame_segment = [coded_bar_jump_frames(jump)-100:coded_bar_jump_frames(jump)+100];
    [~,I_bar_jump_frames(jump)] = max(abs_diff_bar_signal(frame_segment));
    real_bar_jump_frames(jump) = coded_bar_jump_frames(jump) + I_bar_jump_frames(jump) - 100;
end
%correct for the flies for which the method didn't work well
if (fly_ID == 1 | fly_ID == 2)
    real_bar_jump_frames(1) = floor(1489*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2089*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2689*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3289*length(continuous_data.dff_matrix)/continuous_data.time(end));    
end
real_bar_jump_sec = real_bar_jump_frames*continuous_data.time(end)/length(continuous_data.dff_matrix);


%Correct wind jump frames
motor_pos = wrapTo180(rad2deg(continuous_data.motor_pos));
abs_diff_wind_signal = abs(diff(unwrap(deg2rad(motor_pos))));
for jump = 1:length(coded_wind_jump_frames)
    frame_segment = [coded_wind_jump_frames(jump)-100:coded_wind_jump_frames(jump)+100];
    [~,I_wind_jump_frames(jump)] = max(abs_diff_wind_signal(frame_segment));
    real_wind_jump_frames(jump) = coded_wind_jump_frames(jump) + I_wind_jump_frames(jump) - 100;
end
%correct for the flies for which the method didn't work well
if fly_ID == 1
    real_wind_jump_frames(1) = floor(1189*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1789*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2389*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2989*length(continuous_data.dff_matrix)/continuous_data.time(end));    
end

real_wind_jump_sec = real_wind_jump_frames*continuous_data.time(end)/length(continuous_data.dff_matrix);


%% Define block type (1 = bar, 2 = wind, 3 = both)

block = {};

if configuration == 1
    block{1} = 1;
    block{2} = 2;
    block{3} = 3;
else
    block{1} = 2;
    block{2} = 1;
    block{3} = 3;
end




%% Analyze full experiment

figure('Position',[100 100 1800 1000]),
subplot(5,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
hold on
%add jumps
xline(real_wind_jump_frames(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
xline(real_bar_jump_frames(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
legend('Wind jumps','Bar jumps');
for jump = 2:4
   xline(real_wind_jump_frames(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_frames(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
title('EPG activity');
set(gca,'xticklabel',{[]})

subplot(5,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
%add bar position
bar_pos = wrapTo180(continuous_data.panel_angle);
[x_out_bar,bar_to_plot] = removeWrappedLines(continuous_data.time(panels_on),bar_pos(panels_on));
plot(x_out_bar,bar_to_plot,'.','LineWidth',1.5)
%add motor position
[x_out_motor,motor_to_plot] = removeWrappedLines(continuous_data.time(wind_on),motor_pos(wind_on)');
plot(x_out_motor,motor_to_plot,'.','LineWidth',1.5)
%add jumps
xline(real_wind_jump_sec(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
xline(real_bar_jump_sec(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
legend('Wind jumps','Bar jumps');
for jump = 2:4
   xline(real_wind_jump_sec(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_sec(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
title('Bump and stimuli position');
legend('Bump estimate','Bar position','Wind position')
ylim([-180 180]);
xlim([0 x_out_bump(end)]);
set(gca,'xticklabel',{[]})

subplot(5,1,3)
%offset with respect to bar
bar_offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle))));
%offst with respect to wind
wind_offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos,deg2rad(motor_pos))));
[x_out_bar_offset,bar_offset_to_plot] = removeWrappedLines(continuous_data.time,bar_offset);
plot(x_out_bar_offset,bar_offset_to_plot,'LineWidth',1.5,'color',[0.8500 0.3250 0.0980])
hold on
[x_out_wind_offset,wind_offset_to_plot] = removeWrappedLines(continuous_data.time,wind_offset');
plot(x_out_wind_offset,wind_offset_to_plot,'LineWidth',1.5,'color',[0.9290 0.6940 0.1250])
%add jumps
xline(real_wind_jump_sec(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
xline(real_bar_jump_sec(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
legend('Wind jumps','Bar jumps');
for jump = 2:4
   xline(real_wind_jump_sec(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_sec(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
title('Offset');
legend('Bar offset','Wind offset');
ylim([-180 180]);
xlim([0 x_out_bar_offset(end-1)]);
set(gca,'xticklabel',{[]})

subplot(5,1,4)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.','handleVisibility','off')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.','handleVisibility','off')
%add jumps
xline(real_wind_jump_sec(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
xline(real_bar_jump_sec(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
legend('Wind jumps','Bar jumps');
for jump = 2:4
   xline(real_wind_jump_sec(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_sec(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
title('Bump magnitude')
xlim([0 continuous_data.time(end)]);
set(gca,'xticklabel',{[]})

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.','handleVisibility','off')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.','handleVisibility','off')
%add jumps
xline(real_wind_jump_sec(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
xline(real_bar_jump_sec(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
legend('Wind jumps','Bar jumps');
for jump = 2:4
   xline(real_wind_jump_sec(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_sec(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
xlabel('Time (sec)');
title('Bump width');
xlim([0 continuous_data.time(end)]);

saveas(gcf,[path,'\plots\full_experiment.png']);


%% Close-up around the bar jumps

for jump = 1:length(real_bar_jump_frames)
    
   time_zero = continuous_data.time(real_bar_jump_frames(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1200 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_bar_jump_frames(jump)-100:real_bar_jump_frames(jump)+101))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   ax(2) = subplot(8,1,[2 8]);
   time_to_plot = time(real_bar_jump_frames(jump)-100:real_bar_jump_frames(jump)+101);
   phase_to_plot = bump_pos(real_bar_jump_frames(jump)-100:real_bar_jump_frames(jump)+101);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time,bump_pos_to_plot,'linewidth',2)
   hold on
   bar_to_plot = bar_pos(real_bar_jump_frames(jump)-100:real_bar_jump_frames(jump)+101);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'linewidth',2)
   wind_to_plot = motor_pos(real_bar_jump_frames(jump)-100:real_bar_jump_frames(jump)+101);
   [x_out_time,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time,wind_pos_to_plot,'linewidth',2)
   xline(time(real_bar_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(real_bar_jump_frames(jump)-100) time(real_bar_jump_frames(jump)+101)]);
   ylabel('Deg');
   xlabel('Time (sec)');
   legend('Bump estimate','Bar position','Wind position','Bar jump');
   
   saveas(gcf,[path,'\plots\close-up_bar_jump_',num2str(jump),'.png']);
end


%% Close up around the wind jumps

for jump = 1:length(real_wind_jump_frames)
    
   time_zero = continuous_data.time(real_wind_jump_frames(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1200 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_wind_jump_frames(jump)-100:real_wind_jump_frames(jump)+101))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   ax(2) = subplot(8,1,[2 8]);
   time_to_plot = time(real_wind_jump_frames(jump)-100:real_wind_jump_frames(jump)+101);
   phase_to_plot = bump_pos(real_wind_jump_frames(jump)-100:real_wind_jump_frames(jump)+101);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time,bump_pos_to_plot,'linewidth',2)
   hold on
   bar_to_plot = bar_pos(real_wind_jump_frames(jump)-100:real_wind_jump_frames(jump)+101);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'linewidth',2)
   wind_to_plot = motor_pos(real_wind_jump_frames(jump)-100:real_wind_jump_frames(jump)+101);
   [x_out_time,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time,wind_pos_to_plot,'linewidth',2)
   xline(time(real_wind_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(real_wind_jump_frames(jump)-100) time(real_wind_jump_frames(jump)+101)]);
   ylabel('Deg');
   xlabel('Time (sec)');
   legend('Bump estimate','Bar position','Wind position','Wind jump');
   
   saveas(gcf,[path,'\plots\close-up_wind_jump_',num2str(jump),'.png']);
end
