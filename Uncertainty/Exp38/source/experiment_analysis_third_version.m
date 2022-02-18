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

%% Assign fly ID

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
% 
% colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
%     0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
%     0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880;...
%     0 0.4470 0.7410; 0.75, 0.1, 0.75; 0.75, 0.75, 0];
% 
% fly_color = colors_for_plots(fly_ID,:);


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
if fly_ID == 1 
    real_bar_jump_frames(1) = floor(1488.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.5*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.5*length(continuous_data.dff_matrix)/continuous_data.time(end));    
elseif fly_ID == 2  
    real_bar_jump_frames(1) = floor(1488.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.7*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 3
    real_bar_jump_frames(1) = floor(1188.8*length(continuous_data.dff_matrix)/continuous_data.time(end));    
    real_bar_jump_frames(2) = floor(1788.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(3) = floor(2388.5*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(4) = floor(2988.8*length(continuous_data.dff_matrix)/continuous_data.time(end));     
elseif fly_ID == 4
    real_bar_jump_frames(1) = floor(1188.7*length(continuous_data.dff_matrix)/continuous_data.time(end));    
    real_bar_jump_frames(2) = floor(1788.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(3) = floor(2388.7*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(4) = floor(2988.7*length(continuous_data.dff_matrix)/continuous_data.time(end));    
elseif fly_ID == 5
    real_bar_jump_frames(1) = floor(1488.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.6*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 6
    real_bar_jump_frames(1) = floor(1488.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.6*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 7
    real_bar_jump_frames(1) = floor(1189.3*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(1789.3*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(4) = floor(2989.35*length(continuous_data.dff_matrix)/continuous_data.time(end));  
elseif fly_ID == 8
    real_bar_jump_frames(1) = floor(1188.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(1788.9*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(3) = floor(2388.84*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(4) = floor(2988.9*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
elseif fly_ID == 9
    real_bar_jump_frames(1) = floor(1488.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.8*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 10
    real_bar_jump_frames(1) = floor(1188.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(1788.6*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(3) = floor(2388.6*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(4) = floor(2988.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 11
    real_bar_jump_frames(1) = floor(1188.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(1788.6*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(3) = floor(2388.6*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(4) = floor(2988.6*length(continuous_data.dff_matrix)/continuous_data.time(end));   
elseif fly_ID == 12
    real_bar_jump_frames(1) = floor(1488.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.8*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 13
    real_bar_jump_frames(1) = floor(1489*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.9*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.8*length(continuous_data.dff_matrix)/continuous_data.time(end));   
elseif fly_ID == 14
    real_bar_jump_frames(1) = floor(1189.4*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(1789.4*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(3) = floor(2389.4*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
    real_bar_jump_frames(4) = floor(2989.3*length(continuous_data.dff_matrix)/continuous_data.time(end));
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
    real_wind_jump_frames(1) = floor(1188.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.6*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2388.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 2
    real_wind_jump_frames(1) = floor(1188.76*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.88*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2388.76*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.76*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
elseif fly_ID == 3
    real_wind_jump_frames(1) = floor(1489*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2089*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2688.85*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3288.85*length(continuous_data.dff_matrix)/continuous_data.time(end));    
elseif fly_ID == 4
    real_wind_jump_frames(1) = floor(1488.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2088.8*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2688.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3288.7*length(continuous_data.dff_matrix)/continuous_data.time(end));     
elseif fly_ID == 5
    real_wind_jump_frames(1) = floor(1188.7*length(continuous_data.dff_matrix)/continuous_data.time(end));    
    real_wind_jump_frames(2) = floor(1788.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(3) = floor(2388.7*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(4) = floor(2988.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 6
    real_wind_jump_frames(1) = floor(1188.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(3) = floor(2388.4*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 7
    real_wind_jump_frames(1) = floor(1489.39*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2089.37*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2689.22*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3289.44*length(continuous_data.dff_matrix)/continuous_data.time(end));  
elseif fly_ID == 8
    real_wind_jump_frames(1) = floor(1488.95*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2088.95*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2688.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3289*length(continuous_data.dff_matrix)/continuous_data.time(end)); 
elseif fly_ID == 9
    real_wind_jump_frames(1) = floor(1189*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(3) = floor(2388.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 10
    real_wind_jump_frames(1) = floor(1488.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2088.8*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2688.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3288.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 11
    real_wind_jump_frames(1) = floor(1488.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2088.7*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2688.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3288.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 12
    real_wind_jump_frames(1) = floor(1188.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(3) = floor(2388.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.8*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 13
    real_wind_jump_frames(1) = floor(1189*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1789*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(3) = floor(2388.9*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2989*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 14
    real_wind_jump_frames(1) = floor(1489.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2089.5*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2689.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3289.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
end

real_wind_jump_sec = real_wind_jump_frames*continuous_data.time(end)/length(continuous_data.dff_matrix);

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
set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});

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
   
   figure('Position',[100 100 1600 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   colormap(ax(1),flipud(gray));
   pos = get(subplot(8,1,1),'Position');
   %pos2 = get(subplot(12,1,2),'Position');
   h = colorbar('Position', [pos(1)+pos(3)+0.02  0.8054  pos(3)/40  pos(4)+0.055]); 
   
   ax(2) = subplot(8,1,2);
   imagesc(continuous_data.bump_width(:,real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');    
   
   c2 = colormap(ax(2),flipud(bone));
   pos2 = get(subplot(8,1,2),'Position');
   h2 = colorbar('Position', [pos(1)- 0.06  pos2(2)  pos(3)/40  pos(4)+0.055]);
   
   ax(3) = subplot(8,1,[3 6]);
   time_to_plot = time(real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames);
   phase_to_plot = bump_pos(real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time,bump_pos_to_plot,'linewidth',2)
   hold on
   bar_to_plot = bar_pos(real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'linewidth',2)
   wind_to_plot = motor_pos(real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time,wind_pos_to_plot,'linewidth',2)
   xline(time(real_bar_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(real_bar_jump_frames(jump)-floor(10*sec_to_frames)) time(real_bar_jump_frames(jump)+floor(10*sec_to_frames))]);
   ylabel('Deg');
   legend('Bump estimate','Bar position','Wind position','Bar jump');
   set(gca,'xtick',[]);
   
   ax(4) = subplot(8,1,[7:8]);
   bar_offset_to_plot = bar_offset(real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,bar_off_to_plot] = removeWrappedLines(time_to_plot,bar_offset_to_plot);   
   plot(x_out_time,bar_off_to_plot,'linewidth',2,'color',[0.8500 0.3250 0.0980])
   hold on
   wind_offset_to_plot = wind_offset(real_bar_jump_frames(jump)-10*sec_to_frames:real_bar_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,wind_off_to_plot] = removeWrappedLines(time_to_plot,wind_offset_to_plot');   
   plot(x_out_time,wind_off_to_plot,'linewidth',2,'color',[0.9290 0.6940 0.1250])
   xlim([time(real_bar_jump_frames(jump)-floor(10*sec_to_frames)) time(real_bar_jump_frames(jump)+floor(10*sec_to_frames))]);
   xline(time(real_bar_jump_frames(jump)),'k','linestyle','--','linewidth',2)   
   legend('Bar offset','Wind offset','Bar jump');
   xlabel('Time (sec)');ylabel('Offset (deg)');
   ylim([-180 180]);
   
   saveas(gcf,[path,'\plots\closeup_around_bar_jump_',num2str(jump),'.png']);
end


%% Repeat for the wind jumps

for jump = 1:length(real_wind_jump_frames)
    
   time_zero = continuous_data.time(real_wind_jump_frames(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1600 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   colormap(ax(1),flipud(gray));
   pos = get(subplot(8,1,1),'Position');
   %pos2 = get(subplot(12,1,2),'Position');
   h = colorbar('Position', [pos(1)+pos(3)+0.02  0.8054  pos(3)/40  pos(4)+0.055]); 
   
   ax(2) = subplot(8,1,2);
   imagesc(continuous_data.bump_width(:,real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');
         
   c2 = colormap(ax(2),flipud(bone));
   pos2 = get(subplot(8,1,2),'Position');
   h2 = colorbar('Position', [pos(1)- 0.06  pos2(2)  pos(3)/40  pos(4)+0.055]);
   
   ax(3) = subplot(8,1,[3 6]);
   time_to_plot = time(real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames);
   phase_to_plot = bump_pos(real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time,bump_pos_to_plot,'linewidth',2)
   hold on
   bar_to_plot = bar_pos(real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'linewidth',2)
   wind_to_plot = motor_pos(real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time,wind_pos_to_plot,'linewidth',2)
   xline(time(real_wind_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(real_wind_jump_frames(jump)-floor(10*sec_to_frames)) time(real_wind_jump_frames(jump)+floor(10*sec_to_frames))]);
   ylabel('Deg');
   set(gca,'xtick',[]);
   legend('Bump estimate','Bar position','Wind position','Wind jump');
   
   ax(4) = subplot(8,1,[7:8]);
   bar_offset_to_plot = bar_offset(real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,bar_off_to_plot] = removeWrappedLines(time_to_plot,bar_offset_to_plot);   
   plot(x_out_time,bar_off_to_plot,'linewidth',2,'color',[0.8500 0.3250 0.0980])
   hold on
   wind_offset_to_plot = wind_offset(real_wind_jump_frames(jump)-10*sec_to_frames:real_wind_jump_frames(jump)+10*sec_to_frames);
   [x_out_time,wind_off_to_plot] = removeWrappedLines(time_to_plot,wind_offset_to_plot');   
   plot(x_out_time,wind_off_to_plot,'linewidth',2,'color',[0.9290 0.6940 0.1250])
   xlim([time(real_wind_jump_frames(jump)-floor(10*sec_to_frames)) time(real_wind_jump_frames(jump)+floor(10*sec_to_frames))]);
   xline(time(real_wind_jump_frames(jump)),'k','linestyle','--','linewidth',2)   
   legend('Bar offset','Wind offset','Wind jump');
   xlabel('Time (sec)'); ylabel('Offset (deg)');
   ylim([-180 180]);
   
   saveas(gcf,[path,'\plots\closeup_around_wind_jump_',num2str(jump),'.png']);
end


%% Look at a longer around the jump period

for jump = 1:length(real_bar_jump_frames)
    
   time_zero = continuous_data.time(real_bar_jump_frames(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1600 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   colormap(ax(1),flipud(gray));
   pos = get(subplot(8,1,1),'Position');
   %pos2 = get(subplot(12,1,2),'Position');
   h = colorbar('Position', [pos(1)+pos(3)+0.02  0.8054  pos(3)/40  pos(4)+0.055]); 
   
   ax(2) = subplot(8,1,2);
   imagesc(continuous_data.bump_width(:,real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');    
   
   c2 = colormap(ax(2),flipud(bone));
   pos2 = get(subplot(8,1,2),'Position');
   h2 = colorbar('Position', [pos(1)- 0.06  pos2(2)  pos(3)/40  pos(4)+0.055]);
   
   ax(3) = subplot(8,1,[3 6]);
   time_to_plot = time(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
   phase_to_plot = bump_pos(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time,bump_pos_to_plot,'linewidth',2)
   hold on
   bar_to_plot = bar_pos(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'linewidth',2)
   wind_to_plot = motor_pos(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time,wind_pos_to_plot,'linewidth',2)
   xline(time(real_bar_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(real_bar_jump_frames(jump)-floor(120*sec_to_frames)) time(real_bar_jump_frames(jump)+floor(120*sec_to_frames))]);
   ylabel('Deg');
   legend('Bump estimate','Bar position','Wind position','Bar jump');
   set(gca,'xtick',[]);
   
   ax(4) = subplot(8,1,[7:8]);
   bar_offset_to_plot = bar_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,bar_off_to_plot] = removeWrappedLines(time_to_plot,bar_offset_to_plot);   
   plot(x_out_time,bar_off_to_plot,'linewidth',2,'color',[0.8500 0.3250 0.0980])
   hold on
   wind_offset_to_plot = wind_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,wind_off_to_plot] = removeWrappedLines(time_to_plot,wind_offset_to_plot');   
   plot(x_out_time,wind_off_to_plot,'linewidth',2,'color',[0.9290 0.6940 0.1250])
   xlim([time(real_bar_jump_frames(jump)-floor(120*sec_to_frames)) time(real_bar_jump_frames(jump)+floor(120*sec_to_frames))]);
   xline(time(real_bar_jump_frames(jump)),'k','linestyle','--','linewidth',2)   
   legend('Bar offset','Wind offset','Bar jump');
   xlabel('Time (sec)');ylabel('Offset (deg)');
   ylim([-180 180]);
   
   saveas(gcf,[path,'\plots\around_bar_jump_',num2str(jump),'.png']);
end



%% Repeat for wind jumps

for jump = 1:length(real_wind_jump_frames)
    
   time_zero = continuous_data.time(real_wind_jump_frames(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1600 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   colormap(ax(1),flipud(gray));
   pos = get(subplot(8,1,1),'Position');
   %pos2 = get(subplot(12,1,2),'Position');
   h = colorbar('Position', [pos(1)+pos(3)+0.02  0.8054  pos(3)/40  pos(4)+0.055]); 
   
   ax(2) = subplot(8,1,2);
   imagesc(continuous_data.bump_width(:,real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');
         
   c2 = colormap(ax(2),flipud(bone));
   pos2 = get(subplot(8,1,2),'Position');
   h2 = colorbar('Position', [pos(1)- 0.06  pos2(2)  pos(3)/40  pos(4)+0.055]);
   
   ax(3) = subplot(8,1,[3 6]);
   time_to_plot = time(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
   phase_to_plot = bump_pos(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time,bump_pos_to_plot,'linewidth',2)
   hold on
   bar_to_plot = bar_pos(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'linewidth',2)
   wind_to_plot = motor_pos(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time,wind_pos_to_plot,'linewidth',2)
   xline(time(real_wind_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(real_wind_jump_frames(jump)-floor(120*sec_to_frames)) time(real_wind_jump_frames(jump)+floor(120*sec_to_frames))]);
   ylabel('Deg');
   set(gca,'xtick',[]);
   legend('Bump estimate','Bar position','Wind position','Wind jump');
   
   ax(4) = subplot(8,1,[7:8]);
   bar_offset_to_plot = bar_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,bar_off_to_plot] = removeWrappedLines(time_to_plot,bar_offset_to_plot);   
   plot(x_out_time,bar_off_to_plot,'linewidth',2,'color',[0.8500 0.3250 0.0980])
   hold on
   wind_offset_to_plot = wind_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
   [x_out_time,wind_off_to_plot] = removeWrappedLines(time_to_plot,wind_offset_to_plot');   
   plot(x_out_time,wind_off_to_plot,'linewidth',2,'color',[0.9290 0.6940 0.1250])
   xlim([time(real_wind_jump_frames(jump)-floor(120*sec_to_frames)) time(real_wind_jump_frames(jump)+floor(120*sec_to_frames))]);
   xline(time(real_wind_jump_frames(jump)),'k','linestyle','--','linewidth',2)   
   legend('Bar offset','Wind offset','Wind jump');
   xlabel('Time (sec)'); ylabel('Offset (deg)');
   ylim([-180 180]);
   
   saveas(gcf,[path,'\plots\around_wind_jump_',num2str(jump),'.png']);
end


%% Threshold movement data and goodness of fit

% Uncomment below to see the total movement distribution
figure,
histogram(continuous_data.total_mvt_ds)
xlabel('Total movement (deg/s)');
ylabel('Counts');
mvt_thresh = 25;
moving = continuous_data.total_mvt_ds>mvt_thresh;
good_fit = continuous_data.adj_rs>0.5;

%I'm currently not using this

%% Offset and mean bump parameters for initial bouts (last three min of each cue)

if configuration == 1
    
    initial_bar_period = ceil(panels_change_frames(1)-180*sec_to_frames):panels_change_frames(2)-1;
    initial_bar_moving = continuous_data.total_mvt_ds(initial_bar_period) > 25;
    initial_bar_good_fit = continuous_data.adj_rs(initial_bar_period) > 0.5;    
    initial_bar_offset = bar_offset(initial_bar_period);
    initial_bar_offset = initial_bar_offset(initial_bar_moving & initial_bar_good_fit);
    initial_bar_bm = continuous_data.bump_magnitude(initial_bar_period);
    initial_bar_bm = initial_bar_bm(initial_bar_moving & initial_bar_good_fit);
    initial_bar_bw = continuous_data.bump_width(initial_bar_period);
    initial_bar_bw = initial_bar_bw(initial_bar_moving & initial_bar_good_fit);
    
    initial_wind_period = ceil(panels_change_frames(2)-180*sec_to_frames):panels_change_frames(2)-1;
    initial_wind_moving = continuous_data.total_mvt_ds(initial_wind_period) > 25;
    initial_wind_good_fit = continuous_data.adj_rs(initial_wind_period) > 0.5;
    initial_wind_offset = wind_offset(initial_wind_period);
    initial_wind_offset = initial_wind_offset(initial_wind_moving & initial_wind_good_fit);
    initial_wind_bm = continuous_data.bump_magnitude(initial_wind_period);
    initial_wind_bm = initial_wind_bm(initial_wind_moving & initial_wind_good_fit);
    initial_wind_bw = continuous_data.bump_width(initial_wind_period);
    initial_wind_bw = initial_wind_bw(initial_wind_moving & initial_wind_good_fit);
    
     %Offset
    figure,
    subplot(1,2,1)
    polarhistogram(deg2rad(initial_wind_offset),15)
    initial_wind_offset_var = circ_std(deg2rad(initial_wind_offset),[],[],2);
    title({'Initial wind offset';['Offset var = ',num2str(round(initial_wind_offset_var,2))]});
    set(gca,'ThetaZeroLocation','top');
    %set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    subplot(1,2,2)
    polarhistogram(deg2rad(initial_bar_offset),15)
    initial_bar_offset_var = circ_std(deg2rad(initial_bar_offset));
    title({'Initial bar offset';['Bar var = ',num2str(round(initial_bar_offset_var,2))]});
    set(gca,'ThetaZeroLocation','top');
    %set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    saveas(gcf,[path,'\plots\initial_offsets.png']);
    
    
    %Bump parameters
    figure,
    subplot(1,2,1)
    plot([mean(initial_wind_bm), mean(initial_bar_bm)],'-o')
    title('Bump magnitude');
    ylim([0 3]);
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Wind','Bar'});
    
    subplot(1,2,2)
    plot([mean(initial_wind_bw), mean(initial_bar_bw)],'-o')
    title('Bump width');
    ylim([0 3]);
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Wind','Bar'});
    
    saveas(gcf,[path,'\plots\initial_bump_pars.png']);
    
else
    
    initial_bar_period = ceil(wind_change_frames(2)-180*sec_to_frames):1:wind_change_frames(2)-1;
    initial_bar_moving = continuous_data.total_mvt_ds(initial_bar_period) > 25;
    initial_bar_good_fit = continuous_data.adj_rs(initial_bar_period) > 0.5;    
    initial_bar_offset = bar_offset(initial_bar_period);
    initial_bar_offset = initial_bar_offset(initial_bar_moving & initial_bar_good_fit);
    initial_bar_bm = continuous_data.bump_magnitude(initial_bar_period);
    initial_bar_bm = initial_bar_bm(initial_bar_moving & initial_bar_good_fit);
    initial_bar_bw = continuous_data.bump_width(initial_bar_period);
    initial_bar_bw = initial_bar_bw(initial_bar_moving & initial_bar_good_fit);
    
    initial_wind_period = ceil(wind_change_frames(1)-180*sec_to_frames):wind_change_frames(1)-1;
    initial_wind_moving = continuous_data.total_mvt_ds(initial_wind_period) > 25;
    initial_wind_good_fit = continuous_data.adj_rs(initial_wind_period) > 0.5;
    initial_wind_offset = wind_offset(initial_wind_period);
    initial_wind_offset = initial_wind_offset(initial_wind_moving & initial_wind_good_fit);
    initial_wind_bm = continuous_data.bump_magnitude(initial_wind_period);
    initial_wind_bm = initial_wind_bm(initial_wind_moving & initial_wind_good_fit);
    initial_wind_bw = continuous_data.bump_width(initial_wind_period);
    initial_wind_bw = initial_wind_bw(initial_wind_moving & initial_wind_good_fit);
    
    %Offset
    figure,
    subplot(1,2,1)
    polarhistogram(deg2rad(initial_wind_offset),15)
    initial_wind_offset_var = circ_std(deg2rad(initial_wind_offset),[],[],2);
    title({'Initial wind offset';['Offset var = ',num2str(round(initial_wind_offset_var,2))]});
    set(gca,'ThetaZeroLocation','top');
    %set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    subplot(1,2,2)
    polarhistogram(deg2rad(initial_bar_offset),15)
    initial_bar_offset_var = circ_std(deg2rad(initial_bar_offset));
    title({'Initial bar offset';['Bar var = ',num2str(round(initial_bar_offset_var,2))]});
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Rticklabels',{});
    
    saveas(gcf,[path,'\plots\initial_offsets.png']);
    
    
    %Bump parameters
    figure,
    subplot(1,2,1)
    plot([mean(initial_wind_bm), mean(initial_bar_bm)],'-o')
    title('Bump magnitude');
    ylim([0 3]);
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Wind','Bar'});
    
    subplot(1,2,2)
    plot([mean(initial_wind_bw), mean(initial_bar_bw)],'-o')
    title('Bump width');
    ylim([0 3]);
    xlim([0 3]);
    xticks([1 2]);
    xticklabels({'Wind','Bar'});
    
    saveas(gcf,[path,'\plots\initial_bump_pars.png']);
        
end


%% Bump par before and after bar jumps

for jump = 1:length(real_bar_jump_frames)
    
    %short timescale
    pre_bar_jump_period = ceil(real_bar_jump_frames(jump)-1.5*sec_to_frames):real_bar_jump_frames(jump)-1;
    pre_bar_jump_moving = continuous_data.total_mvt_ds(pre_bar_jump_period) > 25;
    pre_bar_jump_good_fit = continuous_data.adj_rs(pre_bar_jump_period) > 0.5;
    short_bm_pre_bar_jump = continuous_data.bump_magnitude(pre_bar_jump_period);
    mean_short_bm_around_bar_jump(jump,1) = nanmean(short_bm_pre_bar_jump(pre_bar_jump_moving & pre_bar_jump_good_fit));
    
    post_bar_jump_period = real_bar_jump_frames(jump)+1:ceil(real_bar_jump_frames(jump)+1.5*sec_to_frames);
    post_bar_jump_moving = continuous_data.total_mvt_ds(post_bar_jump_period) > 25;
    post_bar_jump_good_fit = continuous_data.adj_rs(post_bar_jump_period) > 0.5;
    short_bm_post_bar_jump = continuous_data.bump_magnitude(post_bar_jump_period);
    mean_short_bm_around_bar_jump(jump,2) = nanmean(short_bm_post_bar_jump(post_bar_jump_moving & post_bar_jump_good_fit));
    
    short_bw_pre_bar_jump = continuous_data.bump_width(pre_bar_jump_period);
    mean_short_bw_around_bar_jump(jump,1) = nanmean(short_bw_pre_bar_jump(pre_bar_jump_moving & pre_bar_jump_good_fit));
    
    short_bw_post_bar_jump = continuous_data.bump_width(post_bar_jump_period);
    mean_short_bw_around_bar_jump(jump,2) = nanmean(short_bw_post_bar_jump(post_bar_jump_moving & post_bar_jump_good_fit));
    
    
    %long timescale
    long_pre_bar_jump_period = ceil(real_bar_jump_frames(jump)-120*sec_to_frames):real_bar_jump_frames(jump)-1;
    long_pre_bar_jump_moving = continuous_data.total_mvt_ds(long_pre_bar_jump_period) > 25;
    long_pre_bar_jump_good_fit = continuous_data.adj_rs(long_pre_bar_jump_period) > 0.5;
    long_bm_pre_bar_jump = continuous_data.bump_magnitude(long_pre_bar_jump_period);
    mean_long_bm_around_bar_jump(jump,1) = nanmean(long_bm_pre_bar_jump(long_pre_bar_jump_moving & long_pre_bar_jump_good_fit));
    
    long_post_bar_jump_period = real_bar_jump_frames(jump)+1:ceil(real_bar_jump_frames(jump)+120*sec_to_frames);
    long_post_bar_jump_moving = continuous_data.total_mvt_ds(long_post_bar_jump_period) > 25;
    long_post_bar_jump_good_fit = continuous_data.adj_rs(long_post_bar_jump_period) > 0.5;
    long_bm_post_bar_jump = continuous_data.bump_magnitude(long_post_bar_jump_period);
    mean_long_bm_around_bar_jump(jump,2) = nanmean(long_bm_post_bar_jump(long_post_bar_jump_moving & long_post_bar_jump_good_fit));
    
    long_bw_pre_bar_jump = continuous_data.bump_width(long_pre_bar_jump_period);
    mean_long_bw_around_bar_jump(jump,1) = nanmean(long_bw_pre_bar_jump(long_pre_bar_jump_moving & long_pre_bar_jump_good_fit));
    
    long_bw_post_bar_jump = continuous_data.bump_width(long_post_bar_jump_period);
    mean_long_bw_around_bar_jump(jump,2) = nanmean(long_bw_post_bar_jump(long_post_bar_jump_moving & long_post_bar_jump_good_fit));
    
end

figure('Position',[100 100 800 1000]),
subplot(2,2,1)
plot(mean_short_bm_around_bar_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_short_bm_around_bar_jump),'k','linewidth',2)
ylabel('Bump magnitude');
title('Short timescale');
xlim([0 3]);
ylim([0 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,2)
plot(mean_long_bm_around_bar_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_long_bm_around_bar_jump),'k','linewidth',2)
title('Long timescale');
xlim([0 3]);
ylim([0 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,3)
plot(mean_short_bw_around_bar_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_short_bw_around_bar_jump),'k','linewidth',2)
ylabel('Bump width');
xlim([0 3]);
ylim([0 3.5]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,4)
plot(mean_long_bw_around_bar_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_long_bw_around_bar_jump),'k','linewidth',2)
xlim([0 3]);
ylim([0 3.5]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

suptitle('Bar jumps');

saveas(gcf,[path,'\plots\bump_pars_around_bar_jumps.png']);

%% Repeat for wind

for jump = 1:length(real_wind_jump_frames)
    
    %short timescale
    pre_wind_jump_period = ceil(real_wind_jump_frames(jump)-1.5*sec_to_frames):real_wind_jump_frames(jump)-1;
    pre_wind_jump_moving = continuous_data.total_mvt_ds(pre_wind_jump_period) > 25;
    pre_wind_jump_good_fit = continuous_data.adj_rs(pre_wind_jump_period) > 0.5;
    short_bm_pre_wind_jump = continuous_data.bump_magnitude(pre_wind_jump_period);
    mean_short_bm_around_wind_jump(jump,1) = nanmean(short_bm_pre_wind_jump(pre_wind_jump_moving & pre_wind_jump_good_fit));
    
    post_wind_jump_period = real_wind_jump_frames(jump)+1:ceil(real_wind_jump_frames(jump)+1.5*sec_to_frames);
    post_wind_jump_moving = continuous_data.total_mvt_ds(post_wind_jump_period) > 25;
    post_wind_jump_good_fit = continuous_data.adj_rs(post_wind_jump_period) > 0.5;
    short_bm_post_wind_jump = continuous_data.bump_magnitude(post_wind_jump_period);
    mean_short_bm_around_wind_jump(jump,2) = nanmean(short_bm_post_wind_jump(post_wind_jump_moving & post_wind_jump_good_fit));
    
    short_bw_pre_wind_jump = continuous_data.bump_width(pre_wind_jump_period);
    mean_short_bw_around_wind_jump(jump,1) = nanmean(short_bw_pre_wind_jump(pre_wind_jump_moving & pre_wind_jump_good_fit));
    
    short_bw_post_wind_jump = continuous_data.bump_width(post_wind_jump_period);
    mean_short_bw_around_wind_jump(jump,2) = nanmean(short_bw_post_wind_jump(post_wind_jump_moving & post_wind_jump_good_fit));
    
    
    %long timescale
    long_pre_wind_jump_period = ceil(real_wind_jump_frames(jump)-120*sec_to_frames):real_wind_jump_frames(jump)-1;
    long_pre_wind_jump_moving = continuous_data.total_mvt_ds(long_pre_wind_jump_period) > 25;
    long_pre_wind_jump_good_fit = continuous_data.adj_rs(long_pre_wind_jump_period) > 0.5;
    long_bm_pre_wind_jump = continuous_data.bump_magnitude(long_pre_wind_jump_period);
    mean_long_bm_around_wind_jump(jump,1) = nanmean(long_bm_pre_wind_jump(long_pre_wind_jump_moving & long_pre_wind_jump_good_fit));
    
    long_post_wind_jump_period = real_wind_jump_frames(jump)+1:ceil(real_wind_jump_frames(jump)+120*sec_to_frames);
    long_post_wind_jump_moving = continuous_data.total_mvt_ds(long_post_wind_jump_period) > 25;
    long_post_wind_jump_good_fit = continuous_data.adj_rs(long_post_wind_jump_period) > 0.5;
    long_bm_post_wind_jump = continuous_data.bump_magnitude(long_post_wind_jump_period);
    mean_long_bm_around_wind_jump(jump,2) = nanmean(long_bm_post_wind_jump(long_post_wind_jump_moving & long_post_wind_jump_good_fit));
    
    long_bw_pre_wind_jump = continuous_data.bump_width(long_pre_wind_jump_period);
    mean_long_bw_around_wind_jump(jump,1) = nanmean(long_bw_pre_wind_jump(long_pre_wind_jump_moving & long_pre_wind_jump_good_fit));
    
    long_bw_post_wind_jump = continuous_data.bump_width(long_post_wind_jump_period);
    mean_long_bw_around_wind_jump(jump,2) = nanmean(long_bw_post_wind_jump(long_post_wind_jump_moving & long_post_wind_jump_good_fit));
    
end

figure('Position',[100 100 800 1000]),
subplot(2,2,1)
plot(mean_short_bm_around_wind_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_short_bm_around_wind_jump),'k','linewidth',2)
ylabel('Bump magnitude');
title('Short timescale');
xlim([0 3]);
ylim([0 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,2)
plot(mean_long_bm_around_wind_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_long_bm_around_wind_jump),'k','linewidth',2)
title('Long timescale');
xlim([0 3]);
ylim([0 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,3)
plot(mean_short_bw_around_wind_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_short_bw_around_wind_jump),'k','linewidth',2)
ylabel('Bump width');
xlim([0 3]);
ylim([0 3.5]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,4)
plot(mean_long_bw_around_wind_jump','color',[.5 .5 .5])
hold on
plot(nanmean(mean_long_bw_around_wind_jump),'k','linewidth',2)
xlim([0 3]);
ylim([0 3.5]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

suptitle('Wind jumps');


saveas(gcf,[path,'\plots\bump_pars_around_wind_jumps.png']);

%% Raster plots of bump parameters in the short time scale

%I won't threshold by movement of fit since otherwise I can't plot them all
%together

%Bar jumps
for jump = 1:length(real_bar_jump_frames)
    short_bm_bar_jump(jump,:) = continuous_data.bump_magnitude(real_bar_jump_frames(jump)-2*sec_to_frames:real_bar_jump_frames(jump)+2*sec_to_frames);
    short_bw_bar_jump(jump,:) = continuous_data.bump_width(real_bar_jump_frames(jump)-2*sec_to_frames:real_bar_jump_frames(jump)+2*sec_to_frames);    
end
zscored_short_bm_bar_jump = zscore(short_bm_bar_jump,[],2);
zscored_short_bw_bar_jump = zscore(short_bw_bar_jump,[],2);

%Wind jumps
for jump = 1:length(real_wind_jump_frames)
    short_bm_wind_jump(jump,:) = continuous_data.bump_magnitude(real_wind_jump_frames(jump)-2*sec_to_frames:real_wind_jump_frames(jump)+2*sec_to_frames);
    short_bw_wind_jump(jump,:) = continuous_data.bump_width(real_wind_jump_frames(jump)-2*sec_to_frames:real_wind_jump_frames(jump)+2*sec_to_frames);    
end
zscored_short_bm_wind_jump = zscore(short_bm_wind_jump,[],2);
zscored_short_bw_wind_jump = zscore(short_bw_wind_jump,[],2);


%Plot bump magnitude
figure('Position',[100 100 1600 600]),
subplot(2,1,1)
imagesc(zscored_short_bm_bar_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');
frames = [5 10 15 20 25 30 35]; 
xticks(frames);
timestamps = [-37/2*frames_to_sec:0.11:37/2*frames_to_sec];
xticklabels(num2cell(round(timestamps(frames),1)));

subplot(2,1,2)
imagesc(zscored_short_bm_wind_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');
frames = [5 10 15 20 25 30 35]; 
xticks(frames);
timestamps = [-37/2*frames_to_sec:0.11:37/2*frames_to_sec];
xticklabels(num2cell(round(timestamps(frames),1)));

suptitle('Bump magnitude');

saveas(gcf,[path,'\plots\short_raster_plots_bump_mag.png']);



%Plot bump width
figure('Position',[100 100 1600 600]),
subplot(2,1,1)
imagesc(zscored_short_bw_bar_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(bone))
title('Bar jumps');
frames = [5 10 15 20 25 30 35]; 
xticks(frames);
timestamps = [-37/2*frames_to_sec:0.11:37/2*frames_to_sec];
xticklabels(num2cell(round(timestamps(frames),1)));

subplot(2,1,2)
imagesc(zscored_short_bw_wind_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(bone))
title('Wind jumps');
frames = [5 10 15 20 25 30 35]; 
xticks(frames);
timestamps = [-37/2*frames_to_sec:0.11:37/2*frames_to_sec];
xticklabels(num2cell(round(timestamps(frames),1)));

suptitle('Bump width');

saveas(gcf,[path,'\plots\short_raster_plots_bump_width.png']);

%% Repeat for long time scale

%Bar jumps
for jump = 1:length(real_bar_jump_frames)
    long_bm_bar_jump(jump,:) = continuous_data.bump_magnitude(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
    long_bw_bar_jump(jump,:) = continuous_data.bump_width(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);    
end
zscored_long_bm_bar_jump = zscore(long_bm_bar_jump,[],2);
zscored_long_bw_bar_jump = zscore(long_bw_bar_jump,[],2);

%Wind jumps
for jump = 1:length(real_wind_jump_frames)
    long_bm_wind_jump(jump,:) = continuous_data.bump_magnitude(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
    long_bw_wind_jump(jump,:) = continuous_data.bump_width(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);    
end
zscored_long_bm_wind_jump = zscore(long_bm_wind_jump,[],2);
zscored_long_bw_wind_jump = zscore(long_bw_wind_jump,[],2);


%Plot bump magnitude
figure('Position',[100 100 1600 600]),
subplot(2,1,1)
imagesc(zscored_long_bm_bar_jump)
hold on
xline(120*sec_to_frames,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_long_bm_wind_jump)
hold on
xline(120*sec_to_frames,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump magnitude');

saveas(gcf,[path,'\plots\long_raster_plots_bump_mag.png']);



%Plot bump width
figure('Position',[100 100 1600 600]),
subplot(2,1,1)
imagesc(zscored_long_bw_bar_jump)
hold on
xline(120*sec_to_frames,'r','linewidth',2)
colormap(flipud(bone))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_long_bw_wind_jump)
hold on
xline(120*sec_to_frames,'r','linewidth',2)
colormap(flipud(bone))
title('Wind jumps');

suptitle('Bump width');

saveas(gcf,[path,'\plots\long_raster_plots_bump_width.png']);

%% Averaged bump parameters around jumps

figure('Position',[100 100 1400 1000]),

%bump magnitude around bar jumps
subplot(2,2,1)
plot(mean(short_bm_bar_jump),'k','linewidth',2)
hold on
line([length(short_bm_bar_jump)/2 length(short_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(short_bm_bar_jump)]);
xticks([0 length(short_bm_bar_jump)/4 length(short_bm_bar_jump)/2 length(short_bm_bar_jump)*(3/4) length(short_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
ylabel('Bump mangitude');
title('Bar jumps');

%bump magnitude around wind jumps
subplot(2,2,2)
plot(mean(short_bm_wind_jump),'k','linewidth',2)
hold on
line([length(short_bm_wind_jump)/2 length(short_bm_wind_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(short_bm_bar_jump)]);
xticks([0 length(short_bm_bar_jump)/4 length(short_bm_bar_jump)/2 length(short_bm_bar_jump)*(3/4) length(short_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
title('Wind jumps');

%bump width around bar jumps
subplot(2,2,3)
plot(mean(short_bw_bar_jump),'k','linewidth',2)
hold on
line([length(short_bm_wind_jump)/2 length(short_bm_wind_jump)/2],[0 4],'color','r','linewidth',3);
xlim([0 length(short_bm_bar_jump)]);
xticks([0 length(short_bm_bar_jump)/4 length(short_bm_bar_jump)/2 length(short_bm_bar_jump)*(3/4) length(short_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
xlabel('Time around jumps (s)');
ylabel('Bump width');

%bump width around wind jumps
subplot(2,2,4)
plot(mean(short_bw_wind_jump),'k','linewidth',2)
hold on
line([length(short_bm_wind_jump)/2 length(short_bm_wind_jump)/2],[0 4],'color','r','linewidth',3);
xlim([0 length(short_bm_bar_jump)]);
xticks([0 length(short_bm_bar_jump)/4 length(short_bm_bar_jump)/2 length(short_bm_bar_jump)*(3/4) length(short_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
xlabel('Time around jumps (s)');

saveas(gcf,[path,'\plots\average_bump_pars_aj.png']);

%% Looking at a longer period of time


figure('Position',[100 100 1400 1000]),

%bump magnitude around bar jumps
subplot(2,2,1)
plot(mean(long_bm_bar_jump),'k','linewidth',2)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
ylabel('Bump mangitude');
title('Bar jumps');

%bump magnitude around wind jumps
subplot(2,2,2)
plot(mean(long_bm_wind_jump),'k','linewidth',2)
hold on
line([length(long_bm_wind_jump)/2 length(long_bm_wind_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
title('Wind jumps');

%bump width around bar jumps
subplot(2,2,3)
plot(mean(long_bw_bar_jump),'k','linewidth',2)
hold on
line([length(long_bm_wind_jump)/2 length(long_bm_wind_jump)/2],[0 4],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
xlabel('Time around jumps (s)');
ylabel('Bump width');

%bump width around wind jumps
subplot(2,2,4)
plot(mean(long_bw_wind_jump),'k','linewidth',2)
hold on
line([length(long_bm_wind_jump)/2 length(long_bm_wind_jump)/2],[0 4],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
xlabel('Time around jumps (s)');

saveas(gcf,[path,'\plots\long_average_bump_pars_aj.png']);

%% Average total movement around cue jumps

%Bar jumps
for jump = 1:length(real_bar_jump_frames)
    long_total_mvt_bar_jump(jump,:) = continuous_data.total_mvt_ds(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames);
    long_rot_speed_bar_jump(jump,:) = abs(continuous_data.vel_yaw_ds(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames));    
end

%Wind jumps
for jump = 1:length(real_wind_jump_frames)
    long_total_mvt_wind_jump(jump,:) = continuous_data.total_mvt_ds(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames);
    long_rot_speed_wind_jump(jump,:) = abs(continuous_data.vel_yaw_ds(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames));
end

figure('Position',[100 100 1400 1000]),

%total mvt around bar jumps
subplot(2,2,1)
plot(mean(long_total_mvt_bar_jump),'k','linewidth',2)
hold on
line([length(long_total_mvt_bar_jump)/2 length(long_total_mvt_bar_jump)/2],[0 300],'color','r','linewidth',3);
xlim([0 length(long_total_mvt_bar_jump)]);
xticks([0 length(long_total_mvt_bar_jump)/4 length(long_total_mvt_bar_jump)/2 length(long_total_mvt_bar_jump)*(3/4) length(long_total_mvt_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
ylabel('Total movement');
title('Bar jumps');

%total mvt around wind jumps
subplot(2,2,2)
plot(mean(long_total_mvt_wind_jump),'k','linewidth',2)
hold on
line([length(long_total_mvt_wind_jump)/2 length(long_total_mvt_wind_jump)/2],[0 300],'color','r','linewidth',3);
xlim([0 length(long_total_mvt_bar_jump)]);
xticks([0 length(long_total_mvt_bar_jump)/4 length(long_total_mvt_bar_jump)/2 length(long_total_mvt_bar_jump)*(3/4) length(long_total_mvt_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
title('Wind jumps');

%rot speed around bar jumps
subplot(2,2,3)
plot(mean(long_rot_speed_bar_jump),'k','linewidth',2)
hold on
line([length(long_rot_speed_bar_jump)/2 length(long_rot_speed_bar_jump)/2],[0 60],'color','r','linewidth',3);
xlim([0 length(long_rot_speed_bar_jump)]);
xticks([0 length(long_rot_speed_bar_jump)/4 length(long_rot_speed_bar_jump)/2 length(long_rot_speed_bar_jump)*(3/4) length(long_rot_speed_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
ylabel('Total movement');
title('Bar jumps');

%rot speed around wind jumps
subplot(2,2,4)
plot(mean(long_rot_speed_wind_jump),'k','linewidth',2)
hold on
line([length(long_rot_speed_wind_jump)/2 length(long_rot_speed_wind_jump)/2],[0 60],'color','r','linewidth',3);
xlim([0 length(long_rot_speed_bar_jump)]);
xticks([0 length(long_rot_speed_bar_jump)/4 length(long_rot_speed_bar_jump)/2 length(long_rot_speed_bar_jump)*(3/4) length(long_rot_speed_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
title('Wind jumps');

saveas(gcf,[path,'\plots\long_average_total_mvt_aj.png']);

%% Offset mean with respect to both stimuli pre and post-jump

% Compute means around jumps

%bump offsets
bar_offset = circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle));
wind_offset = circ_dist(continuous_data.bump_pos,continuous_data.motor_pos);

%behavioral offsets
bar_heading_offset = circ_dist(continuous_data.heading,deg2rad(continuous_data.panel_angle));
wind_heading_offset = circ_dist(continuous_data.heading,continuous_data.motor_pos');
wind_heading_offset = wind_heading_offset';
                
%For bar jumps
for jump = 1:length(real_bar_jump_frames)
    long_bar_offset_mean_around_bar_jump(jump,1) = circ_mean(bar_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1));
    long_bar_offset_mean_around_bar_jump(jump,2) = circ_mean(bar_offset(real_bar_jump_frames(jump)+40*sec_to_frames:real_bar_jump_frames(jump)+160*sec_to_frames));
    long_wind_offset_mean_around_bar_jump(jump,1) = circ_mean(wind_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1),[],2);
    long_wind_offset_mean_around_bar_jump(jump,2) = circ_mean(wind_offset(real_bar_jump_frames(jump)+40*sec_to_frames:real_bar_jump_frames(jump)+160*sec_to_frames),[],2);
    
    long_bar_heading_offset_mean_around_bar_jump(jump,1) = circ_mean(bar_heading_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1));
    long_bar_heading_offset_mean_around_bar_jump(jump,2) = circ_mean(bar_heading_offset(real_bar_jump_frames(jump)+40*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames));
    long_wind_heading_offset_mean_around_bar_jump(jump,1) = circ_mean(wind_heading_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1),[],2);
    long_wind_heading_offset_mean_around_bar_jump(jump,2) = circ_mean(wind_heading_offset(real_bar_jump_frames(jump)+40*sec_to_frames:real_bar_jump_frames(jump)+160*sec_to_frames),[],2);
end

%Repeat for wind jumps
for jump = 1:length(real_wind_jump_frames)
    long_bar_offset_mean_around_wind_jump(jump,1) = circ_mean(bar_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1));
    long_bar_offset_mean_around_wind_jump(jump,2) = circ_mean(bar_offset(real_wind_jump_frames(jump)+40*sec_to_frames:real_wind_jump_frames(jump)+160*sec_to_frames));
    long_wind_offset_mean_around_wind_jump(jump,1) = circ_mean(wind_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1),[],2);
    long_wind_offset_mean_around_wind_jump(jump,2) = circ_mean(wind_offset(real_wind_jump_frames(jump)+40*sec_to_frames:real_wind_jump_frames(jump)+160*sec_to_frames),[],2);
    
    long_bar_heading_offset_mean_around_wind_jump(jump,1) = circ_mean(bar_heading_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1));
    long_bar_heading_offset_mean_around_wind_jump(jump,2) = circ_mean(bar_heading_offset(real_wind_jump_frames(jump)+40*sec_to_frames:real_wind_jump_frames(jump)+160*sec_to_frames));
    long_wind_heading_offset_mean_around_wind_jump(jump,1) = circ_mean(wind_heading_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1),[],2);
    long_wind_heading_offset_mean_around_wind_jump(jump,2) = circ_mean(wind_heading_offset(real_wind_jump_frames(jump)+40*sec_to_frames:real_wind_jump_frames(jump)+160*sec_to_frames),[],2);
end


% Get all_PI_values

%1) For EPG PI
%Bar jumps
mean_bar_offset_diff_bj = abs(circ_dist(long_bar_offset_mean_around_bar_jump(:,2), long_bar_offset_mean_around_bar_jump(:,1)));
mean_wind_offset_diff_bj = abs(circ_dist(long_wind_offset_mean_around_bar_jump(:,2), long_wind_offset_mean_around_bar_jump(:,1)));
%Plot differences in mean offset
figure('Position',[100 100 800 800]),
plot([rad2deg(mean_bar_offset_diff_bj),rad2deg(mean_wind_offset_diff_bj)]','-ro')

%Wind jumps
mean_bar_offset_diff_wj = abs(circ_dist(long_bar_offset_mean_around_wind_jump(:,2), long_bar_offset_mean_around_wind_jump(:,1)));
mean_wind_offset_diff_wj = abs(circ_dist(long_wind_offset_mean_around_wind_jump(:,2), long_wind_offset_mean_around_wind_jump(:,1)));
%add to plot
hold on
plot([rad2deg(mean_bar_offset_diff_wj),rad2deg(mean_wind_offset_diff_wj)]','-ko')
xlim([0 3]);
ylim([0 180]);
xticks([1 2]);
xticklabels({'Bar offset diff','Wind offset diff'});
ylabel('Post-pre mean offset difference (deg)','fontsize',12);
saveas(gcf,[path,'\plots\change_in_mean_offset.png']);
%compute preference index
mean_bar_offset_diff = [mean_bar_offset_diff_bj;mean_bar_offset_diff_wj];
mean_wind_offset_diff = [mean_wind_offset_diff_bj;mean_wind_offset_diff_wj];
PI = (mean_bar_offset_diff-mean_wind_offset_diff)./(mean_bar_offset_diff+mean_wind_offset_diff);


%1) For heading PI
%Bar jumps
mean_bar_heading_offset_diff_bj = abs(circ_dist(long_bar_heading_offset_mean_around_bar_jump(:,2),long_bar_heading_offset_mean_around_bar_jump(:,1)));
mean_wind_heading_offset_diff_bj = abs(circ_dist(long_wind_heading_offset_mean_around_bar_jump(:,2), long_wind_heading_offset_mean_around_bar_jump(:,1)));
%Wind jumps
mean_bar_heading_offset_diff_wj = abs(circ_dist(long_bar_heading_offset_mean_around_wind_jump(:,2), long_bar_heading_offset_mean_around_wind_jump(:,1)));
mean_wind_heading_offset_diff_wj = abs(circ_dist(long_wind_heading_offset_mean_around_wind_jump(:,2), long_wind_heading_offset_mean_around_wind_jump(:,1)));
%compute heading preference index
mean_bar_heading_offset_diff = [mean_bar_heading_offset_diff_bj;mean_bar_heading_offset_diff_wj];
mean_wind_heading_offset_diff = [mean_wind_heading_offset_diff_bj;mean_wind_heading_offset_diff_wj];
PI_heading = (mean_bar_heading_offset_diff-mean_wind_heading_offset_diff)./(mean_bar_heading_offset_diff+mean_wind_heading_offset_diff);

% Compute the stickiness index
%Bar moves
mean_move_cue_offset_diff_bj = abs( circ_dist(long_bar_offset_mean_around_bar_jump(:,2), long_bar_offset_mean_around_bar_jump(:,1)));
mean_stay_cue_offset_diff_bj = abs(circ_dist(long_wind_offset_mean_around_bar_jump(:,2), long_wind_offset_mean_around_bar_jump(:,1)));
%Wind moves
mean_move_cue_offset_diff_wj = abs(circ_dist(long_wind_offset_mean_around_wind_jump(:,2), long_wind_offset_mean_around_wind_jump(:,1)));
mean_stay_cue_offset_diff_wj = abs(circ_dist(long_bar_offset_mean_around_wind_jump(:,2), long_bar_offset_mean_around_wind_jump(:,1)));
%compute stickiness index
mean_move_cue_offset_diff = [mean_move_cue_offset_diff_bj;mean_move_cue_offset_diff_wj];
mean_stay_cue_offset_diff = [mean_stay_cue_offset_diff_bj;mean_stay_cue_offset_diff_wj];
SI = (mean_move_cue_offset_diff - mean_stay_cue_offset_diff)./(mean_move_cue_offset_diff + mean_stay_cue_offset_diff);

%Plot preference indices                
figure,
boxplot([PI,PI_heading],'color','k')
hold on
yline(0);
scatter(repmat([1:2],8,1),[PI,PI_heading],[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Preference index');
ylim([-1 1]);
xticks([1 2]);
xticklabels({'Bump','Behavior'});
saveas(gcf,[path,'\plots\preference_indices.png']);

%Compute the mean preference index (proxy for which cue the bump and heading are focusing on)
    %the closer to -1, the more the bar is preferred
    %the closer to 1, the more the wind is preferred
pref_index = mean(PI); %preference index as the mean of ratios
pref_index_2 = mean(mean_bar_offset_diff-mean_wind_offset_diff)/mean(mean_bar_offset_diff+mean_wind_offset_diff); %preference index as the ratio of means 
heading_pref_index = mean(PI_heading);
heading_pref_index_2 = mean(mean_bar_heading_offset_diff-mean_wind_heading_offset_diff)/mean(mean_bar_heading_offset_diff+mean_wind_heading_offset_diff);;

%Plot and save stickiness index
figure,
boxplot(SI,'color','k')
hold on
yline(0);
scatter(repelem(1,8,1),SI,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Stickiness index');
ylim([-1 1]);
saveas(gcf,[path,'\plots\stickiness_index.png']);

stick_index = mean(SI);
stick_index_2 = mean(mean_move_cue_offset_diff-mean_stay_cue_offset_diff)/mean(mean_move_cue_offset_diff+mean_stay_cue_offset_diff);

%% Sort bump parameters around jumps into 'preferred' and 'non-preferred'

%1) Get the bump parameters around all jumps, sorted by jump order
short_bm_around_jumps = [];
short_bw_around_jumps = [];

if configuration == 1
    
   short_bm_around_jumps([1,3,5,7],:) = short_bm_bar_jump([1,2,3,4],:);
   short_bw_around_jumps([1,3,5,7],:) = short_bw_bar_jump([1,2,3,4],:);
   short_bm_around_jumps([2,4,6,8],:) = short_bm_wind_jump([1,2,3,4],:);
   short_bw_around_jumps([2,4,6,8],:) = short_bm_wind_jump([1,2,3,4],:);

else
   
   short_bm_around_jumps([1,3,5,7],:) = short_bm_wind_jump([1,2,3,4],:);
   short_bw_around_jumps([1,3,5,7],:) = short_bw_wind_jump([1,2,3,4],:);
   short_bm_around_jumps([2,4,6,8],:) = short_bm_bar_jump([1,2,3,4],:);
   short_bw_around_jumps([2,4,6,8],:) = short_bm_bar_jump([1,2,3,4],:);
    
end


%2) Determine for each jump if it was her preferred cue that jumped or not
pref_cue = [];

for jump = 1:8
    
    if configuration == 1
        
        if mod(jump,2) == 0 % if this is an even jump number (which in this case is a wind jump)
            if PI(jump) > 0 %if she prefers the wind in this jump
                pref_cue(jump) = 1;
            else
                pref_cue(jump) = 0;
            end
        else %if this is a bar jump
            if PI(jump) > 0 %if she prefers the wind in this jump
                pref_cue(jump) = 0;
            else
                pref_cue(jump) = 1;
            end            
        end
        
    else
        
        if mod(jump,2) == 0 % if this is a bar jump
            if PI(jump) > 0 %if she prefers the wind in this jump
                pref_cue(jump) = 0;
            else
                pref_cue(jump) = 1;
            end
        else %if this is a wind jump
            if PI(jump) > 0 %if she prefers the wind in this jump
                pref_cue(jump) = 1;
            else
                pref_cue(jump) = 0;
            end
        end
        
    end
      
end

%3) Sort bump parameters by preference, using the pref_cue vector
short_bm_pref_cue = short_bm_around_jumps(pref_cue==1,:);
short_bw_pref_cue = short_bw_around_jumps(pref_cue==1,:);
short_bm_non_pref_cue = short_bm_around_jumps(pref_cue==0,:);
short_bw_non_pref_cue = short_bw_around_jumps(pref_cue==0,:);

%% Look at offset variance vs time: how does average variance change over time after a jump?

%I will analyze offset variability by computing a rolling window of circular
%standard deviation of offset
fcn = @(x) circ_std(x);
%Compute the offset variability over different window sizes, from ~1s to
%100s
window_sizes = [50,75,100,150];
for window = 1:length(window_sizes)
    bar_offset_variability(:,window) = matlab.tall.movingWindow(fcn,window_sizes(window),bar_offset);
    wind_offset_variability(:,window) = matlab.tall.movingWindow(fcn,window_sizes(window),wind_offset');    
end

%% Plotting it in various ways

%Regular line plots
for window = 1:length(window_sizes)
    
    figure,
    subplot(2,1,1)
    plot(bar_offset_variability(:,window),'k')
    hold on
    for jump = 1:length(real_bar_jump_frames)
        line([real_bar_jump_frames(jump) real_bar_jump_frames(jump)],[0 1.5],'color','r')
    end
    for jump = 1:length(real_wind_jump_frames)
        line([real_wind_jump_frames(jump) real_wind_jump_frames(jump)],[0 1.5],'color','b')
    end
    title('Bar offset variability');
    set(gca,'xticklabel',{[]})
        
    subplot(2,1,2)
    plot(wind_offset_variability(:,window),'k')
    hold on
    for jump = 1:length(real_bar_jump_frames)
        line([real_bar_jump_frames(jump) real_bar_jump_frames(jump)],[0 1.5],'color','r')
    end
    for jump = 1:length(real_wind_jump_frames)
        line([real_wind_jump_frames(jump) real_wind_jump_frames(jump)],[0 1.5],'color','b')
    end
    title('Wind offset variability');
    set(gca,'xticklabel',{[]})
    xlabel('Time');    
    
    suptitle(['Window size = ',num2str(window_sizes(window)),' frames']);

end

%% Averaged offset variability around jumps

for window = 1:length(window_sizes)
    figure('Position',[100 100 1400 1000]),
    
    %bar offset variability around bar jumps
    subplot(2,2,1)
    for jump = 1:length(real_bar_jump_frames)
        bar_offset_variability_abj(jump,:) = bar_offset_variability(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames,window);
    end
    plot(mean(bar_offset_variability_abj),'k','linewidth',2)
    title('Bar offset variability around bar jumps');
    hold on
    line([length(bar_offset_variability_abj)/2 length(bar_offset_variability_abj)/2],[0 1.5],'color','r','linewidth',3);
    xlim([0 length(bar_offset_variability_abj)]);
    xticks([0 length(bar_offset_variability_abj)/4 length(bar_offset_variability_abj)/2 length(bar_offset_variability_abj)*(3/4) length(bar_offset_variability_abj)]);
    xticklabels({'-120','-60','0','60','120'});
    xlabel('Time around jumps (s)');
    
    %bar offset variability around wind jumps
    subplot(2,2,2)
    for jump = 1:length(real_wind_jump_frames)
        bar_offset_variability_awj(jump,:) = bar_offset_variability(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames,window);
    end
    plot(mean(bar_offset_variability_awj),'k','linewidth',2)
    title('Bar offset variability around wind jumps');
    hold on
    line([length(bar_offset_variability_abj)/2 length(bar_offset_variability_abj)/2],[0 1.5],'color','b','linewidth',3);
    xlim([0 length(bar_offset_variability_abj)]);
    xticks([0 length(bar_offset_variability_abj)/4 length(bar_offset_variability_abj)/2 length(bar_offset_variability_abj)*(3/4) length(bar_offset_variability_abj)]);
    xticklabels({'-120','-60','0','60','120'});
    xlabel('Time around jumps (s)');
    
    %wind offset variability around bar jumps
    subplot(2,2,3)
    for jump = 1:length(real_bar_jump_frames)
        wind_offset_variability_abj(jump,:) = wind_offset_variability(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames,window);
    end
    plot(mean(wind_offset_variability_abj),'k','linewidth',2)
    title('Wind offset variability around bar jumps');
    hold on
    line([length(bar_offset_variability_abj)/2 length(bar_offset_variability_abj)/2],[0 1.5],'color','r','linewidth',3);
    xlim([0 length(bar_offset_variability_abj)]);
    xticks([0 length(bar_offset_variability_abj)/4 length(bar_offset_variability_abj)/2 length(bar_offset_variability_abj)*(3/4) length(bar_offset_variability_abj)]);
    xticklabels({'-120','-60','0','60','120'});
    xlabel('Time around jumps (s)');
    
    %wind offset variability around wind jumps
    subplot(2,2,4)
    for jump = 1:length(real_wind_jump_frames)
        wind_offset_variability_awj(jump,:) = wind_offset_variability(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames,window);
    end
    plot(mean(wind_offset_variability_awj),'k','linewidth',2)
    title('Wind offset variability around wind jumps');
    hold on
    line([length(wind_offset_variability_awj)/2 length(wind_offset_variability_awj)/2],[0 1.5],'color','b','linewidth',3);
    xlim([0 length(wind_offset_variability_awj)]);
    xticks([0 length(wind_offset_variability_awj)/4 length(wind_offset_variability_awj)/2 length(wind_offset_variability_awj)*(3/4) length(wind_offset_variability_awj)])
    xticklabels({'-120','-60','0','60','120'});
    xlabel('Time around jumps (s)');
    
    suptitle(['Window size = ',num2str(window_sizes(window)),' frames']);
    
    saveas(gcf,[path,'\plots\average_offset_var_window',num2str(window),'.png']);
    
    
    %store the offset variability values I want
    if window == 3
        bar_offset_variability_abj_75 = bar_offset_variability_abj;
        wind_offset_variability_abj_75 = wind_offset_variability_abj;
        bar_offset_variability_awj_75 = bar_offset_variability_awj;
        wind_offset_variability_awj_75 = wind_offset_variability_awj;
    elseif window == 4
        bar_offset_variability_abj_100 = bar_offset_variability_abj;
        wind_offset_variability_abj_100 = wind_offset_variability_abj;
        bar_offset_variability_awj_100 = bar_offset_variability_awj;
        wind_offset_variability_awj_100 = wind_offset_variability_awj;
    end
    
end


%% Raster plots 

figure('Position',[100 100 1600 1000])
for window = 1:length(window_sizes)
   subplot(length(window_sizes),1,window)
   imagesc(bar_offset_variability(:,window)')   
   colormap(flipud(gray))
   title(['Window size = ',num2str(window_sizes(window)),' frames']);
   hold on
    for jump = 1:length(real_bar_jump_frames)
        line([real_bar_jump_frames(jump) real_bar_jump_frames(jump)],[0.5 1.5],'color','r','linewidth',3) 
    end
    for jump = 1:length(real_wind_jump_frames)
        line([real_wind_jump_frames(jump) real_wind_jump_frames(jump)],[0.5 1.5],'color','b','linewidth',3) 
    end
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
end
suptitle('Bar offset variability');
saveas(gcf,[path,'\plots\raster_plot_bar_offset_var.png']);


figure('Position',[100 100 1600 1000])
for window = 1:length(window_sizes)
   subplot(length(window_sizes),1,window)
   imagesc(wind_offset_variability(:,window)')   
   colormap(flipud(gray))
   title(['Window size = ',num2str(window_sizes(window)),' frames']);
   hold on
    for jump = 1:length(real_bar_jump_frames)
        line([real_bar_jump_frames(jump) real_bar_jump_frames(jump)],[0.5 1.5],'color','r','linewidth',3) 
    end
    for jump = 1:length(real_wind_jump_frames)
        line([real_wind_jump_frames(jump) real_wind_jump_frames(jump)],[0.5 1.5],'color','b','linewidth',3) 
    end
    set(gca,'xticklabel',{[]})
    set(gca,'yticklabel',{[]})
end
suptitle('Wind offset variability');
saveas(gcf,[path,'\plots\raster_plot_bar_offset_var.png']);

    
%%% Scatterplot bump amplitude vs offset variance: is there a negative relationship? how strong/significant?
   
% for window = 1:length(window_sizes)
%     
%     figure,
%     
%     %Define bins
%     nbins = 15;
%     max_bin = prctile(bar_offset_variability(:,window),99,'all');
%     min_bin = prctile(bar_offset_variability(:,window),5,'all');
%     binWidth = (max_bin-min_bin)/nbins;
%     Bins = [min_bin:binWidth:max_bin];
%     
%     %Create axes for plot
%     mvtAxes = Bins - binWidth;
%     mvtAxes = mvtAxes(2:end);
%     
%     subplot(1,2,1)
%     for bin = 1:length(Bins)-1
%         meanBin(bin) = nanmean(continuous_data.bump_magnitude(bar_offset_variability(moving & good_fit,window) > Bins(bin) & (bar_offset_variability(moving & good_fit,window)) < Bins(bin+1)));
%     end
%     plot(mvtAxes,meanBin,'color',[.5 .5 .5])
%     hold on
%     ylabel('Bump magnitude'); xlabel('Bar offset variability');
%     ylim([0 max(max(meanBin))+0.5]);
%     xlim([mvtAxes(1) mvtAxes(end)]);
%     title('Bump magnitude');
%     
%     subplot(1,2,2)
%     for bin = 1:length(Bins)-1
%         meanBinw(bin) = nanmean(continuous_data.bump_width(bar_offset_variability(moving & good_fit,window) > Bins(bin) & (bar_offset_variability(moving & good_fit,window)) < Bins(bin+1)));
%     end
%     plot(mvtAxes,meanBinw,'color',[.5 .5 .5])
%     hold on
%     ylabel('Bump half width'); xlabel('Bar offset variability');
%     ylim([0 max(max(meanBinw))+0.5]);
%     xlim([mvtAxes(1) mvtAxes(end)]);
%     title('Bump half width');
%     
%     suptitle(['Window size = ',num2str(window_sizes(window)),' frames']);
%     saveas(gcf,[path,'\plots\bumps_pars_vs_bar_offset_var_',num2str(window),'.png']);
%      
% end

%% Excluding initial single cue timepoints

for window = 1:length(window_sizes)
    
    figure,
    
    %Define bins
    nbins = 15;
    max_bin = prctile(bar_offset_variability(:,window),99,'all');
    min_bin = prctile(bar_offset_variability(:,window),5,'all');
    binWidth = (max_bin-min_bin)/nbins;
    Bins = [min_bin:binWidth:max_bin];
    
    %Create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);
    
    %Define frames to include
    both_cues_on = panels_on & wind_on;
    
    subplot(1,2,1)
    for bin = 1:length(Bins)-1
        meanBin(bin) = nanmean(continuous_data.bump_magnitude(bar_offset_variability(moving & good_fit & both_cues_on',window) > Bins(bin) & (bar_offset_variability(moving & good_fit & both_cues_on',window)) < Bins(bin+1)));
    end
    plot(mvtAxes,meanBin,'color',[.5 .5 .5])
    hold on
    ylabel('Bump magnitude'); xlabel('Bar offset variability');
    ylim([0 max(max(meanBin))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump magnitude');
    
    subplot(1,2,2)
    for bin = 1:length(Bins)-1
        meanBinw(bin) = nanmean(continuous_data.bump_width(bar_offset_variability(moving & good_fit & both_cues_on',window) > Bins(bin) & (bar_offset_variability(moving & good_fit & both_cues_on',window)) < Bins(bin+1)));
    end
    plot(mvtAxes,meanBinw,'color',[.5 .5 .5])
    hold on
    ylabel('Bump half width'); xlabel('Bar offset variability');
    ylim([0 max(max(meanBinw))+0.5]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    title('Bump half width');
    
    suptitle(['Window size = ',num2str(window_sizes(window)),' frames']);
    saveas(gcf,[path,'\plots\bumps_pars_vs_bar_offset_var_only_cc_',num2str(window),'.png']);
     
end


%% Look at cross-correlation between bump parameters and offset variance: what is the time of the peak?

for window = 1:length(window_sizes)
    
    figure('Position',[100 100 1600 600]),
    subplot(1,2,1)
    [c,lags] = xcorr(bar_offset_variability(moving & good_fit & both_cues_on',window),continuous_data.bump_magnitude(moving & good_fit & both_cues_on'));
    stem(lags,c)
    corr_val = corrcoef(bar_offset_variability(moving & good_fit & both_cues_on',window),continuous_data.bump_magnitude(moving & good_fit & both_cues_on'));
    title(['Bump magnitude, corr = ',num2str(corr_val(2,1))]);
    
    subplot(1,2,2)
    [c,lags] = xcorr(bar_offset_variability(moving & good_fit & both_cues_on',window),continuous_data.bump_width(moving & good_fit & both_cues_on'));
    stem(lags,c)
    corr_val2 =  corrcoef(bar_offset_variability(moving & good_fit & both_cues_on',window),continuous_data.bump_width(moving & good_fit & both_cues_on'));
    title(['Bump width, corr = ',num2str(corr_val2(2,1))]);
    
    suptitle(['Window size = ',num2str(window_sizes(window)),' frames']);
    saveas(gcf,[path,'\plots\cross_corr_',num2str(window),'.png']);
    
end

%% Save data

save([path,'\data.mat'],'short_bm_bar_jump','short_bw_bar_jump','short_bm_wind_jump','short_bw_wind_jump','long_bm_bar_jump','long_bw_bar_jump','long_bm_wind_jump','long_bw_wind_jump','pref_index','pref_index_2','heading_pref_index','heading_pref_index_2','stick_index','stick_index_2','initial_bar_bm','initial_wind_bm','initial_bar_bw','initial_wind_bw','initial_bar_offset','initial_wind_offset','PI','PI_heading','SI','wind_offset_variability_awj_75','wind_offset_variability_abj_75','bar_offset_variability_abj_75','bar_offset_variability_awj_75','wind_offset_variability_awj_100','wind_offset_variability_abj_100','bar_offset_variability_abj_100','bar_offset_variability_awj_100','configuration','long_total_mvt_wind_jump','long_total_mvt_bar_jump','long_rot_speed_wind_jump','long_rot_speed_bar_jump','short_bm_pref_cue','short_bw_pref_cue','short_bm_non_pref_cue','short_bw_non_pref_cue')

%% Clear space

clear all; close all;
