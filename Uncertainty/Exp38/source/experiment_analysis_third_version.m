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
if fly_ID == 1 
    real_bar_jump_frames(1) = floor(1488.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.5*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.5*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.5*length(continuous_data.dff_matrix)/continuous_data.time(end));    
elseif (fly_ID == 2 | fly_ID == 6)  
    real_bar_jump_frames(1) = floor(1488.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2088.7*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(3) = floor(2688.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(4) = floor(3288.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif (fly_ID == 3 | fly_ID == 4 | fly_ID == 5)
    real_bar_jump_frames(1) = floor(1188.7*length(continuous_data.dff_matrix)/continuous_data.time(end));    
    real_bar_jump_frames(2) = floor(1788.7*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(3) = floor(2388.5*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_bar_jump_frames(4) = floor(2988.7*length(continuous_data.dff_matrix)/continuous_data.time(end));   
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
    real_wind_jump_frames(1) = floor(1188.56*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.78*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2388.66*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.76*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 4
    real_wind_jump_frames(1) = floor(1488.74*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2088.7*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2688.54*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3288.72*length(continuous_data.dff_matrix)/continuous_data.time(end));   
elseif(fly_ID == 3 | fly_ID == 4 | fly_ID == 5)
    real_wind_jump_frames(1) = floor(1489*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2089*length(continuous_data.dff_matrix)/continuous_data.time(end));  
    real_wind_jump_frames(3) = floor(2689*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(3289*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 6
    real_wind_jump_frames(1) = floor(1188.96*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(1788.6*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(3) = floor(2388.26*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(4) = floor(2988.46*length(continuous_data.dff_matrix)/continuous_data.time(end));
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

%% Bump par before and after bar jumps

for jump = 1:length(real_bar_jump_frames)
    
    short_bm_around_bar_jump(jump,1) = mean(continuous_data.bump_magnitude(real_bar_jump_frames(jump)-1*sec_to_frames:real_bar_jump_frames(jump)-1));
    short_bm_around_bar_jump(jump,2) = mean(continuous_data.bump_magnitude(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+1*sec_to_frames));
    
    short_bw_around_bar_jump(jump,1) = mean(continuous_data.bump_width(real_bar_jump_frames(jump)-1*sec_to_frames:real_bar_jump_frames(jump)-1));
    short_bw_around_bar_jump(jump,2) = mean(continuous_data.bump_width(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+1*sec_to_frames));
    
    long_bm_around_bar_jump(jump,1) = mean(continuous_data.bump_magnitude(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1));
    long_bm_around_bar_jump(jump,2) = mean(continuous_data.bump_magnitude(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+120*sec_to_frames));
    
    long_bw_around_bar_jump(jump,1) = mean(continuous_data.bump_width(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1));
    long_bw_around_bar_jump(jump,2) = mean(continuous_data.bump_width(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+120*sec_to_frames));
    
end

figure('Position',[100 100 800 1000]),
subplot(2,2,1)
plot(short_bm_around_bar_jump','color',[.5 .5 .5])
hold on
plot(median(short_bm_around_bar_jump),'k','linewidth',2)
ylabel('Bump magnitude');
title('Short timescale');
xlim([0 3]);
ylim([0.5 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,2)
plot(long_bm_around_bar_jump','color',[.5 .5 .5])
hold on
plot(median(long_bm_around_bar_jump),'k','linewidth',2)
title('Long timescale');
xlim([0 3]);
ylim([0.5 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,3)
plot(short_bw_around_bar_jump','color',[.5 .5 .5])
hold on
plot(median(short_bw_around_bar_jump),'k','linewidth',2)
ylabel('Bump width');
xlim([0 3]);
ylim([1 4]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,4)
plot(long_bw_around_bar_jump','color',[.5 .5 .5])
hold on
plot(median(long_bw_around_bar_jump),'k','linewidth',2)
xlim([0 3]);
ylim([1 4]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

suptitle('Bar jumps');

saveas(gcf,[path,'\plots\bump_pars_around_bar_jumps.png']);

%% Repeat for wind

for jump = 1:length(real_wind_jump_frames)
    
    short_bm_around_wind_jump(jump,1) = mean(continuous_data.bump_magnitude(real_wind_jump_frames(jump)-1*sec_to_frames:real_wind_jump_frames(jump)-1));
    short_bm_around_wind_jump(jump,2) = mean(continuous_data.bump_magnitude(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+1*sec_to_frames));
    
    short_bw_around_wind_jump(jump,1) = mean(continuous_data.bump_width(real_wind_jump_frames(jump)-1*sec_to_frames:real_wind_jump_frames(jump)-1));
    short_bw_around_wind_jump(jump,2) = mean(continuous_data.bump_width(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+1*sec_to_frames));
    
    long_bm_around_wind_jump(jump,1) = mean(continuous_data.bump_magnitude(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1));
    long_bm_around_wind_jump(jump,2) = mean(continuous_data.bump_magnitude(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+120*sec_to_frames));
    
    long_bw_around_wind_jump(jump,1) = mean(continuous_data.bump_width(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1));
    long_bw_around_wind_jump(jump,2) = mean(continuous_data.bump_width(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+120*sec_to_frames));
    
end

figure('Position',[100 100 800 1000]),
subplot(2,2,1)
plot(short_bm_around_wind_jump','color',[.5 .5 .5])
hold on
plot(median(short_bm_around_wind_jump),'k','linewidth',2)
ylabel('Bump magnitude');
title('Short timescale');
xlim([0 3]);
ylim([0.5 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,2)
plot(long_bm_around_wind_jump','color',[.5 .5 .5])
hold on
plot(median(long_bm_around_wind_jump),'k','linewidth',2)
title('Long timescale');
xlim([0 3]);
ylim([0.5 3]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,3)
plot(short_bw_around_bar_jump','color',[.5 .5 .5])
hold on
plot(median(short_bw_around_bar_jump),'k','linewidth',2)
ylabel('Bump width');
xlim([0 3]);
ylim([1 4]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

subplot(2,2,4)
plot(long_bw_around_wind_jump','color',[.5 .5 .5])
hold on
plot(median(long_bw_around_wind_jump),'k','linewidth',2)
xlim([0 3]);
ylim([1 4]);
xticks([1 2]);
xticklabels({'pre jump','post jump'});

suptitle('Wind jumps');

saveas(gcf,[path,'\plots\bump_pars_around_wind_jumps.png']);

%% Raster plots of bump parameters in the short time scale

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


%I'm pretty sure we need to zscore, but is this zscoring correct?

%% Repeat for long time scale

%change timescale for seconds

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


%% Compute offset variability pre and post-jump



%compare to changes in bump parameters.



%% Offset mean with respect to both stimuli pre and post-jump

%Compute offset mean around bar jumps
for jump = 1:length(real_bar_jump_frames)
    
    short_bar_offset_mean_around_bar_jump(jump,1) = circ_mean(deg2rad(bar_offset(real_bar_jump_frames(jump)-2*sec_to_frames:real_bar_jump_frames(jump)-1)));
    short_bar_offset_mean_around_bar_jump(jump,2) = circ_mean(deg2rad(bar_offset(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+2*sec_to_frames)));
    
    short_wind_offset_mean_around_bar_jump(jump,1) = circ_mean(deg2rad(wind_offset(real_bar_jump_frames(jump)-2*sec_to_frames:real_bar_jump_frames(jump)-1)),[],2);
    short_wind_offset_mean_around_bar_jump(jump,2) = circ_mean(deg2rad(wind_offset(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+2*sec_to_frames)),[],2);
    
    long_bar_offset_mean_around_bar_jump(jump,1) = circ_mean(deg2rad(bar_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1)));
    long_bar_offset_mean_around_bar_jump(jump,2) = circ_mean(deg2rad(bar_offset(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+120*sec_to_frames)));
    
    long_wind_offset_mean_around_bar_jump(jump,1) = circ_mean(deg2rad(wind_offset(real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)-1)),[],2);
    long_wind_offset_mean_around_bar_jump(jump,2) = circ_mean(deg2rad(wind_offset(real_bar_jump_frames(jump)+1:real_bar_jump_frames(jump)+120*sec_to_frames)),[],2);
    
end

%Compute differences in the offset mean pre and post jump
short_bar_offset_mean_diff_around_bar_jumps = abs(rad2deg(circ_dist(short_bar_offset_mean_around_bar_jump(:,1),short_bar_offset_mean_around_bar_jump(:,2))));
short_wind_offset_mean_diff_around_bar_jumps = abs(rad2deg(circ_dist(short_wind_offset_mean_around_bar_jump(:,1),short_wind_offset_mean_around_bar_jump(:,2))));
long_bar_offset_mean_diff_around_bar_jumps = abs(rad2deg(circ_dist(long_bar_offset_mean_around_bar_jump(:,1),long_bar_offset_mean_around_bar_jump(:,2))));
long_wind_offset_mean_diff_around_bar_jumps = abs(rad2deg(circ_dist(long_wind_offset_mean_around_bar_jump(:,1),long_wind_offset_mean_around_bar_jump(:,2))));

%The short bar offset mean will change quite a bit by design because of the
%bar jumps, so not sure if it's the best metric.


%Repeat for wind jumps
for jump = 1:length(real_wind_jump_frames)
    
    short_bar_offset_mean_around_wind_jump(jump,1) = circ_mean(deg2rad(bar_offset(real_wind_jump_frames(jump)-2*sec_to_frames:real_wind_jump_frames(jump)-1)));
    short_bar_offset_mean_around_wind_jump(jump,2) = circ_mean(deg2rad(bar_offset(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+2*sec_to_frames)));
    
    short_wind_offset_mean_around_wind_jump(jump,1) = circ_mean(deg2rad(wind_offset(real_wind_jump_frames(jump)-2*sec_to_frames:real_wind_jump_frames(jump)-1)),[],2);
    short_wind_offset_mean_around_wind_jump(jump,2) = circ_mean(deg2rad(wind_offset(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+2*sec_to_frames)),[],2);
    
    long_bar_offset_mean_around_wind_jump(jump,1) = circ_mean(deg2rad(bar_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1)));
    long_bar_offset_mean_around_wind_jump(jump,2) = circ_mean(deg2rad(bar_offset(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+120*sec_to_frames)));
    
    long_wind_offset_mean_around_wind_jump(jump,1) = circ_mean(deg2rad(wind_offset(real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)-1)),[],2);
    long_wind_offset_mean_around_wind_jump(jump,2) = circ_mean(deg2rad(wind_offset(real_wind_jump_frames(jump)+1:real_wind_jump_frames(jump)+120*sec_to_frames)),[],2);
    
end

%Compute differences in the offset mean pre and post jump
short_bar_offset_mean_diff_around_wind_jumps = abs(rad2deg(circ_dist(short_bar_offset_mean_around_wind_jump(:,1),short_wind_offset_mean_around_wind_jump(:,2))));
short_wind_offset_mean_diff_around_wind_jumps = abs(rad2deg(circ_dist(short_wind_offset_mean_around_wind_jump(:,1),short_wind_offset_mean_around_wind_jump(:,2))));
long_bar_offset_mean_diff_around_wind_jumps = abs(rad2deg(circ_dist(long_bar_offset_mean_around_wind_jump(:,1),long_wind_offset_mean_around_wind_jump(:,2))));
long_wind_offset_mean_diff_around_wind_jumps = abs(rad2deg(circ_dist(long_wind_offset_mean_around_wind_jump(:,1),long_wind_offset_mean_around_wind_jump(:,2))));


%Plot differences in mean offset
figure('Position',[100 100 800 800]),
plot([1,1,1,1],long_bar_offset_mean_diff_around_bar_jumps,'ro')
hold on
plot([1,1,1,1],long_bar_offset_mean_diff_around_wind_jumps,'ko')
plot([2,2,2,2],long_wind_offset_mean_diff_around_bar_jumps,'ro')
plot([2,2,2,2],long_wind_offset_mean_diff_around_wind_jumps,'ko')
xlim([0 3]);
ylim([0 180]);
xticks([1 2]);
xticklabels({'Bar offset diff','Wind offset diff'});
ylabel('Post-pre mean offset difference (deg)','fontsize',12);
legend('Around bar jumps','Around wind jumps');


%Compute the preference index (proxy for which cue the bump is focusing on)
    %the closer to -1, the more the bar is preferred
    %the closer to 1, the more the wind is preferred
pref_index = (mean([long_bar_offset_mean_diff_around_bar_jumps;long_bar_offset_mean_diff_around_wind_jumps])-mean([long_wind_offset_mean_diff_around_bar_jumps;long_wind_offset_mean_diff_around_wind_jumps]))/(mean([long_bar_offset_mean_diff_around_bar_jumps;long_bar_offset_mean_diff_around_wind_jumps]) + mean([long_wind_offset_mean_diff_around_bar_jumps;long_wind_offset_mean_diff_around_wind_jumps]));

title(['Preference index = ',num2str(round(pref_index,2))]);

saveas(gcf,[path,'\plots\change_in_mean_offset.png']);

%% Save data

%save PI, and bump parameters behavior around the jumps
save([path,'\data.mat'],'short_bm_bar_jump','short_bw_bar_jump','short_bm_wind_jump','short_bw_wind_jump','long_bm_bar_jump','long_bw_bar_jump','long_bm_wind_jump','long_bw_wind_jump','pref_index')

%% Clear space

clear all; close all;
