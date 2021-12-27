%Full experiment analysis
%Code to analyze the full experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp38\data\second_version\');
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

folderNames = dir(path(1:68));
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
    if strcmp(flyNames(fly).name,path(69:end-10))
        fly_ID = fly;
    end
end

colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
    0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880;...
    0 0.4470 0.7410; 0.75, 0.1, 0.75; 0.75, 0.75, 0];

fly_color = colors_for_plots(fly_ID,:);

%% Determine when the stimuli are on

%Uncomment the below figure to check when the panels and wind are on
% figure,
% plot(continuous_data.fr_y_ds)
% hold on
% plot(continuous_data.wind_valve)
% legend('Panels on','Wind on');

panels_on = continuous_data.fr_y_ds>7;
wind_on = continuous_data.wind_valve>2;

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


%% Define block type (1 = bar, 2 = wind, 3 = both pre jump, 4 = both post jump)

block = {};

if configuration == 1
    block{1} = 1;
    block{2} = 2;
    block{3} = 3;
    block{4} = 1;
    block{5} = 2;
    block{6} = 3;
    block{7} = 3;
    block{8} = 3;
    block{9} = 3;
    block{10} = 3;
    block{11} = 3;
    block{12} = 3;
else
    block{1} = 2;
    block{2} = 1;
    block{3} = 3;
    block{4} = 2;
    block{5} = 1;
    block{6} = 3;
    block{7} = 3;
    block{8} = 3;
    block{9} = 3;
    block{10} = 3;
    block{11} = 3;
    block{12} = 3;
end

%% Find the jump frames

if configuration == 1
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1380*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(1980*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2580*length(continuous_data.dff_matrix)/continuous_data.time(end))];
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1680*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2280*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2880*length(continuous_data.dff_matrix)/continuous_data.time(end))];   
else
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1380*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(1980*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2580*length(continuous_data.dff_matrix)/continuous_data.time(end))];
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1680*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2280*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2880*length(continuous_data.dff_matrix)/continuous_data.time(end))];   
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
    real_bar_jump_frames(2) = floor(2270*length(continuous_data.dff_matrix)/continuous_data.time(end));
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
    real_wind_jump_frames(2) = floor(1970*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID ==2
    real_wind_jump_frames(3) = floor(2870*length(continuous_data.dff_matrix)/continuous_data.time(end));
end
real_wind_jump_sec = real_wind_jump_frames*continuous_data.time(end)/length(continuous_data.dff_matrix);


%% Define the boundaries for the different blocks

boundaries = {};

if configuration == 1
    
    boundaries{1} = [480*sec_to_frames,panels_change_frames(1)];
    boundaries{2} = [panels_change_frames(2)-120*sec_to_frames,panels_change_frames(2)];
    boundaries{3} = [panels_change_frames(2),wind_change_frames(2)];
    boundaries{4} = [wind_change_frames(2),panels_change_frames(3)];
    boundaries{5} = [panels_change_frames(3),panels_change_frames(4)];
    boundaries{6} = [panels_change_frames(4),real_bar_jump_frames(1)];
    boundaries{7} = [real_bar_jump_frames(1),real_wind_jump_frames(1)];
    boundaries{8} = [real_wind_jump_frames(1),real_bar_jump_frames(2)];
    boundaries{9} = [real_bar_jump_frames(2),real_wind_jump_frames(2)];
    boundaries{10} = [real_wind_jump_frames(2),real_bar_jump_frames(3)];
    boundaries{11} = [real_bar_jump_frames(3),real_wind_jump_frames(3)];
    boundaries{12} = [real_wind_jump_frames(3),wind_change_frames(4)];  
   
else

    boundaries{1} = [480*sec_to_frames,wind_change_frames(1)];
    boundaries{2} = [wind_change_frames(1)-120*sec_to_frames,wind_change_frames(2)];
    boundaries{3} = [wind_change_frames(2),panels_change_frames(2)];
    boundaries{4} = [panels_change_frames(2),panels_change_frames(3)];
    boundaries{5} = [panels_change_frames(3),wind_change_frames(4)];
    boundaries{6} = [wind_change_frames(4),real_wind_jump_frames(1)];
    boundaries{7} = [real_wind_jump_frames(1),real_bar_jump_frames(1)];
    boundaries{8} = [real_bar_jump_frames(1),real_wind_jump_frames(2)];
    boundaries{9} = [real_wind_jump_frames(2),real_bar_jump_frames(2)];
    boundaries{10} = [real_bar_jump_frames(2),real_wind_jump_frames(3)];
    boundaries{11} = [real_wind_jump_frames(3),real_bar_jump_frames(3)];
    boundaries{12} = [real_bar_jump_frames(3),panels_change_frames(4)];
        
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
for jump = 2:3
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
for jump = 2:3
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
for jump = 2:3
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
for jump = 2:3
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
for jump = 2:3
   xline(real_wind_jump_sec(jump),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
   xline(real_bar_jump_sec(jump),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
end
xlabel('Time (sec)');
title('Bump width');
xlim([0 continuous_data.time(end)]);

saveas(gcf,[path,'\plots\full_experiment.png']);


%% Close-up around the jumps

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


%% Bump magnitude around jumps

all_jumps = [real_bar_jump_frames,real_wind_jump_frames];
jump_nature = zeros(1,length(all_jumps));
jump_nature(1:length(real_bar_jump_frames)) = 1;

for jump = 1:length(all_jumps)
    
    prej_moving = continuous_data.total_mvt_ds(all_jumps(jump)-floor(3*sec_to_frames):all_jumps(jump))>25;
    prej_fit = continuous_data.adj_rs(all_jumps(jump)-floor(3*sec_to_frames):all_jumps(jump))>0.5;
    
    pre_jump_BM = continuous_data.bump_magnitude(:,all_jumps(jump)-floor(3*sec_to_frames):all_jumps(jump));
    pre_jump_bm(jump) = nanmean(pre_jump_BM(prej_fit));
    pre_jump_BW = continuous_data.bump_width(:,all_jumps(jump)-floor(3*sec_to_frames):all_jumps(jump));
    pre_jump_bw(jump) = nanmean(pre_jump_BW(prej_fit));
    
    
    postj_moving = continuous_data.total_mvt_ds(all_jumps(jump)+1:all_jumps(jump)+floor(3*sec_to_frames))>25;
    postj_fit = continuous_data.adj_rs(all_jumps(jump)+1:all_jumps(jump)+floor(3*sec_to_frames))>0.5;
    
    post_jump_BM = continuous_data.bump_magnitude(:,all_jumps(jump)+1:all_jumps(jump)+floor(3*sec_to_frames));
    post_jump_bm(jump) = nanmean(post_jump_BM(postj_fit));
    post_jump_BW = continuous_data.bump_width(:,all_jumps(jump)+1:all_jumps(jump)+floor(3*sec_to_frames));
    post_jump_bw(jump) = nanmean(post_jump_BW(postj_fit));
    
end

around_jump_bm = [pre_jump_bm;post_jump_bm];
around_jump_bw = [pre_jump_bw;post_jump_bw];

figure,
subplot(1,2,1)
for jump = 1:length(jump_nature)
    if jump_nature(jump) == 1 %if a bar
        plot(around_jump_bm(:,jump),'-o','color',[0.8500 0.3250 0.0980])
        hold on
    else
        plot(around_jump_bm(:,jump),'-o','color',[0.9290 0.6940 0.1250])
        hold on
    end
end
xticks([1 2]); xlim([0 3]);
ylim([0 2.5]);
ylabel('Bump magnitude');
xticklabels({'pre jump','post jump'});

subplot(1,2,2)
for jump = 1:length(jump_nature)
    if jump_nature(jump) == 1 %if a bar
        plot(around_jump_bw(:,jump),'-o','color',[0.8500 0.3250 0.0980])
        hold on
    else
        plot(around_jump_bw(:,jump),'-o','color',[0.9290 0.6940 0.1250])
        hold on
    end
end
xticks([1 2]); xlim([0 3]);
ylim([0 3.5]);
ylabel('Bump width');
xticklabels({'pre jump','post jump'});

saveas(gcf,[path,'\plots\aj_bump_parameters.png']);


%% Look at a longer around the jump period

for jump = 1:length(real_bar_jump_frames)
    
   time_zero = continuous_data.time(real_bar_jump_frames(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1300 800]),
   ax(1) = subplot(12,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   colormap(ax(1),flipud(gray));
   pos = get(subplot(12,1,1),'Position');
   %pos2 = get(subplot(12,1,2),'Position');
   h = colorbar('Position', [pos(1)+pos(3)+0.02  0.8054  pos(3)/40  pos(4)+0.055]); 
   
   ax(2) = subplot(12,1,2);
   imagesc(continuous_data.bump_width(:,real_bar_jump_frames(jump)-120*sec_to_frames:real_bar_jump_frames(jump)+120*sec_to_frames))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');    
   
   c2 = colormap(ax(2),flipud(bone));
   pos2 = get(subplot(12,1,2),'Position');
   h2 = colorbar('Position', [pos(1)- 0.06  pos2(2)  pos(3)/40  pos(4)+0.055]);
   
   ax(3) = subplot(12,1,[3 9]);
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
   
   ax(4) = subplot(12,1,[10:11]);
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
   
   figure('Position',[100 100 1200 800]),
   ax(1) = subplot(12,1,1);
   imagesc(continuous_data.bump_magnitude(:,real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   colormap(ax(1),flipud(gray));
   pos = get(subplot(12,1,1),'Position');
   %pos2 = get(subplot(12,1,2),'Position');
   h = colorbar('Position', [pos(1)+pos(3)+0.02  0.8054  pos(3)/40  pos(4)+0.055]); 
   
   ax(2) = subplot(12,1,2);
   imagesc(continuous_data.bump_width(:,real_wind_jump_frames(jump)-120*sec_to_frames:real_wind_jump_frames(jump)+120*sec_to_frames))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');
         
   c2 = colormap(ax(2),flipud(bone));
   pos2 = get(subplot(12,1,2),'Position');
   h2 = colorbar('Position', [pos(1)- 0.06  pos2(2)  pos(3)/40  pos(4)+0.055]);
   
   ax(3) = subplot(12,1,[3 9]);
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
   
   ax(4) = subplot(12,1,[10:12]);
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


%% Look at offset variability of last darkness bout

figure,
polarhistogram(deg2rad(bar_offset(~panels_on & ~wind_on)),15,'FaceColor','k')
set(gca,'ThetaZeroLocation','top');
title('No stimulus');
saveas(gcf,[path,'\plots\no_stim_offset.png']);


%% Threshold movement data and goodness of fit

% Uncomment below to see the total movement distribution
% figure,
% histogram(continuous_data.total_mvt_ds)
% xlabel('Total movement (deg/s)');
% ylabel('Counts');
mvt_thresh = 25;
moving = continuous_data.total_mvt_ds>mvt_thresh;
good_fit = continuous_data.adj_rs>0.5;

%% Offset and bump parameter analysis thresholding the movement and the fit

figure('Position',[0 300 1800 800]),

subplot(6,12,[5:8])
set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});
text(0.5, 0.5,'Bar offset','fontsize',18)
axis off;

subplot(6,12,[29:32])
set(gca,'xticklabel',{[]});
set(gca,'yticklabel',{[]});
text(0.48, 0.5,'Wind offset','fontsize',18)
axis off;

for i = 1:length(boundaries)
    subplot(6,length(boundaries),i+length(boundaries))
    baroff = bar_offset(boundaries{i}(1):boundaries{i}(2));
    if block{i} == 1
        polarhistogram(deg2rad(baroff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),15,'FaceColor',[0.8500 0.3250 0.0980])
        title('Bar');
    elseif block{i} ==2
        polarhistogram(deg2rad(baroff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),15,'FaceColor',[0.9290 0.6940 0.1250])
        title('Wind');
    else
        polarhistogram(deg2rad(baroff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),15,'FaceColor',[0.4660 0.6740 0.1880])
        title('Both');
    end
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    
    subplot(6,length(boundaries),i+length(boundaries)*3)
    windoff = wind_offset(boundaries{i}(1):boundaries{i}(2));
    if block{i} == 1
        polarhistogram(deg2rad(windoff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),15,'FaceColor',[0.8500 0.3250 0.0980])
        title('Bar');
    elseif block{i} ==2
        polarhistogram(deg2rad(windoff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),15,'FaceColor',[0.9290 0.6940 0.1250])
        title('Wind');
    else
        polarhistogram(deg2rad(windoff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),15,'FaceColor',[0.4660 0.6740 0.1880])
        title('Both');
    end
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
end


%add the bump parameters
subplot(6,length(boundaries),[49:72])
meanBM_perblock_thresh = [];
meanBW_perblock_thresh = [];
for bout = 1:length(boundaries)
    bump_mag = continuous_data.bump_magnitude(boundaries{bout}(1):boundaries{bout}(2));
    bump_width = continuous_data.bump_width(boundaries{bout}(1):boundaries{bout}(2));
    meanBM_perblock_thresh(bout) = nanmean(bump_mag(continuous_data.total_mvt_ds(boundaries{bout}(1):boundaries{bout}(2))>25 & continuous_data.adj_rs(boundaries{bout}(1):boundaries{bout}(2))>0.5));
    meanBW_perblock_thresh(bout) = nanmean(bump_width(continuous_data.total_mvt_ds(boundaries{bout}(1):boundaries{bout}(2))>25 & continuous_data.adj_rs(boundaries{bout}(1):boundaries{bout}(2))>0.5));
end
plot(meanBM_perblock_thresh,'-ko')
hold on
plot(meanBW_perblock_thresh,'-ro')
xlim([0.6 12.4]); ylim([0 3]);
xlabel('Bout #');
legend('Bump mag','Bump width');

saveas(gcf,[path,'\plots\bump_parameters_and_offset_per_block_mvt_and_fit_thresh.png']);


%% Offset means and variability

for i = 1:length(boundaries)
    baroff = bar_offset(boundaries{i}(1):boundaries{i}(2));
    bar_offset_mean(i) = rad2deg(circ_mean(deg2rad(baroff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5))));
    [~, bar_offset_var(i)] = circ_std(deg2rad(baroff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)));
    
    windoff = wind_offset(boundaries{i}(1):boundaries{i}(2));
    wind_offset_mean(i) = rad2deg(circ_mean(deg2rad(windoff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)),[],2));    
    [~, wind_offset_var(i)] = circ_std(deg2rad(baroff(continuous_data.total_mvt_ds(boundaries{i}(1):boundaries{i}(2))>25 & continuous_data.adj_rs(boundaries{i}(1):boundaries{i}(2))>0.5)));
end




%% Bump parameters vs total movement

nbins = 20;
maxBin = max(continuous_data.total_mvt_ds);
minBin = min(continuous_data.total_mvt_ds);
binWidth = (maxBin-minBin)/nbins;
mvtBins = [minBin:binWidth:maxBin]; 

%getting binned means 
for bin = 1:length(mvtBins)-1
    meanBinnedBMBar(bin) = mean(continuous_data.bump_magnitude((continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) > mvtBins(bin)) & panels_on' & (continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) < mvtBins(bin+1))));
    meanBinnedBMWind(bin) = mean(continuous_data.bump_magnitude((continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) > mvtBins(bin)) & wind_on' & (continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) < mvtBins(bin+1))));

    meanBinnedBWBar(bin) = mean(continuous_data.bump_width((continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) > mvtBins(bin)) & panels_on' & (continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) < mvtBins(bin+1))));
    meanBinnedBWWind(bin) = mean(continuous_data.bump_width((continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) > mvtBins(bin)) & wind_on' & (continuous_data.total_mvt_ds(1:length(continuous_data.bump_magnitude)) < mvtBins(bin+1))));
end

%create axes for plot
mvtAxes = mvtBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
figure('Position',[200 200 1400 600]),
subplot(1,2,1)
plot(mvtAxes,meanBinnedBMBar,'-o','color',[0.8500 0.3250 0.0980])
hold on
plot(mvtAxes,meanBinnedBMWind,'-o','color',[0.9290 0.6940 0.1250])
ylabel('Mean bump magnitude'); xlabel('Total movement (deg/s)');
ylim([0 2]);
legend('Bar','Wind');

subplot(1,2,2)
plot(mvtAxes,meanBinnedBWBar,'-o','color',[0.8500 0.3250 0.0980])
hold on
plot(mvtAxes,meanBinnedBWWind,'-o','color',[0.9290 0.6940 0.1250])
ylabel('Mean bump width'); xlabel('Total movement (deg/s)');
ylim([0 3]);
legend('Bar','Wind');

saveas(gcf,[path,'\plots\bump_par_vs_mvt.png']);


%% Get mean total movement for the different blocks

mean_total_mvt = [];

for bout = 1:length(boundaries)
    total_mvt = continuous_data.total_mvt_ds(boundaries{bout}(1):boundaries{bout}(2));
    mean_total_mvt(bout) = nanmean(total_mvt(continuous_data.total_mvt_ds(boundaries{bout}(1):boundaries{bout}(2))>25));
end


%% Reproduce Exp35 analysis (focus only on first 5 blocks)

figure('Position',[100 100 1000 600]),

subplot(1,2,1)
yyaxis left
plot(meanBM_perblock_thresh(1:5),'-o')
ylim([0 3]);
ylabel('Bump magnitude');
yyaxis right
plot(mean_total_mvt(1:5),'-o')
xlim([0 6]);
ylim([0 300]);
ylabel('Total movement (deg/s)');
xlabel('Block #');

subplot(1,2,2)
yyaxis left
plot(meanBW_perblock_thresh(1:5),'-o')
ylim([0 3.5]);
ylabel('Bump width');
yyaxis right
plot(mean_total_mvt(1:5),'-o')
xlim([0 6]);
ylim([0 300]);
ylabel('Total movement (deg/s)');
xlabel('Block #');

saveas(gcf,[path,'\plots\Exp35_rep.png']);


%% Save data

%save([path,'data.mat'],'meanBM_perblock_thresh','mean_total_mvt','meanBW_perblock_thresh','alldata_BM_thresh','alldata_BW_thresh','stim_ID_thresh','pre_BM_prejump','pre_BM_postjump','pre_BW_prejump','bar_offset_diff','wind_offset_diff','initial_cue_diff')
save([path,'data.mat'],'meanBM_perblock_thresh','mean_total_mvt','meanBW_perblock_thresh','bar_offset_mean','wind_offset_mean','wind_offset_var','bar_offset_var','around_jump_bm','around_jump_bw','jump_nature','configuration')


%% Close

clear all; close all;