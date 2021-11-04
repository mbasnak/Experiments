%Full experiment analysis
%Code to analyze the full experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp38\data\');
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
    if strcmp(flyNames(fly).name,path(54:end-10))
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
    block{7} = 4;
    block{8} = 1;
    block{9} = 2;
    block{10} = 3;
    block{11} = 4;
    block{12} = 1;
    block{13} = 2;
    block{14} = 3;
    block{15} = 4;
    block{16} = 1;
    block{17} = 2;
else
    block{1} = 2;
    block{2} = 1;
    block{3} = 3;
    block{4} = 2;
    block{5} = 1;
    block{6} = 3;
    block{7} = 4;
    block{8} = 2;
    block{9} = 1;
    block{10} = 3;
    block{11} = 4;
    block{12} = 2;
    block{13} = 1;
    block{14} = 3;
    block{15} = 4;
    block{16} = 2;
    block{17} = 1;
end

%% Define the boundaries for the different blocks

boundaries = {};

if configuration == 1
    
    boundaries{1} = [480*sec_to_frames,panels_change_frames(1)];
    boundaries{2} = [panels_change_frames(2)-120*sec_to_frames,panels_change_frames(2)];
    boundaries{3} = [panels_change_frames(2),wind_change_frames(2)];
    boundaries{4} = [wind_change_frames(2),panels_change_frames(3)];
    boundaries{5} = [panels_change_frames(3),panels_change_frames(4)];
    both_cues_2_full = wind_change_frames(4) - panels_change_frames(4);
    boundaries{6} = [panels_change_frames(4),floor(panels_change_frames(4)+0.5*both_cues_2_full)];
    boundaries{7} = [floor(panels_change_frames(4)+0.5*both_cues_2_full), wind_change_frames(4)];
    boundaries{8} = [wind_change_frames(4),wind_change_frames(5)];
    boundaries{9} = [wind_change_frames(5),panels_change_frames(6)];
    both_cues_3_full = wind_change_frames(6) - panels_change_frames(6);
    boundaries{10} = [panels_change_frames(6),floor(panels_change_frames(6)+0.5*both_cues_3_full)];
    boundaries{11} = [floor(panels_change_frames(6)+0.5*both_cues_3_full),wind_change_frames(6)];
    boundaries{12} = [wind_change_frames(6),wind_change_frames(7)];
    boundaries{13} = [wind_change_frames(7),panels_change_frames(8)];
    both_cues_4_full = wind_change_frames(8) - panels_change_frames(8);
    boundaries{14} = [panels_change_frames(8),floor(panels_change_frames(8)+0.5*both_cues_4_full)];
    boundaries{15} = [floor(panels_change_frames(8)+0.5*both_cues_4_full),wind_change_frames(8)];
    boundaries{16} = [wind_change_frames(8),wind_change_frames(9)];
    if fly_ID > 3
        boundaries{17} = [wind_change_frames(9),wind_change_frames(10)];
    else
        boundaries{17} = [wind_change_frames(9),length(wind_on)];
    end
else

    boundaries{1} = [480*sec_to_frames,wind_change_frames(1)];
    boundaries{2} = [wind_change_frames(1)-120*sec_to_frames,wind_change_frames(2)];
    boundaries{3} = [wind_change_frames(2),panels_change_frames(2)];
    boundaries{4} = [panels_change_frames(2),panels_change_frames(3)];
    boundaries{5} = [panels_change_frames(3),wind_change_frames(4)];
    both_cues_2_full = panels_change_frames(4) - wind_change_frames(4);
    boundaries{6} = [wind_change_frames(4),floor(wind_change_frames(4)+0.5*both_cues_2_full)];
    boundaries{7} = [floor(wind_change_frames(4)+0.5*both_cues_2_full),panels_change_frames(4)];
    boundaries{8} = [panels_change_frames(4),panels_change_frames(5)];
    boundaries{9} = [panels_change_frames(5),wind_change_frames(6)];
    both_cues_3_full = panels_change_frames(6) - wind_change_frames(6);
    boundaries{10} = [wind_change_frames(6),floor(wind_change_frames(6)+0.5*both_cues_3_full)];
    boundaries{11} = [floor(wind_change_frames(6)+0.5*both_cues_3_full),panels_change_frames(6)];
    boundaries{12} = [panels_change_frames(6),panels_change_frames(7)];
    boundaries{13} = [panels_change_frames(7),wind_change_frames(8)];
    both_cues_4_full = panels_change_frames(8) - wind_change_frames(8);
    boundaries{14} = [wind_change_frames(8),floor(wind_change_frames(8)+0.5*both_cues_4_full)];
    boundaries{15} = [floor(wind_change_frames(8)+0.5*both_cues_4_full),panels_change_frames(8)];
    boundaries{16} = [panels_change_frames(8),panels_change_frames(9)];
    if fly_ID > 3
        boundaries{17} = [panels_change_frames(9),panels_change_frames(10)];
    else
        boundaries{17} = [panels_change_frames(9),length(panels_on)];
    end
        
end

%% Import wind data for the first 3 flies

%import hdf5 data
if fly_ID < 4
    if fly_ID ==1
        motorpos = h5read([path(1:end-10),'\ball\hdf5_Bar_wind_jump_panels_Closed_Loop_X_Closed_Loop_Y_wind_Closed_Loop_20211026_120722_sid_1_tid_1_arduino000001.hdf5'],'/motor');
    elseif fly_ID == 2
        motorpos = h5read([path(1:end-10),'\ball\hdf5_Bar_wind_jump_panels_Closed_Loop_X_Closed_Loop_Y_wind_Closed_Loop_20211026_134625_sid_0_tid_1_arduino000001.hdf5'],'/motor');
    elseif fly_ID == 3
        motorpos = h5read([path(1:end-10),'\ball\hdf5_Bar_wind_jump_panels_Closed_Loop_X_Closed_Loop_Y_wind_Closed_Loop_20211026_145901_sid_0_tid_1_arduino000001.hdf5'],'/motor');
    end
    motor_pos_un = motorpos(round(linspace(1, length(motorpos), length(continuous_data.motor_pos))));
    motor_pos_un = wrapTo180(rad2deg(motor_pos_un));
    
    %compute the cross-correlation with respect to the bar, to adjust the time
    %a bit for it to match at its best
    bar_pos = wrapTo180(continuous_data.panel_angle);
    
    %Uncomment to see what the uncorrected wind data looks like
%     figure,
%     plot(bar_pos)
%     hold on
%     plot(motor_pos_un)
%     legend('Bar position', 'Uncorrected wind position');
    
    [c,lags] = xcorr(motor_pos_un,bar_pos);
    %Uncomment to see the correlation at different lags.
%     figure,
%     stem(lags,c)
%     xlabel('Lag (frames)'); ylabel('Correlation');
    %find max correlation lag
    [~, opt_corr] = max(c);
    opt_lag = lags(opt_corr);
    
    motor_pos = [motor_pos_un(opt_lag+1:end);repelem(NaN,opt_lag,1)];
    motor_pos = motor_pos';
    
    %figure to check if the lag is correct
%     figure,
%     plot(bar_pos)
%     hold on
%     plot(motor_pos)
%     legend('bar position', 'wind position');
    
else
    motor_pos = wrapTo180(rad2deg(continuous_data.motor_pos));
    
end


%% Find the jump frames

if configuration == 1
   %Find the bar jump frames
   coded_bar_jump_frames = [floor(1380*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2340*length(continuous_data.dff_matrix)/continuous_data.time(end))];
   %Find the wind jump frames
   coded_wind_jump_frames = floor(1860*length(continuous_data.dff_matrix)/continuous_data.time(end));   
else
   %Find the wind jump frames
   coded_wind_jump_frames = [floor(1380*length(continuous_data.dff_matrix)/continuous_data.time(end)),floor(2340*length(continuous_data.dff_matrix)/continuous_data.time(end))];
   %Find the bar jump frames
   coded_bar_jump_frames = floor(1860*length(continuous_data.dff_matrix)/continuous_data.time(end));    
end

%Correct bar jump frames
abs_diff_bar_signal = abs(diff(unwrap(deg2rad(continuous_data.panel_angle))));
for jump = 1:length(coded_bar_jump_frames)
    frame_segment = [coded_bar_jump_frames(jump)-100:coded_bar_jump_frames(jump)+100];
    [~,I_bar_jump_frames(jump)] = max(abs_diff_bar_signal(frame_segment));
    real_bar_jump_frames(jump) = coded_bar_jump_frames(jump) + I_bar_jump_frames(jump) - 100;
end
real_bar_jump_sec = real_bar_jump_frames*continuous_data.time(end)/length(continuous_data.dff_matrix);

%Correct wind jump frames
abs_diff_wind_signal = abs(diff(unwrap(deg2rad(motor_pos))));
for jump = 1:length(coded_wind_jump_frames)
    frame_segment = [coded_wind_jump_frames(jump)-100:coded_wind_jump_frames(jump)+100];
    [~,I_wind_jump_frames(jump)] = max(abs_diff_wind_signal(frame_segment));
    real_wind_jump_frames(jump) = coded_wind_jump_frames(jump) + I_wind_jump_frames(jump) - 100;
end
%correct for the flies for which the method didn't work well
if fly_ID == 4
    real_wind_jump_frames(1) = 12592;
elseif fly_ID ==5
    real_wind_jump_frames(1) = floor(1851*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_bar_jump_frames(2) = floor(2331*length(continuous_data.dff_matrix)/continuous_data.time(end));
elseif fly_ID == 6
    real_wind_jump_frames(1) = floor(1372*length(continuous_data.dff_matrix)/continuous_data.time(end));
    real_wind_jump_frames(2) = floor(2332*length(continuous_data.dff_matrix)/continuous_data.time(end));
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
if configuration == 1
    xline(real_bar_jump_frames(1),'color',[0.4940 0.1840 0.5560],'linewidth',4)
    xline(real_bar_jump_frames(2),'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
    xline(real_wind_jump_frames,'color',[0.4660 0.6740 0.1880],'linewidth',4)
    legend('Bar jumps','Wind jump');
else
    xline(real_wind_jump_frames(1),'color',[0.4660 0.6740 0.1880],'linewidth',4)
    xline(real_wind_jump_frames(2),'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
    xline(real_bar_jump_frames,'color',[0.4940 0.1840 0.5560],'linewidth',4)
    legend('Wind jumps','Bar jump');
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
if configuration == 1
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4)
else
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4)
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
if configuration == 1
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4)
else
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4)
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
if configuration == 1
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4)
else
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4)
end
title('Bump magnitude')
xlim([0 continuous_data.time(end)]);
set(gca,'xticklabel',{[]})

subplot(5,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.','handleVisibility','off')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.','handleVisibility','off')
%add jumps
if configuration == 1
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4,'handlevisibility','off')
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4)
else
    xline(real_wind_jump_sec,'color',[0.4660 0.6740 0.1880],'linewidth',4,'handlevisibility','off')
    xline(real_bar_jump_sec,'color',[0.4940 0.1840 0.5560],'linewidth',4)
end
xlabel('Time (sec)');
title('Bump width');
xlim([0 continuous_data.time(end)]);

saveas(gcf,[path,'\plots\full_experiment.png']);
% 
% %% Check the boundaries
% 
% figure('Position',[100 100 1400 800]),
% subplot(2,1,1)
% plot(x_out_bar,bar_to_plot,'.','LineWidth',1.5)
% hold on
% plot(x_out_motor,motor_to_plot,'.','LineWidth',1.5)
% %add boundaries
% for bout = 1:17
%    xline(frames_to_sec*boundaries{bout}(1),'b','linewidth',2);
% end
% title('Start of bout');
% legend('Bar position','Wind position')
% ylim([-180 180]);
% xlim([0 x_out_bump(end)]);
% set(gca,'xticklabel',{[]})
% 
% subplot(2,1,2)
% plot(x_out_bar,bar_to_plot,'.','LineWidth',1.5)
% hold on
% plot(x_out_motor,motor_to_plot,'.','LineWidth',1.5)
% %add boundaries
% for bout = 1:17
%    xline(frames_to_sec*boundaries{bout}(2),'r','linewidth',2);
% end
% title('End of bout');
% legend('Bar position','Wind position')
% ylim([-180 180]);
% xlim([0 x_out_bump(end)]);
% set(gca,'xticklabel',{[]})

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


%% Look at offset variability of last darkness bout

if fly_ID >3
    
    figure,
    polarhistogram(deg2rad(bar_offset(~panels_on & ~wind_on)),15,'FaceColor','k')
    set(gca,'ThetaZeroLocation','top');
    title('No stimulus');
    saveas(gcf,[path,'\plots\no_stim_offset.png']);

end


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

for i = 1:length(boundaries)
    subplot(4,length(boundaries),i)
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
    
    
    subplot(4,length(boundaries),i+17)
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
subplot(4,length(boundaries),[35:68])
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
xlim([0.6 17.4]);
legend('Bump mag','Bump width');

saveas(gcf,[path,'\plots\bump_parameters_and_offset_per_block_mvt_and_fit_thresh.png']);


%% Bump parameter distribution thresholded with movement and fit data

%Bump magnitude per section
bar_only_BM_mvt_thresh = continuous_data.bump_magnitude(panels_on & ~wind_on & moving' & good_fit');
wind_only_BM_mvt_thresh = continuous_data.bump_magnitude(~panels_on & wind_on & moving' & good_fit');
both_cues_BM_mvt_thresh = continuous_data.bump_magnitude(panels_on & wind_on & moving' & good_fit');
allData_mvt_thresh = [bar_only_BM_mvt_thresh,wind_only_BM_mvt_thresh,both_cues_BM_mvt_thresh];
stim_ID_mvt_thresh = [repelem(1,1,length(bar_only_BM_mvt_thresh)),repelem(2,1,length(wind_only_BM_mvt_thresh)),repelem(3,1,length(both_cues_BM_mvt_thresh))];
alldata_BM_mvt_thresh = table(allData_mvt_thresh,stim_ID_mvt_thresh);

figure('Position',[100 100 1200 800]),
subplot(1,2,1)
boxplot(alldata_BM_mvt_thresh.allData_mvt_thresh,alldata_BM_mvt_thresh.stim_ID_mvt_thresh,'color','k')
set(findobj(gca,'type','line'),'linew',2)
ylim([0 4]);
xticks([1 2 3])
xticklabels({'Bar only','Wind only','Both cues'})
ylabel('Bump magnitude');


% Bump width per section
bar_only_BW_mvt_thresh = continuous_data.bump_width(panels_on & ~wind_on & moving' & good_fit');
wind_only_BW_mvt_thresh = continuous_data.bump_width(~panels_on & wind_on & moving' & good_fit');
both_cues_BW_mvt_thresh = continuous_data.bump_width(panels_on & wind_on & moving' & good_fit');
allData_BW_mvt_thresh = [bar_only_BW_mvt_thresh,wind_only_BW_mvt_thresh,both_cues_BW_mvt_thresh];
alldata_BW_mvt_thresh = table(allData_BW_mvt_thresh,stim_ID_mvt_thresh);

subplot(1,2,2)
boxplot(alldata_BW_mvt_thresh.allData_BW_mvt_thresh,alldata_BW_mvt_thresh.stim_ID_mvt_thresh,'color','k')
set(findobj(gca,'type','line'),'linew',2)
ylim([0 5]);
xticks([1 2 3])
xticklabels({'Bar only','Wind only','Both cues'})
ylabel('Bump width');

saveas(gcf,[path,'\plots\bump_parameters_mvt_and_fit_thresh.png']);

%% Bump parameter distribution, dividing the cue combination bouts into pre and post jump

%pre jump frames
both_cues_pre_jump_frames = boundaries([block{:}]==3);
both_cues_pre_jump = [both_cues_pre_jump_frames{1}(1):both_cues_pre_jump_frames{1}(2),both_cues_pre_jump_frames{2}(1):both_cues_pre_jump_frames{2}(2),both_cues_pre_jump_frames{3}(1):both_cues_pre_jump_frames{3}(2),both_cues_pre_jump_frames{4}(1):both_cues_pre_jump_frames{4}(2)];
both_cues_pre = zeros(length(panels_on),1);
both_cues_pre(both_cues_pre_jump) = 1;

%post jump frames
both_cues_post_jump_frames = boundaries([block{:}]==4);
both_cues_post_jump = [both_cues_post_jump_frames{1}(1):both_cues_post_jump_frames{1}(2),both_cues_post_jump_frames{2}(1):both_cues_post_jump_frames{2}(2),both_cues_post_jump_frames{3}(1):both_cues_post_jump_frames{3}(2)];
both_cues_post = zeros(length(panels_on),1);
both_cues_post(both_cues_post_jump) = 1;

both_cues_BM_mvt_thresh_pre = continuous_data.bump_magnitude(panels_on & wind_on & moving' & good_fit' & both_cues_pre);
both_cues_BM_mvt_thresh_post = continuous_data.bump_magnitude(panels_on & wind_on & moving' & good_fit' & both_cues_post);
allData_thresh = [bar_only_BM_mvt_thresh,wind_only_BM_mvt_thresh,both_cues_BM_mvt_thresh_pre,both_cues_BM_mvt_thresh_post];
stim_ID_thresh = [repelem(1,1,length(bar_only_BM_mvt_thresh)),repelem(2,1,length(wind_only_BM_mvt_thresh)),repelem(3,1,length(both_cues_BM_mvt_thresh_pre)),repelem(4,1,length(both_cues_BM_mvt_thresh_post))];
alldata_BM_thresh = table(allData_thresh,stim_ID_thresh);

both_cues_BW_mvt_thresh_pre = continuous_data.bump_width(panels_on & wind_on & moving' & good_fit' & both_cues_pre);
both_cues_BW_mvt_thresh_post = continuous_data.bump_width(panels_on & wind_on & moving' & good_fit' & both_cues_post);
allDataBW_thresh = [bar_only_BW_mvt_thresh,wind_only_BW_mvt_thresh,both_cues_BW_mvt_thresh_pre,both_cues_BW_mvt_thresh_post];
alldata_BW_thresh = table(allDataBW_thresh,stim_ID_thresh);


figure('Position',[100 100 1200 800]),
subplot(1,2,1)
boxplot(alldata_BM_thresh.allData_thresh,alldata_BM_thresh.stim_ID_thresh,'color','k')
set(findobj(gca,'type','line'),'linew',2)
ylim([0 5]);
xticks([1 2 3 4])
xticklabels({'Bar only','Wind only','Both cues pre','Both cues post'})
ylabel('Bump magnitude');

subplot(1,2,2)
boxplot(alldata_BW_thresh.allDataBW_thresh,alldata_BW_thresh.stim_ID_thresh,'color','k')
set(findobj(gca,'type','line'),'linew',2)
ylim([0 5]);
xticks([1 2 3 4])
xticklabels({'Bar only','Wind only','Both cues pre','Both cues post'})
ylabel('Bump width');

saveas(gcf,[path,'\plots\bump_parameters_mvt_and_fit_thresh_four_blocks.png']);

%% Compute the wind plasticity 
%We will quantify the plasticity here as the difference in the mean offset
%between two blocks with only that cue

wind_blocks = find([block{:}] == 2);

for bout = 1:length(wind_blocks)-1
    
    initial_windoff = deg2rad(wind_offset(boundaries{wind_blocks(bout)}(1):boundaries{wind_blocks(bout)}(2)));
    i_thresh_wind_off = initial_windoff(continuous_data.total_mvt_ds(boundaries{wind_blocks(bout)}(1):boundaries{wind_blocks(bout)}(2))>25);
    final_windoff = deg2rad(wind_offset(boundaries{wind_blocks(bout+1)}(1):boundaries{wind_blocks(bout+1)}(2)));
    f_thresh_wind_off = final_windoff(continuous_data.total_mvt_ds(boundaries{wind_blocks(bout+1)}(1):boundaries{wind_blocks(bout+1)}(2))>25);
    
    wind_offset_diff(bout) = rad2deg(circ_dist(circ_mean(f_thresh_wind_off(~isnan(f_thresh_wind_off)),[],2),circ_mean(i_thresh_wind_off(~isnan(i_thresh_wind_off)),[],2)));

    figure,
    subplot(2,2,1)
    plot(rad2deg(i_thresh_wind_off),'color',[0.9290 0.6940 0.1250])
    title(['Initial wind offset = ',num2str(round(rad2deg(circ_mean(i_thresh_wind_off(~isnan(i_thresh_wind_off)),[],2)))),'deg']);
    
    subplot(2,2,2)
    polarhistogram(i_thresh_wind_off,15,'FaceColor',[0.9290 0.6940 0.1250])
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    subplot(2,2,3)
    plot(rad2deg(f_thresh_wind_off),'color',[0.9290 0.6940 0.1250])
    title(['Final wind offset = ',num2str(round(rad2deg(circ_mean(f_thresh_wind_off(~isnan(f_thresh_wind_off)),[],2)))),'deg']);
    
    subplot(2,2,4)
    polarhistogram(f_thresh_wind_off,15,'FaceColor',[0.9290 0.6940 0.1250])
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    suptitle(['Wind offset difference = ',num2str(wind_offset_diff(bout)),' deg']);

    saveas(gcf,[path,'\plots\wind_offset_diff_block',num2str(bout),'.png']);

end


%% Link wind plasticity to bump parameters

%problem: which period to choose, pre-jump and/or post-jump?
%second problem: what about bouts in which the offset toggles between two?
for bout = 1:length(wind_blocks)-1
    if configuration == 1
        %Find the bump parameters in the preceding bout (post jump cue combination)
        pre_BM_postjump(bout) = meanBM_perblock_thresh(wind_blocks(bout+1)-2);
        pre_BW_postjump(bout) = meanBW_perblock_thresh(wind_blocks(bout+1)-2);        
    else
        %Find the bump parameters in the preceding bout (post jump cue combination)
        pre_BM_postjump(bout) = meanBM_perblock_thresh(wind_blocks(bout+1)-1);
        pre_BW_postjump(bout) = meanBW_perblock_thresh(wind_blocks(bout+1)-1);
    end
end

%Scatter plot of wind plasticity and bump parameters
figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(pre_BM_postjump,abs(wind_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'post jump'}); ylabel('Wind offset difference');
ylim([0 180]); xlim([0 3]);

subplot(1,2,2)
plot(pre_BW_postjump,abs(wind_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'post jump'});
ylim([0 180]); xlim([0 3]);

saveas(gcf,[path,'\plots\wind_offset_diff_and_bump_par_postjump.png']);

%% Link the plasticity to the bump parameters in the pre-jump cue combination portion

%for the first bout, take the cue combination portion
if configuration == 1
    pre_BM_prejump(1) = meanBM_perblock_thresh(wind_blocks(bout+1)-2);
    pre_BW_prejump(1) = meanBW_perblock_thresh(wind_blocks(bout+1)-2);            
else
    pre_BM_prejump(1) = meanBM_perblock_thresh(wind_blocks(bout+1)-1);
    pre_BW_prejump(1) = meanBW_perblock_thresh(wind_blocks(bout+1)-1);            
end
for bout = 2:length(wind_blocks)-1
    if configuration == 1
        %Find the bump parameters in the preceding bout (pre jump cue combination)
        pre_BM_prejump(bout) = meanBM_perblock_thresh(wind_blocks(bout+1)-3);
        pre_BW_prejump(bout) = meanBW_perblock_thresh(wind_blocks(bout+1)-3);        
    else
        %Find the bump parameters in the preceding bout (pre jump cue combination)
        pre_BM_prejump(bout) = meanBM_perblock_thresh(wind_blocks(bout+1)-2);
        pre_BW_prejump(bout) = meanBW_perblock_thresh(wind_blocks(bout+1)-2);
    end
end

%Scatter plot of wind plasticity and bump parameters
figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(pre_BM_prejump,abs(wind_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'pre jump'}); ylabel('Wind offset difference');
ylim([0 180]); xlim([0 3]);

subplot(1,2,2)
plot(pre_BW_prejump,abs(wind_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'pre jump'});
ylim([0 180]); xlim([0 3]);

saveas(gcf,[path,'\plots\wind_offset_diff_and_bump_par_prejump.png']);

%% Compute the bar plasticity 
%We will quantify the plasticity here as the difference in the mean offset
%between two blocks with only that cue

bar_blocks = find([block{:}] == 1);

for bout = 1:length(bar_blocks)-1
    
    initial_baroff = deg2rad(bar_offset(boundaries{bar_blocks(bout)}(1):boundaries{bar_blocks(bout)}(2)));
    i_thresh_bar_off = initial_baroff(continuous_data.total_mvt_ds(boundaries{bar_blocks(bout)}(1):boundaries{bar_blocks(bout)}(2))>25);
    final_baroff = deg2rad(bar_offset(boundaries{bar_blocks(bout+1)}(1):boundaries{bar_blocks(bout+1)}(2)));
    f_thresh_bar_off = final_baroff(continuous_data.total_mvt_ds(boundaries{bar_blocks(bout+1)}(1):boundaries{bar_blocks(bout+1)}(2))>25);
    
    bar_offset_diff(bout) = rad2deg(circ_dist(circ_mean(f_thresh_bar_off(~isnan(f_thresh_bar_off))),circ_mean(i_thresh_bar_off(~isnan(i_thresh_bar_off)))));

    figure,
    subplot(2,2,1)
    plot(rad2deg(i_thresh_bar_off),'color',[0.8500 0.3250 0.0980])
    title(['Initial bar offset = ',num2str(round(rad2deg(circ_mean(i_thresh_bar_off(~isnan(i_thresh_bar_off)))))),'deg']);
    
    subplot(2,2,2)
    polarhistogram(i_thresh_bar_off,15,'FaceColor',[0.8500 0.3250 0.0980])
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    subplot(2,2,3)
    plot(rad2deg(f_thresh_bar_off),'color',[0.8500 0.3250 0.0980])
    title(['Final wind offset = ',num2str(round(rad2deg(circ_mean(f_thresh_bar_off(~isnan(f_thresh_bar_off)))))),'deg']);
    
    subplot(2,2,4)
    polarhistogram(f_thresh_bar_off,15,'FaceColor',[0.8500 0.3250 0.0980])
    set(gca,'ThetaZeroLocation','top');
    set(gca,'Thetaticklabels',{});
    set(gca,'Rticklabels',{});
    
    suptitle(['Bar offset difference = ',num2str(bar_offset_diff(bout)),' deg']);

    saveas(gcf,[path,'\plots\bar_offset_diff_block',num2str(bout),'.png']);

end


%% Link bar plasticity to bump parameters

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(pre_BM_postjump,abs(bar_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'post jump'}); ylabel('Bar offset difference');
ylim([0 180]); xlim([0 3]);

subplot(1,2,2)
plot(pre_BW_postjump,abs(bar_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'post jump'});
ylim([0 180]); xlim([0 3]);

saveas(gcf,[path,'\plots\bar_offset_diff_and_bump_par_postjump.png']);

%% Repeat focusing on the pre jump parameters

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(pre_BM_prejump,abs(bar_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'pre jump'}); ylabel('Bar offset difference');
ylim([0 180]); xlim([0 3]);

subplot(1,2,2)
plot(pre_BW_prejump,abs(bar_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'pre jump'});
ylim([0 180]); xlim([0 3]);

saveas(gcf,[path,'\plots\bar_offset_diff_and_bump_par_prejump.png']);


%% Plot the wind plasticity as a function of the difference between the bar and wind offset

%% Plot the bar plasticity as a function of the difference between the bar and wind offset


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

saveas(gcf,[path,'\plots\Exp35_rep.png']);

%% Save data

save([path,'data.mat'],'meanBM_perblock_thresh','mean_total_mvt','meanBW_perblock_thresh','alldata_BM_thresh','alldata_BW_thresh','stim_ID_thresh','pre_BM_prejump','pre_BM_postjump','pre_BW_prejump','bar_offset_diff','wind_offset_diff')

%% Close

clear all; close all;