%analysis of the cue conflict session

%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp37\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])
load([path,'\analysis\continuous_analysis_sid_',num2str(sessions.cue_conflict),'_tid_0.mat'])


%% Make directory to save plots

%Move to the analysis folder
cd([path,'\analysis\'])
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(path,'\analysis\plots'); 
end
%List the contents of the 'plots' folder
cd([path,'\analysis\plots\'])

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
    if strcmp(flyNames(fly).name,path(54:end))
        fly_ID = fly;
    end
end

colors_for_plots = [0.2 0.8 0.8 ; 1 0.5 0; 0 0.5 1;...
    0 0.6 0.3;  1 0.2 0.2; 0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880;...
    0 0.4470 0.7410; 0.75, 0.1, 0.75; 0.75, 0.75, 0];

fly_color = colors_for_plots(fly_ID,:);

%% Get the jumps

%Use the y dimension of the panels to find the jumps, taking its derivative
jumps = diff(continuous_data.fr_y_ds);
jumps(abs(jumps)>0.4 & abs(jumps)<1)=1;
jumps = round(jumps);

figure,
suptitle('Bar jumps');
subplot(3,1,1)
plot(continuous_data.fr_y_ds)
ylabel('Voltage (V)');xlabel('Time');
subplot(3,1,2)
plot(jumps);
ylabel('Voltage difference (V)');xlabel('Time');

j = find(jumps); %indices of the actual bar jumps, taken from the y signal

%plot the data from the yPanels and add lines of the previously determined
%bar jumps
subplot(3,1,3),
plot(continuous_data.fr_y_ds)
title('Bar jumps');
xlabel('Time (frames)'); ylabel('Voltage (V)');
hold on
%add the bar jumps
for i = 1:length(j)
    plot([j(i) j(i)],[0 10],'r');
end

%% Plot activity heatmap with bump parameters

figure('Position',[100 100 1200 800]),
subplot(6,1,1)
dff = continuous_data.dff_matrix';
imagesc(flip(dff))
colormap(flipud(gray))
title('EPG activity');
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
%add bar jumps
for jump = 1:length(j)
    xline(j(jump),'color',[0.4660, 0.6740, 0.1880],'linewidth',3)
end

subplot(6,1,2)
bump_pos = wrapTo180(rad2deg(continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
%bar position
bar_pos = wrapTo180(-continuous_data.panel_angle);
[x_out_bar,bar_to_plot] = removeWrappedLines(continuous_data.time,bar_pos);
plot(x_out_bar,bar_to_plot,'LineWidth',1.5)
%wind position
wind_pos = wrapTo180(rad2deg(continuous_data.motor_pos));
[x_out_wind,wind_to_plot] = removeWrappedLines(continuous_data.time,wind_pos');
plot(x_out_wind,wind_to_plot,'LineWidth',1.5)
title('Bump and stim position');
legend('Bump estimate','Bar position','Wind position','Location','best')
ylim([-180 180]);
set(gca,'xticklabel',{[]})
%add bar jumps
for jump = 1:length(j)
    xline(continuous_data.time(j(jump)),'color',[0.4660, 0.6740, 0.1880],'linewidth',3,'HandleVisibility','off')
end

%offset with respect to the bar
subplot(6,1,3)
bar_offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',deg2rad(-continuous_data.panel_angle))));
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,bar_offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Bar offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})
%add bar jumps
for jump = 1:length(j)
    xline(continuous_data.time(j(jump)),'color',[0.4660, 0.6740, 0.1880],'linewidth',3)
end

%offset with respect to the wind
subplot(6,1,4)
wind_offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',continuous_data.motor_pos')));
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,wind_offset);
plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
title('Wind offset')
ylim([-180 180]);
set(gca,'xticklabel',{[]})
%add bar jumps
for jump = 1:length(j)
    xline(continuous_data.time(j(jump)),'color',[0.4660, 0.6740, 0.1880],'linewidth',3)
end

subplot(6,1,5)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'k.')
title('Bump magnitude')
%add bar jumps
for jump = 1:length(j)
    xline(continuous_data.time(j(jump)),'color',[0.4660, 0.6740, 0.1880],'linewidth',3)
end

subplot(6,1,6)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'r.')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'k.')
title('Bump width')
%add bar jumps
for jump = 1:length(j)
    xline(continuous_data.time(j(jump)),'color',[0.4660, 0.6740, 0.1880],'linewidth',3)
end

%save plot
saveas(gcf,[path,'\analysis\plots\summary_results.png']);

%% Compare bar and wind offset variability

figure,
subplot(121)
polarhistogram(deg2rad(bar_offset),20,'FaceColor',fly_color)
set(gca,'ThetaZeroLocation','top');
Ax = gca; 
Ax.RGrid = 'off';
Ax.RTickLabel = [];
title('Bar offset')

subplot(122)
polarhistogram(deg2rad(wind_offset),20,'FaceColor',fly_color)
set(gca,'ThetaZeroLocation','top');
Ax = gca; 
Ax.RGrid = 'off';
Ax.RTickLabel = [];
title('Wind offset')

%save plot
saveas(gcf,[path,'\analysis\plots\offsets.png']);

%% Close up around the jumps

for jump = 1:length(j)
   
   time_zero = continuous_data.time(j(jump));
   time = continuous_data.time-time_zero;
   
   figure('Position',[100 100 1200 500]),
   ax(1) = subplot(8,1,1);
   imagesc(continuous_data.bump_magnitude(:,j(jump)-25:j(jump)+26))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   ax(2) = subplot(8,1,[2 8]);
   time_to_plot = time(j(jump)-25:j(jump)+26);
   phase_to_plot = bump_pos(j(jump)-25:j(jump)+26);
   [x_out_time_bump,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot');   
   plot(x_out_time_bump,bump_pos_to_plot,'color',[0 0.6 0.3],'linewidth',2)
   hold on
   bar_to_plot = bar_pos(j(jump)-25:j(jump)+26);
   [x_out_time_bar,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time_bar,bar_pos_to_plot,'color',[1 0.2 0.2],'linewidth',2)
   wind_to_plot = wind_pos(j(jump)-25:j(jump)+26);
   [x_out_time_wind,wind_pos_to_plot] = removeWrappedLines(time_to_plot,wind_to_plot');
   plot(x_out_time_wind,wind_pos_to_plot,'color',[0.2 0.2 0.7],'linewidth',2)
   xline(time(j(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(j(jump)-25) time(j(jump)+26)]);
   ylabel('Heading (deg)');
   xlabel('Time (sec)');
   legend('Bump estimate','Bar position','Wind position');
   
end
   

%% Bump magnitude around the jumps

for jump = 1:length(j)
    
   figure('Position',[100 100 1200 800]),
     
   subplot(3,1,1)
   dff_matrix = continuous_data.dff_matrix';
   imagesc(flip(dff_matrix(:,j(jump)-25:j(jump)+25)))
   colormap(flipud(gray))
   hold on
   xline(26,'linewidth',2,'color','r')
   title('EPG activity');
     
   subplot(3,1,2)
   plot(bar_offset(j(jump)-25:j(jump)+25))
   hold on
   plot(wind_offset(j(jump)-25:j(jump)+25))
   legend('Bar offset','Wind offset')
   xline(26,'linewidth',2,'color','r','HandleVisibility','off')
   title(['Jump #',num2str(jump)]); 
   xlim([1 51]);
   ylim([-180 180]);
   title('Offset (deg)');
   
   subplot(3,1,3)
   %compute the mean bump magnitude around the jumps
   mean_bm_pre_jump(jump) = mean(continuous_data.bump_magnitude(j(jump)-25:j(jump)));
   mean_bm_post_jump(jump) = mean(continuous_data.bump_magnitude(j(jump)+1:j(jump)+26));
   bm_around_jump{jump} = [mean_bm_pre_jump(jump),mean_bm_post_jump(jump)];
   plot(bm_around_jump{jump},'-ko')
   xlim([0 3]);
   set(gca,'xticklabel',{[]})
   ylim([min(bm_around_jump{jump})-0.5, max(bm_around_jump{jump})+0.5]);  
   title('Bump magnitude');
   
   saveas(gcf,[path,'\analysis\plots\AJ_bump_magnitude_',num2str(jump),'.png']);
   
end


%% Bump width around the jumps

for jump = 1:length(j)
    
   figure('Position',[100 100 1200 800]),
     
   subplot(3,1,1)
   dff_matrix = continuous_data.dff_matrix';
   imagesc(flip(dff_matrix(:,j(jump)-25:j(jump)+25)))
   colormap(flipud(gray))
   hold on
   xline(26,'linewidth',2,'color','r')
   title('EPG activity');
     
   subplot(3,1,2)
   plot(bar_offset(j(jump)-25:j(jump)+25))
   hold on
   plot(wind_offset(j(jump)-25:j(jump)+25))
   legend('Bar offset','Wind offset')
   xline(26,'linewidth',2,'color','r','HandleVisibility','off')
   title(['Jump #',num2str(jump)]); 
   xlim([1 51]);
   ylim([-180 180]);
   title('Offset (deg)');
   
   subplot(3,1,3)
   %compute the mean bump width around the jump
   mean_bw_pre_jump(jump) = mean(continuous_data.bump_width(j(jump)-25:j(jump)));
   mean_bw_post_jump(jump) = mean(continuous_data.bump_width(j(jump)+1:j(jump)+26));
   bw_around_jump{jump} = [mean_bw_pre_jump(jump),mean_bw_post_jump(jump)];
   plot(bw_around_jump{jump},'-ko')
   xlim([0 3]);
   set(gca,'xticklabel',{[]})
   ylim([min(bw_around_jump{jump})-0.5, max(bw_around_jump{jump})+0.5]);  
   title('Bump width');
   
   saveas(gcf,[path,'\analysis\plots\AJ_bump_width_',num2str(jump),'.png']);
   
end

%% Looking at both bump parameters combined

bm_aj = cell2mat(bm_around_jump);
bm_aj = reshape(bm_aj,length(bm_aj)/length(bm_around_jump),length(bm_around_jump));
bw_aj = cell2mat(bw_around_jump);
bw_aj = reshape(bw_aj,length(bw_aj)/length(bw_around_jump),length(bw_around_jump));

figure,
subplot(121)
plot(bm_aj,'-o','color',[.5 .5 .5])
title('Bump magnitude')
xlim([0 3]); ylim([0 3]);
hold on
plot(median(bm_aj,2),'-ko','linewidth',2)
xticks([1 2]);
xticklabels({'pre-jump','post-jump'});

subplot(1,2,2)
plot(bw_aj,'-o','color',[.5 .5 .5])
title('Bump width')
xlim([0 3]); ylim([0 4]);
hold on
plot(median(bw_aj,2),'-ko','linewidth',2)
xticks([1 2]);
xticklabels({'pre-jump','post-jump'});

saveas(gcf,[path,'\analysis\plots\AJ_bump_parameters.png']);

%% Look at angular velocity around the jumps

for jump = 1:length(j)
    
   figure('Position',[100 100 1200 800]),
     
   plot(bar_pos(j(jump)-25:j(jump)+25))
   hold on
   vel_yaw{jump} = -continuous_data.vel_yaw_ds(j(jump)-25:j(jump)+25);
   plot(vel_yaw{jump})
   legend('Bar position','Angular velocity')
   xline(26,'linewidth',2,'color','r','HandleVisibility','off')
   title(['Jump #',num2str(jump)]); 
   xlim([1 51]);
   ylim([-180 180]);  
      
   %saveas(gcf,[path,'\analysis\plots\AJ_yaw_vel_',num2str(jump),'.png']);
   
end


jumps = [-120,120,-120,-120,120,120,-120,120,120,-120];

for jump = 1:length(jumps)
    if jumps(jump) == 120
        vel_yaw_120{jump} = vel_yaw{jump};
    else
        vel_yaw_neg120{jump} = vel_yaw{jump};
    end
end
vel_yaw_120 =  vel_yaw_120(~cellfun('isempty',vel_yaw_120));
vel_yaw120 = cell2mat(vel_yaw_120);
vel_yaw120 = reshape(vel_yaw120,length(vel_yaw_120{1}),length(vel_yaw_120));
vel_yaw_neg120 =  vel_yaw_neg120(~cellfun('isempty',vel_yaw_neg120));
vel_yawneg120 = cell2mat(vel_yaw_neg120);
vel_yawneg120 = reshape(vel_yawneg120,length(vel_yaw_neg120{1}),length(vel_yaw_neg120));

figure,
plot(mean(vel_yaw120,2),'r')
hold on
plot(mean(vel_yawneg120,2),'b')
xline(26,'linewidth',2,'color','k','HandleVisibility','off')
yline(0,'linewidth',2,'color','k','linestyle','--','HandleVisibility','off')
xlim([1 51]);
ylim([-100 100]);
legend('Positive bar jumps','Negative bar jumps');
set(gca,'xticklabels',{[]});
title('Angular velocity around the jumps');
ylabel('Angular velocity (deg/s)');
saveas(gcf,[path,'\analysis\plots\AJ_mean_yaw_vel.png']);

%% Clear 
clear all; close all; clc;