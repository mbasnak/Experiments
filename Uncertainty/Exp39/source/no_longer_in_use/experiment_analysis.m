%% Analysis for the EB-DAN experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[path] = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp39\data\');

% Import sessions information
load([path,'\analysis\sessions_info.mat'])

%% Set colormap

folderNames = dir(path(1:53));
flyNames = struct();
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'75B10') & ~contains(folderNames(folder).name,'txt'))
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

%% Analyze initial closed-loop panels

load([path,'\analysis\analysis_sid_',num2str(sessions.initial_cl_bar),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff{1};
plot(dff);
title('EB-DAN activity');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,2)
plot(data.vel_for_ds);
title('Forward velocity (mm/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
plot(abs(data.vel_yaw_ds));
title('Yaw speed (deg/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.total_mvt_ds);
xlabel('Time (sec)');
title('Total velocity (deg/s)');
xlim([0 data.time(end)]);
set(gca,'xticklabel',{[]})

suptitle('Initial trial with just panels');

saveas(gcf,[path,'\analysis\plots\initial_panels.png']);

%Get mean parameters
mean_EBDAN_act(1) = nanmean(dff);
mean_total_mvt(1) = nanmean(data.total_mvt_ds);

%% Analyze initial closed-loop wind

load([path,'\analysis\analysis_sid_',num2str(sessions.initial_cl_wind),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff{1};
plot(dff);
title('EB-DAN activity');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,2)
plot(data.vel_for_ds);
title('Forward velocity (mm/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
plot(abs(data.vel_yaw_ds));
title('Yaw speed (deg/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.total_mvt_ds);
xlabel('Time (sec)');
title('Total velocity (deg/s)');
xlim([0 data.time(end)]);
set(gca,'xticklabel',{[]})

suptitle('Initial trial with just wind');

saveas(gcf,[path,'\analysis\plots\initial_wind.png']);

%Get mean parameters
mean_EBDAN_act(2) = nanmean(dff);
mean_total_mvt(2) = nanmean(data.total_mvt_ds);

%% Analyze cue combination trial

load([path,'\analysis\analysis_sid_',num2str(sessions.cue_combination),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff{1};
plot(dff);
title('EB-DAN activity');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,2)
plot(data.vel_for_ds);
title('Forward velocity (mm/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
plot(abs(data.vel_yaw_ds));
title('Yaw speed (deg/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.total_mvt_ds);
xlabel('Time (sec)');
title('Total velocity (deg/s)');
xlim([0 data.time(end)]);
set(gca,'xticklabel',{[]})

suptitle('Cue combination trial');

saveas(gcf,[path,'\analysis\plots\cue_combination.png']);

%Get mean parameters
mean_EBDAN_act(3) = nanmean(dff);
mean_total_mvt(3) = nanmean(data.total_mvt_ds);

%% Analyze final closed-loop panels

load([path,'\analysis\analysis_sid_',num2str(sessions.final_cl_bar),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff{1};
plot(dff);
title('EB-DAN activity');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,2)
plot(data.vel_for_ds);
title('Forward velocity (mm/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
plot(abs(data.vel_yaw_ds));
title('Yaw speed (deg/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.total_mvt_ds);
xlabel('Time (sec)');
title('Total velocity (deg/s)');
xlim([0 data.time(end)]);
set(gca,'xticklabel',{[]})

suptitle('Final trial with just panels');

saveas(gcf,[path,'\analysis\plots\final_panels.png']);

%Get mean parameters
mean_EBDAN_act(4) = nanmean(dff);
mean_total_mvt(5) = nanmean(data.total_mvt_ds);

%% Analyze final closed-loop wind

load([path,'\analysis\analysis_sid_',num2str(sessions.final_cl_wind),'_tid_0.mat'])

figure('Position',[100 100 1200 800]),
subplot(4,1,1)
dff = data.dff{1};
plot(dff);
title('EB-DAN activity');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,2)
plot(data.vel_for_ds);
title('Forward velocity (mm/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,3)
plot(abs(data.vel_yaw_ds));
title('Yaw speed (deg/s)');
xlim([0 length(dff)]);
set(gca,'xticklabel',{[]})

subplot(4,1,4)
plot(data.time,data.total_mvt_ds);
xlabel('Time (sec)');
title('Total velocity (deg/s)');
xlim([0 data.time(end)]);
set(gca,'xticklabel',{[]})

suptitle('Final trial with just wind');

saveas(gcf,[path,'\analysis\plots\final_wind.png']);

%Get mean parameters
mean_EBDAN_act(5) = nanmean(dff);
mean_total_mvt(5) = nanmean(data.total_mvt_ds);

%% Analyze average EB-DAN and movement

figure,
yyaxis left
plot(mean_EBDAN_act,'-o')
ylabel('Mean DFF');
yyaxis right
plot(mean_total_mvt,'-o')
ylabel('Mean total movement (deg/s)');
xlim([0 6]);
xlabel('Block number');

saveas(gcf,[path,'\analysis\plots\activity_and_mvt_per_block.png']);

