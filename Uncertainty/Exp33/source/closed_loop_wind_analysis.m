%Analysis code for the closed-loop bouts of the change in contrast experiment


%% Load data

clear all; close all;

%Get the pre-processed data
[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp33\data\');
load([path,file])

%Determine the sid we're working with, to save the plots to specific folder
%later
sid = file(24:26);

%% Make directory to save plots

%Move to the analysis folder
cd(path)
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(path,'plots'); 
end
%List the contents of the 'plots' folder
cd([path,'plots\'])


%% Recover PMT gain and zoom factor
% 
% imagingDir = [path(1:end-10), '\2p\'];
% expression = ['*sid' sid 'tid_0_*'];
% imagingFile = dir(fullfile(imagingDir, expression));
% 
% import ScanImageTiffReader.ScanImageTiffReader; %imports functinos necessary to open the scanimage tiff files
% reader = ScanImageTiffReader(fullfile(imagingDir, imagingFile(1).name));
% rawFile_original = reader.data();
% metadata_char = reader.metadata();
% %descriptions = reader.descriptions();
% list = strsplit(metadata_char,'RoiGroups');%,'CollapseDelimiters',true);
% list = char(list{1});
% list = list(1:end-7);
% eval(list);
% zoom = SI.hRoiManager.scanZoomFactor;
% PMT_gain = SI.hPmts.gains(1);

%% Determine airflow

%import run_obj
run_obj_dir = dir([path(1:end-9),'ball\runobj\']);
for run_obj_file = 1:length(run_obj_dir)
    if contains(run_obj_dir(run_obj_file).name,sid)
        load(fullfile(run_obj_dir(run_obj_file).folder,run_obj_dir(run_obj_file).name))
    end
end

airflow = run_obj.airflow.Value;

%% Plot the activity heatmap, phase and heading positions, and offset in time

figure,
subplot(3,1,1)
dff_data = continuous_data.dff_matrix';
imagesc(dff_data)
colormap(flipud(gray))
title('EPG activity in the PB');

subplot(3,1,2)
bump_pos = wrapTo180(rad2deg(-continuous_data.bump_pos));
%Remove wrapped lines to plot
[x_out_bump,bump_to_plot] = removeWrappedLines(continuous_data.time,bump_pos');
plot(x_out_bump,bump_to_plot,'LineWidth',1.5)
hold on
motor_pos = wrapTo180(rad2deg(-continuous_data.motor_pos));
%Remove wrapped lines to plot
[x_out_motor,motor_to_plot] = removeWrappedLines(continuous_data.time,motor_pos');
plot(x_out_motor,motor_to_plot,'LineWidth',1.5)
ylim([-180 180]);
title('Bump and motor position');
legend('bump position','motor position')

subplot(3,1,3)
offset = wrapTo180(rad2deg(circ_dist(-continuous_data.bump_pos,-continuous_data.motor_pos)));
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset');
plot(x_out_offset,offset_to_plot,'LineWidth',1.5)
ylim([-180 180]);
title('Offset');
xlabel('Time (sec)');

suptitle(['Airflow = ', num2str(airflow),' L/min']);

%save figure
saveas(gcf,[path,'plots\closedLoopWind_',num2str(airflow),'.png']);

%% Calculate and plot bump magnitude and bump width in time

figure('Position',[200 200 1600 600]),
%Plot EPG activity
subplot(3,1,1)
imagesc(flip(dff_data))
colormap(flipud(gray))
title('EPG activity in the PB','fontweight','bold','fontsize',12);
set(gca,'xtick',[]);
set(gca,'XTickLabel',[]);
legend('Change in stimulus');

%Plot bump magnitude
subplot(3,1,2)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5),'.r')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_magnitude(continuous_data.adj_rs<0.5),'.k')
title('Bump magnitude','fontweight','bold','fontsize',12);
ylabel({'Bump magnitude';'(from von Mises fit)'},'fontweight','bold','fontsize',10);
xlim([0 continuous_data.time(end)]);

%Plot bump width
subplot(3,1,3)
plot(continuous_data.time(continuous_data.adj_rs>=0.5),continuous_data.bump_width(continuous_data.adj_rs>=0.5),'.r')
hold on
plot(continuous_data.time(continuous_data.adj_rs<0.5),continuous_data.bump_width(continuous_data.adj_rs<0.5),'.k')
title('Bump width','fontweight','bold','fontsize',12);
ylabel({'Bump width';'(from von Mises fit)'},'fontweight','bold','fontsize',10);
xlabel('Time (sec)','fontweight','bold','fontsize',10);
xlim([0 continuous_data.time(end)]);


%save figure
saveas(gcf,[path,'plots\closedLoopBMandBWinTime_',num2str(airflow),'.png']);


%% Compute some important variables

meanBM = mean(continuous_data.bump_magnitude(continuous_data.adj_rs>=0.5));
meanBW = mean(continuous_data.bump_width(continuous_data.adj_rs>=0.5));
[~, offset_var] = circ_std(deg2rad(offset(continuous_data.adj_rs>=0.5)),[],[],2);
[~, heading_var] = circ_std(deg2rad(motor_pos(continuous_data.adj_rs>=0.5)),[],[],2);

%% Save relevant data

save([path,'summary_data_sid',sid,'.mat'],'meanBM','meanBW','offset_var','heading_var','airflow')

close all; clear all;