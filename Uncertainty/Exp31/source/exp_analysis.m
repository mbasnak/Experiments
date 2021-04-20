%% Load data

clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp31\data');
load([path,file])

%determine the sid we're working with, to save the plots to specific folder
%later
sid = file(10:14);
session_num = str2num(sid(end));

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

%% Load airflow data for this fly

load([path,'airflow.mat'])
sid_airflow = airflow(session_num+1);

%% Plot heatmap with wind position, phase and offset

figure('Position',[200 200 1000 600]),
subplot(3,6,[1 5])
%Plot heatmap of EPG activity
imagesc(data.dff_matrix)
colormap(flipud(gray))
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
title('EPG activity in the PB');
set(gca,'xtick',[]);
set(gca,'xticklabel',{[]});
title('EPG activity in the PB');

subplot(3,6,[7 11])
%Get wind position to plot
wind_position = wrapTo180(data.motor_angle);
[x_out_wind,wind_position_to_plot] = removeWrappedLines(data.time,wind_position); 
plot(x_out_wind,wind_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
hold on
%Get EPG phase to plot
phase = wrapTo180(rad2deg(data.dff_pva));
[x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase); 
plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
xlim([0 x_out_phase(end-1)]);
ylim([-180 180]);
axP = get(gca,'Position');
legend('Wind position','EPG phase','Location','EastOutside')
set(gca, 'Position', axP);
ylabel('Deg');
title('wind and bump position');

subplot(3,6,[13 17])
%Get offset to plot
offset = wrapTo180(data.wind_offset);
change_offset = abs([0;diff(smooth(offset))]);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset); 
plot(x_out_offset,offset_to_plot,'k','LineWidth',1.5)
xlim([0 x_out_offset(end-1)]);
ylim([-180 180]);
xlabel('Time');
ylabel('Deg');
title('Offset');

subplot(3,6,18)
polarhistogram(deg2rad(offset),'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title('Offset distribution');
pax = gca;
pax.RTickLabels = [];

suptitle(['Airflow = ',num2str(sid_airflow),' L/min']);

%Save figure
saveas(gcf,[path,'plots\',sid,'_heatmap.png']);
saveas(gcf,[path,'plots\',sid,'_heatmap.svg']);
saveas(gcf,[path,'plots\',sid,'_heatmap.pdf']);


%% Binned plots of bump magnitude and movement

% bump_mag = data.bump_magnitude;
% 
% figure('Position',[200 200 1400 600]),
% 
% nbins = 20;
% max_fwd_bin = max(data.vel_for_deg_ds)-10;
% min_fwd_bin = 0;
% binWidth = (max_fwd_bin-min_fwd_bin)/nbins;
% fwdBins = [min_fwd_bin:binWidth:max_fwd_bin]; 
% 
% %getting binned means 
% for bin = 1:length(fwdBins)-1
%     meanBin(bin) = mean(bump_mag((data.vel_for_deg_ds > fwdBins(bin)) & (data.vel_for_deg_ds < fwdBins(bin+1))));
%     stdBin(bin) = std(bump_mag((data.vel_for_deg_ds > fwdBins(bin)) & (data.vel_for_deg_ds < fwdBins(bin+1))));
%     errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((data.vel_for_deg_ds > fwdBins(bin)) & (data.vel_for_deg_ds < fwdBins(bin+1))))));
% end
% 
% %create axes for plot
% mvtAxes = fwdBins - binWidth;
% mvtAxes = mvtAxes(2:end);
% mvtAxes(end) = mvtAxes(end-1)+binWidth;
% 
% %Plot
% subplot(1,4,1)
% boundedline(mvtAxes,meanBin,errBin)
% ylabel('Bump magnitude (max-min)'); xlabel('Fwd velocity (deg/s)');
% ylim([0 max(bump_mag)]);
% title('Forward velocity');
% 
% 
% %side speed
% max_side_bin = max(abs(data.vel_side_deg_ds))-10;
% binWidth = max_side_bin/nbins;
% sideBins = [0:binWidth:max_side_bin]; 
% 
% %getting binned means 
% for bin = 1:length(sideBins)-1
%     meanBin(bin) = mean(bump_mag((abs(data.vel_side_deg_ds) > sideBins(bin)) & (abs(data.vel_side_deg_ds) < sideBins(bin+1))));
%     stdBin(bin) = std(bump_mag((abs(data.vel_side_deg_ds)> sideBins(bin)) & (abs(data.vel_side_deg_ds) < sideBins(bin+1))));
%     errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((abs(data.vel_side_deg_ds) > sideBins(bin)) & (abs(data.vel_side_deg_ds) < sideBins(bin+1))))));
% end
% 
% %create axes for plot
% mvtAxes = sideBins - binWidth;
% mvtAxes = mvtAxes(2:end);
% mvtAxes(end) = mvtAxes(end-1)+binWidth;
% 
% %Plot
% subplot(1,4,2)
% boundedline(mvtAxes,meanBin,errBin,'m')
% xlabel('Side speed (deg/s)');
% ylim([0 max(bump_mag)]);
% title('Side speed');
% 
% 
% %yaw speed
% max_yaw_bin = max(abs(data.vel_yaw_ds))-10;
% binWidth = max_yaw_bin/nbins;
% yawBins = [0:binWidth:max_yaw_bin]; 
% 
% %getting binned means 
% for bin = 1:length(yawBins)-1
%     meanBin(bin) = mean(bump_mag((abs(data.vel_yaw_ds) > yawBins(bin)) & (abs(data.vel_yaw_ds) < yawBins(bin+1))));
%     stdBin(bin) = std(bump_mag((abs(data.vel_yaw_ds)> yawBins(bin)) & (abs(data.vel_yaw_ds) < yawBins(bin+1))));
%     errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((abs(data.vel_yaw_ds) > yawBins(bin)) & (abs(data.vel_yaw_ds) < yawBins(bin+1))))));
% end
% 
% %create axes for plot
% mvtAxes = yawBins - binWidth;
% mvtAxes = mvtAxes(2:end);
% mvtAxes(end) = mvtAxes(end-1)+binWidth;
% 
% %Plot
% subplot(1,4,3)
% boundedline(mvtAxes,meanBin,errBin,'r')
% xlabel('Yaw speed (deg/s)');
% ylim([0 max(bump_mag)]);
% title('Yaw speed');
% 
% 
% 
% %total movement
% max_mvt_bin = max(data.total_mvt_ds)-10;
% binWidth = max_mvt_bin/nbins;
% mvtBins = [0:binWidth:max_mvt_bin]; 
% 
% %getting binned means 
% for bin = 1:length(mvtBins)-1
%     meanBin(bin) = mean(bump_mag((data.total_mvt_ds > mvtBins(bin)) & (data.total_mvt_ds < mvtBins(bin+1))));
%     stdBin(bin) = std(bump_mag((data.total_mvt_ds > mvtBins(bin)) & (data.total_mvt_ds < mvtBins(bin+1))));
%     errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((data.total_mvt_ds > mvtBins(bin)) & (data.total_mvt_ds < mvtBins(bin+1))))));
% end
% 
% %create axes for plot
% mvtAxes = mvtBins - binWidth;
% mvtAxes = mvtAxes(2:end);
% mvtAxes(end) = mvtAxes(end-1)+binWidth;
% 
% %Plot
% subplot(1,4,4)
% boundedline(mvtAxes,meanBin,errBin,'k')
% xlabel('Total movement (deg/s)');
% ylim([0 max(bump_mag)]);
% title('Total movement');
% 
% %Save figure
% saveas(gcf,[path,'plots\',sid,'_mvt_trend.png']);

%% Bump mag

mean_bm = mean(data.bump_magnitude);
median_bm = median(data.bump_magnitude);

%% Compute heading and offset variation

[~, offset_var] = circ_std(deg2rad(offset),[],[],1); 
[~, heading_var] = circ_std(data.heading,[],[],1); 

%% save data and close

save([path,'\summary_data',sid,'.mat'],'offset_var','heading_var','sid_airflow','mean_bm','median_bm');
clear all; close all
