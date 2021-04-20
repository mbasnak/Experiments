%analysis for the stabilizing offset bout


%% Load data

clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');
load([path,file])

%determine the sid we're working with, to save the plots to specific folder
%later
sid = file(10:14);

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

%% Plot heatmap with bar position, phase and offset

figure('Position',[200 200 1000 600]),
subplot(3,6,[1 5])
%Plot heatmap of EPG activity
imagesc(data.dff_matrix)
colormap(flipud(gray))
ylabel('PB glomerulus');
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
title('EPG activity in the PB');
set(gca,'xtick',[]);
set(gca,'xticklabel',{[]});

subplot(3,6,[7 11])
%Get bar position to plot
bar_position = wrapTo180(data.panel_angle);
change_bar_position = abs([0;diff(smooth(bar_position))]);
bar_position_to_plot = smooth(bar_position);
bar_position_to_plot(change_bar_position>40==1) = NaN;
plot(bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
hold on
%Get EPG phase to plot
phase = wrapTo180(rad2deg(data.dff_pva));
change_phase = abs([0;diff(smooth(phase))]);
phase_to_plot = smooth(phase);
phase_to_plot(change_phase>40==1) = NaN;
plot(phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
xlim([0 length(data.panel_angle)]);
ylim([-180 180]);
axP = get(gca,'Position');
legend('Bar position','EPG phase','Location','EastOutside')
set(gca, 'Position', axP);
ylabel('Deg');
title('Bar and bump position');
set(gca,'xticklabel',{[]})

subplot(3,6,[13 17])
%Get offset to plot
offset = wrapTo180(data.bar_offset);
change_offset = abs([0;diff(smooth(offset))]);
offset_to_plot = smooth(offset);
offset_to_plot(change_offset>30==1) = NaN;
plot(data.time,offset_to_plot,'k','LineWidth',1.5)
xlim([0 data.time(end)]);
ylim([-180 180]);
xlabel('Time (sec)');
ylabel('Deg');
title('Offset');

subplot(3,6,18)
polarhistogram(deg2rad(offset),'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title('Offset distribution');
rticklabels({}); %remove rticklabels

%Save figure
saveas(gcf,[path,'plots\Offset_stabilizer_block_heatmap.png']);

% %% Compute bump magnitude and plot in time
% 
% bump_mag = data.bump_magnitude;
% 
% %plot in time
% figure('Position',[200 200 1400 600]),
% subplot(2,1,1)
% %Plot heatmap of EPG activity
% imagesc(data.dff_matrix)
% colormap(flipud(gray))
% ylabel('PB glomerulus');
% title('EPG activity in the PB');
% set(gca,'XTickLabel',[]);
% 
% subplot(2,1,2)
% plot(data.time,bump_mag)
% title('Bump magnitude');
% xlabel('Time (sec)');
% ylabel('Bump magnitude (max-min)');
% 
% %Save figure
% saveas(gcf,[path,'plots\Offset_stabilizer_block_bm.png']);
% %% Scatter plots of bump magnitude and movement
% 
% figure('Position',[100 100 1600 1000]),
% subplot(1,4,1)
% binscatter(data.vel_for_deg_ds',bump_mag)
% title('Forward velocity');
% 
% subplot(1,4,2)
% binscatter(abs(data.vel_side_deg_ds'),bump_mag)
% title('Side speed');
% 
% subplot(1,4,3)
% binscatter(abs(data.vel_yaw_ds),bump_mag)
% title('Yaw speed')
% 
% subplot(1,4,4)
% binscatter(data.total_mvt_ds,bump_mag)
% title('Total movement');
% 
% %Save figure
% saveas(gcf,[path,'plots\Offset_stabilizer_block_mvt_scatter.png']);

%% Binned plots of bump magnitude and movement

bump_mag = data.bump_magnitude;

figure('Position',[200 200 1400 600]),

nbins = 20;
max_fwd_bin = max(data.vel_for_deg_ds)-10;
min_mvt_bin = 0;
binWidth = max_fwd_bin/nbins;
fwdBins = [0:binWidth:max_fwd_bin]; 

%getting binned means 
for bin = 1:length(fwdBins)-1
    meanBin(bin) = mean(bump_mag((data.vel_for_deg_ds > fwdBins(bin)) & (data.vel_for_deg_ds < fwdBins(bin+1))));
    stdBin(bin) = std(bump_mag((data.vel_for_deg_ds > fwdBins(bin)) & (data.vel_for_deg_ds < fwdBins(bin+1))));
    errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((data.vel_for_deg_ds > fwdBins(bin)) & (data.vel_for_deg_ds < fwdBins(bin+1))))));
end

%create axes for plot
mvtAxes = fwdBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,1)
boundedline(mvtAxes,meanBin,errBin)
ylabel('Bump magnitude (max-min)'); xlabel('Fwd velocity (deg/s)');
ylim([0 max(bump_mag)]);
title('Forward velocity');


%side speed
max_side_bin = max(abs(data.vel_side_deg_ds))-10;
binWidth = max_side_bin/nbins;
sideBins = [0:binWidth:max_side_bin]; 

%getting binned means 
for bin = 1:length(sideBins)-1
    meanBin(bin) = mean(bump_mag((abs(data.vel_side_deg_ds) > sideBins(bin)) & (abs(data.vel_side_deg_ds) < sideBins(bin+1))));
    stdBin(bin) = std(bump_mag((abs(data.vel_side_deg_ds)> sideBins(bin)) & (abs(data.vel_side_deg_ds) < sideBins(bin+1))));
    errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((abs(data.vel_side_deg_ds) > sideBins(bin)) & (abs(data.vel_side_deg_ds) < sideBins(bin+1))))));
end

%create axes for plot
mvtAxes = sideBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,2)
boundedline(mvtAxes,meanBin,errBin,'m')
xlabel('Side speed (deg/s)');
ylim([0 max(bump_mag)]);
title('Side speed');


%yaw speed
max_yaw_bin = max(abs(data.vel_yaw_ds))-10;
binWidth = max_yaw_bin/nbins;
yawBins = [0:binWidth:max_yaw_bin]; 

%getting binned means 
for bin = 1:length(yawBins)-1
    meanBin(bin) = mean(bump_mag((abs(data.vel_yaw_ds) > yawBins(bin)) & (abs(data.vel_yaw_ds) < yawBins(bin+1))));
    stdBin(bin) = std(bump_mag((abs(data.vel_yaw_ds)> yawBins(bin)) & (abs(data.vel_yaw_ds) < yawBins(bin+1))));
    errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((abs(data.vel_yaw_ds) > yawBins(bin)) & (abs(data.vel_yaw_ds) < yawBins(bin+1))))));
end

%create axes for plot
mvtAxes = yawBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,3)
boundedline(mvtAxes,meanBin,errBin,'r')
xlabel('Yaw speed (deg/s)');
ylim([0 max(bump_mag)]);
title('Yaw speed');



%total movement
max_mvt_bin = max(data.total_mvt_ds)-10;
binWidth = max_mvt_bin/nbins;
mvtBins = [0:binWidth:max_mvt_bin]; 

%getting binned means 
for bin = 1:length(mvtBins)-1
    meanBin(bin) = mean(bump_mag((data.total_mvt_ds > mvtBins(bin)) & (data.total_mvt_ds < mvtBins(bin+1))));
    stdBin(bin) = std(bump_mag((data.total_mvt_ds > mvtBins(bin)) & (data.total_mvt_ds < mvtBins(bin+1))));
    errBin(bin) = stdBin(bin)./sqrt(length((bump_mag((data.total_mvt_ds > mvtBins(bin)) & (data.total_mvt_ds < mvtBins(bin+1))))));
end

%create axes for plot
mvtAxes = mvtBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
subplot(1,4,4)
boundedline(mvtAxes,meanBin,errBin,'k')
xlabel('Total movement (deg/s)');
ylim([0 max(bump_mag)]);
title('Total movement');

%Save figure
saveas(gcf,[path,'plots\Offset_stabilizer_block_mvt_trend.png']);

%% Get reference offset from the last part of the block

visual_offset = circ_dist(data.dff_pva,deg2rad(wrapTo180(data.panel_angle)));
mean_reference_offset = rad2deg(circ_mean(visual_offset(600:900),[],1));

save([path,'\reference_offset.mat'],'mean_reference_offset')

close all; clc;