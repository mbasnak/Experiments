%analysis for the stabilizing offset bout


%% Load data

clear all; close all;

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');

load([path,'\sessions_info.mat'])

sid = sessions_info.offset_stabilizer;
load([path,'\analysis\continuous_analysis_sid_',num2str(sid),'_tid_0.mat'])

%% Make directory to save plots

%Move to the analysis folder
cd([path,'\analysis'])
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'continuous_plots') == 0)
   mkdir(path,'continuous_plots'); 
end
%List the contents of the 'plots' folder
cd([path,'\analysis\continuous_plots\'])

%% Plot heatmap with bar position, phase and offset

figure('Position',[200 200 1000 600]),
subplot(3,6,[1 5])
%Plot heatmap of EPG activity
dff_matrix = continuous_data.dff_matrix';
imagesc(flip(dff_matrix))
colormap(flipud(gray))
ylabel('PB glomerulus','fontweight','bold','fontsize',10);
set(gca,'ytick',[]);
set(gca,'yticklabel',{[]});
title('EPG activity in the PB','fontweight','bold','fontsize',12);
set(gca,'xtick',[]);
set(gca,'xticklabel',{[]});

subplot(3,6,[7 11])
%Get bar position to plot
bar_position = wrapTo180(continuous_data.panel_angle);
%Remove lines that wrap around using auxiliary function
[x_bar_pos,bar_position_to_plot] = removeWrappedLines(continuous_data.time,bar_position);
plot(x_bar_pos,bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
hold on
%Get EPG phase to plot
phase = wrapTo180(rad2deg(continuous_data.bump_pos'));
%Remove lines that wrap around using auxiliary function
[x_out_phase,phase_to_plot] = removeWrappedLines(continuous_data.time,phase);
plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
xlim([0 x_out_phase(end-1)]);
ylim([-180 180]);
legend('Bar position','EPG phase','Location','EastOutside')
%Set legend outside of the plot so that it doesn't occlude the traces
axP = get(gca,'Position');
set(gca, 'Position', axP);
ylabel('Deg','fontweight','bold','fontsize',10);
title('Bar and bump position','fontweight','bold','fontsize',12);
set(gca,'xticklabel',{[]})

subplot(3,6,[13 17])
%Get offset to plot
offset = wrapTo180(continuous_data.offset);
%Remove lines that wrap around using auxiliary function
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'k','LineWidth',1.5)
xlim([0 x_out_offset(end)]);
ylim([-180 180]);
xlabel('Time (sec)','fontweight','bold','fontsize',10);
ylabel('Deg','fontweight','bold','fontsize',10);
title('Offset','fontweight','bold','fontsize',12);

%Plot circular offset distribution
subplot(3,6,18)
polarhistogram(deg2rad(offset),'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title('Offset distribution','fontweight','bold','fontsize',12);
rticklabels({}); %remove rticklabels

%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\Offset_stabilizer_block_heatmap.png']);


%% Analyze relationship between yaw speed data and BM

nbins = 20;
yaw_speed = abs(continuous_data.vel_yaw_ds);
maxBinYS = max(yaw_speed);
binWidthYS = maxBinYS/nbins;
YSBins = [0:binWidthYS:maxBinYS]; 

%getting binned means 
for bin = 1:length(YSBins)-1
    meanBin(bin) = mean(continuous_data.bump_magnitude((yaw_speed(1:length(continuous_data.bump_magnitude)) > YSBins(bin)) & (yaw_speed(1:length(continuous_data.bump_magnitude)) < YSBins(bin+1))));
end

%create axes for plot
YSAxes = YSBins - binWidthYS;
YSAxes = YSAxes(2:end);
YSAxes(end) = YSAxes(end-1)+binWidthYS;

%Plot
figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),continuous_data.bump_magnitude, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot total movement 
subplot(2,4,[5 7])
plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),yaw_speed(1:length(continuous_data.bump_magnitude)),'k')
xlabel('Time (sec)');
ylabel('Yaw speed (deg/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(YSAxes,meanBin,'-ko')
ylabel('Mean bump magnitude'); xlabel('Yaw speed (deg/s)');
ylim([0 (max(meanBin)+0.5)]);


%% Bin for vel data

for_vel = abs(continuous_data.vel_for_ds');
maxBinFV = max(for_vel);
binWidthFV = maxBinFV/nbins;
FVBins = [0:binWidthFV:maxBinFV]; 

%getting binned means 
for bin = 1:length(FVBins)-1
    meanBin(bin) = mean(continuous_data.bump_magnitude((for_vel(1:length(continuous_data.bump_magnitude)) > FVBins(bin)) & (for_vel(1:length(continuous_data.bump_magnitude)) < FVBins(bin+1))));
end

%create axes for plot
fvAxes = FVBins - binWidthFV;
fvAxes = fvAxes(2:end);
fvAxes(end) = fvAxes(end-1)+binWidthFV;

%Plot
figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),continuous_data.bump_magnitude, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot total movement 
subplot(2,4,[5 7])
plot(continuous_data.time(1:length(continuous_data.bump_magnitude)),for_vel(1:length(continuous_data.bump_magnitude)),'k')
xlabel('Time (sec)');
ylabel('Forward velocity (mm/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(fvAxes,meanBin,'-ko')
ylabel('Mean bump magnitude'); xlabel('Forward velocity (mm/s)');
ylim([0 (max(meanBin)+0.5)]);

%% Bin both parameters to build heatmap

%getting binned means 
for ys_bin = 1:length(YSBins)-1
    for fv_bin = 1:length(FVBins)-1
        doubleBin(ys_bin,fv_bin) = mean(continuous_data.bump_magnitude((yaw_speed(1:length(continuous_data.bump_magnitude)) > YSBins(ys_bin)) & (yaw_speed(1:length(continuous_data.bump_magnitude)) < YSBins(ys_bin+1)) & (for_vel(1:length(continuous_data.bump_magnitude)) > FVBins(fv_bin)) & (for_vel(1:length(continuous_data.bump_magnitude)) < FVBins(fv_bin+1))));
    end
end

%flip the data such that high forward velocity values are at the top
binned_data = flip(doubleBin);

figure,
imagesc(binned_data)
xticks([1:4:20])
set(gca, 'XTickLabel', round(FVBins(1:4:20)))
xlabel('Forward velocity (mm/s)','fontsize',12,'fontweight','bold');
ylabel('Yaw speed (deg/sec)','fontsize',12,'fontweight','bold');
yticks([1:4:20])
set(gca, 'YTickLabel', round(YSBins(20:-4:1)))
c = colorbar;
ylabel(c, 'Bump magnitude')


%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\bm_heatmap.png']);

%%
clear all; close all;