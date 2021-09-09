%analysis for the stabilizing offset bout


%% Load data

clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
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
offset = wrapTo180(data.offset);
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
