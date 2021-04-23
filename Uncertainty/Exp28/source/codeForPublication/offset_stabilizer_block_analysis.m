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
ylabel('PB glomerulus','fontweight','bold','fontsize',10);
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
title('EPG activity in the PB','fontweight','bold','fontsize',12);
set(gca,'xtick',[]);
set(gca,'xticklabel',{[]});

subplot(3,6,[7 11])
%Get bar position to plot
bar_position = wrapTo180(data.panel_angle);
%Remove lines that wrap around using auxiliary function
[x_bar_pos,bar_position_to_plot] = removeWrappedLines(data.time,bar_position);
plot(x_bar_pos,bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
hold on
%Get EPG phase to plot
phase = wrapTo180(rad2deg(data.dff_pva));
%Remove lines that wrap around using auxiliary function
[x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
plot(x_out_phase,phase_to_plot,'LineWidth',1.5,'color',[0.9 0.3 0.4])
xlim([0 x_out_phase(end)]);
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
offset = wrapTo180(data.bar_offset);
%Remove lines that wrap around using auxiliary function
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
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
saveas(gcf,[path,'plots\Offset_stabilizer_block_heatmap.png']);