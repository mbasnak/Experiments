%analysis for the verifying offset bout


%% Load data

clear all; close all;

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');

load([path,'\sessions_info.mat'])

sid = sessions_info.offset_verifier;
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

%% Plot heatmap with bar position, phase and offset

figure('Position',[200 200 1000 600]),
subplot(3,6,[1 5])
%Plot heatmap of EPG activity
dff_matrix = continuous_data.dff_matrix';
imagesc(flip(dff_matrix))
colormap(flipud(gray))
ylabel('PB glomerulus','fontweight','bold','fontsize',12);
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
ylabel('Deg','fontweight','bold','fontsize',12);
title('Bar and bump position','fontweight','bold','fontsize',12);
set(gca,'xticklabel',{[]})

subplot(3,6,[13 17])
%Get offset to plot
offset = wrapTo180(rad2deg(circ_dist(continuous_data.bump_pos',-continuous_data.heading)));
%Remove lines that wrap around using auxiliary function
[x_out_offset,offset_to_plot] = removeWrappedLines(continuous_data.time,offset);
plot(x_out_offset,offset_to_plot,'k','LineWidth',1.5)
xlim([0 x_out_offset(end)]);
ylim([-180 180]);
xlabel('Time (sec)','fontweight','bold','fontsize',12);
ylabel('Deg','fontweight','bold','fontsize',12);
title('Offset','fontweight','bold','fontsize',12);

%Plot circular offset distribution
subplot(3,6,18)
polarhistogram(deg2rad(offset),'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title('Offset distribution','fontweight','bold','fontsize',12);
rticklabels({}); %remove rticklabels

%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\offset_verifyier_block_heatmap.png']);

%% Get vefying offset

visual_offset = circ_dist(continuous_data.bump_pos',deg2rad(continuous_data.panel_angle));
mean_verifying_offset = rad2deg(circ_mean(visual_offset,[],1));

%save 
save([path,'\analysis\continuous_verifying_offset.mat'],'mean_verifying_offset');

%% 
clear all; close all; clc;