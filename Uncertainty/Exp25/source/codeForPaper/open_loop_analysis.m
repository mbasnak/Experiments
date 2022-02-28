%Code to analyze the open-loop bouts

clear all; close all;

%% Load data

%Get directory you're interested in
path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts');

%I have a .txt file for each fly with the open-loop session numbers
%Read the .txt file containing the ID of the open loop sessions.
fileID = fopen([path,'\2p\open_loop_sessions.txt'],'r');
formatSpec = '%d';
open_loop_sessions = fscanf(fileID,formatSpec);

%Load the data
for session = 1:length(open_loop_sessions)
   data{session} = load([path,'\analysis\continuous_analysis_sid_',num2str(open_loop_sessions(session)),'_tid_0.mat']);    
end

%% For each open-loop bout, plot the heatmap of PB activity, the trajectory of the stimulus and the 'stimulus offset'
%defined as the circular distance between the phase of EPG activity and the stimulus
%position.

for session = 1:length(open_loop_sessions)
    
    figure('Position',[100, 100, 1000, 800]),
    
    %PB heatmap
    subplot(3,5,[1 4])
    dff_matrix = data{1,session}.continuous_data.dff_matrix';
    imagesc(flip(dff_matrix))
    colormap(flipud(gray))
    ylabel('PB glomerulus','fontsize',12)
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    set(gca,'yticklabel',{[]});
    set(gca,'YTick',[]);
    title('EPG activity in the PB','fontsize',12);
    
    %Stimulus position
    subplot(3,5,[6 9])
    bar_position = wrapTo180(data{1,session}.continuous_data.panel_angle);
    [x_out_bar, bar_position_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time, bar_position);
    %Get EPG phase to plot
    phase = wrapTo180(rad2deg(data{1,session}.continuous_data.bump_pos'));
    [x_out_phase, phase_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time, phase);
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_bar,bar_position_to_plot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
    else
        plot(x_out_bar,bar_position_to_plot,'color',[0 0 0.6],'LineWidth',1.5)
    end  
    hold on
    plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
    axP = get(gca,'Position');
    legend('Bar position','EPG phase','Location','EastOutside')
    set(gca, 'Position', axP);
    if ~isnan(x_out_phase(end))
        xlim([0, x_out_phase(end)]);
    else
        xlim([0, x_out_phase(end-1)]);
    end
    ylim([-180 180]);
    ylabel('Angular position (deg)','fontsize',12);
    set(gca,'xticklabel',{[]})
    title('Stimulus position','fontsize',12);
    
    %Stim offset
    %Calculate stim offset as the circular distance between EPG activity
    %phase and the stimulus posotion
    stim_offset{session} = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos',deg2rad(data{1,session}.continuous_data.panel_angle)));
    stimoffset = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos',deg2rad(data{1,session}.continuous_data.panel_angle)));    
    [x_out_stimoffset, stimOffsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time, stimoffset);
    subplot(3,5,[11 14])
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_stimoffset,stimOffsetToPlot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
    else
        plot(x_out_stimoffset,stimOffsetToPlot,'color',[0 0 0.6],'LineWidth',1.5)
    end  
    xlim([0, data{1,session}.continuous_data.time(end)]);
    ylim([-180 180]);
    xlabel('Time (sec)','fontsize',12); ylabel('Degrees','fontsize',12);
    title('Stimulus offset','fontsize',12);
    
    %Plot stimulus offset distribution
    subplot(3,5,15)
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        polarhistogram(deg2rad(stim_offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','High contrast'});
        Ax = gca;
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    else
        polarhistogram(deg2rad(stim_offset{session}),20,'FaceColor',[ 0 0 0.6])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','Low contrast'});
        Ax = gca; 
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    end
    
    %Save figures
    saveas(gcf,[path,'\analysis\continuous_plots\open_loop_heatmap_and_stim_offset_sid',num2str(open_loop_sessions(session)),'.png']);
end


%% Summarize relevant data in table

close all;

%Calculate offset variation as the circular standard deviation of the
%'stimulus offset'
for session = 1:length(open_loop_sessions)
   [~,stim_offset_var(session)] = circ_std(deg2rad(stim_offset{session}),[],[],1);  
end

%Define stimulus velocity conditions according to which panel functions were used
low_stim_vel = [50,51]; %20 deg/s
medium_stim_vel = [52,53]; %30 deg/s
high_stim_vel = [54,55]; %60 deg/s

%Define stimulus contrasts according to which panel pattern was used
low_stim_contrast = 56; %1/15
high_stim_contrast = 57; %15/15

%Get bump parameters and compute the mean
for session = 1:length(open_loop_sessions)
    
    moving = data{1,session}.continuous_data.total_mvt_ds > 20;
    good_fit = data{1,session}.continuous_data.adj_rs >= 0.5;    
    
    BumpMagnitude{session} = data{1,session}.continuous_data.bump_magnitude;
    meanBM(session) = mean(BumpMagnitude{session});  
    meanBM_thresh(session) = mean(BumpMagnitude{session}(moving & good_fit));
    
    bump_width{session} = data{1,session}.continuous_data.bump_width;        
    meanBW(session) = mean(bump_width{session}); 
    meanBW_thresh(session) = mean(bump_width{session}(moving & good_fit));    
end

%Create a table with relevant variables for each session
summarydata = array2table(zeros(0,7), 'VariableNames',{'stim_offset_var','bump_mag','bump_width','bump_mag_thresh','bump_width_thresh','contrast_level','stim_vel'});
warning('off');
for session = 1:length(open_loop_sessions)
    summarydata{session,'stim_offset_var'} = stim_offset_var(session);
    summarydata{session,'bump_mag'} = meanBM(session);
    summarydata{session,'bump_width'} = meanBW(session);
    summarydata{session,'bump_mag_thresh'} = meanBM_thresh(session);
    summarydata{session,'bump_width_thresh'} = meanBW_thresh(session);
    summarydata{session,'contrast_level'} = data{1,session}.continuous_data.run_obj.pattern_number;
    if (data{1,session}.continuous_data.run_obj.function_number == 50 | data{1,session}.continuous_data.run_obj.function_number == 51)
        summarydata{session,'stim_vel'} = 1; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 52 | data{1,session}.continuous_data.run_obj.function_number == 53)
        summarydata{session,'stim_vel'} = 2; 
    else
        summarydata{session,'stim_vel'} = 3; 
    end
end

%% Compare the offset variability for each speed between contrasts

%Average offset variability, grouped by stimulus contrast and velocity
mean_stim_offset_var_data = varfun(@mean,summarydata,'InputVariables','stim_offset_var',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_stim_offset_var_data,1),1); 
for contrast = 1:length(mean_stim_offset_var_data.contrast_level)
    if mean_stim_offset_var_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_stim_offset_var_data.stim_vel,mean_stim_offset_var_data.mean_stim_offset_var,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_stim_offset_var_data.mean_stim_offset_var) - 0.2, max(mean_stim_offset_var_data.mean_stim_offset_var) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20','30','60'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[~,h2] = suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Offset variability (rad)','FontSize',12);

%save
saveas(gcf,[path,'\analysis\continuous_plots\open_loop_offset_var.png']);

%% Compare the bump magnitude for each speed between contrasts

%Average relevant information
mean_bump_data = varfun(@mean,summarydata,'InputVariables','bump_mag',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_bump_data,1),1); 
for contrast = 1:length(mean_bump_data.contrast_level)
    if mean_bump_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_bump_data.stim_vel,mean_bump_data.mean_bump_mag,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_bump_data.mean_bump_mag) - 0.2, max(mean_bump_data.mean_bump_mag) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20','30','60'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[~,h2]=suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Bump magnitude (DF/F)');

%save
saveas(gcf,[path,'\analysis\continuous_plots\mean_bump_mag.png']);

%% Compare the width at half max for each speed between contrasts

%Average relevant information
mean_bw_data = varfun(@mean,summarydata,'InputVariables','bump_width',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_bw_data,1),1); 
for contrast = 1:length(mean_bw_data.contrast_level)
    if mean_bw_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_bw_data.stim_vel,mean_bw_data.mean_bump_width,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_bw_data.mean_bump_width) - 0.2, max(mean_bw_data.mean_bump_width) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20','30','60'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[~,h2]=suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Bump width at half max (rad)');

%save
saveas(gcf,[path,'\analysis\continuous_plots\mean_bump_bw.png']);


%% Heatmap with offset vs movement

for session = 1:length(open_loop_sessions)
    
    figure('Position',[100, 100, 1000, 800]),
    
    %PB heatmap
    subplot(3,5,[1 4])
    dff_matrix = data{1,session}.continuous_data.dff_matrix';
    %I will now flip the matrix to match the fly's heading
    imagesc(flip(dff_matrix))
    colormap(flipud(gray))
    ylabel('PB glomerulus','fontsize',12)
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    set(gca,'yticklabel',{[]});
    set(gca,'YTick',[]);
    title('EPG activity in the PB','fontsize',12);
    
    %Fly heading
    subplot(3,5,[6 9])
    %Heading
    heading = wrapTo180(data{1,session}.continuous_data.heading_deg);
    [x_out_heading,heading_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,heading);  
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_heading,heading_to_plot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
    else
        plot(x_out_heading,heading_to_plot,'color',[0 0 0.6],'LineWidth',1.5)
    end  
    hold on
    %Phase
    %I'm negating the phase because it should go against the heading
    phase = -wrapTo180(rad2deg(data{1,session}.continuous_data.bump_pos'));
    [x_out_phase,phase_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,phase);
    plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
    axP = get(gca,'Position');
    legend('Fly heading','EPG phase','Location','EastOutside');
    set(gca, 'Position', axP);
    if ~isnan(x_out_phase(end))
        xlim([0, x_out_phase(end)]);
    else
        xlim([0, x_out_phase(end-1)]);
    end
    ylim([-180 180]);
    ylabel('Angular position (deg)','fontsize',12);
    set(gca,'xticklabel',{[]})
    title('Fly heading','fontsize',12);
    
    %Mvt offset
    mvt_offset{session} = wrapTo180(rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos',-data{1,session}.continuous_data.heading)));
    [x_out_offset,offsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,mvt_offset{session});
    subplot(3,5,[11 14])
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_offset,offsetToPlot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
    else
        plot(x_out_offset,offsetToPlot,'color',[0 0 0.6],'LineWidth',1.5)
    end
    if ~isnan(x_out_offset(end))
        xlim([0, x_out_offset(end)]);
    else
        xlim([0, x_out_offset(end-1)]);
    end
    ylim([-180 180]);
    xlabel('Time (sec)','fontsize',12); ylabel('Degrees','fontsize',12);
    title('Movement offset','fontsize',12);
    
    %Plot stimulus offset distribution
    subplot(3,5,15)
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        polarhistogram(deg2rad(mvt_offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','High contrast'});
        Ax = gca;
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    else
        polarhistogram(deg2rad(mvt_offset{session}),20,'FaceColor',[ 0 0 0.6])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','Low contrast'});
        Ax = gca; 
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    end  

    
   %save figures
   saveas(gcf,[path,'\analysis\continuous_plots\open_loop_heatmap_and_mvt_offset_sid',num2str(open_loop_sessions(session)),'.png']);
end

%% Get movement offset variation 

close all;
warning('off');

mvt_offset_var = zeros(1,length(open_loop_sessions));
for session = 1:length(open_loop_sessions)
    [~, mvt_offset_var(session)] = circ_std(mvt_offset{session},[],[],1);
end

mvt_offset_var = mvt_offset_var';
%add variable to existing table
if sum(strcmpi(summarydata.Properties.VariableNames,'mvt_offset_var')) == 0
    summarydata = addvars(summarydata,mvt_offset_var); 
end

%% Compare the movement offset variation for each speed between contrasts

%Average relevant information
mean_mvt_offset_var_data = varfun(@mean,summarydata,'InputVariables','mvt_offset_var',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_mvt_offset_var_data,1),1); 
for contrast = 1:length(mean_mvt_offset_var_data.contrast_level)
    if mean_mvt_offset_var_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_mvt_offset_var_data.stim_vel,mean_mvt_offset_var_data.mean_mvt_offset_var,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_mvt_offset_var_data.mean_mvt_offset_var) - 0.2, max(mean_mvt_offset_var_data.mean_mvt_offset_var) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20','30','60'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[ax,h2] = suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Movement offset variability (rad)');

%Save
saveas(gcf,[path,'\analysis\continuous_plots\open_loop_mvt_offset_var.png']);

%% Plots with heatmap, stimulus position, fly heading, and offset

for session = 1:length(open_loop_sessions)
    
    figure('Position',[100, 100, 1000, 800]),
    
    %PB heatmap
    subplot(5,1,1)
    dff_matrix = data{1,session}.continuous_data.dff_matrix';
    imagesc(flip(dff_matrix))
    colormap(flipud(gray))
    ylabel('PB glomerulus')
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    set(gca,'yticklabel',{[]});
    set(gca,'YTick',[]);
    title('EPG activity in the PB');
    
    %Stimulus position
    subplot(5,1,2)
    bar_position = wrapTo180(data{1,session}.continuous_data.panel_angle);
    [x_out_bar,bar_position_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,bar_position);
    %Get EPG phase to plot
    phase = wrapTo180(rad2deg(data{1,session}.continuous_data.bump_pos'));
    [x_out_phase,phase_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,phase);
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_bar,bar_position_to_plot,'color',[0.5 0.8 0.9],'lineWidth',1.5)
    else
        plot(x_out_bar,bar_position_to_plot,'color',[0 0 0.6],'lineWidth',1.5)
    end
    hold on
    plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'lineWidth',1.5)
    legend('Stimulus position','EPG phase');
    if ~isnan(x_out_phase(end))
        xlim([0, x_out_phase(end)]);
    else
        xlim([0, x_out_phase(end-1)]);
    end
    ylim([-180 180]);
    ylabel('Angular position (deg)');
    set(gca,'xticklabel',{[]})
    title('Stimulus position');
    
    %Stim offset
    %Calculate stim offset as the circular distance between EPG activity
    %phase and the stimulus posotion
    stim_offset{session} = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos',deg2rad(data{1,session}.continuous_data.panel_angle)));
    [x_out_stimoffset,stimOffsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,stim_offset{session});
    subplot(5,1,3)
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_stimoffset,stimOffsetToPlot,'color',[ 0.5 0.8 0.9],'lineWidth',1.5)
    else
        plot(x_out_stimoffset,stimOffsetToPlot,'color',[0 0 0.6],'lineWidth',1.5)
    end  
    if ~isnan(x_out_stimoffset(end))
        xlim([0, x_out_stimoffset(end)]);
    else
        xlim([0, x_out_stimoffset(end-1)]);
    end
    ylim([-180 180]);
    set(gca,'xticklabel',{[]})
    ylabel('Degrees');
    title('Stimulus offset');
    
     
    %add movement plots
    %Fly heading
    subplot(5,1,4)
    %Heading
    heading = wrapTo180(-data{1,session}.continuous_data.heading_deg);
    [x_out_heading,heading_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,heading); 
    %Plot using different colors if the stimulus was low or high contrast
    plot(x_out_heading,heading_to_plot,'color',[ 0.6 0.8 0.2],'lineWidth',1.5)
    hold on
    %Phase
    plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'lineWidth',1.5)
    legend('-Fly heading','EPG phase');
    if ~isnan(x_out_heading(end))
        xlim([0, x_out_heading(end)]);
    else
        xlim([0, x_out_heading(end-1)]);
    end
    ylim([-180 180]);
    ylabel('Angular position (deg)');
    set(gca,'xticklabel',{[]})
    title('Fly heading');
    
    %Mvt offset
    mvt_offset{session} = wrapTo180(rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos',-data{1,session}.continuous_data.heading)));
    [x_out_offset,offsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,mvt_offset{session});
    subplot(5,1,5)
    %Plot using different colors if the stimulus was low or high contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        plot(x_out_offset,offsetToPlot,'color',[ 0.5 0.8 0.9],'lineWidth',1.5)
    else
        plot(x_out_offset,offsetToPlot,'color',[0 0 0.6],'lineWidth',1.5)
    end
    if ~isnan(x_out_offset(end))
        xlim([0, x_out_offset(end)]);
    else
        xlim([0, x_out_offset(end-1)]);
    end
    ylim([-180 180]);
    xlabel('Time (sec)'); ylabel('Degrees');
    title('Movement offset');
       
    %save figures
    saveas(gcf,[path,'\analysis\continuous_plots\open_loop_combined_heatmap_sid',num2str(open_loop_sessions(session)),'.png']);
end

%% Plot BM as function of contrast

%Average relevant information
mean_bm_data = varfun(@mean,summarydata,'InputVariables','bump_mag',...
       'GroupingVariables',{'contrast_level'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_bm_data,1),1); 
for contrast = 1:length(mean_bm_data.contrast_level)
    if mean_bm_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_bm_data.contrast_level,mean_bm_data.mean_bump_mag,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([55 58]);
ylim([min(mean_bm_data.mean_bump_mag) - 0.2, max(mean_bm_data.mean_bump_mag) + 0.2]);
xticks(56:57);
xticklabels({'Low contrast','High contrast'});
ylabel('Bump magntiude (DF/F)');

%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\open_loop_mean_bm_vs_contrast.png']);

%% Plot bw as function of contrast

%Average relevant information
mean_bw_data = varfun(@mean,summarydata,'InputVariables','bump_width',...
       'GroupingVariables',{'contrast_level'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_bw_data,1),1); 
for contrast = 1:length(mean_bw_data.contrast_level)
    if mean_bw_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_bw_data.contrast_level,mean_bw_data.mean_bump_width,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([55 58]);
ylim([min(mean_bw_data.mean_bump_width) - 0.2, max(mean_bw_data.mean_bump_width) + 0.2]);
xticks(56:57);
xticklabels({'Low contrast','High contrast'});
ylabel('Bump width (rad)');

%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\open_loop_mean_bw_vs_contrast.png']);

%% Model bump velocity as a function of visual and proprioceptive cues

all_bump_vel = [];
all_fly_vel = [];
all_bump_estimate = [];
all_stim_position = [];
all_stim_contrast = [];

for session = 1:length(open_loop_sessions)
    
    %Compute bump velocity
    bump_pos = data{1,session}.continuous_data.bump_pos';
    unwrapped_bump_pos = unwrap(bump_pos);
    smooth_bump_pos = smooth(rad2deg(unwrapped_bump_pos));
    bump_vel = diff(smooth_bump_pos).*(length(data{1,session}.continuous_data.time)/data{1,session}.continuous_data.time(end));
    smooth_bump_vel = smooth(bump_vel);
    all_bump_vel = [all_bump_vel;smooth_bump_vel];
 
    %Compute fly angular velocity
    fly_heading = -data{1,session}.continuous_data.heading;
    unwrapped_fly_heading = unwrap(fly_heading);
    smooth_fly_heading = smooth(rad2deg(unwrapped_fly_heading));
    fly_vel = diff(smooth_fly_heading).*(length(data{1,session}.continuous_data.time)/data{1,session}.continuous_data.time(end));
    smooth_fly_vel = smooth(fly_vel);
    all_fly_vel = [all_fly_vel;smooth_fly_vel];

    %Compute bump estimate of the stimulus' position (this will be a
    %differnce between the bump's position and this fly's offset, which we
    %have stored from the closed-loop session
        %1) Import reference offset
    load([path,'\analysis\continuous_summary_data.mat']);       
        %2) Compute bump estimate using the offset
    bump_estimate = wrapTo180(rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos', deg2rad(mean_reference_offset))));
    all_bump_estimate = [all_bump_estimate;bump_estimate(2:end)];

    %Get stimulus position
    stim_position = wrapTo180(data{1,session}.continuous_data.panel_angle);
    all_stim_position = [all_stim_position;stim_position(2:end)];
    
    %Get contrast
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        stim_contrast = repelem(2,length(bump_estimate),1);
    else
        stim_contrast = repelem(1,length(bump_estimate),1);
    end
    all_stim_contrast = [all_stim_contrast;stim_contrast(2:end)];
end

%Make stim contrast a categorical variable
all_stim_contrast = categorical(all_stim_contrast);

%Compute the visual cue drive (this is the difference between the
%stimulus position and the bump's estimate of it
all_visual_cue_drive = [wrapTo180(rad2deg(circ_dist(deg2rad(all_stim_position),deg2rad(all_bump_estimate))))];

%Combine variables in table
bump_vel_model_data = table(all_fly_vel,all_visual_cue_drive,all_stim_contrast,all_bump_vel,'VariableNames',{'FlyAngVel','VisualCueDrive','contrast','BumpAngVel'});

%Fit models
bump_vel_model = fitlm(bump_vel_model_data,'BumpAngVel ~ VisualCueDrive:contrast + FlyAngVel:contrast','Intercept',false)
bump_vel_model2 = fitlm(bump_vel_model_data,'BumpAngVel ~ VisualCueDrive:contrast + FlyAngVel','Intercept',false)
bump_vel_model3 = fitlm(bump_vel_model_data,'BumpAngVel ~ VisualCueDrive*contrast + FlyAngVel*contrast','Intercept',false)
bump_vel_model4 = fitlm(bump_vel_model_data,'BumpAngVel ~ VisualCueDrive*contrast + FlyAngVel','Intercept',false)
bump_vel_model5 = fitlm(bump_vel_model_data,'BumpAngVel ~ VisualCueDrive:contrast + FlyAngVel:contrast + VisualCueDrive + FlyAngVel','Intercept',false)

%% Save relevant data

save([path,'\analysis\continuous_open_loop_data.mat'],'summarydata','bump_vel_model_data');

clear all; close all; clc