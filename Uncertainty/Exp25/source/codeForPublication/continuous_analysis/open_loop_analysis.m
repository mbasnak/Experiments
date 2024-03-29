%Code to analyze the open-loop bouts

clear all; close all;

%% Load data

%Get directory you're interested in
path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts');

%I will create a .txt file for each fly with the session numbers
%Read the .txt file containing the ID of the open loop sessions.
fileID = fopen([path,'\2p\open_loop_sessions.txt'],'r');
formatSpec = '%d';
open_loop_sessions = fscanf(fileID,formatSpec);

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
    ylabel('PB glomerulus')
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    set(gca,'yticklabel',{[]});
    set(gca,'YTick',[]);
    title('EPG activity in the PB');
    
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
    ylabel('Degrees');
    set(gca,'xticklabel',{[]})
    title('Stimulus position');
    
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
    xlabel('Time (sec)'); ylabel('Degrees');
    title('Stimulus offset');
    
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
    saveas(gcf,[path,'\analysis\continuous_plots\openLoopHeatmapAndOffset_session',num2str(session),'.png']);
end


%% Summarize relevant data in table

close all;
%Calculate offset variation as the circular standard deviation of the
%'stimulus offset'
for session = 1:length(open_loop_sessions)
   [~,stim_offset_var(session)] = circ_std(deg2rad(stim_offset{session}),[],[],1);  
end

%Define stimulus velocity conditions-which panel functions were used?
low_stim_vel = [50,51];
medium_stim_vel = [52,53];
high_stim_vel = [54,55];

%Define stimulus contrasts-which panel patterns were used?
low_stim_contrast = 56;
high_stim_contrast = 57;

%Compute bump magnitude as max-min
for session = 1:length(open_loop_sessions)
    BumpMagnitude{session} = data{1,session}.continuous_data.bump_magnitude;
    meanBM(session) = mean(BumpMagnitude{session});   
end

%compute bump width at half max
for session = 1:length(open_loop_sessions)    
    half_max_width{session} = data{1,session}.continuous_data.bump_width;        
    meanHW(session) = mean(half_max_width{session});  
end


%Obtain for each session a table with their offset var, their median bump magnitude their stim
%velocity, and their contrast
summarydata = array2table(zeros(0,5), 'VariableNames',{'offset_var','bump_mag','half_width','contrast_level','stim_vel'});
warning('off');
for session = 1:length(open_loop_sessions)
    summarydata{session,'offset_var'} = stim_offset_var(session);
    summarydata{session,'bump_mag'} = meanBM(session);
    summarydata{session,'half_width'} = meanHW(session);
    summarydata{session,'contrast_level'} = data{1,session}.continuous_data.run_obj.pattern_number;
    if (data{1,session}.continuous_data.run_obj.function_number == 50 | data{1,session}.continuous_data.run_obj.function_number == 51)
        summarydata{session,'stim_vel'} = 1; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 52 | data{1,session}.continuous_data.run_obj.function_number == 53)
        summarydata{session,'stim_vel'} = 2; 
    else
        summarydata{session,'stim_vel'} = 3; 
    end
end

%% Compare the offset variation for each speed between contrasts

%Average relevant information
mean_data = varfun(@mean,summarydata,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_data,1),1); 
for contrast = 1:length(mean_data.contrast_level)
    if mean_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_data.stim_vel,mean_data.mean_offset_var,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_data.mean_offset_var) - 0.2, max(mean_data.mean_offset_var) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20 deg/s','30 deg/s','60 deg/s'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Offset variation (circ std)','FontSize',12);

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
set(gca,'xticklabel',{'20 deg/s','30 deg/s','60 deg/s'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel({'Bump magnitude';'(amplitude of Fourier component)'});

%save
saveas(gcf,[path,'\analysis\continuous_plots\mean_bump_mag.png']);


%% Compare the width at half max for each speed between contrasts

%Average relevant information
mean_hw_data = varfun(@mean,summarydata,'InputVariables','half_width',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_hw_data,1),1); 
for contrast = 1:length(mean_hw_data.contrast_level)
    if mean_hw_data.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_hw_data.stim_vel,mean_hw_data.mean_half_width,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_hw_data.mean_half_width) - 0.2, max(mean_hw_data.mean_half_width) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20 deg/s','30 deg/s','60 deg/s'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Bump width at half max');

%save
saveas(gcf,[path,'\analysis\continuous_plots\mean_bump_hw.png']);


%% Heatmap with offset vs movement

for session = 1:length(open_loop_sessions)
    
    figure('Position',[100, 100, 1000, 800]),
    
    %PB heatmap
    subplot(3,5,[1 4])
    %I will now flip the matrix to match the fly's heading
    imagesc(flip(dff_matrix))
    colormap(flipud(gray))
    ylabel('PB glomerulus')
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    set(gca,'yticklabel',{[]});
    set(gca,'YTick',[]);
    title('EPG activity in the PB');
    
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
    legend('Fly heading','EPG phase');
    if ~isnan(x_out_phase(end))
        xlim([0, x_out_phase(end)]);
    else
        xlim([0, x_out_phase(end-1)]);
    end
    ylim([-180 180]);
    ylabel('Degrees');
    set(gca,'xticklabel',{[]})
    title('Fly heading');
    
    %Mvt offset
    offset{session} = wrapTo180(data{1,session}.continuous_data.offset);
    [x_out_offset,offsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,offset{session});
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
    xlabel('Time (sec)'); ylabel('Degrees');
    title('Movement offset');
    
    %Plot stimulus offset distribution
    subplot(3,5,15)
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        polarhistogram(deg2rad(offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','High contrast'});
        Ax = gca;
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    else
        polarhistogram(deg2rad(offset{session}),20,'FaceColor',[ 0 0 0.6])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','Low contrast'});
        Ax = gca; 
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    end  

    
   %save figures
   saveas(gcf,[path,'\analysis\continuous_plots\openLoopHeatmapAndMovement_session',num2str(session),'.png']);
end

%% Get movement offset variation 

close all;
warning('off');

for session = 1:length(open_loop_sessions)
    [~, mvt_offset_var(session)] = circ_std(offset{session},[],[],1);
end

mvt_offset_var = mvt_offset_var';
summarydata = addvars(summarydata,mvt_offset_var);

%% Compare the movement offset variation for each speed between contrasts

%Average relevant information
mean_mvt_offset = varfun(@mean,summarydata,'InputVariables','mvt_offset_var',...
       'GroupingVariables',{'contrast_level','stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
colorID = cell(size(mean_mvt_offset,1),1); 
for contrast = 1:length(mean_mvt_offset.contrast_level)
    if mean_mvt_offset.contrast_level(contrast) == 57
        colorID{contrast} = [ 0.5 0.8 0.9]; 
    else
        colorID{contrast} = [ 0 0 0.6]; 
    end
end
scatter(mean_mvt_offset.stim_vel,mean_mvt_offset.mean_mvt_offset_var,60,cell2mat(colorID),'filled')
%Add custom legend
hold on
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
xlim([0, 4]);
ylim([min(mean_mvt_offset.mean_mvt_offset_var) - 0.2, max(mean_mvt_offset.mean_mvt_offset_var) + 0.2]);
xticks(1:3);
set(gca,'xticklabel',{'20 deg/s','30 deg/s','60 deg/s'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',10,'FontWeight','bold');
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Movement offset variation (circ std)');

%Save
saveas(gcf,[path,'\analysis\continuous_plots\open_loop_mvt_offset_var.png']);

%% Plots with heatmap, stimulus position, fly heading, and offset

for session = 1:length(open_loop_sessions)
    
    figure('Position',[100, 100, 1000, 800]),
    
    %PB heatmap
    subplot(5,4,[1 3])
    imagesc(flip(dff_matrix))
    colormap(flipud(gray))
    ylabel('PB glomerulus')
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    set(gca,'yticklabel',{[]});
    set(gca,'YTick',[]);
    title('EPG activity in the PB');
    
    %Stimulus position
    subplot(5,4,[5 7])
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
    ylabel('Degrees');
    set(gca,'xticklabel',{[]})
    title('Stimulus position');
    
    %Stim offset
    %Calculate stim offset as the circular distance between EPG activity
    %phase and the stimulus posotion
    stim_offset{session} = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos',deg2rad(data{1,session}.continuous_data.panel_angle)));
    [x_out_stimoffset,stimOffsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,stim_offset{session});
    subplot(5,4,[9 11])
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
    
    %Plot stimulus offset distribution
    subplot(5,4,12)
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
    
    %add movement plots
    %Fly heading
    subplot(5,4,[13 15])
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
    ylabel('Degrees');
    set(gca,'xticklabel',{[]})
    title('Fly heading');
    
    %Mvt offset
    offset{session} = wrapTo180(data{1,session}.continuous_data.offset);
    [x_out_offset,offsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,offset{session});
    subplot(5,4,[17 19])
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
    
    %Plot stimulus offset distribution
    subplot(5,4,20)
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        polarhistogram(deg2rad(offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','High contrast'});
        Ax = gca;
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    else
        polarhistogram(deg2rad(offset{session}),20,'FaceColor',[ 0 0 0.6])
        set(gca,'ThetaZeroLocation','top');
        title({'Offset distribution','Low contrast'});
        Ax = gca; 
        Ax.RGrid = 'off';
        Ax.RTickLabel = [];
    end  
    
    %save figures
    saveas(gcf,[path,'\analysis\continuous_plots\openLoopCombinedHeatmap_session',num2str(session),'.png']);
end

%% Model BM and HW as a function of total movement, time and contrast level

close all;
%create table with relevant variables
allBumpMag = [];
allHalfWidth = [];
allTotalMvt = [];
allZscoredMvt = [];
allTime = [];
contrastLevel = [];
for session = 1:length(open_loop_sessions)
   allBumpMag = [allBumpMag,BumpMagnitude{1,session}];
   allHalfWidth = [allHalfWidth,half_max_width{1,session}]; 
   allTotalMvt = [allTotalMvt,fillmissing(data{1,session}.continuous_data.total_mvt_ds,'linear')];
   allZscoredMvt = [allZscoredMvt,zscore(fillmissing(data{1,session}.continuous_data.total_mvt_ds,'linear'))];
   allTime = [allTime,data{1,session}.continuous_data.time'];
   if data{1,session}.continuous_data.run_obj.pattern_number == 57
       contrastLevel = [contrastLevel,repelem(2,length(data{1,session}.continuous_data.time))];
   else
       contrastLevel = [contrastLevel,repelem(1,length(data{1,session}.continuous_data.time))];
   end
end
bump_mag_data = table(nominal(contrastLevel)',allTime',allTotalMvt',allZscoredMvt',allBumpMag','VariableNames',{'ContrastLevel','Time','TotalMovement','ZscoredMvt','BumpMagnitude'});
half_width_data = table(nominal(contrastLevel)',allTime',allTotalMvt',allZscoredMvt',allHalfWidth','VariableNames',{'ContrastLevel','Time','TotalMovement','ZscoredMvt','HalfWidth'});

mdl_BM = fitlm(bump_mag_data,'BumpMagnitude ~ ContrastLevel+ZscoredMvt')
mdl_HW = fitlm(half_width_data,'HalfWidth ~ ContrastLevel+ZscoredMvt')


%% Plot BM as function of contrast

%Average relevant information
mean_BM = varfun(@mean,bump_mag_data,'InputVariables','BumpMagnitude',...
       'GroupingVariables',{'ContrastLevel'});
   
figure,
plot(1,mean_BM.mean_BumpMagnitude(1),'o','color',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',8)
hold on
plot(2,mean_BM.mean_BumpMagnitude(2),'o','color',[ 0.5 0.8 0.9],'MarkerFaceColor',[ 0.5 0.8 0.9],'MarkerSize',8)
xticks(1:2);
xticklabels({'Low contrast','High contrast'});
xlim([0 3]);
ylabel({'Bump magnitude';'(amplitude of Fourier component)'});
ylim([0 max(mean_BM.mean_BumpMagnitude)+0.4]);

%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\openLoopMeanBMvsContrast.png']);

%% Plot HW as function of contrast

%Average relevant information
mean_HW = varfun(@mean,half_width_data,'InputVariables','HalfWidth',...
       'GroupingVariables',{'ContrastLevel'});
   
figure,
plot(1,mean_HW.mean_HalfWidth(1),'o','color',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',8)
hold on
plot(2,mean_HW.mean_HalfWidth(2),'o','color',[ 0.5 0.8 0.9],'MarkerFaceColor',[ 0.5 0.8 0.9],'MarkerSize',8)
xticks(1:2);
xticklabels({'Low contrast','High contrast'});
xlim([0 3]);
ylabel('Mean bump width at half max');
ylim([0 max(mean_HW.mean_HalfWidth)+0.4]);

%Save figure
saveas(gcf,[path,'\analysis\continuous_plots\openLoopMeanHWvsContrast.png']);

%% Model bump velocity as a function of visual and proprioceptive cues

all_bump_vel_HC = [];
all_fly_vel_HC = [];
all_bump_estimate_HC = [];
all_stim_position_HC = [];

all_bump_vel_LC = [];
all_fly_vel_LC = [];
all_bump_estimate_LC = [];
all_stim_position_LC = [];

for session = 1:length(open_loop_sessions)
    
    %Get bump velocity
    bump_pos = data{1,session}.continuous_data.bump_pos';
    unwrapped_bump_pos = unwrap(bump_pos);
    smooth_bump_pos = smooth(rad2deg(unwrapped_bump_pos));
    bump_vel = diff(smooth_bump_pos).*(length(data{1,session}.continuous_data.time)/data{1,session}.continuous_data.time(end));
    smooth_bump_vel = smooth(bump_vel);
    %Save in the corresponding array
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        all_bump_vel_HC = [all_bump_vel_HC;smooth_bump_vel];
    else
        all_bump_vel_LC = [all_bump_vel_LC;smooth_bump_vel];
    end
    
    %Get fly angular velocity
    fly_heading = -data{1,session}.continuous_data.heading;
    unwrapped_fly_heading = unwrap(fly_heading);
    smooth_fly_heading = smooth(rad2deg(unwrapped_fly_heading));
    fly_vel = diff(smooth_fly_heading).*(length(data{1,session}.continuous_data.time)/data{1,session}.continuous_data.time(end));
    smooth_fly_vel = smooth(fly_vel);
    %Save in the corresponding array
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        all_fly_vel_HC = [all_fly_vel_HC;smooth_fly_vel];
    else
        all_fly_vel_LC = [all_fly_vel_LC;smooth_fly_vel];
    end
    
    %Import reference offset
    load([path,'\analysis\continuous_summary_data.mat']);
        
    %Compute bump estimate using the offset
    bump_estimate = wrapTo180(rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos', deg2rad(mean_reference_offset))));
    
    %Save in the corresponding array
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        all_bump_estimate_HC = [all_bump_estimate_HC;bump_estimate(2:end)];
    else
        all_bump_estimate_LC = [all_bump_estimate_LC;bump_estimate(2:end)];
    end

    %Get stimulus position
    stim_position = wrapTo180(data{1,session}.continuous_data.panel_angle);
        %Save in the corresponding array
    if data{1,session}.continuous_data.run_obj.pattern_number == 57
        all_stim_position_HC = [all_stim_position_HC;stim_position(2:end)];
    else
        all_stim_position_LC = [all_stim_position_LC;stim_position(2:end)];
    end
    
end

all_stim_difference_HC = [wrapTo180(rad2deg(circ_dist(deg2rad(all_stim_position_HC),deg2rad(all_bump_estimate_HC))))];
all_stim_difference_LC = [wrapTo180(rad2deg(circ_dist(deg2rad(all_stim_position_LC),deg2rad(all_bump_estimate_LC))))];

%Combine variables in table
bump_vel_data_HC = table(all_fly_vel_HC,all_stim_difference_HC,all_bump_vel_HC,'VariableNames',{'FlyAngVel','VisualCueDrive','BumpAngVel'});
bump_vel_data_LC = table(all_fly_vel_LC,all_stim_difference_LC,all_bump_vel_LC,'VariableNames',{'FlyAngVel','VisualCueDrive','BumpAngVel'});

%Fit models
bump_vel_model_HC = fitlm(bump_vel_data_HC,'Intercept',false)
bump_vel_model_LC = fitlm(bump_vel_data_LC,'Intercept',false)


%% Fit single model for both

high_contrast = repelem(2,size(bump_vel_data_HC,1),1);
bump_vel_data_HC = addvars(bump_vel_data_HC,nominal(high_contrast),'NewVariableName','contrast');

low_contrast = repelem(1,size(bump_vel_data_LC,1),1);
bump_vel_data_LC = addvars(bump_vel_data_LC,nominal(low_contrast),'NewVariableName','contrast');

%combine both tables
bump_vel_data = [bump_vel_data_LC;bump_vel_data_HC];

%Fit models
bump_vel_model = fitlm(bump_vel_data,'BumpAngVel ~ VisualCueDrive:contrast + FlyAngVel:contrast','Intercept',false)
bump_vel_model2 = fitlm(bump_vel_data,'BumpAngVel ~ VisualCueDrive:contrast + FlyAngVel','Intercept',false)
bump_vel_model3 = fitlm(bump_vel_data,'BumpAngVel ~ VisualCueDrive*contrast + FlyAngVel*contrast','Intercept',false)
bump_vel_model4 = fitlm(bump_vel_data,'BumpAngVel ~ VisualCueDrive*contrast + FlyAngVel','Intercept',false)
bump_vel_model5 = fitlm(bump_vel_data,'BumpAngVel ~ VisualCueDrive:contrast + FlyAngVel:contrast + VisualCueDrive + FlyAngVel','Intercept',false)


%% Save relevant data

save([path,'\analysis\continuous_open_loop_data.mat'],'bump_mag_data','half_width_data','mean_BM','summarydata','mean_data','bump_vel_model_HC','bump_vel_model_LC','mean_bump_data','bump_vel_data_HC','bump_vel_data_LC');

clear all; close all; clc