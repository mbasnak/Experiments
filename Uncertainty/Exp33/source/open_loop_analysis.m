%Code to analyze the open-loop bouts

clear all; close all;

%% Load data

%Get directory you're interested in
path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp33\data');

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
    imagesc(data{1,session}.continuous_data.dff_matrix')
    colormap(flipud(gray))
    ylabel('PB distance')
    set(gca,'xticklabel',{[]});
    set(gca,'XTick',[]);
    title('EPG activity in the PB');
    
    %Stimulus position
    subplot(3,5,[6 9])
    bar_position = wrapTo180(data{1,session}.continuous_data.panel_angle);
    [x_out_bar, bar_position_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time, bar_position);
    %Get EPG phase to plot
    phase = wrapTo180(rad2deg(data{1,session}.continuous_data.bump_pos));
    [x_out_phase, phase_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time, phase');
    %Plot using different colors if the stimulus was low or high contrast
    plot(x_out_bar,bar_position_to_plot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
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
    stim_offset{session} = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos,deg2rad(data{1,session}.continuous_data.panel_angle')));
    stimoffset = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos,deg2rad(data{1,session}.continuous_data.panel_angle')));    
    [x_out_stimoffset, stimOffsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time, stimoffset');
    subplot(3,5,[11 14])
    %Plot using different colors if the stimulus was low or high contrast
    plot(x_out_stimoffset,stimOffsetToPlot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
    xlim([0, data{1,session}.continuous_data.time(end)]);
    ylim([-180 180]);
    xlabel('Time (sec)'); ylabel('Degrees');
    title('Stimulus offset');
    
    %Plot stimulus offset distribution
    subplot(3,5,15)
    polarhistogram(deg2rad(stim_offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
    set(gca,'ThetaZeroLocation','top');
    title({'Offset distribution','High contrast'});
    Ax = gca;
    Ax.RGrid = 'off';
    Ax.RTickLabel = [];
    
    %Save figures
    saveas(gcf,[path,'\analysis\plots\openLoopHeatmapAndOffset_session',num2str(session),'.png']);
end

%% Summarize relevant data in table

close all;
%Calculate offset variation as the circular standard deviation of the
%'stimulus offset'
for session = 1:length(open_loop_sessions)
   [~,stim_offset_var(session)] = circ_std(deg2rad(stim_offset{session}),[],[],2);  
end

%Define stimulus velocity conditions-which panel functions were used?
stim_vel_30 = [52,53];
stim_vel_60 = [54,55];
stim_vel_90 = [206,207];
stim_vel_120 = [208,209];
stim_vel_150 = [210,211];
stim_vel_180 = [212,213];
stim_vel_210 = [214,215];


%Obtain bump parameters
for session = 1:length(open_loop_sessions)
    BumpMagnitude{session} = data{1,session}.continuous_data.bump_magnitude;
    meanBM(session) = mean(BumpMagnitude{session});   

    BumpWidth{session} = data{1,session}.continuous_data.bump_width;
    meanBW(session) = mean(BumpWidth{session});   
end


%Obtain for each session a table with their offset var, their mean bump magnitude and bump width and their stim
%velocity
summarydata = array2table(zeros(0,4), 'VariableNames',{'offset_var','bump_mag','bump_width','stim_vel'});
warning('off');
for session = 1:length(open_loop_sessions)
    summarydata{session,'offset_var'} = stim_offset_var(session);
    summarydata{session,'bump_mag'} = meanBM(session);
    summarydata{session,'bump_width'} = meanBW(session);
    if (data{1,session}.continuous_data.run_obj.function_number == 52 | data{1,session}.continuous_data.run_obj.function_number == 53)
        summarydata{session,'stim_vel'} = 30; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 54 | data{1,session}.continuous_data.run_obj.function_number == 55)
        summarydata{session,'stim_vel'} = 60; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 206 | data{1,session}.continuous_data.run_obj.function_number == 207)
        summarydata{session,'stim_vel'} = 90; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 208 | data{1,session}.continuous_data.run_obj.function_number == 209)
        summarydata{session,'stim_vel'} = 120; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 210 | data{1,session}.continuous_data.run_obj.function_number == 211)
        summarydata{session,'stim_vel'} = 150; 
    elseif (data{1,session}.continuous_data.run_obj.function_number == 212 | data{1,session}.continuous_data.run_obj.function_number == 213)
        summarydata{session,'stim_vel'} = 180; 
    else
        summarydata{session,'stim_vel'} = 210; 
    end
end

%% Compare the offset variation across speeds

%Average relevant information
mean_data = varfun(@mean,summarydata,'InputVariables','offset_var',...
       'GroupingVariables',{'stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
scatter(mean_data.stim_vel,mean_data.mean_offset_var,60)
ylim([min(mean_data.mean_offset_var) - 0.2, max(mean_data.mean_offset_var) + 0.2]);
[ax,h2]=suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Offset variation (circ std)','FontSize',12);

%save
saveas(gcf,[path,'\analysis\plots\open_loop_offset_var.png']);


%% Compare the bump magnitude for each speed between contrasts

%Average relevant information
mean_bump_data = varfun(@mean,summarydata,'InputVariables','bump_mag',...
       'GroupingVariables',{'stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
scatter(mean_bump_data.stim_vel,mean_bump_data.mean_bump_mag,60)
ylim([min(mean_bump_data.mean_bump_mag) - 0.2, max(mean_bump_data.mean_bump_mag) + 0.2]);
[ax,h2]=suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel({'Bump magnitude';'(from von Mises fit'});

%save
saveas(gcf,[path,'\analysis\plots\mean_bump_mag.png']);


%% Compare the width at half max for each speed between contrasts

%Average relevant information
mean_hw_data = varfun(@mean,summarydata,'InputVariables','bump_width',...
       'GroupingVariables',{'stim_vel'});
   
%Plot
figure('Position',[200 200 800 800]),
scatter(mean_hw_data.stim_vel,mean_hw_data.mean_bump_width,60)
ylim([min(mean_hw_data.mean_bump_width) - 0.2, max(mean_hw_data.mean_bump_width) + 0.2]);
[ax,h2]=suplabel('Stimulus angular velocity (deg/s)','x');
set(h2,'FontSize',12,'FontWeight','bold')
ylabel('Bump width at half max');

%save
saveas(gcf,[path,'\analysis\plots\mean_bump_hw.png']);


%% Plots with heatmap, stimulus position, fly heading, and offset

for session = 1:length(open_loop_sessions)
    
    figure('Position',[100, 100, 1000, 800]),
    
    %PB heatmap
    subplot(5,4,[1 3])
    imagesc(data{1,session}.continuous_data.dff_matrix')
    colormap(flipud(gray))
    ylabel('PB glomerulus')
    set(gca,'xticklabel',{[]});
    yticks(1:2:16);
    yticklabels({'8R','6R','4R','2R','1L','3L','5L','7L'});
    title('EPG activity in the PB');
    
    %Stimulus position
    subplot(5,4,[5 7])
    bar_position = wrapTo180(data{1,session}.continuous_data.panel_angle);
    [x_out_bar,bar_position_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,bar_position);
    %Get EPG phase to plot
    phase = wrapTo180(rad2deg(data{1,session}.continuous_data.bump_pos));
    [x_out_phase,phase_to_plot] = removeWrappedLines(data{1,session}.continuous_data.time,phase');
    %Plot using different colors if the stimulus was low or high contrast
    plot(x_out_bar,bar_position_to_plot,'color',[0.5 0.8 0.9],'lineWidth',1.5)
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
    stim_offset{session} = rad2deg(circ_dist(data{1,session}.continuous_data.bump_pos,deg2rad(data{1,session}.continuous_data.panel_angle')));
    [x_out_stimoffset,stimOffsetToPlot] = removeWrappedLines(data{1,session}.continuous_data.time,stim_offset{session}');
    subplot(5,4,[9 11])
    %Plot using different colors if the stimulus was low or high contrast
    plot(x_out_stimoffset,stimOffsetToPlot,'color',[ 0.5 0.8 0.9],'lineWidth',1.5)
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
    polarhistogram(deg2rad(stim_offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
    set(gca,'ThetaZeroLocation','top');
    title({'Offset distribution','High contrast'});
    Ax = gca;
    Ax.RGrid = 'off';
    Ax.RTickLabel = [];
    
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
    plot(x_out_offset,offsetToPlot,'color',[ 0.5 0.8 0.9],'lineWidth',1.5)
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
    polarhistogram(deg2rad(offset{session}),20,'FaceColor',[ 0.5 0.8 0.9])
    set(gca,'ThetaZeroLocation','top');
    title({'Offset distribution','High contrast'});
    Ax = gca;
    Ax.RGrid = 'off';
    Ax.RTickLabel = [];
    
    %save figures
    saveas(gcf,[path,'\analysis\plots\openLoopCombinedHeatmap_session',num2str(session),'.png']);
end

%% Model BM and BW as a function of total movement, and stim velocity

close all;
%create table with relevant variables
allBumpMag = [];
allBumpWidth = [];
allTotalMvt = [];
allZscoredMvt = [];
allTime = [];
stimVel = [];
for session = 1:length(open_loop_sessions)
   allBumpMag = [allBumpMag,BumpMagnitude{1,session}];
   allBumpWidth = [allBumpWidth,BumpWidth{1,session}]; 
   allTotalMvt = [allTotalMvt,fillmissing(data{1,session}.continuous_data.total_mvt_ds,'linear')];
   allZscoredMvt = [allZscoredMvt,zscore(fillmissing(data{1,session}.continuous_data.total_mvt_ds,'linear'))];
   allTime = [allTime,data{1,session}.continuous_data.time'];
   if (data{1,session}.continuous_data.run_obj.function_number == 52 | data{1,session}.continuous_data.run_obj.function_number == 53)
       stimVel = [stimVel,repelem(30,length(data{1,session}.continuous_data.time))];
   elseif (data{1,session}.continuous_data.run_obj.function_number == 54 | data{1,session}.continuous_data.run_obj.function_number == 55)
       stimVel = [stimVel,repelem(60,length(data{1,session}.continuous_data.time))];
   elseif (data{1,session}.continuous_data.run_obj.function_number == 206 | data{1,session}.continuous_data.run_obj.function_number == 207)
       stimVel = [stimVel,repelem(90,length(data{1,session}.continuous_data.time))];
   elseif (data{1,session}.continuous_data.run_obj.function_number == 208 | data{1,session}.continuous_data.run_obj.function_number == 209)
       stimVel = [stimVel,repelem(120,length(data{1,session}.continuous_data.time))];
   elseif (data{1,session}.continuous_data.run_obj.function_number == 210 | data{1,session}.continuous_data.run_obj.function_number == 211)
       stimVel = [stimVel,repelem(150,length(data{1,session}.continuous_data.time))];
   elseif (data{1,session}.continuous_data.run_obj.function_number == 212 | data{1,session}.continuous_data.run_obj.function_number == 213)
       stimVel = [stimVel,repelem(180,length(data{1,session}.continuous_data.time))];
   elseif (data{1,session}.continuous_data.run_obj.function_number == 214 | data{1,session}.continuous_data.run_obj.function_number == 215)
       stimVel = [stimVel,repelem(210,length(data{1,session}.continuous_data.time))];
   end
end
bump_mag_data = table(stimVel',allTime',allTotalMvt',allZscoredMvt',allBumpMag','VariableNames',{'stimVel','Time','TotalMovement','ZscoredMvt','BumpMagnitude'});
bump_width_data = table(stimVel',allTime',allTotalMvt',allZscoredMvt',allBumpWidth','VariableNames',{'stimVel','Time','TotalMovement','ZscoredMvt','BumpWidth'});

mdl_BM = fitlm(bump_mag_data,'BumpMagnitude ~ stimVel+ZscoredMvt')
mdl_BW = fitlm(bump_width_data,'BumpWidth ~ stimVel+ZscoredMvt')


%% Save relevant data

save([path,'\analysis\open_loop_data.mat'],'bump_mag_data','bump_width_data','summarydata','mean_data');

clear all; close all;