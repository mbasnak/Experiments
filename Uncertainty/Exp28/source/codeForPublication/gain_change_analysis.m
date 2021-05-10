%analysis for the closed-loop bouts of the change in contrast experiment


%% Load data

clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');
load([path,file])

%determine the sid we're working with, to save the plots to specific folder
%later
sid = file(10:end-10);


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

%% Determine changes in gain

%Get hdf5 file from this session
file_names = dir(fullfile([path(1:end-9),'ball\'],'*hdf5'));
for file = 1:length(file_names)
    if contains(file_names(file).name,[sid,'_'])
        hdf5_file_to_read = fullfile(file_names(file).folder,file_names(file).name);
    end
end
%Read the gain yaw variable
gain_yaw = double(h5read(hdf5_file_to_read,'/gain_yaw'));
%Downsample to match data length
gain_yaw_ds = resample(gain_yaw,length(data.time),length(gain_yaw));
gain_yaw_ds(gain_yaw_ds<0) = -1;
gain_yaw_ds(gain_yaw_ds>0) = 1;

%Determine frames where the gain changes
gain_changes = find(abs(diff(gain_yaw_ds))>0.5);
gain_changes = gain_changes(1:2);

%% Set block limits

blockLimits{1} = [1,gain_changes(1)-1];
blockLimits{2} = [gain_changes(1),gain_changes(2)];
blockLimits{3} = [gain_changes(2)+1,length(data.time)];

%Define gain per block
gain_per_block = [1,-1,1];

%Set color palette based on gain
for block = 1:3
    if gain_per_block(block) == 1
        color_gradient{block} = [.4 0.1 0.1];
    else
        color_gradient{block} = [.1 0.4 0.2];
    end
end

%% Compute heading offset variability

heading_offset = wrapTo180(data.heading_offset);
[x_out_heading_offset, heading_offset_to_plot] = removeWrappedLines(data.time,heading_offset);

%I will analyze offset variability by computing a rolling window of circular
%standard deviation of offset and taking the inverse
fcn = @(x) circ_std(x);
%Compute the offset variability over different window sizes, from 1 s to
%10 s
window_sizes = [10,30,50,100,500,1000];
for window = 1:length(window_sizes)
    heading_offset_variability(:,window) = matlab.tall.movingWindow(fcn,window_sizes(window),deg2rad(heading_offset));
end
smoothed_heading_offset_variability = smooth(heading_offset_variability(:,3),150,'rloess');

figure('Position',[100 100 1600 400]),
subplot(3,1,1)
plot(x_out_heading_offset,heading_offset_to_plot,'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0, 0.5, 0]);
end
ylim([-180 180]);
title('Heading offset');
set(gca,'xticklabel',{[]});

subplot(3,1,2)
plot(data.time,heading_offset_variability(:,3),'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
end
title('Heading offset variability');
set(gca,'xticklabel',{[]});

subplot(3,1,3)
plot(data.time,smoothed_heading_offset_variability,'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
end
title('Smoothed heading offset variability');
xlabel('Time (sec)');

%save figure
saveas(gcf,[path,'plots\HeadingOffsetVariability.png']);

%% Repeat for bar offset

offset = wrapTo180(data.bar_offset);
[x_out_offset, offset_to_plot] = removeWrappedLines(data.time,offset);

for window = 1:length(window_sizes)
    bar_offset_variability(:,window) = matlab.tall.movingWindow(fcn,window_sizes(window),deg2rad(offset));
end
smoothed_bar_offset_variability = smooth(bar_offset_variability(:,3),150,'rloess');


figure('Position',[100 100 1600 400]),
subplot(3,1,1)
plot(x_out_offset,offset_to_plot,'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0, 0.5, 0]);
end
ylim([-180 180]);
title('Bar offset');
set(gca,'xticklabel',{[]});

subplot(3,1,2)
plot(data.time,bar_offset_variability(:,3),'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
end
title('Bar offset variability');
set(gca,'xticklabel',{[]});

subplot(3,1,3)
plot(data.time,smoothed_bar_offset_variability,'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0, 0.5, 0]);
end
title('Smoothed bar offset variability');
xlabel('Time (sec)');

%save figure
saveas(gcf,[path,'plots\BarOffsetVariability.png']);

%% Determine the type of fly based on the ratio between bar and offset variability

ratio_means = mean(bar_offset_variability(:,3))/mean(heading_offset_variability(:,3));
ratio_medians = median(bar_offset_variability(:,3))/median(heading_offset_variability(:,3));

all_ratios = [ratio_means,ratio_medians];

figure,
plot(all_ratios,'-ko','linewidth',2)
hold on
yline(1,'r','linewidth',2)
text(0.2,1.2,'Type 1');
text(0.2,0.8,'Type 2');
ylim([0 2]);
xlim([0 3]);
xticks([1 2])
xticklabels({'ratio of means','ratio of medians'})
ylabel('Bar offset variability / heading offset variability');
%Save figure
saveas(gcf,[path,'plots\FlyClassification.png']);

%Classify fly
if mean(all_ratios) > 1
    type_of_fly = 1; %fly that learns the new mapping (and where heading offset is more stable)
else
    type_of_fly = 2; %fly that ignores proprioceptive cues (and where bar offset is more stable)
end
    
%% Obtain bump parameters

%Bump magnitude
allBumpMag = [];
for block = 1:length(blockLimits)
   bump_mag{block} = data.bump_magnitude(:,blockLimits{block}(1):blockLimits{block}(2)); 
   allBumpMag = [allBumpMag,bump_mag{block}];
end

%Bump width at half max
half_max_width = compute_bump_width(data.mean_dff_EB);    
allHalfWidth = [];
for block = 1:length(blockLimits)
   width_half_max{block} = half_max_width(blockLimits{block}(1):blockLimits{block}(2)); 
   allHalfWidth = [allHalfWidth,width_half_max{block}];
end

%% Heatmap plot including bump magnitude

if type_of_fly == 1
    
    % Plot the heatmap of EPG activity
    figure('Position',[100 100 1400 800]),
    subplot(5,1,1)
    %I'm flipping the dff matrix for it to make sense along with the fly's
    %heading
    imagesc(flip(data.dff_matrix))
    colormap(flipud(gray))
    hold on
    %add the changes in stim
    for change = 1:length(gain_changes)
        line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', [0, 0.5, 0]);
    end
    yticks(1:2:16);
    yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
    ylabel('PB glomerulus','fontweight','bold','fontsize',10);
    title('EPG activity in the PB','fontweight','bold','fontsize',12);
    set(gca,'XTickLabel',[]);
    legend('Change in stimulus');
    
    % Plot the bar position and the EPG phase
    subplot(5,1,2)
    %Get heading position to plot
    heading = wrapTo180(data.heading_deg);
    [x_out_heading, heading_to_plot] = removeWrappedLines(data.time,heading);
    plot(x_out_heading,heading_to_plot,'color',[0.6 0.3 0.8],'LineWidth',1.5)
    hold on
    %I'm now going to negate the phase, since I'm plotting heading instead of
    %bar position, to the bump moves in the other direction
    phase = wrapTo180(rad2deg(-data.dff_pva));
    [x_out_phase, phase_to_plot] = removeWrappedLines(data.time,phase);
    plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
    %add the changes in stim
    for change = 1:length(gain_changes)
        line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color',  [0, 0.5, 0]);
    end
    legend('Fly heading', 'EPG phase');
    title('EPG phase and fly heading','fontweight','bold','fontsize',12);
    ylim([-180, 180]);
    xlim([0,data.time(end)]);
    ylabel('Deg','fontweight','bold','fontsize',10);
    set(gca,'XTickLabel',[]);
    
    % Plot the offset
    subplot(5,1,3)
    plot(x_out_heading_offset,heading_offset_to_plot,'LineWidth',1.5,'color','k')
    %Add the changes in stim
    for change = 1:length(gain_changes)
        line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color',  [0, 0.5, 0]);
    end
    ylim([-180 180]);
    ylabel('Deg','fontweight','bold','fontsize',10);
    set(gca,'XTickLabel',[]);
    title('Heading offset','fontweight','bold','fontsize',12);
    
    %Plot the bump magnitude
    subplot(5,1,4)
    for block = 1:length(blockLimits)
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(bump_mag{block}),'color',color_gradient{block},'linewidth',1.5)
        hold on
    end
    title('Bump magnitude','fontweight','bold','fontsize',12);
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'},'fontweight','bold','fontsize',10);
    set(gca,'xticklabel',{[]});
    xlim([0 data.time(end)]);
    
    subplot(5,1,5)
    for block = 1:length(blockLimits)
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(width_half_max{block}),'color',color_gradient{block},'linewidth',1.5)
        hold on
    end
    title('Bump half width','fontweight','bold','fontsize',12);
    ylabel('Bump half width','fontweight','bold','fontsize',10);
    xlabel('Time (sec)','fontweight','bold','fontsize',10);
    xlim([0 data.time(end)]);
    
    %save figure
    saveas(gcf,[path,'plots\HeadingOffsetWithBumpParameters.png']);
    
    
else
    
    % Plot the heatmap of EPG activity
    figure('Position',[100 100 1400 800]),
    subplot(5,1,1)
    imagesc(data.dff_matrix)
    colormap(flipud(gray))
    hold on
    %Add the changes in stim
    for change = 1:length(gain_changes)
        line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', [0,0.5,0]);
    end
    yticks(1:2:16);
    yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
    ylabel('PB glomerulus','fontweight','bold','fontsize',10);
    title('EPG activity in the PB','fontweight','bold','fontsize',12);
    set(gca,'XTickLabel',[]);
    legend('Change in stimulus');
    
    % Plot the bar position and the EPG phase
    subplot(5,1,2)
    heading = wrapTo180(data.heading_deg);
    [x_out_heading, heading_to_plot] = removeWrappedLines(data.time,heading);
    bar_position = wrapTo180(data.panel_angle);
    [x_out_bar, bar_pos_to_plot] = removeWrappedLines(data.time,bar_position);
    plot(x_out_bar,bar_pos_to_plot,'color',[0.2 0.6 0.7],'LineWidth',1.5)
    hold on
    %Get EPG phase to plot
    phase = wrapTo180(rad2deg(data.dff_pva));
    [x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
    plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
    %Add the changes in stim
    for change = 1:length(gain_changes)
        line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0,0.5,0]);
    end
    legend('Bar position', 'EPG phase');
    title('EPG phase and bar position','fontweight','bold','fontsize',12);
    ylim([-180, 180]);
    xlim([0,data.time(end)]);
    ylabel('Deg','fontweight','bold','fontsize',10);
    set(gca,'XTickLabel',[]);
    
    % Plot the offset
    subplot(5,1,3)
    offset = wrapTo180(data.bar_offset);
    plot(x_out_offset,offset_to_plot,'LineWidth',1.5,'color','k')
    %Add the changes in stim
    for change = 1:length(gain_changes)
        line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0,0.5,0]);
    end
    ylim([-180 180]);
    title('Bar offset','fontweight','bold','fontsize',12);
    ylabel('Deg','fontweight','bold','fontsize',10); 
    
    %Plot the bump magnitude
    subplot(5,1,4)
    for block = 1:length(blockLimits)
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(bump_mag{block}),'color',color_gradient{block},'linewidth',1.5)
        hold on
    end
    title('Bump magnitude','fontweight','bold','fontsize',12);
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'},'fontweight','bold','fontsize',10);
    set(gca,'xticklabel',{[]});
    xlim([0 data.time(end)]);
    
    subplot(5,1,5)
    for block = 1:length(blockLimits)
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),smooth(width_half_max{block}),'color',color_gradient{block},'linewidth',1.5)
        hold on
    end
    title('Bump half width','fontweight','bold','fontsize',12);
    ylabel('Bump half width','fontweight','bold','fontsize',10);
    xlabel('Time (sec)','fontweight','bold','fontsize',10);
    xlim([0 data.time(end)]);
    
    %Save figure
    saveas(gcf,[path,'plots\BarOffsetWithBumpParameters.png']);   
end


%% Heading stability

heading_variability = matlab.tall.movingWindow(fcn,50,deg2rad(heading));
smoothed_heading_variability = smooth(heading_variability,150,'rloess');


figure('Position',[100 100 1600 400]),
subplot(3,1,1)
plot(x_out_heading,heading_to_plot,'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', [0,0.5,0]);
end
ylim([-180 180]);
title('Heading');

subplot(3,1,2)
plot(data.time,smooth(heading_variability),'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0,0.5,0]);
end
title('Heading variability');

subplot(3,1,3)
plot(data.time,smoothed_heading_variability,'k')
hold on
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [0 2], 'LineWidth', 2, 'color', [0,0.5,0]);
end
title('Smoothed heading variability');

%save figure
saveas(gcf,[path,'plots\HeadingVariability.png']);

%% Correlate BM and HW with heading or bar offset variability during the inverted gain portion, depending on the type of fly

%Set a mvt threshold
mvt_thresh = 50;

BumpMagIG = bump_mag{2};
HalfWidthIG = width_half_max{2};
heading_offset_variabilityIG = heading_offset_variability(gain_changes(1):gain_changes(2),:);
bar_offset_variabilityIG = bar_offset_variability(gain_changes(1):gain_changes(2),:);
total_mvtIG = data.total_mvt_ds(gain_changes(1):gain_changes(2));
yaw_speedIG = abs(data.vel_yaw_ds(gain_changes(1):gain_changes(2)));


if type_of_fly == 1
    
    figure,
    
    nbins = 10;
    
    %Define bins
    max_bin = prctile(heading_offset_variabilityIG,95,'all');
    min_bin = prctile(heading_offset_variabilityIG,5,'all');
    binWidth = (max_bin-min_bin)/nbins;
    Bins = [min_bin:binWidth:max_bin];
    
    %Create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);      
    
    %Getting binned means for the different windowSizes
    for bin = 1:length(Bins)-1
        for window = 1:length(window_sizes)
            meanBin(bin,window) = nanmean(BumpMagIG((heading_offset_variabilityIG(:,window) > Bins(bin)) & (heading_offset_variabilityIG(:,window) < Bins(bin+1)) & (total_mvtIG' > mvt_thresh)));
        end
    end
       
    %Plot
    subplot(1,2,1)
    plot(mvtAxes,meanBin)
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
    ylim([0 1]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    legend({'10 frames','30','50','100','500','1000'},'location','best');
    title('Bump magnitude');
    
    %Getting binned means
    for bin = 1:length(Bins)-1
        for window = 1:length(window_sizes)
            meanBin(bin,window) = nanmean(HalfWidthIG((heading_offset_variabilityIG(:,window) > Bins(bin)) & (total_mvtIG' > mvt_thresh) &(heading_offset_variabilityIG(:,window) < Bins(bin+1))));
        end
    end
    
    %Plot
    subplot(1,2,2)
    plot(mvtAxes,meanBin)
    ylabel('Bump half width'); xlabel('Heading offset variability');
    ylim([0 4]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    legend({'10 frames','30','50','100','500','1000'},'location','best');
    title('Bump half width');

    
    %save figure
    saveas(gcf,[path,'plots\headingOffsetVsMedianBumpParametersIG.png']);
    
else
    
    figure,
    
    nbins = 10;
    
    %Define bins
    max_bin = prctile(heading_offset_variabilityIG,95,'all');
    min_bin = prctile(heading_offset_variabilityIG,5,'all');
    binWidth = (max_bin-min_bin)/nbins;
    Bins = [min_bin:binWidth:max_bin];
    
    %Create axes for plot
    mvtAxes = Bins - binWidth;
    mvtAxes = mvtAxes(2:end);      
    
    %Getting binned means for the different windowSizes
    for bin = 1:length(Bins)-1
        for window = 1:length(window_sizes)
            meanBin(bin,window) = nanmean(BumpMagIG((bar_offset_variabilityIG(:,window) > Bins(bin)) & (bar_offset_variabilityIG(:,window) < Bins(bin+1)) & (total_mvtIG' > mvt_thresh)));
        end
    end
       
    %Plot
    subplot(1,2,1)
    plot(mvtAxes,meanBin)
    ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Bar offset variability');
    ylim([0 1]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    legend({'10 frames','30','50','100','500','1000'},'location','best');
    title('Bump magnitude');
    
    %Getting binned means
    for bin = 1:length(Bins)-1
        for window = 1:length(window_sizes)
            meanBin(bin,window) = nanmean(HalfWidthIG((bar_offset_variabilityIG(:,window) > Bins(bin)) & (total_mvtIG' > mvt_thresh) &(bar_offset_variabilityIG(:,window) < Bins(bin+1))));
        end
    end
    
    %Plot
    subplot(1,2,2)
    plot(mvtAxes,meanBin)
    ylabel('Bump half width'); xlabel('Bar offset variability');
    ylim([0 4]);
    xlim([mvtAxes(1) mvtAxes(end)]);
    legend({'10 frames','30','50','100','500','1000'},'location','best');
    title('Bump half width');
    
    %save figure
    saveas(gcf,[path,'plots\barOffsetVsMedianBumpParametersIG.png']);
end

%% Models of bump parameters

%We will model the two bump parameters as a function of total movement and
%offset variability 

%Create table with the model's variables
for window = 1:length(window_sizes)
    modelTable{window} = table(bar_offset_variabilityIG(:,window),heading_offset_variabilityIG(:,window),total_mvtIG',yaw_speedIG',BumpMagIG',HalfWidthIG','VariableNames',{'BarOffsetVariability','HeadingOffsetVariability','TotalMovement','YawSpeed','BumpMagnitude','BumpWidth'});
    if type_of_fly == 1
        mdl_BM{window} = fitlm(modelTable{window},'BumpMagnitude~HeadingOffsetVariability+TotalMovement');
        mdl_HW{window} = fitlm(modelTable{window},'BumpWidth~HeadingOffsetVariability+TotalMovement');
        %Model Rsquared
        Rsquared_BM(window) = mdl_BM{window}.Rsquared.Adjusted;
        Rsquared_HW(window) = mdl_HW{window}.Rsquared.Adjusted;
    else
        mdl_BM{window} = fitlm(modelTable{window},'BumpMagnitude~BarOffsetVariability+TotalMovement');
        mdl_HW{window} = fitlm(modelTable{window},'BumpWidth~BarOffsetVariability+TotalMovement');
        %Model Rsquared
        Rsquared_BM(window) = mdl_BM{window}.Rsquared.Adjusted;
        Rsquared_HW(window) = mdl_HW{window}.Rsquared.Adjusted;        
    end
end

%Plot model fit with the different time windows
figure,
subplot(1,2,1)
plot(Rsquared_BM,'-o')
title('Bump magnitude');
ylabel('Rsquared');
xlabel('window #');

subplot(1,2,2)
plot(Rsquared_HW,'-o')
title('Bump width');
ylabel('Rsquared');
xlabel('window #');

saveas(gcf,[path,'plots\model_fit.png']);

%% Correlate BM and HW with offset variability during the normal gain portion

BumpMagNG = allBumpMag([1:gain_changes(1),gain_changes(2):end]);
HalfWidthNG = allHalfWidth([1:gain_changes(1),gain_changes(2):end]);
heading_offset_variabilityNG = heading_offset_variability([1:gain_changes(1),gain_changes(2):end],:);
total_mvtNG = data.total_mvt_ds([1:gain_changes(1),gain_changes(2):end]);
heading_variabilityNG = heading_variability([1:gain_changes(1),gain_changes(2):end]);

figure,

nbins = 10;

%Define bins
max_bin = prctile(heading_offset_variabilityNG,95,'all');
min_bin = prctile(heading_offset_variabilityNG,5,'all');
binWidth = (max_bin-min_bin)/nbins;
Bins = [min_bin:binWidth:max_bin];

%Create axes for plot
mvtAxes = Bins - binWidth;
mvtAxes = mvtAxes(2:end);

%Getting binned means for the different windowSizes
for bin = 1:length(Bins)-1
    for window = 1:length(window_sizes)
        meanBin(bin,window) = nanmean(BumpMagNG((heading_offset_variabilityNG(:,window) > Bins(bin)) & (heading_offset_variabilityNG(:,window) < Bins(bin+1)) & (total_mvtNG' > mvt_thresh)));
    end
end

%Plot
subplot(1,2,1)
plot(mvtAxes,meanBin)
ylabel({'Bump magnitude';'(amplitude of Fourier component)'}); xlabel('Heading offset variability');
ylim([0 max(max(meanBin))+0.5]);
xlim([mvtAxes(1) mvtAxes(end)]);
legend({'10 frames','30','50','100','500','1000'},'location','best');
title('Bump magnitude');

%Getting binned means
for bin = 1:length(Bins)-1
    for window = 1:length(window_sizes)
        meanBin(bin,window) = nanmean(HalfWidthNG((heading_offset_variabilityNG(:,window) > Bins(bin)) & (total_mvtNG' > mvt_thresh) &(heading_offset_variabilityNG(:,window) < Bins(bin+1))));
    end
end

%Plot
subplot(1,2,2)
plot(mvtAxes,meanBin)
ylabel('Bump half width'); xlabel('Heading offset variability');
ylim([0 max(max(meanBin))+0.5]);
xlim([mvtAxes(1) mvtAxes(end)]);
legend({'10 frames','30','50','100','500','1000'},'location','best');
title('Bump half width');

%Save figure
saveas(gcf,[path,'plots\OffsetVsMedianBumpParametersNG.png']);

%% Run models for the normal gain portion

%Create table with the model's variables
for window = 1:length(window_sizes)
    modelTableNG{window} = table(heading_offset_variabilityNG(:,window),total_mvtNG',BumpMagNG',HalfWidthNG',heading_variabilityNG,'VariableNames',{'HeadingOffsetVariability','TotalMovement','BumpMagnitude','BumpWidth','HeadingVariability'});
    mdl_BM_NG{window} = fitlm(modelTableNG{window},'BumpMagnitude~HeadingOffsetVariability+TotalMovement');
    mdl_HW_NG{window} = fitlm(modelTableNG{window},'BumpWidth~HeadingOffsetVariability+TotalMovement');
    %Model Rsquared
    Rsquared_BM_NG(window) = mdl_BM_NG{window}.Rsquared.Adjusted;
    Rsquared_HW_NG(window) = mdl_HW_NG{window}.Rsquared.Adjusted;
end

%Plot model fit with the different time windows
figure,
subplot(1,2,1)
plot(Rsquared_BM_NG,'-o')
title('Bump magnitude');
ylabel('Rsquared');
xlabel('window #');

subplot(1,2,2)
plot(Rsquared_HW_NG,'-o')
title('Bump width');
ylabel('Rsquared');
xlabel('window #');

saveas(gcf,[path,'plots\model_fit_NG.png']);

%% Divide the inverted gain portion into 4 parts and compare bump parameters

quartile_changes = [gain_changes(1):floor(gain_changes(2)-gain_changes(1))/4:gain_changes(2)];
warning('off','all')

for quartile = 1:4
   BM = data.bump_magnitude(quartile_changes(quartile):quartile_changes(quartile+1));
   mvt_q = data.total_mvt_ds(quartile_changes(quartile):quartile_changes(quartile+1));
   BMq(quartile) = nanmean(BM(mvt_q > mvt_thresh));
   
   BW = half_max_width(quartile_changes(quartile):quartile_changes(quartile+1));
   BWq(quartile) = nanmean(BW(mvt_q > mvt_thresh));   
end

figure,
subplot(1,6,1)
BMi = data.bump_magnitude(1:gain_changes(1));
mvt_i = data.total_mvt_ds(1:gain_changes(1));
plot(nanmean(BMi(mvt_i > mvt_thresh)),'o')
ylabel('Bump magnitude');
title('N gain');
set(gca,'xticklabel',{[]});
ylim([0 1]);

subplot(1,6,[2:5])
plot(BMq,'o')
set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Inverted gain');
xlim([0 5]);
ylim([0 1]);

subplot(1,6,6)
BMe = data.bump_magnitude(gain_changes(2):end);
mvt_e = data.total_mvt_ds(gain_changes(2):end);
plot(nanmean(BMe(mvt_e > mvt_thresh)),'o')
title('N gain');
set(gca,'yticklabel',{[]});
set(gca,'xticklabel',{[]});
ylim([0 1]);

saveas(gcf,[path,'plots\BM_quartiles.png']);


%Repeat for bump width
figure,
subplot(1,6,1)
BWi = half_max_width(1:gain_changes(1));
plot(nanmean(BWi(mvt_i > mvt_thresh)),'o')
ylabel('Bump width');
title('N gain');
set(gca,'xticklabel',{[]});
ylim([0 4]);

subplot(1,6,[2:5])
plot(BWq,'o')
set(gca,'yticklabel',{[]});
xticks([1:4])
xlabel('Quartile');
title('Inverted gain');
xlim([0 5]);
ylim([0 4]);

subplot(1,6,6)
BWe = half_max_width(gain_changes(2):end);
plot(nanmean(BWe(mvt_e > mvt_thresh)),'o')
title('N gain');
set(gca,'yticklabel',{[]});
set(gca,'xticklabel',{[]});
ylim([0 4]);

saveas(gcf,[path,'plots\BW_quartiles.png']);

%% Save relevant data

save([path,'\gain_change_data.mat'],'modelTable','modelTableNG','type_of_fly');

close all; clear all;
