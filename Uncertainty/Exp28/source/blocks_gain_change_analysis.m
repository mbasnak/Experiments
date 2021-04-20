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

%determine initial gain
%get hdf5 files in directory
file_names = dir(fullfile([path(1:end-9),'ball\'],'*hdf5'));
for file = 1:length(file_names)
    if contains(file_names(file).name,[sid,'_'])
        hdf5_file_to_read = fullfile(file_names(file).folder,file_names(file).name);
    end
end

gain_yaw = double(h5read(hdf5_file_to_read,'/gain_yaw'));
%downsample to match data length
gain_yaw_ds = resample(gain_yaw,length(data.time),length(gain_yaw));
gain_yaw_ds(gain_yaw_ds<0) = -1;
gain_yaw_ds(gain_yaw_ds>0) = 1;

%determine gain changes
gain_changes = find(abs(diff(gain_yaw_ds))>0.5);
gain_changes = gain_changes(1:5);



%% Set block limits

blockLimits{1} = [1,gain_changes(1)-1];
for block = 2:5
    blockLimits{block} = [gain_changes(block-1),gain_changes(block)-1];
end
blockLimits{6} = [gain_changes(5),length(data.time)];

%determine gain per block
if gain_yaw(1) == 1
    gain_per_block(1:2:5) = 1;
    gain_per_block(2:2:6) = -1;
else
    gain_per_block(1:2:5) = -1;
    gain_per_block(2:2:6) = 1;
end

%set color palette based on gain
for block = 1:6
    if gain_per_block(block) == 1
        color_gradient{block} = [.8 0.2 0.2];
    else
        color_gradient{block} = [.3 0.8 0.2];
    end
end

%% Plot heatmap and bar offset

%Plot
% Plot the heatmap of EPG activity
figure('Position',[100 100 1400 800]),
subplot(4,6,[1 6])
%I'm flipping the dff matrix for it to make sense along with the fly's
%heading
imagesc(flip(data.dff_matrix))
colormap(gray)
hold on
%add the changes in stim
for change = 1:length(gain_changes)
    line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', 'r');
end
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
ylabel('PB glomerulus');
title('EPG activity in the PB');
set(gca,'XTickLabel',[]);
legend('Change in stimulus');

% Plot the heading and the EPG phase
subplot(4,6,[7 12])
%Get heading to plot
heading = wrapTo180(data.heading_deg);
change_heading = abs([0;diff(smooth(heading))]);
heading_to_plot = smooth(heading);
heading_to_plot(change_heading>40==1) = NaN;
plot(heading_to_plot,'LineWidth',1.5)
plot(data.time,heading_to_plot,'color',[0.6 0.2 0.4],'LineWidth',1.5)
hold on
%Get EPG phase to plot
%I'm now going to negate the phase, since I'm plotting heading instead of
%bar position, to the bump moves in the other direction
phase = wrapTo180(rad2deg(-data.dff_pva));
change_phase = abs([0;diff(smooth(phase))]);
phase_to_plot = smooth(phase);
phase_to_plot(change_phase>40==1) = NaN;
plot(data.time,phase_to_plot,'color',[0.2 0.6 0.4],'LineWidth',1.5)
%add the changes in stim
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', 'r');
end
legend('Bar position', 'EPG phase');
ylim([-180, 180]);
xlim([0,data.time(end)]);
ylabel('Deg');
set(gca,'XTickLabel',[]);

% Plot the offset
subplot(4,6,[13 18])
%Get offset to plot
offset = wrapTo180(data.offset);
change_offset = abs([0;diff(smooth(offset))]);
offset_to_plot = smooth(offset);
offset_to_plot(change_offset>30==1) = NaN;
plot(data.time,offset_to_plot,'LineWidth',1.5,'color','k')
%add the changes in stim
for change = 1:length(gain_changes)
    line([data.time(gain_changes(change)) data.time(gain_changes(change))], [-180 180], 'LineWidth', 2, 'color', 'r');
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
legend('Offset');

% Polar histograms of offset
% Color histograms acording to the intensity level of the bar, using the
% vector Intensities
for block = 1:6
    subplot(4,6,18+block)
    polarhistogram(deg2rad(offset(blockLimits{block}(1):blockLimits{block}(2))),15,'FaceColor',color_gradient{block})
    set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
end

%save figure
saveas(gcf,[path,'plots\gainChangeHeatmapAndOffset.png']);

%% Offset distribution

figure('Position',[200 200 1000 800]),
subplot(2,1,1)
for block = 1:length(blockLimits)
   [~, offset_var(block)] = circ_std(deg2rad(offset(blockLimits{block}(1):blockLimits{block}(2))),[],[],1); 
   plot(block,offset_var(block),'ko','MarkerFaceColor',color_gradient{block},'MarkerSize',8)
   hold on
end
xlim([0 7]);
xticks(1:6);
title('Offset variation per block');
ylabel({'Circular standard deviation','of the offset'});
%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
if gain_yaw(1) == 1
    legend(h,'Normal gain','Inverted gain');
else
    legend(h,'Inverted gain','Normal gain');
end

subplot(2,1,2)
%Create table with contrast level and offset variation
summary_data = table(gain_per_block',offset_var','VariableNames',{'gain','offset_var'});
%Get mean offset var by contrast
mean_data = varfun(@mean,summary_data,'InputVariables','offset_var',...
       'GroupingVariables',{'gain'});

%plot
plot(mean_data.mean_offset_var,'-ko','MarkerSize',8)
xlim([0 3]);
ylim([min(mean_data.mean_offset_var)-0.5 max(mean_data.mean_offset_var)+0.5]); 
xticks(1:3);
if gain_yaw(1) == 1
    xticklabels({'Normal gain','Inverted gain'});
else
    xticklabels({'Inverted gain','Normal gain'});
end
title('Mean offset variation per gain');
ylabel({'Circular standard deviation','of the offset'});

%save figure
saveas(gcf,[path,'plots\GainChangeOffsetVariation.png']);

%% Plot heatmap and heading offset



%% Get offset value from last bout of gain = 1

%Find limits of last high contrast bout
normal_gain_bouts = find(gain_per_block==1);
last_normal_gain = normal_gain_bouts(end);

%This will be the 'visual' offset, so we need to compute it as the circular
%distance between the stimulus and the phase
visual_offset = circ_dist(data.dff_pva,deg2rad(wrapTo180(data.panel_angle)));
mean_reference_offset = rad2deg(circ_mean(visual_offset(blockLimits{last_normal_gain}(1):blockLimits{last_normal_gain}(2)),[],1));

%% Calculate and plot bump magnitude in time

%compute bump magnitude as max-min
for block = 1:length(blockLimits)
   bump_mag{block} = max(data.mean_dff_EB(:,blockLimits{block}(1):blockLimits{block}(2)))-min(data.mean_dff_EB(:,blockLimits{block}(1):blockLimits{block}(2))); 
end

figure('Position',[200 200 1600 600]),
%plot EPG activity
subplot(2,1,1)
imagesc(flip(data.dff_matrix))
colormap(gray)
hold on
%add the changes in stim
for change = 1:length(gain_changes)
   line([gain_changes(change) gain_changes(change)], [0 17], 'LineWidth', 2, 'color', 'r'); 
end
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
ylabel('PB glomerulus');
title('EPG activity in the PB');
set(gca,'XTickLabel',[]);
legend('Change in stimulus');

%plot bump magnitude
subplot(2,1,2)
for block = 1:length(blockLimits)
    plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{block})
    hold on
end
title('Bump magnitude');
ylabel('Bump magnitude (max-min)');
xlabel('Time (sec)');
xlim([0 data.time(end)]);

%save figure
saveas(gcf,[path,'plots\gainChangeBMinTime.png']);

%% Compute and plot mean bump magnitude per block

figure('Position',[200 200 1000 800]),
subplot(2,1,1)
for block = 1:length(blockLimits)
   mean_bump_mag(block) = mean(bump_mag{block}); 
   plot(block,mean_bump_mag(block),'ko','MarkerFaceColor',color_gradient{block},'MarkerSize',8)
   hold on
end
xlim([0 7]);
xticks(1:6);
title('Mean bump magnitude per block');
ylabel('Bump magnitude (max-min)');
%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
if gain_yaw(1) == 1
    legend(h,'Normal gain','Inverted gain');
else
    legend(h,'Inverted gain','Normal gain');
end

subplot(2,1,2)
%add bump mag data to table
if size(summary_data,2) == 2
    summary_data = addvars(summary_data, mean_bump_mag','NewVariableNames','mean_bump_mag');
end
%Get mean bump mag by contrast
mean_bump_data = varfun(@mean,summary_data,'InputVariables','mean_bump_mag',...
       'GroupingVariables',{'gain'});
%Plot
plot(mean_bump_data.mean_mean_bump_mag,'-ko','MarkerSize',8)
xlim([0 3]);
ylim([min(mean_bump_data.mean_mean_bump_mag)-0.5 max(mean_bump_data.mean_mean_bump_mag)+0.5]); 
xticks(1:3);
if gain_yaw(1) == 1
    xticklabels({'Normal gain','Inverted gain'});
else
    xticklabels({'Inverted gain','Normal gain'});
end
title('Mean bump magnitude per gain');
ylabel('Bump magnitude (max-min)');

%save figure
saveas(gcf,[path,'plots\gainChangeMeanBM.png']);

%% Bump width at half max

%Plot example data point

for timepoint = [1,400,1000,5000]
    
    figure,
    %add bump mag
    ex_bump_mag = max(data.mean_dff_EB(:,timepoint)) - min(data.mean_dff_EB(:,timepoint));
    %linearly interpolate to have 1000 datapoints instead of 8
    interp_ex_data = interp1([1:8],data.mean_dff_EB(:,timepoint),[1:7/1000:8]);
    plot(interp_ex_data)
    hold on
    xlabel('EB tile');
    ylabel('DFF');
    %Find the half max point
    half_max = (max(data.mean_dff_EB(:,timepoint))-min(data.mean_dff_EB(:,timepoint)))/2 + min(data.mean_dff_EB(:,timepoint));
    [ex_bump_mag_interp I_interp] = max(interp_ex_data);
    plot([I_interp I_interp],[ex_bump_mag_interp ex_bump_mag_interp],'ko')
    %Find in each half the index closest to the half max
    diff_data = abs(interp_ex_data-half_max);
    [sortedVals,indexes] = sort(diff_data);
    %remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
    %index(1)
    diff_indexes = abs(indexes-indexes(1));
    indexes(diff_indexes<175 & diff_indexes>0)=NaN;
    indexes = indexes(~isnan(indexes));
    two_indexes = [indexes(1), indexes(2)];
    I1 = min(two_indexes);
    I2 = max(two_indexes);
    plot(I1,interp_ex_data(I1),'ro')
    plot(I2,interp_ex_data(I2),'ro')
    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
        line([1 I1],[half_max half_max],'color','r','LineWidth',2);
        line([I2 1000],[half_max half_max],'color','r','LineWidth',2);
        half_max_w = I1+1000-I2;  
    else
        line([I1 I2],[half_max half_max],'color','r','LineWidth',2);
        half_max_w = I2-I1;
    end  
    text(500,ex_bump_mag,num2str(half_max_w))

    %convert to EB tiles
    half_max_width = half_max_w*8/1001;    
end


%% Repeat with all the data points

for timepoint = 1:size(data.mean_dff_EB,2)
    
    %linearly interpolate to have 1000 datapoints instead of 8
    interp_ex_data = interp1([1:8],data.mean_dff_EB(:,timepoint),[1:7/1000:8]);
    %Find the half max point
    half_max = (max(data.mean_dff_EB(:,timepoint))-min(data.mean_dff_EB(:,timepoint)))/2 + min(data.mean_dff_EB(:,timepoint));
    [ex_bump_mag_interp I_interp] = max(interp_ex_data);
    %Find in each half the index closest to the half max
    diff_data = abs(interp_ex_data-half_max);
    [sortedVals,indexes] = sort(diff_data);
    %remove all the indexes that are less than 175 datapoints away (i.e., a little more than 1 glomerulus) from
    %index(1)
    diff_indexes = abs(indexes-indexes(1));
    indexes(diff_indexes<175 & diff_indexes>0)=NaN;
    indexes = indexes(~isnan(indexes));
    two_indexes = [indexes(1), indexes(2)];
    I1 = min(two_indexes);
    I2 = max(two_indexes);
    if (all(two_indexes>I_interp) | all(two_indexes<I_interp))
        half_max_w = I1+1000-I2;  
    else
        half_max_w = I2-I1;
    end  
    %convert to EB tiles
    half_max_width(timepoint) = half_max_w*8/1001; 
    
end

%% Analyze per block

for block = 1:length(blockLimits)
   width_half_max{block} = half_max_width(blockLimits{block}(1):blockLimits{block}(2)); 
end

%plot
figure('Position',[200 200 1000 800]),
subplot(2,1,1)
for block = 1:length(blockLimits)
   mean_half_w(block) = nanmean(width_half_max{block}); 
   plot(block,mean_half_w(block),'ko','MarkerFaceColor',color_gradient{block},'MarkerSize',8)
   hold on
end
xlim([0 7]);
xticks(1:6);
title('Mean half max width per block');
ylabel('Half max width');
%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
if gain_yaw(1) == 1
    legend(h,'Normal gain','Inverted gain');
else
    legend(h,'Inverted gain','Normal gain');
end

subplot(2,1,2)
%add bump mag data to table
if size(summary_data,2) == 3
    summary_data = addvars(summary_data, mean_half_w','NewVariableNames','mean_half_width');
end
%Get mean bump mag by contrast
mean_half_width = varfun(@mean,summary_data,'InputVariables','mean_half_width',...
       'GroupingVariables',{'gain'});
%Plot
plot(mean_half_width.mean_mean_half_width,'-ko','MarkerSize',8)
xlim([0 3]);
ylim([min(mean_half_width.mean_mean_half_width)-0.5 max(mean_half_width.mean_mean_half_width)+0.5]); 
xticks(1:3);
if gain_yaw(1) == 1
    xticklabels({'Normal gain','Inverted gain'});
else
    xticklabels({'Inverted gain','Normal gain'});
end
title('Width at half max per gain');
ylabel('Full width at half max');

%save figure
saveas(gcf,[path,'plots\gainChangeMeanHW.png']);

% %% Plot fly's velocity in all 3 axes  with all vel in deg/s
% 
% %correct forward velocity units!
% 
% figure('Position',[100, 100, 1600, 1000]),
% subplot(4,5,[1 4])
% plot(data.time,data.vel_for_deg_ds)
% title('Forward velocity');
% ylabel('Forward velocity (deg/s)');
% 
% subplot(4,5,5)
% histogram(data.vel_for_deg_ds)
% title('Velocity distributions');
% xlabel('Forward velocity (deg/s)');
% ylabel('Counts');
% 
% subplot(4,5,[6 9])
% plot(data.time,data.vel_side_deg_ds,'color',[0.8 0.2 0.6])
% title('Side velocity');
% ylabel('Side velocity (deg/s)');
% 
% subplot(4,5,10)
% histogram(data.vel_side_deg_ds,'FaceColor',[0.8 0.2 0.6])
% xlabel('Side velocity (deg/s)');
% ylabel('Counts');
% 
% subplot(4,5,[11 14])
% plot(data.time,data.vel_yaw_ds,'color',[0.4 0.2 0.8])
% title('Angular velocity');
% ylabel('Angular velocity (deg/s)');
% 
% subplot(4,5,15)
% histogram(data.vel_yaw_ds,'FaceColor',[0.4 0.2 0.8])
% xlabel('Angular velocity (deg/s)');
% ylabel('Counts');
% 
% subplot(4,5,[16 19])
% plot(data.time,data.total_mvt_ds,'k')
% title('Total movement');
% ylabel('Total movement (deg/s)');
% xlabel('Time (sec)');
% 
% subplot(4,5,20)
% histogram(data.total_mvt_ds,'FaceColor','k')
% xlabel('Total movement (deg/s)');
% ylabel('Counts');
% 
% %save figure
% saveas(gcf,[path,'plots\gainChangeVelDistributions.png']);
% 
% %% Compute and plot median total movement per block
% 
% figure('Position',[200 200 1400 800]),
% %Plot total movement in time per block, colored appropriately
% subplot(2,5,[1 3])
% for block = 1:length(blockLimits)
%     if contains(contrasts(block),'Dark')
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{1})
%     elseif contains(contrasts(block),'Low')
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{2})
%     else
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{3})
%     end
%     hold on
% end
% title('Total movement in time');
% ylabel('Total movement (deg/s)');
% xlabel('Time (sec)');
% xlim([0 data.time(end)]);
% ylim([0 max(data.total_mvt_ds)+10]);
% 
% %Plot median total movement per block
% subplot(2,5,[6 8])
% for block = 1:length(blockLimits)
%    mean_total_mvt(block) = nanmean(data.total_mvt_ds(1,blockLimits{block}(1):blockLimits{block}(2))); 
%    plot(block,mean_total_mvt(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
%    hold on
% end
% if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
%     xlim([0 6]);
%     xticks(1:5);    
% else
%     xlim([0 7]);
%     xticks(1:6);
% end
% title('Mean total movement per block');
% ylabel('Total movement (deg/s)');
% ylim([0 max(mean_total_mvt)+10]);
% xlabel('Block #');
% %Add custom legend
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
% h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
% h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
% legend(h, 'Darkness','Low contrast','High Contrast');
% 
% %Plot median total movement per block type
% subplot(2,5,[4,5,9,10])
% if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
%     Intensities = Intensities(1:5);
% end
% for contrast = 1:3
%    mean_mean_total_mvt(contrast) = mean(mean_total_mvt(find(Intensities==contrast)));
% end
% plot(mean_mean_total_mvt,'-ko','LineWidth',2)
% ylabel('Total movement (deg/s)');
% title('Mean total movement');
% ylim([0 max(mean_mean_total_mvt+10)]);
% xlim([0 4]);
% xticks(1:3);
% xticklabels({'Darkness','Low contrast','High contrast'});
% 
% %save figure
% saveas(gcf,[path,'plots\closedLoopTotalMvtPerBlock.png']);
% 
% %% Repeat for fwd vel
% 
% figure('Position',[200 200 1400 800]),
% %Plot fwd vel in time per block, colored appropriately
% subplot(2,5,[1 3])
% for block = 1:length(blockLimits)
%     if contains(contrasts(block),'Dark')
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{1})
%     elseif contains(contrasts(block),'Low')
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{2})
%     else
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2)),'color',color_gradient{3})
%     end
%     hold on
% end
% title('Forward velocity in time');
% ylabel('Foward velocity (mm/s)');
% xlabel('Time (sec)');
% xlim([0 data.time(end)]);
% ylim([0 max(data.vel_for_ds)+1]);
% 
% %Plot median total movement per block
% subplot(2,5,[6 8])
% for block = 1:length(blockLimits)
%    mean_fwd_vel(block) = nanmean(data.vel_for_ds(blockLimits{block}(1):blockLimits{block}(2))); 
%    plot(block,mean_fwd_vel(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
%    hold on
% end
% if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
%     xlim([0 6]);
%     xticks(1:5);    
% else
%     xlim([0 7]);
%     xticks(1:6);
% end
% title('Mean forward velocity per block');
% ylabel('Forward velocity (mm/s)');
% ylim([0 max(mean_fwd_vel)+1]);
% xlabel('Block #');
% %Add custom legend
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
% h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
% h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
% legend(h, 'Darkness','Low contrast','High Contrast');
% 
% %Plot median total movement per block type
% subplot(2,5,[4,5,9,10])
% for contrast = 1:3
%    mean_mean_fwd_vel(contrast) = mean(mean_fwd_vel(find(Intensities==contrast)));
% end
% plot(mean_mean_fwd_vel,'-ko','LineWidth',2)
% ylabel('Forward velocity (mm/s)');
% title('Mean forward velocity');
% ylim([0 max(mean_mean_fwd_vel)+1]);
% xlim([0 4]);
% xticks(1:3);
% xticklabels({'Darkness','Low contrast','High contrast'});
% 
% %save figure
% saveas(gcf,[path,'plots\closedLoopFwdVelPerBlock.png']);
% 
% %% Repeat for angular speed
% 
% figure('Position',[200 200 1400 800]),
% %Plot fwd vel in time per block, colored appropriately
% subplot(2,5,[1 3])
% for block = 1:length(blockLimits)
%     if contains(contrasts(block),'Dark')
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))),'color',color_gradient{1})
%     elseif contains(contrasts(block),'Low')
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))),'color',color_gradient{2})
%     else
%         plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2))),'color',color_gradient{3})
%     end
%     hold on
% end
% title('Angular speed in time');
% ylabel('Angular speed (deg/s)');
% xlabel('Time (sec)');
% xlim([0 data.time(end)]);
% ylim([0 max(abs(data.vel_yaw_ds))+10]);
% 
% %Plot median total movement per block
% subplot(2,5,[6 8])
% for block = 1:length(blockLimits)
%    mean_ang_speed(block) = nanmean(abs(data.vel_yaw_ds(blockLimits{block}(1):blockLimits{block}(2)))); 
%    plot(block,mean_ang_speed(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
%    hold on
% end
% if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
%     xlim([0 6]);
%     xticks(1:5);    
% else
%     xlim([0 7]);
%     xticks(1:6);
% end
% title('Mean angular speed per block');
% ylabel('Angular speed (deg/s)');
% ylim([0 max(mean_ang_speed)+10]);
% xlabel('Block #');
% %Add custom legend
% h = zeros(3, 1);
% h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
% h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
% h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
% legend(h, 'Darkness','Low contrast','High Contrast');
% 
% %Plot median total movement per block type
% subplot(2,5,[4,5,9,10])
% for contrast = 1:3
%    mean_mean_ang_speed(contrast) = mean(mean_ang_speed(find(Intensities==contrast)));
% end
% plot(mean_mean_ang_speed,'-ko','LineWidth',2)
% ylabel('Angular speed (deg/s)');
% title('Mean angular speed');
% ylim([0 max(mean_mean_ang_speed)+10]);
% xlim([0 4]);
% xticks(1:3);
% xticklabels({'Darkness','Low contrast','High contrast'});
% 
% %save figure
% saveas(gcf,[path,'plots\closedLoopAngSpeedPerBlock.png']);
% 
% %% Relationship between bump magnitude and velocity
% 
% allBumpMag = [];
% for block = 1:length(blockLimits)
%     allBumpMag = [allBumpMag,bump_mag{block}];
% end
% 
% figure('Position',[200 200 1600 800]),
% subplot(1,4,1)
% scatter(data.vel_for_deg_ds(1:length(allBumpMag)),allBumpMag)
% xlim([0 max(data.vel_for_deg_ds)]);
% xlabel('Forward velocity (deg/s)');
% ylabel('Bump magnitude (max-min)');
% 
% subplot(1,4,2)
% scatter(abs(data.vel_side_deg_ds(1:length(allBumpMag))),allBumpMag,[],[0.8 0.2 0.6])
% xlim([0 max(abs(data.vel_side_deg_ds))]);
% xlabel('Side speed (deg/s)');
% 
% subplot(1,4,3)
% scatter(abs(data.vel_yaw_ds(1:length(allBumpMag))),allBumpMag,[],[0.4 0.2 0.8])
% xlim([0 max(abs(data.vel_yaw_ds))]);
% xlabel('Angular speed (deg/s)');
% 
% subplot(1,4,4)
% scatter(data.total_mvt_ds(1:length(allBumpMag)),allBumpMag,[],'k')
% xlim([0 max(data.total_mvt_ds)]);
% xlabel('Total movement (deg/s)');
% 
% %save figure
% saveas(gcf,[path,'plots\closedLoopBMvsVel.png']);
% 
% %% Binning it
% 
% %Let's always define a fixed number of vel bins
% nbins = 20;
% maxBin = max(data.total_mvt_ds);
% binWidth = maxBin/nbins;
% mvtBins = [0:binWidth:maxBin]; 
% 
% %getting binned means 
% for bin = 1:length(mvtBins)-1
%     meanBin(bin) = mean(allBumpMag((data.total_mvt_ds(1:length(allBumpMag)) > mvtBins(bin)) & (data.total_mvt_ds(1:length(allBumpMag)) < mvtBins(bin+1))));
% end
% 
% %create axes for plot
% mvtAxes = mvtBins - binWidth;
% mvtAxes = mvtAxes(2:end);
% mvtAxes(end) = mvtAxes(end-1)+binWidth;
% 
% %Plot
% figure('Position',[200 200 1400 600]),
% %Plot bump magnitude in time
% subplot(2,4,[1 3])
% plot(data.time(1:length(allBumpMag)),allBumpMag, 'k')
% ylabel('Bump magnitude');
% ylim([0 5]);
% set(gca,'xticklabel',{[]});
% 
% %Plot total movement 
% subplot(2,4,[5 7])
% plot(data.time(1:length(allBumpMag)),data.total_mvt_ds(1:length(allBumpMag)),'k')
% xlabel('Time (sec)');
% ylabel('Total movement (deg/s)');
% 
% %Plot relationship between both parameters
% subplot(2,4,[4,8]);
% plot(mvtAxes,meanBin,'-ko')
% ylabel('Mean bump magnitude'); xlabel('Total movement (deg/s)');
% ylim([0 (max(meanBin)+0.5)]);
% 
% %save figure
% saveas(gcf,[path,'plots\closedLoopBMvsVelBinned.png']);
% 
% %% Model bump magnitude as a function of contrastLevel and total movement
% 
% %create vector with the contrast level for each timepoint
% all_contrast_levels = [];
% for block = 1:length(blockLimits)
%     contrast_level{block} = repelem(Intensities(block),blockLimits{block}(2)+1-blockLimits{block}(1));
%     all_contrast_levels = [all_contrast_levels,contrast_level{block}];
% end
% 
% %crate table with the model's variables
% modelTable = table(all_contrast_levels',data.total_mvt_ds(1:length(allBumpMag))',data.time(1:length(allBumpMag)),allBumpMag','VariableNames',{'ContrastLevel','TotalMovement','Time','BumpMagnitude'});
% 
% %fit linear model using contrast level as a categorical variable
% mdl_BM = fitlm(modelTable,'BumpMagnitude~ContrastLevel+TotalMovement','CategoricalVars',1)
% %This results in very small weigths being assigned to total movement, but
% %that probably is because being in deg/s, total movement has large units.
% 
% %fit linear model using contrast level as a numerical variable
% mdl_BM2 = fitlm(modelTable,'BumpMagnitude~ContrastLevel+TotalMovement+Time');

%% Heading variation per stim

figure('Position',[200 200 1000 800]),
subplot(2,1,1)
for block = 1:length(blockLimits)
   [~, heading_var(block)] = circ_std(deg2rad(heading(blockLimits{block}(1):blockLimits{block}(2))),[],[],1); 
   plot(block,heading_var(block),'ko','MarkerFaceColor',color_gradient{block},'MarkerSize',8)
   hold on
end
xlim([0 7]);
xticks(1:6);
title('Heading variation per block');
ylabel({'Circular standard deviation','of the heading'});
%Add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
if gain_yaw(1) == 1
    legend(h,'Normal gain','Inverted gain');
else
    legend(h,'Inverted gain','Normal gain');
end

subplot(2,1,2)
%Create table with contrast level and offset variation
if size(summary_data,2) == 4
    summary_data = addvars(summary_data, heading_var','NewVariableNames','heading_var');
end
%Get mean heading var by contrast
mean_heading_data = varfun(@mean,summary_data,'InputVariables','heading_var',...
       'GroupingVariables',{'gain'});
%Plot
plot(mean_heading_data.mean_heading_var,'-ko','MarkerSize',8)
xlim([0 3]);
ylim([min(mean_heading_data.mean_heading_var)-0.5 max(mean_heading_data.mean_heading_var)+0.5]); 
xticks(1:3);
if gain_yaw(1) == 1
    xticklabels({'Normal gain','Inverted gain'});
else
    xticklabels({'Inverted gain','Normal gain'});
end
title('Mean heading variation per gain');
ylabel({'Circular standard deviation','of the heading'});

%save figure
saveas(gcf,[path,'plots\gainChangeHeadingVariation.png']);



%% Save useful data

save([path,'\summary_data.mat'],'summary_data','mean_reference_offset')

close all; clc;