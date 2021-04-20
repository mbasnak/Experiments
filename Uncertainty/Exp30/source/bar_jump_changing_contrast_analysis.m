%analysis for the stabilizing offset bout


%% Load data

clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data');
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


%% Get jump frames

jump_frames = uint32(data.jump_frames_ds);

%for flies that had fictrac or panels errors, remove the corresponding
%jumps
if (contains(path,'20210413_60D05_7f') & ~contains(path,'fly'))
   jump_frames = jump_frames(1:end-4); 
end

%% Calibrate jump frames

bar_position = wrapTo180(data.panel_angle);

recalibrated_jump_frames = jump_frames;
for jump = 1:length(jump_frames)
    figure,
    subplot(1,2,1)
    time_around_jump_cu = data.time(jump_frames(jump)-5:jump_frames(jump)+5) - data.time(jump_frames(jump));
    bar_position_AJ_cu(jump,:) = bar_position(jump_frames(jump)-5:jump_frames(jump)+5);
    plot(bar_position_AJ_cu(jump,:),'-o','linewidth',2,'color',[0.4 0.8 0.3])
    %determine jump size
    bar_jump_size(jump) = wrapTo180(bar_position_AJ_cu(jump,7)-bar_position_AJ_cu(jump,6));
    text(-0.4,100,['Bar jump size = ',num2str(bar_jump_size(jump))],'fontsize',14);
    xline(6,'linewidth',2,'color','r');
    ylim([-180 180]);
    %xlim([-1 1]);
    title(['Bar position, Jump #',num2str(jump)]);
    
    %if the bar jump is incorrect (one or two datapoints away from truth),
    %recalibrate
    for frame = 1:size(bar_position_AJ_cu,2)-1
        diff_bar_pos(frame) = abs(rad2deg(circ_dist(deg2rad(bar_position_AJ_cu(jump,frame+1)),deg2rad(bar_position_AJ_cu(jump,frame)))));
    end
    pos = find(diff_bar_pos>85 & diff_bar_pos<150);
    recalibrated_jump_frames(jump) = recalibrated_jump_frames(jump) + pos(end) - 6;
    
    %Plot again to check
    subplot(1,2,2)
    bar_position_AJ_cu(jump,:) = bar_position(recalibrated_jump_frames(jump)-5:recalibrated_jump_frames(jump)+5);
    plot(bar_position_AJ_cu(jump,:),'-o','linewidth',2,'color',[0.4 0.8 0.3])
    %determine jump size
    bar_jump_size(jump) = wrapTo180(bar_position_AJ_cu(jump,7)-bar_position_AJ_cu(jump,6));
    text(-0.4,100,['Bar jump size = ',num2str(bar_jump_size(jump))],'fontsize',14);
    xline(6,'linewidth',2,'color','r');
    ylim([-180 180]);
    %xlim([-1 1]);
    title(['Bar position, Jump #',num2str(jump)]);
    
end

% Plot the bar jump sizes

figure,
plot(bar_jump_size,'-ro')
ylim([-180 180]);
xlabel('Jump #');
ylabel('Bar jump size');

close all;
%% Plot heatmap with bar position, phase and offset

figure('Position',[200 200 1000 600]),
subplot(3,6,[1 5])
%Plot heatmap of EPG activity
imagesc(data.dff_matrix)
colormap(flipud(gray))
ylabel('PB glomerulus');
title('EPG activity in the PB');

subplot(3,6,[7 11])
%Get bar position to plot
[x_out_bar,bar_position_to_plot] = removeWrappedLines(data.time,bar_position);
plot(x_out_bar,bar_position_to_plot,'LineWidth',1.5)
hold on
%Get EPG phase to plot
phase = wrapTo180(rad2deg(data.dff_pva));
[x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
plot(x_out_phase,phase_to_plot,'LineWidth',1.5)
xlim([0 x_out_phase(end-1)]);
ylim([-180 180]);
legend('Bar position','EPG phase');
ylabel('Deg');
title('Bar and bump position');

subplot(3,6,[13 17])
%Get offset to plot
offset = wrapTo180(data.bar_offset);
offset_hc = offset(data.contrast_ds == 1);
time_hc = data.time(data.contrast_ds == 1);
offset_hc(diff(time_hc)>100) = NaN; %add nans when the contrast changes to remove spurious plot lines
[x_out_offset_hc,offset_hc_to_plot] = removeWrappedLines(time_hc,offset_hc);
offset_lc = offset(data.contrast_ds == 0);
time_lc = data.time(data.contrast_ds == 0);
offset_lc(diff(time_lc)>30) = NaN; %add nans when the contrast changes to remove spurious plot lines
[x_out_offset_lc,offset_lc_to_plot] = removeWrappedLines(time_lc,offset_lc);
plot(x_out_offset_hc,offset_hc_to_plot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
hold on
plot(x_out_offset_lc,offset_lc_to_plot,'color',[0 0 0.6],'LineWidth',1.5)
jump_times = recalibrated_jump_frames.*data.time(end)/length(data.time); %convert jump frames to time units
for jump = 1:length(recalibrated_jump_frames)
    line([jump_times(jump) jump_times(jump)],[-180 180],'color','r','linewidth',2)
end
%xlim([0 length(data.bar_offset)]);
ylim([-180 180]);
xlabel('Time');
ylabel('Deg');
title('Offset');

subplot(3,6,18)
polarhistogram(deg2rad(offset),'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title('Offset distribution');

saveas(gcf,[path,'plots\Heatmap.png']);

%% Compute bump magnitude and bump width at half max

bump_mag = data.bump_magnitude;

% Bump width at half max
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


%% Close up around the jumps

for jump = 1:length(recalibrated_jump_frames)
    
   figure,
   subplot(2,2,[1 2])
   imagesc(data.dff_matrix(:,recalibrated_jump_frames(jump)-20:recalibrated_jump_frames(jump)+20))
   colormap(flipud(gray))
   hold on
   xline(21,'linewidth',2,'color','r')
   title(['Jump #',num2str(jump)]);
   
   %compute the mean bump parameters right before the jump
   mean_bm_pre_jump(jump) = mean(bump_mag(recalibrated_jump_frames(jump)-10:recalibrated_jump_frames(jump)));
   mean_bw_pre_jump(jump) = mean(half_max_width(recalibrated_jump_frames(jump)-10:recalibrated_jump_frames(jump)));
   
   %add as text
   subplot(2,2,3)
   plot(0,0)
   ylim([-1 1]);
   xlim([-1 1]);
   text(-0.8,0.5,['Mean bump mag = ',num2str(round(mean_bm_pre_jump(jump),2))]);
   text(-0.8,-0.5,['Mean bump width = ',num2str(round(mean_bw_pre_jump(jump),2))]);
   set(gca,'xticklabel',{[]});
   set(gca,'yticklabel',{[]});
   
   %get bump jump magnitude
   bump_jump_mag(jump) = abs(rad2deg(circ_dist(circ_mean(data.dff_pva(recalibrated_jump_frames(jump)-5:recalibrated_jump_frames(jump))),circ_mean(data.dff_pva(recalibrated_jump_frames(jump)+5:recalibrated_jump_frames(jump)+10)))));
   subplot(2,2,4)
   plot(0,0)
   ylim([-1 1]);
   xlim([-1 1]);
   text(-0.8,0,['Bump jump mag = ',num2str(round(bump_jump_mag(jump),2))]);
   set(gca,'xticklabel',{[]});
   set(gca,'yticklabel',{[]});
      
end

close all;
%% Plot the offsets around the jumps

offset_AJ = [];
for jump = 1:length(recalibrated_jump_frames)
    offset_AJ(jump,:) = offset(recalibrated_jump_frames(jump)-90:recalibrated_jump_frames(jump)+90); 
end

figure('Position',[100 100 1200 800]),
for jump = 1:length(recalibrated_jump_frames)
    subplot(length(recalibrated_jump_frames),1,jump)
    plot(offset_AJ(jump,:),'-ko','linewidth',2)
    hold on
    xline(91,'linewidth',2,'color','r');
    ylim([-180 180]);
    if jump ~= (length(recalibrated_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Offset around the jumps');


%% Shift the offset right before the jump

shifted_offset_AJ_short = wrapTo180(offset_AJ-offset_AJ(:,91));

% Sort by the value of the mean BM and BW pre-jump
[sorted_BM,I] = sort(mean_bm_pre_jump);
[sorted_BW,Iw] = sort(mean_bw_pre_jump);

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
imagesc(shifted_offset_AJ_short(I,:))
colorbar
gray_wrap = vertcat(gray,flipud(gray));
colormap(gray_wrap);
hold on
xline(91.5,'linewidth',2,'color','r');
title('Offset around the jumps, sorted by BM');
ylabel('Sorted jump #');

subplot(1,2,2)
imagesc(shifted_offset_AJ_short(Iw,:))
colorbar
gray_wrap = vertcat(gray,flipud(gray));
colormap(gray_wrap);
hold on
xline(91.5,'linewidth',2,'color','r');
title('Offset around the jumps, sorted by BW');
ylabel('Sorted jump #');


%% Sort by the contrast pre jump

pre_jump_hc = find(data.contrast_ds(uint32(recalibrated_jump_frames)-5)==1);
pre_jump_lc = find(data.contrast_ds(uint32(recalibrated_jump_frames)-5)==0);

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
imagesc(shifted_offset_AJ_short(pre_jump_hc,:))
colorbar
gray_wrap = vertcat(gray,flipud(gray));
colormap(gray_wrap);
hold on
xline(91.5,'linewidth',2,'color','r');
title('Offset around the jumps for high contrast');
ylabel('Jump #');

subplot(1,2,2)
imagesc(shifted_offset_AJ_short(pre_jump_lc,:))
colorbar
gray_wrap = vertcat(gray,flipud(gray));
colormap(gray_wrap);
hold on
xline(91.5,'linewidth',2,'color','r');
title('Offset around the jumps for low contrast');
ylabel('Jump #');

%% Get offset around the jumps for more frames, to compute the time it takes to go back to its pre-jump value

longer_offset_AJ = [];
for jump = 1:length(recalibrated_jump_frames)
    longer_offset_AJ(jump,:) = offset(recalibrated_jump_frames(jump)-250:recalibrated_jump_frames(jump)+250); 
end

figure('Position',[100 100 1200 800]),
for jump = 1:length(recalibrated_jump_frames)
    recalibrated_jump_frames = uint32(recalibrated_jump_frames);
    time_around_jump = data.time(recalibrated_jump_frames(jump)-250:recalibrated_jump_frames(jump)+250) - data.time(recalibrated_jump_frames(jump));
    subplot(length(recalibrated_jump_frames),1,jump)
    plot(longer_offset_AJ(jump,:),'-k','linewidth',2)
    hold on
    xline(251,'linewidth',2,'color','r');
    ylim([-180 180]);
    %xlim([-26 26]);
    if jump ~= (length(recalibrated_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Longer offset around the jumps');


%% Sort by BM value pre jump

figure('Position',[100 100 1200 800]),
for jump = 1:length(recalibrated_jump_frames)
    subplot(length(recalibrated_jump_frames),1,jump)
    plot(longer_offset_AJ(I(jump),:),'-k','linewidth',2)
    hold on
    xline(251,'linewidth',2,'color','r');
    ylim([-180 180]);
    if jump ~= (length(recalibrated_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Longer offset around the jumps,sorted');


%% Find the number of frames it takes for the offset to return to its pre-jump value

close all;
%Define offset value pre-jump as the circular mean in the 10 sec preceding the jump
pre_jump_offset = rad2deg(circ_mean(deg2rad(longer_offset_AJ(:,151:251)),[],2));

%find frames after the jump that are within 10 degrees of that offset
return_frames = {};
for jump = 1:length(recalibrated_jump_frames)
    offset_difference = rad2deg(circ_dist(deg2rad(longer_offset_AJ(jump,252:end)),deg2rad(pre_jump_offset(jump))));
    return_frames{jump} = find(abs(offset_difference)<=20);
    if ~isempty(return_frames{jump})
        first_return_frame(jump) = return_frames{jump}(1);
    else
        first_return_frame(jump) = NaN;
    end
end

%Plot
for jump = 1:length(recalibrated_jump_frames)
    figure('Position',[100 100 1600 400]),
    plot(time_around_jump,longer_offset_AJ(I(jump),:),'-k','linewidth',2)
    hold on
    xline(0,'linewidth',2,'color','r');
    yline(pre_jump_offset(I(jump)),'-.b','linewidth',2);
    legend('offset','cue jump','pre jump mean offset');
    if ~isnan(first_return_frame(I(jump)))
        plot(time_around_jump(first_return_frame(I(jump))+251),longer_offset_AJ(I(jump),first_return_frame(I(jump))+251),'ro','MarkerSize',12,'MarkerFaceColor','r')
        legend('offset','cue jump','pre jump mean offset','return frame')
    end
    ylim([-180 180]);
    xlim([-6 26]);
end

%Plot the offset return frames sorted by BM
figure,
subplot(1,2,1)
plot(first_return_frame(I),'o')
ylabel('Frame that returned to pre-jump offset');
xlabel('Jumps sorted by BM');

% Plot the offset return frames sorted by BW
[sorted_BW,Iw] = sort(mean_bw_pre_jump);
subplot(1,2,2)
plot(first_return_frame(Iw),'o')
ylabel('Frame that returned to pre-jump offset');
xlabel('Jumps sorted by BW');
%saveas(gcf,[path,'plots\offset_return_sorted_frames.png']);


%% Classify these return frames by the pre-jump cue they got

figure,
plot(first_return_frame(pre_jump_hc),'o')
hold on
plot(first_return_frame(pre_jump_lc),'ro')
ylabel('Frame that returned to pre-jump offset');
xlabel('Jumps');
legend({'high contrast','low contrast'})


%% Close up plot of the epg activity, the fly and bump position, and the offset, around the jumps, coloring the pre-cue bar contrast, and adding the return frame for the offset and the mean BM and BW

% 
% for jump = 1:length(recalibrated_jump_frames)
%     
%     figure('Position',[200 200 1000 600]),
%     subplot(3,1,1)
%     %Plot heatmap of EPG activity
%     imagesc(data.dff_matrix)
%     colormap(flipud(gray))
%     hold on
%     line([recalibrated_jump_frames(jump) recalibrated_jump_frames(jump)],[0 17],'color','r','linewidth',2);
%     ylabel('PB glomerulus');
%     title('EPG activity in the PB');
%     xlim([recalibrated_jump_frames(jump)-90 recalibrated_jump_frames(jump)+180]);
%     
%     subplot(3,1,2)
%     %Get bar position to plot
%     [x_out_bar,bar_position_to_plot] = removeWrappedLines(data.time,bar_position);
%     plot(x_out_bar,bar_position_to_plot,'LineWidth',1.5)
%     hold on
%     %Get EPG phase to plot
%     phase = wrapTo180(rad2deg(data.dff_pva));
%     [x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
%     plot(x_out_phase,phase_to_plot,'LineWidth',1.5)
%     %add cue jump
%     frames_to_sec = data.time(end)/length(data.time);
%     line([(recalibrated_jump_frames(jump))*frames_to_sec (recalibrated_jump_frames(jump))*frames_to_sec],[-180 180],'color','r','linewidth',2);   
%     xlim([(recalibrated_jump_frames(jump)-90)*frames_to_sec (recalibrated_jump_frames(jump)+180)*frames_to_sec]);
%     ylim([-180 180]);
%     legend('Bar position','EPG phase');
%     ylabel('Deg');
%     title('Bar and bump position');
%     
%     subplot(3,1,3)
%     %Get offset to plot
%     offset = wrapTo180(data.bar_offset);
%     offset_hc = offset(data.contrast_ds == 1);
%     time_hc = data.time(data.contrast_ds == 1);
%     offset_hc(diff(time_hc)>100) = NaN; %add nans when the contrast changes to remove spurious plot lines
%     [x_out_offset_hc,offset_hc_to_plot] = removeWrappedLines(time_hc,offset_hc);
%     offset_lc = offset(data.contrast_ds == 0);
%     time_lc = data.time(data.contrast_ds == 0);
%     offset_lc(diff(time_lc)>30) = NaN; %add nans when the contrast changes to remove spurious plot lines
%     [x_out_offset_lc,offset_lc_to_plot] = removeWrappedLines(time_lc,offset_lc);
%     plot(x_out_offset_hc,offset_hc_to_plot,'color',[ 0.5 0.8 0.9],'LineWidth',1.5)
%     hold on
%     plot(x_out_offset_lc,offset_lc_to_plot,'color',[0 0 0.6],'LineWidth',1.5)
%     line([(recalibrated_jump_frames(jump))*frames_to_sec (recalibrated_jump_frames(jump))*frames_to_sec],[-180 180],'color','r','linewidth',2);
%     %add bump parameters
%     text(((recalibrated_jump_frames(jump))*frames_to_sec)-3,100,['BM = ',num2str(round(mean_bm_pre_jump(jump),2))])
%     text(((recalibrated_jump_frames(jump))*frames_to_sec)-3,50,['BW = ',num2str(round(mean_bw_pre_jump(jump),2))])
%     %add offset return
%     if ~isnan(first_return_frame(jump))
%         return_time = recalibrated_jump_frames(jump)*frames_to_sec + first_return_frame(jump)*frames_to_sec;
%         plot(return_time,offset(return_time),'ro')
%         plot(return_time,offset(return_time),'ro')
%     end
%     xlim([(recalibrated_jump_frames(jump)-90)*frames_to_sec (recalibrated_jump_frames(jump)+180)*frames_to_sec]);
%     ylim([-180 180]);
%     xlabel('Time');
%     ylabel('Deg');
%     title('Offset');
% 
% end

%% Save relevant data

save([path,'\offset_change.mat'],'time_around_jump','mean_bm_pre_jump','mean_bw_pre_jump','offset_AJ','longer_offset_AJ','first_return_frame','pre_jump_hc','pre_jump_lc');

close all; clear all;