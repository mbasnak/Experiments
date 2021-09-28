%bar jump analysis


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


%% Find the jump frames

%use the panels signal to detect the jumps
figure,
subplot(2,1,1)
bar_position = data.panel_angle;
plot(bar_position)
title('Bar position')
ylim([0 360]);
xlim([0 length(bar_position)]);

subplot(2,1,2)
abs_diff_signal = abs(diff(unwrap(deg2rad(bar_position))));
plot(abs_diff_signal)
ds_jump_frames = find(abs_diff_signal>1.5);
hold on
plot(ds_jump_frames,abs_diff_signal(ds_jump_frames),'ro')
title('Absolute of the bar position derivative')
ylim([0 3.5]);
xlim([0 length(bar_position)]);

%% Remove some of the frames depending on the fly id (due to the fly not moving anymore or other errors)

if contains(path,'20210305_60D05_7f_fly4')==1
   ds_jump_frames = ds_jump_frames([1:11,13:15],:); %remove the jumps after that because she stopped walking
elseif (contains(path,'20210311_60D05_7f')==1) & (contains(path,'20210311_60D05_7f_fly3')==0)
   ds_jump_frames(20,:) = [];
elseif (contains(path,'20210312_60D05_7f')==1) & (contains(path,'20210312_60D05_7f_fly2')==0) & (contains(path,'20210312_60D05_7f_fly3')==0)
   ds_jump_frames = ds_jump_frames(1:15);
elseif contains(path,'20210312_60D05_7f_fly2')==1
   ds_jump_frames(10,:) = []; 
elseif contains(path,'20210312_60D05_7f_fly3')==1
   ds_jump_frames([13:20],:) = []; 
elseif contains(path,'20210323_60D05_7f_fly2')==1
   ds_jump_frames(20,:) = []; 
elseif contains(path,'20210324_60D05_7f')==1
   ds_jump_frames(20,:) = [];   
elseif contains(path,'20210325_60D05_7f')==1
   ds_jump_frames(20,:) = [];
end

%convert to sec
ds_jump_sec = ds_jump_frames.*(data.time(end)/length(data.time));
%% Plot heatmap with bar position, phase and offset

figure('Position',[200 200 1000 600]),
%Plot heatmap of EPG activity
subplot(3,1,1)
imagesc(data.dff_matrix)
colormap(flipud(gray))
ylabel('PB glomerulus');
title('EPG activity in the PB');
hold on
for jump = 1:length(ds_jump_frames)
   xline(ds_jump_frames(jump),'color',[0, 0.5, 0],'lineWidth',3); 
end

subplot(3,1,2)
%Get bar position to plot
bar_position = wrapTo180(data.panel_angle);
[x_out_bar,bar_position_to_plot] = removeWrappedLines(data.time,bar_position);
plot(x_out_bar,bar_position_to_plot,'LineWidth',1.5,'color',[0.2 0.6 0.7])
hold on
%Get EPG phase to plot
phase = wrapTo180(rad2deg(data.dff_pva));
[x_out_phase,phase_to_plot] = removeWrappedLines(data.time,phase);
plot(x_out_phase,phase_to_plot,'color',[0.9 0.3 0.4],'LineWidth',1.5)
xlim([0 x_out_phase(end)]);
ylim([-180 180]);
legend('Bar position','EPG phase')
ylabel('Deg');
title('Bar and bump position');
for jump = 1:length(ds_jump_frames)
   xline(ds_jump_sec(jump),'color',[0, 0.5, 0],'lineWidth',3,"HandleVisibility", 'off'); 
end

subplot(3,1,3)
%Get offset to plot
offset = wrapTo180(data.bar_offset);
[x_out_offset,offset_to_plot] = removeWrappedLines(data.time,offset);
plot(x_out_offset,offset_to_plot,'k','LineWidth',1.5)
xlim([0 x_out_offset(end)]);
ylim([-180 180]);
xlabel('Time');
ylabel('Deg');
title('Offset');
hold on
for jump = 1:length(ds_jump_frames)
   xline(ds_jump_sec(jump),'color',[0, 0.5, 0],'lineWidth',3); 
end

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

for jump = 1:length(ds_jump_frames)
    
   figure,
   subplot(2,2,[1 2])
   imagesc(data.dff_matrix(:,ds_jump_frames(jump)-20:ds_jump_frames(jump)+20))
   colormap(flipud(gray))
   hold on
   xline(21,'linewidth',2,'color','r')
   title(['Jump #',num2str(jump)]);
   
   %compute the mean bump parameters right before the jump
   mean_bm_pre_jump(jump) = mean(bump_mag(ds_jump_frames(jump)-10:ds_jump_frames(jump)));
   mean_bw_pre_jump(jump) = mean(half_max_width(ds_jump_frames(jump)-10:ds_jump_frames(jump)));
   
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
   bump_jump_mag(jump) = abs(rad2deg(circ_dist(circ_mean(data.dff_pva(ds_jump_frames(jump)-5:ds_jump_frames(jump))),circ_mean(data.dff_pva(ds_jump_frames(jump)+5:ds_jump_frames(jump)+10)))));
   subplot(2,2,4)
   plot(0,0)
   ylim([-1 1]);
   xlim([-1 1]);
   text(-0.8,0,['Bump jump mag = ',num2str(round(bump_jump_mag(jump),2))]);
   set(gca,'xticklabel',{[]});
   set(gca,'yticklabel',{[]});
   
   saveas(gcf,[path,'plots\AJ',num2str(jump),'.png']);
   
end

%% Reproduce the Kim et al plot

for jump = 1:length(ds_jump_frames)
   
   time_zero = data.time(ds_jump_frames(jump));
   time = data.time-time_zero;
   
   figure('Position',[100 100 1200 500]),
   ax(1) = subplot(8,1,1);
   imagesc(bump_mag(:,ds_jump_frames(jump)-25:ds_jump_frames(jump)+26))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump magnitude');
   
   ax(2) = subplot(8,1,[2 8]);
   time_to_plot = time(ds_jump_frames(jump)-25:ds_jump_frames(jump)+26);
   bar_to_plot = bar_position(ds_jump_frames(jump)-25:ds_jump_frames(jump)+26);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'color',[1 0.2 0.2],'linewidth',2)
   hold on
   phase_to_plot = phase(ds_jump_frames(jump)-25:ds_jump_frames(jump)+26);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot);   
   plot(x_out_time,bump_pos_to_plot,'color',[0 0.6 0.3],'linewidth',2)
   xline(time(ds_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(ds_jump_frames(jump)-25) time(ds_jump_frames(jump)+26)]);
   ylabel('Heading (deg)');
   xlabel('Time (sec)');
   legend('Bar position','Bump estimate');
   
   %cb = colorbar(ax(2),'orientation','horizontal','Location','NorthOutside');
   
   saveas(gcf,[path,'plots\bm_evolution',num2str(jump),'.png'])
end

%% Repeat with bump width

for jump = 1:length(ds_jump_frames)
   
   time_zero = data.time(ds_jump_frames(jump));
   time = data.time-time_zero;
   
   figure('Position',[100 100 1200 500]),
   ax(1) = subplot(8,1,1);
   imagesc(half_max_width(:,ds_jump_frames(jump)-25:ds_jump_frames(jump)+26))
   colormap(flipud(gray))
   set(gca,'xtick',[])
   set(gca,'ytick',[])
   title('Bump width');
   
   ax(2) = subplot(8,1,[2 8]);
   time_to_plot = time(ds_jump_frames(jump)-25:ds_jump_frames(jump)+26);
   bar_to_plot = bar_position(ds_jump_frames(jump)-25:ds_jump_frames(jump)+26);
   [x_out_time,bar_pos_to_plot] = removeWrappedLines(time_to_plot,bar_to_plot);
   plot(x_out_time,bar_pos_to_plot,'color',[1 0.2 0.2],'linewidth',2)
   hold on
   phase_to_plot = phase(ds_jump_frames(jump)-25:ds_jump_frames(jump)+26);
   [x_out_time,bump_pos_to_plot] = removeWrappedLines(time_to_plot,phase_to_plot);   
   plot(x_out_time,bump_pos_to_plot,'color',[0 0.6 0.3],'linewidth',2)
   xline(time(ds_jump_frames(jump)),'k','linestyle','--','linewidth',2)
   ylim([-180 180]);
   xlim([time(ds_jump_frames(jump)-25) time(ds_jump_frames(jump)+26)]);
   ylabel('Heading (deg)');
   xlabel('Time (sec)');
   legend('Bar position','Bump estimate');
   
   %cb = colorbar(ax(2),'orientation','horizontal','Location','NorthOutside');
   
   saveas(gcf,[path,'plots\bw_evolution',num2str(jump),'.png'])
end

%% Plot the bump parameters vs the bump jump magnitude

figure,
subplot(1,2,1)
plot(mean_bm_pre_jump,bump_jump_mag,'ro')
bm_corr = corrcoef(mean_bm_pre_jump,bump_jump_mag);
xlabel('Bump magnitude (max-min)');
ylabel('Bump jump magnitude (deg)');
text(1,120,['corr = ',num2str(round(bm_corr(2,1),2))])

subplot(1,2,2)
plot(mean_bw_pre_jump,bump_jump_mag,'ko')
bw_corr = corrcoef(mean_bw_pre_jump,bump_jump_mag);
xlabel('Bump width (EB wedges)');
text(3,120,['corr = ',num2str(round(bw_corr(2,1),2))])

saveas(gcf,[path,'plots\BumpJumpCorrelation.png']);


%% Look at bar position around the jumps

figure('Position',[100 100 1200 800]),
for jump = 1:length(ds_jump_frames)
    subplot(length(ds_jump_frames),1,jump)
    time_around_jump = data.time(ds_jump_frames(jump)-25:ds_jump_frames(jump)+25) - data.time(ds_jump_frames(jump));
    bar_position_AJ(jump,:) = bar_position(ds_jump_frames(jump)-25:ds_jump_frames(jump)+25);
    plot(time_around_jump,bar_position_AJ(jump,:),'-o','linewidth',2,'color',[0.4 0.8 0.3])
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-3 3]);
    if jump ~= (length(ds_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Bar position around the jumps');

saveas(gcf,[path,'plots\AJ_bar_pos.png']);


%% Close up around the jumps

for jump = 1:length(ds_jump_frames)
    figure,
    time_around_jump_cu = data.time(ds_jump_frames(jump)-5:ds_jump_frames(jump)+5) - data.time(ds_jump_frames(jump));
    bar_position_AJ_cu(jump,:) = bar_position(ds_jump_frames(jump)-5:ds_jump_frames(jump)+5);
    plot(time_around_jump_cu,bar_position_AJ_cu(jump,:),'-o','linewidth',2,'color',[0.4 0.8 0.3])
    %determine jump size
    bar_jump_size(jump) = wrapTo180(bar_position_AJ_cu(jump,7)-bar_position_AJ_cu(jump,6));
    text(-0.4,100,['Bar jump size = ',num2str(bar_jump_size(jump))],'fontsize',14);
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-1 1]);
    title(['Bar position, Jump #',num2str(jump)]);
    saveas(gcf,[path,'plots\AJ_bar_pos_closeup_',num2str(jump),'.png']);
end


%% Plot the bar jump sizes

figure,
plot(bar_jump_size,'-ro')
ylim([-180 180]);
xlabel('Jump #');
ylabel('Bar jump size');

saveas(gcf,[path,'plots\bar_jump_size.png']);

%% Plot the offsets around the jumps

offset_AJ = [];
for jump = 1:length(ds_jump_frames)
    offset_AJ(jump,:) = offset(ds_jump_frames(jump)-25:ds_jump_frames(jump)+25); 
end

figure('Position',[100 100 1200 800]),
for jump = 1:length(ds_jump_frames)
    time_around_jump = data.time(ds_jump_frames(jump)-25:ds_jump_frames(jump)+25) - data.time(ds_jump_frames(jump));
    subplot(length(ds_jump_frames),1,jump)
    plot(time_around_jump,offset_AJ(jump,:),'-ko','linewidth',2)
    hold on
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-3 3]);
    if jump ~= (length(ds_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Offset around the jumps');

saveas(gcf,[path,'plots\AJ_offset.png']);

%% Shift all the offsets to start at 0 for comparison

shifted_offset_AJ = wrapTo180(offset_AJ-offset_AJ(:,1));

figure('Position',[100 100 1200 800]),
for jump = 1:length(ds_jump_frames)
    subplot(length(ds_jump_frames),1,jump)
    plot(time_around_jump,offset_AJ(jump,:),'-ko','linewidth',2)
    hold on
    plot(time_around_jump,shifted_offset_AJ(jump,:),'-bo','linewidth',2)
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-3 3]);
    if jump ~= (length(ds_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
    legend('offset','shifted offset');
end
suptitle('Shifted offset around the jumps');

saveas(gcf,[path,'plots\AJ_offset_r_and_s.png']);


%% Plot a heatmap of the offset around the jump to see the change dynamics

figure,
imagesc(shifted_offset_AJ)
colormap(flipud(gray))
hold on
colorbar
xline(26.5,'linewidth',2,'color','r');
title('Offset around the jumps');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(time_around_jump), max(time_around_jump), numel(xt));                    
set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
xlabel('Time around the jumps (sec)');
ylabel('Jump #');

saveas(gcf,[path,'plots\AJ_shifted_offset_heatmap.png']);

%% Sort by the value of the mean BM pre-jump

gray_wrap = vertcat(gray,flipud(gray));

%get indexes for sorted BM
[sorted_BM,I] = sort(mean_bm_pre_jump);

figure('Position',[100 100 1000 800]),
imagesc(shifted_offset_AJ(I,:))
colorbar
%colormap(flipud(gray))
colormap(gray_wrap);
hold on
xline(26.5,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(time_around_jump), max(time_around_jump), numel(xt));                    
set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Sorted offset around the jumps');
ylabel('Jumps number, sorted by BM');

saveas(gcf,[path,'plots\AJ_shifted_offset_heatmap_sorted.png']);

%% Shift the offset right before the jump

shifted_offset_AJ_short = wrapTo180(offset_AJ-offset_AJ(:,26));

figure('Position',[100 100 1000 800]),
imagesc(shifted_offset_AJ_short(I,:))
colorbar
%colormap(flipud(gray))
colormap(gray_wrap);
hold on
xline(26.5,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(time_around_jump), max(time_around_jump), numel(xt));                    
set(gca,'XTick',xt-2.5, 'XTickLabel',round(xtlbl,2))  
title('Sorted offset around the jumps'); %says sorted but doesn't seem sorted?
ylabel('Sorted jump #');



%% change in offset in lower vs higher bump magnitude pre-jump

%midle jump
mid_value = median(mean_bm_pre_jump);

%divide into low and high values
low_pre_jump_bm = mean_bm_pre_jump < median(mean_bm_pre_jump);
high_pre_jump_bm = mean_bm_pre_jump >= median(mean_bm_pre_jump);

%compute change in offset
one_sec_offset_change = abs(shifted_offset_AJ_short(:,36));
one_sec_offset_change_l = mean(one_sec_offset_change(low_pre_jump_bm));
one_sec_offset_change_h = mean(one_sec_offset_change(high_pre_jump_bm));

figure,
plot([one_sec_offset_change_l,one_sec_offset_change_h],'-ko')
xlim([0 3]);
ylim([0 180]);
ylabel('abs(Offset one sec after jump - offset right before jump)');
xticks([1 2])
xticklabels({'Low pre-jump BM','High pre-jump BM'})

saveas(gcf,[path,'plots\one_sec_offset_change_by_pre_jump_bm.png']);


%Combine the data in table
category = [repelem(1,sum(low_pre_jump_bm)),repelem(2,sum(high_pre_jump_bm))];
all_pre_offset_change = [one_sec_offset_change(low_pre_jump_bm);one_sec_offset_change(high_pre_jump_bm)];

figure,
boxplot(all_pre_offset_change,category,'color','k')
hold on
scatter(category,all_pre_offset_change,[],[.5 .5 .5],'filled','jitter', 'on', 'jitterAmount', 0.05)
ylabel('Offset one sec after jump - offset right before jump');
set(findobj(gca,'type','line'),'linew',2)
xticks([1 2])
xticklabels({'Low pre-jump BM','High pre-jump BM'})
ylim([0 180]);

saveas(gcf,[path,'plots\one_sec_offset_change_by_pre_jump_bm_boxplots.png']);


%% change in offset in lower vs higher bump magnitude pre-jump

%midle jump
mid_value = median(mean_bw_pre_jump);

%divide into low and high values
low_pre_jump_bw = mean_bw_pre_jump < median(mean_bw_pre_jump);
high_pre_jump_bw = mean_bw_pre_jump >= median(mean_bw_pre_jump);

%compute change in offset
one_sec_offset_change = abs(shifted_offset_AJ_short(:,36));
one_sec_offset_change_l = mean(one_sec_offset_change(low_pre_jump_bw));
one_sec_offset_change_h = mean(one_sec_offset_change(high_pre_jump_bw));

figure,
plot([one_sec_offset_change_l,one_sec_offset_change_h],'-ko')
xlim([0 3]);
ylim([0 180]);
ylabel('abs(Offset one sec after jump - offset right before jump)');
xticks([1 2])
xticklabels({'Low pre-jump BW','High pre-jump BW'})

saveas(gcf,[path,'plots\one_sec_offset_change_by_pre_jump_bw.png']);


%Combine the data in table
category = [repelem(1,sum(low_pre_jump_bw)),repelem(2,sum(high_pre_jump_bw))];
all_pre_offset_change = [one_sec_offset_change(low_pre_jump_bw);one_sec_offset_change(high_pre_jump_bw)];

figure,
boxplot(all_pre_offset_change,category,'color','k')
hold on
scatter(category,all_pre_offset_change,[],[.5 .5 .5],'filled','jitter', 'on', 'jitterAmount', 0.05)
ylabel('Offset one sec after jump - offset right before jump');
set(findobj(gca,'type','line'),'linew',2)
xticks([1 2])
xticklabels({'Low pre-jump BW','High pre-jump BW'})
ylim([0 180]);

saveas(gcf,[path,'plots\one_sec_offset_change_by_pre_jump_bw_boxplots.png']);

%% Look at the change in offset

change_in_offset = wrapTo180(rad2deg(circ_dist(deg2rad(offset(2:end)),deg2rad(offset(1:end-1))))); 

figure('Position',[100 100 1200 800]),
for jump = 1:length(ds_jump_frames)
    time_around_jump = data.time(ds_jump_frames(jump)-45:ds_jump_frames(jump)+45) - data.time(ds_jump_frames(jump));
    subplot(length(ds_jump_frames),1,jump)
    plot(time_around_jump,change_in_offset(ds_jump_frames(jump)-45:ds_jump_frames(jump)+45),'k','linewidth',2)
    hold on
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-5 5]);
    if jump ~= (length(ds_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Change in offset around the jumps');

saveas(gcf,[path,'plots\AJ_offset_change.png']);

%% Plot it as a heatmap

figure,
for jump = 1:length(ds_jump_frames)
    change_in_offset_AJ(jump,:) = abs(change_in_offset(ds_jump_frames(jump)-45:ds_jump_frames(jump)+45)); 
end
imagesc(change_in_offset_AJ)
colormap(flipud(gray))
hold on
colorbar
xline(45,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(time_around_jump), max(time_around_jump), numel(xt));                    
set(gca,'XTick',xt-4.5, 'XTickLabel',round(xtlbl,2))  
title('Change in offset around the jumps');
xlabel('Time around the jumps (sec)');
ylabel('Jump #');

saveas(gcf,[path,'plots\AJ_offset_change_heatmap.png']);


%% Sort by the value of the mean BM pre-jump


figure,
imagesc(change_in_offset_AJ(I,:))
colormap(flipud(gray))
%colormap(gray_wrap)
hold on
colorbar
xline(45.5,'linewidth',2,'color','r');
xt = get(gca, 'XTick');                                            
xtlbl = linspace(min(time_around_jump), max(time_around_jump), numel(xt));                    
set(gca,'XTick',xt-4.5, 'XTickLabel',round(xtlbl,2))  
title('Change in offset around the jumps');
xlabel('Time around the jumps (sec)');
ylabel('Sorted jumps');

saveas(gcf,[path,'plots\AJ_offset_change_heatmap_sorted.png']);

%% Get offset around the jumps for more frames, to compute the time it takes to go back to its pre-jump value

longer_offset_AJ = [];
for jump = 1:length(ds_jump_frames)
    longer_offset_AJ(jump,:) = offset(ds_jump_frames(jump)-250:ds_jump_frames(jump)+250); 
end

figure('Position',[100 100 1200 800]),
for jump = 1:length(ds_jump_frames)
    time_around_jump = data.time(ds_jump_frames(jump)-250:ds_jump_frames(jump)+250) - data.time(ds_jump_frames(jump));
    subplot(length(ds_jump_frames),1,jump)
    plot(time_around_jump,longer_offset_AJ(jump,:),'-k','linewidth',2)
    hold on
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-26 26]);
    if jump ~= (length(ds_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Longer offset around the jumps');

saveas(gcf,[path,'plots\AJ_longer_offset.png']);

%% Sort by BM value pre jump

figure('Position',[100 100 1200 800]),
for jump = 1:length(ds_jump_frames)
    subplot(length(ds_jump_frames),1,jump)
    plot(time_around_jump,longer_offset_AJ(I(jump),:),'-k','linewidth',2)
    hold on
    xline(0,'linewidth',2,'color','r');
    ylim([-180 180]);
    xlim([-26 26]);
    if jump ~= (length(ds_jump_frames))
        set(gca,'xticklabel',{[]})
    else
        xlabel('Time around the jump (sec)');
    end
end
suptitle('Longer offset around the jumps,sorted');

saveas(gcf,[path,'plots\AJ_longer_offset_sorted.png']);

%% Find the number of frames it takes for the offset to return to its pre-jump value

close all;
%Define offset value pre-jump as the circular mean in the 5 sec preceding the jump
pre_jump_offset = rad2deg(circ_mean(deg2rad(longer_offset_AJ(:,206:251)),[],2));

%find frames after the jump that are within 10 degrees of that offset
return_frames = {};
for jump = 1:length(ds_jump_frames)
    offset_difference = rad2deg(circ_dist(deg2rad(longer_offset_AJ(jump,252:end)),deg2rad(pre_jump_offset(jump))));
    return_frames{jump} = find(abs(offset_difference)<=6);
    if ~isempty(return_frames{jump})
        first_return_frame(jump) = return_frames{jump}(1);
    else
        first_return_frame(jump) = NaN;
    end
end

%Plot
for jump = 1:length(ds_jump_frames)
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
    saveas(gcf,[path,'plots\offset_return_sortedj_',num2str(jump),'.png']);
end

%Plot the offset return frames sorted by BM
figure,
subplot(1,2,1)
plot(first_return_frame(I),'o')
ylabel('Frame that returned to pre-jump offset');
xlabel('Jumps sorted by BM');

%% Plot the offset return frames sorted by BM

[sorted_BW,Iw] = sort(mean_bw_pre_jump);
subplot(1,2,2)
plot(first_return_frame(Iw),'o')
ylabel('Frame that returned to pre-jump offset');
xlabel('Jumps sorted by BW');
saveas(gcf,[path,'plots\offset_return_sorted_frames.png']);

%% Bump magnitude around the jumps

for jump = 1:length(ds_jump_frames)
    
   figure('Position',[100 100 1200 800]),
     
   subplot(3,1,1)
   imagesc(data.dff_matrix(:,ds_jump_frames(jump)-25:ds_jump_frames(jump)+25))
   colormap(flipud(gray))
   hold on
   xline(26,'linewidth',2,'color','r')
   title('EPG activity');
     
   subplot(3,1,2)
   plot(offset(ds_jump_frames(jump)-25:ds_jump_frames(jump)+25))
   xline(26,'linewidth',2,'color','r')
   title(['Jump #',num2str(jump)]); 
   xlim([1 51]);
   ylim([-180 180]);
   title('Offset (deg)');
   
   subplot(3,1,3)
   %compute the mean bump parameters right before the jump
   mean_bm_pre_jump(jump) = mean(bump_mag(ds_jump_frames(jump)-25:ds_jump_frames(jump)));
   mean_bm_post_jump(jump) = mean(bump_mag(ds_jump_frames(jump)+1:ds_jump_frames(jump)+26));
   bm_around_jump{jump} = [mean_bm_pre_jump(jump),mean_bm_post_jump(jump)];
   plot(bm_around_jump{jump},'-ko')
   xlim([0 3]);
   set(gca,'xticklabel',{[]})
   ylim([min(bm_around_jump{jump})-0.5, max(bm_around_jump{jump})+0.5]);  
   title('Bump magnitude (max-min)');
   
   saveas(gcf,[path,'plots\AJ_bump_magnitude_',num2str(jump),'.png']);
   
end

%% Combining all the data on bump magnitude around the jumps

figure,
bm_aj = cell2mat(bm_around_jump);
bm_aj = reshape(bm_aj,2,length(bm_around_jump));
plot(bm_aj,'-o','color',[.5 .5 .5])
hold on
plot(mean(bm_aj,2),'-ko','linewidth',2)
xlim([0 3]);
xticks([1 2])
xticklabels({'pre-jump','post-jump'})
ylabel('Mean bump magnitude (max-min)');
saveas(gcf,[path,'plots\AJ_bump_magnitude.png']);

%% Bump width aroung the jumps

for jump = 1:length(ds_jump_frames)
    
   figure('Position',[100 100 1200 800]),
     
   subplot(3,1,1)
   imagesc(data.dff_matrix(:,ds_jump_frames(jump)-25:ds_jump_frames(jump)+25))
   colormap(flipud(gray))
   hold on
   xline(26,'linewidth',2,'color','r')
   title('EPG activity');
     
   subplot(3,1,2)
   plot(offset(ds_jump_frames(jump)-25:ds_jump_frames(jump)+25))
   xline(26,'linewidth',2,'color','r')
   title(['Jump #',num2str(jump)]); 
   xlim([1 51]);
   ylim([-180 180]);
   title('Offset (deg)');
   
   subplot(3,1,3)
   %compute the mean bump parameters right before the jump
   mean_bw_pre_jump(jump) = mean(half_max_width(ds_jump_frames(jump)-25:ds_jump_frames(jump)));
   mean_bw_post_jump(jump) = mean(half_max_width(ds_jump_frames(jump)+1:ds_jump_frames(jump)+26));
   bw_around_jump{jump} = [mean_bw_pre_jump(jump),mean_bw_post_jump(jump)];
   plot(bm_around_jump{jump},'-ko')
   xlim([0 3]);
   set(gca,'xticklabel',{[]})
   %ylim([min(bw_around_jump{jump})-0.5, max(bw_around_jump{jump})+0.5]);  
   title('Bump width');
   
   saveas(gcf,[path,'plots\AJ_bump_width_',num2str(jump),'.png']);
   
end


%% Combining all the data on bump width around the jumps

figure,
bw_aj = cell2mat(bw_around_jump);
bw_aj = reshape(bw_aj,2,length(bw_around_jump));
plot(bw_aj,'-o','color',[.5 .5 .5])
hold on
plot(mean(bw_aj,2),'-ko','linewidth',2)
xlim([0 3]);
xticks([1 2])
xticklabels({'pre-jump','post-jump'})
ylabel('Mean bump width');
saveas(gcf,[path,'plots\AJ_bump_width.png']);

%% Save data

save([path,'\offset_change.mat'],'change_in_offset_AJ','time_around_jump','mean_bm_pre_jump','mean_bw_pre_jump','offset_AJ','longer_offset_AJ','first_return_frame');
save([path,'\bump_parameter_change.mat'],'bm_aj','bw_aj')

close all; clear all;