%Group analysis code for the cue combination experiment

clear all; close all;

%% Load data

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\high_reliability';

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\data.mat']);
   end
end

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%% Combine all data

offset_var_r = [];
offset_var_pre_panels = {};
offset_var_pre_wind = {};
offset_var_combination = {};
offset_var_post_panels = {};
offset_var_post_wind = {};
offset_mean = [];
offset_mean_final = [];
heading_mean = [];
thresh_BM_mean = [];
thresh_BW_mean = [];
final_BM_mean = [];
final_BW_mean = [];
total_mvt = [];
mean_bout_BM_pre_panels = {};
mean_bout_BW_pre_panels = {};
mean_bout_BM_pre_wind = {};
mean_bout_BW_pre_wind = {};
mean_bout_BM_cue_combination = {};
mean_bout_BW_cue_combination = {};
mean_bout_BM_post_panels = {};
mean_bout_BW_post_panels = {};
mean_bout_BM_post_wind = {};
mean_bout_BW_post_wind = {};
initial_cue_diff = [];
bar_offset_diff = [];
wind_offset_diff = [];
initial_cue_diff_last_part = [];
bar_offset_diff_last_part = [];
wind_offset_diff_last_part = [];
configuration = [];
p_stop_pre_panels = {};
p_stop_pre_wind = {};
p_stop_cue_combination = {};
p_stop_post_panels = {};
p_stop_post_wind = {};
p_stopped = [];
allSummaryData = array2table(zeros(0,4),'VariableNames', {'BumpPos','gof','TotalMvt','BlockType'});
flyNumber = [];

for fly = 1:length(data)
    
    offset_var_r = [offset_var_r;data(fly).offset_var_r];
    offset_var_pre_panels{fly} = data(fly).offset_var_pre_panels;
    offset_var_pre_wind{fly} = data(fly).offset_var_pre_wind;
    offset_var_cue_combination{fly} = data(fly).offset_var_cue_combination;
    offset_var_post_panels{fly} = data(fly).offset_var_post_panels;
    offset_var_post_wind{fly} = data(fly).offset_var_post_wind;
    offset_mean = [offset_mean;data(fly).offset_mean];
    offset_mean_final = [offset_mean_final;data(fly).offset_mean_final];
    heading_mean = [heading_mean;data(fly).heading_mean];
    thresh_BM_mean = [thresh_BM_mean;data(fly).allBM_thresh'];
    thresh_BW_mean = [thresh_BW_mean;data(fly).allBW_thresh'];  
    final_BM_mean = [final_BM_mean;data(fly).allBM_thresh_final'];
    final_BW_mean = [final_BW_mean;data(fly).allBW_thresh_final'];  
    total_mvt = [total_mvt;data(fly).all_total_mvt_thresh'];
    mean_bout_BM_pre_panels{fly} = data(fly).mean_bout_BM_pre_panels;
    mean_bout_BW_pre_panels{fly} = data(fly).mean_bout_BW_pre_panels;
    mean_bout_BM_pre_wind{fly} = data(fly).mean_bout_BM_pre_wind;
    mean_bout_BW_pre_wind{fly} = data(fly).mean_bout_BW_pre_wind;
    mean_bout_BM_cue_combination{fly} = data(fly).mean_bout_BM_cue_combination;
    mean_bout_BW_cue_combination{fly} = data(fly).mean_bout_BW_cue_combination;
    mean_bout_BM_post_panels{fly} = data(fly).mean_bout_BM_post_panels;
    mean_bout_BW_post_panels{fly} = data(fly).mean_bout_BW_post_panels;
    mean_bout_BM_post_wind{fly} = data(fly).mean_bout_BM_post_wind;
    mean_bout_BW_post_wind{fly} = data(fly).mean_bout_BW_post_wind;
    initial_cue_diff(fly) = data(fly).initial_cue_diff;
    bar_offset_diff(fly) = data(fly).bar_offset_diff;
    wind_offset_diff(fly) = data(fly).wind_offset_diff;
    initial_cue_diff_last_part(fly) = data(fly).initial_cue_diff_last_part;
    bar_offset_diff_last_part(fly) = data(fly).bar_offset_diff_last_part;
    wind_offset_diff_last_part(fly) = data(fly).wind_offset_diff_last_part;
    configuration(fly) = data(fly).configuration;
    p_stop_pre_panels{fly} = data(fly).p_stop_pre_panels;
    p_stop_pre_wind{fly} = data(fly).p_stop_pre_wind;
    p_stop_cue_combination{fly} = data(fly).p_stop_cue_combination;
    p_stop_post_panels{fly} = data(fly).p_stop_post_panels;
    p_stop_post_wind{fly} = data(fly).p_stop_post_wind;
    p_stopped = [p_stopped;data(fly).p_stopped];  
    allSummaryData = [allSummaryData;data(fly).summaryData];
    flyNumber = [flyNumber,repelem(fly,length(data(fly).summaryData.BumpPos))];
end
allSummaryData = addvars(allSummaryData,flyNumber');
allSummaryData.Properties.VariableNames{'Var5'} = 'Fly';


%% Effects of bump position on plasticity
%To determine whether plasticity is influenced by where the bump lives during the training and testing periods:
%divide test period into two categories: timepoints when the bump is in a place that was 'well-sampled'
%during the cue combination period (i.e., that the bump spent more than x% of the session within an x window size of 
%that value) and timepoints when that location was poorly sampled.

%For each fly
for fly = 1:length(data)
    %Extract bump pos from cue combination period
    cc_bump_pos =  allSummaryData.BumpPos(allSummaryData.Fly == fly & allSummaryData.BlockType == 3 & allSummaryData.TotalMvt > 25 & allSummaryData.gof >= 0.5);
    %Focus on the test periods
    for cue = [4,5]
        %Define vector for bump positions
        sampling_of_bump_pos{fly,cue-3} = zeros(1,sum(allSummaryData.Fly == fly & allSummaryData.BlockType == cue & allSummaryData.TotalMvt > 25 & allSummaryData.gof >= 0.5));
        all_bump_pos = allSummaryData.BumpPos(allSummaryData.Fly == fly & allSummaryData.BlockType == cue & allSummaryData.TotalMvt > 25 & allSummaryData.gof >= 0.5);    
        %For each timepoint
        for timepoint = 1:length(all_bump_pos)
            %Collect the bump position
            bump_pos = all_bump_pos(timepoint);
            %Focus on the cue combination bout and
            %determine if that bump position -+ 15 deg was sampled for more than
            %10% of the frames in which the fly was walking
            if (sum(abs(rad2deg(circ_dist(bump_pos,cc_bump_pos)))< 20) / length(cc_bump_pos)) > 0.1
                %If it was, assign a value of 1
                sampling_of_bump_pos{fly,cue-3}(timepoint) = 1;     
            end
        end
    end
end


%Compute for each fly a ratio of 'well sampled':'poor sampled' bump pos
for fly = 1:length(data)
    for cue = 1:2
        goodness_of_sampling(fly,cue) = sum(sampling_of_bump_pos{fly,cue})/length(sampling_of_bump_pos{fly,cue});
    end
end

%Correlate the ratio to offset variability for those bouts
figure,
plot(goodness_of_sampling(:,1),offset_var_r(:,4),'ko')
hold on
plot(goodness_of_sampling(:,2),offset_var_r(:,5),'ko')
xlabel({'Goodness of bump position sampling','%bump pos well sampled'});
ylabel('Offset variability in test block');

saveas(gcf,[path,'\groupPlots\gos_bump_pos_vs_offset_var.png']);


%% Plot offset variation per block

set(0,'DefaultTextInterpreter','none')

figure('Position',[100 100 800 800]),
plot(offset_var_r','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmedian(offset_var_r),'-ko','linewidth',3,'MarkerFaceColor','k')
ylabel('Offset variability (circular std)','fontweight','bold','fontsize',16);
xticks([1 2 3 4 5]);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16,'FontWeight','bold')

%we downsampled the blocks to compute the offset variability so that we
%used the number of frames of the shorter block

saveas(gcf,[path,'\groupPlots\offset_var.png']);


%% Offset variability in 60 sec bouts across sessions

figure('Position',[100 100 800 1000]),
for fly = 1:length(data)
    
    %Plot initial block
    subplot(5,1,1)
    max_length_of_first_bout = max([cellfun(@length,offset_var_pre_panels),cellfun(@length,offset_var_pre_wind)]);
    all_offset_var_1 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_offset_var_1(fly,1:length(offset_var_pre_panels{fly})) = offset_var_pre_panels{fly};
        else
            all_offset_var_1(fly,1:length(offset_var_pre_wind{fly})) = offset_var_pre_wind{fly};            
        end
    end
    plot(all_offset_var_1','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_offset_var_1),'-ko','linewidth',2)
    ylim([0 2]);
    title('First block'); ylabel('Bump magnitude');
    
    subplot(5,1,2)
    all_offset_var_2 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_offset_var_2(fly,1:length(offset_var_pre_panels{fly})) = offset_var_pre_panels{fly};
        else
            all_offset_var_2(fly,1:length(offset_var_pre_wind{fly})) = offset_var_pre_wind{fly};            
        end
    end
    plot(all_offset_var_2','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_offset_var_2),'-ko','linewidth',2)
    ylim([0 2]);
    title('Second block'); ylabel('Bump magnitude');
    
    subplot(5,1,3)
    max_length_of_cc_bouts = max(cellfun(@length,offset_var_cue_combination));
    all_offset_var_cc = NaN(length(data),max_length_of_cc_bouts);
    for fly = 1:length(data)
        all_offset_var_cc(fly,1:length(offset_var_cue_combination{fly})) = offset_var_cue_combination{fly};
    end
    plot(all_offset_var_cc','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_offset_var_cc),'-ko','linewidth',2)
    ylim([0 2]); xlim([1 13]);
    title('Third block'); ylabel('Bump magnitude');
    
    subplot(5,1,4)
    max_length_of_fourth_bout = max([cellfun(@length,offset_var_post_panels),cellfun(@length,offset_var_post_wind)]);
    all_offset_var_4 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_offset_var_4(fly,1:length(offset_var_post_panels{fly})) = offset_var_post_panels{fly};
        else
            all_offset_var_4(fly,1:length(offset_var_post_wind{fly})) = offset_var_post_wind{fly};            
        end
    end
    plot(all_offset_var_4','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_offset_var_4),'-ko','linewidth',2)
    ylim([0 2]); xlim([1 5]);
    title('Fourth block'); ylabel('Bump magnitude');
    
    subplot(5,1,5)
    all_offset_var_5 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_offset_var_5(fly,1:length(offset_var_post_panels{fly})) = offset_var_post_panels{fly};
        else
            all_offset_var_5(fly,1:length(offset_var_post_wind{fly})) = offset_var_post_wind{fly};            
        end
    end
    plot(all_offset_var_5','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_offset_var_5),'-ko','linewidth',2)
    ylim([0 2]); xlim([1 5]);
    title('Fifth block'); ylabel('Bump magnitude');
    
end

saveas(gcf,[path,'\groupPlots\offset_var_across_session_number.png']);

%% Plot change in offset mean per block

diff_offset_mean = [];
diff_offet_mean_deg = [];

for fly = 1:length(data)
   for block = 2:5
      diff_offset_mean(fly,block-1) = circ_dist(offset_mean(fly,block),offset_mean(fly,block-1));
   end
end
diff_offset_mean = [repelem(nan,fly,1),diff_offset_mean];
diff_offset_mean_deg = rad2deg(diff_offset_mean);

figure('Position',[100 100 800 800]),
plot(abs(diff_offset_mean_deg'),'-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmedian(abs(diff_offset_mean_deg)),'-ko','linewidth',2)
xlabel('Block #','fontweight','bold');
xticks([1 2 3 4 5]);
ylabel([{'Absolute difference in offset mean'};{'with respect to previous block (deg)'}],'fontweight','bold');

saveas(gcf,[path,'\groupPlots\diff_offset.png']);

%% Focusing on the last 120 sec of each bout to determine the offset mean

diff_offset_mean_final = [];
diff_offet_mean_deg_final = [];

for fly = 1:length(data)
   for block = 2:5
      diff_offset_mean_final(fly,block-1) = circ_dist(offset_mean_final(fly,block),offset_mean_final(fly,block-1));
   end
end
diff_offset_mean_final = [repelem(nan,fly,1),diff_offset_mean_final];
diff_offset_mean_deg_final = rad2deg(diff_offset_mean_final);

figure('Position',[100 100 800 800]),
plot(abs(diff_offset_mean_deg_final'),'-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmedian(abs(diff_offset_mean_deg_final)),'-ko','linewidth',2)
xlabel('Block #','fontweight','bold');
xticks([1 2 3 4 5]);
ylabel([{'Absolute difference in offset mean'};{'with respect to previous block (deg)'}],'fontweight','bold');

saveas(gcf,[path,'\groupPlots\diff_offset_final.png']);

%% Matrix of offset differences

for fly = 1:length(data)
    my_vct =  [offset_mean(fly,:)];
    L = numel(my_vct);
    tmp_vct = [my_vct(2:end) my_vct(end)];
    x = circ_dist(my_vct,tmp_vct);
    sol = zeros(L-1,L);
    for k=2:L
        tmp_vct = [my_vct(k:end) my_vct((L-k+2):end)];
        offset_differences(k-1,:,fly) = round(abs(rad2deg(circ_dist(my_vct,tmp_vct))));
    end
end

mean_offset_differences = mean(offset_differences,3);

figure,
h = heatmap(mean_offset_differences);
h.Colormap = gray;
h.ColorbarVisible = 'off';
ax = gca;
ax.YData = ['5'; '4'; '3'; '2'];
h.Title = 'Mean offset differences';
h.XLabel = 'Block number';
h.YLabel = 'Block number';

saveas(gcf,[path,'\groupPlots\all_offset_differences.png']);

%% Matrix of offset differences focusing on the last 120 sec

for fly = 1:length(data)
    my_vct =  [offset_mean_final(fly,:)];
    L = numel(my_vct);
    tmp_vct = [my_vct(2:end) my_vct(end)];
    x = circ_dist(my_vct,tmp_vct);
    sol = zeros(L-1,L);
    for k=2:L
        tmp_vct = [my_vct(k:end) my_vct((L-k+2):end)];
        offset_differences_final(k-1,:,fly) = round(abs(rad2deg(circ_dist(my_vct,tmp_vct))));
    end
end

mean_offset_differences_final = mean(offset_differences_final,3);

figure,
h = heatmap(mean_offset_differences_final);
h.Colormap = gray;
h.ColorbarVisible = 'off';
ax = gca;
ax.YData = ['5'; '4'; '3'; '2'];
h.Title = 'Mean offset differences';
h.XLabel = 'Block number';
h.YLabel = 'Block number';

saveas(gcf,[path,'\groupPlots\all_offset_differences_final.png']);

%% Plot thresholded bump parameters

figure,
subplot(1,2,1)
plot(thresh_BM_mean','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmean(thresh_BM_mean),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold'); xticks([1 2 3 4 5]);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean bump magnitude','fontweight','bold');

subplot(1,2,2)
plot(thresh_BW_mean','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmean(thresh_BW_mean),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold'); xticks([1 2 3 4 5]);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean bump width','fontweight','bold');

saveas(gcf,[path,'\groupPlots\thresh_bump_parameters_evolution.png']);

%% Plot thresholded bump parameters focusing on last part of bouts

figure,
subplot(1,2,1)
plot(final_BM_mean','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(final_BM_mean),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold');
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean bump magnitude','fontweight','bold');

subplot(1,2,2)
plot(final_BW_mean','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(final_BW_mean),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold');
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean bump width','fontweight','bold');

saveas(gcf,[path,'\groupPlots\final_bump_parameters_evolution.png']);

%% Plot zscored-bump parameters

zscored_BM = zscore(thresh_BM_mean,[],2);
zscored_BW = zscore(thresh_BW_mean,[],2);

figure,
subplot(1,2,1)
plot(zscored_BM','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_BM),'-ko','linewidth',3,'MarkerFaceColor','k')
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean zscore','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
title('Bump magnitude','fontsize',16);

subplot(1,2,2)
plot(zscored_BW','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_BW),'-ko','linewidth',3,'MarkerFaceColor','k')
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean zscore','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')
title('Bump width','fontsize',16);

saveas(gcf,[path,'\groupPlots\zscored_bump_parameters_evolution.png']);

%% Repeat zscoring data for last 120 sec of each bout

zscored_BM_final = zscore(final_BM_mean,[],2);
zscored_BW_final = zscore(final_BW_mean,[],2);

figure,
subplot(1,2,1)
plot(zscored_BM_final','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_BM_final),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold');
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean zscored bump magnitude','fontweight','bold');

subplot(1,2,2)
plot(zscored_BW_final','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_BW_final),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold');
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Mean zscored bump width','fontweight','bold');

saveas(gcf,[path,'\groupPlots\zscored_bump_parameters_evolution_final.png']);

%% Plot mean total movement to check that bump parameter evolution is not driven by that

figure,
plot(total_mvt','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(total_mvt),'-ko','linewidth',3,'MarkerFaceColor','k')
xlabel('Block type','fontweight','bold','fontsize',14);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Total movement (deg/s)','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')

saveas(gcf,[path,'\groupPlots\total_mvt.png']);

%% Repeat zscoring

zscored_mvt = zscore(total_mvt,[],2);

figure,
plot(zscored_mvt','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_mvt),'-ko','linewidth',3,'MarkerFaceColor','k')
xlabel('Block type','fontweight','bold','fontsize',14);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Zscored total movement','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')

saveas(gcf,[path,'\groupPlots\zscored_total_mvt.png']);


%% Plot change in thresholded heading mean per block

diff_heading_mean = [];
diff_heading_mean_deg = [];

for fly = 1:length(data)
   for block = 2:5
      diff_heading_mean(fly,block-1) = circ_dist(heading_mean(fly,block),heading_mean(fly,block-1));
   end
end
diff_heading_mean = [repelem(nan,fly,1),diff_heading_mean];
diff_heading_mean_deg = rad2deg(diff_heading_mean);

figure,
plot(abs(diff_heading_mean_deg'),'-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(mean(abs(diff_heading_mean_deg)),'-ko','linewidth',2)
xlabel('Block #','fontweight','bold');
ylabel([{'Absolute difference in heading mean'};{'with respect to previous block (deg)'}],'fontweight','bold');

saveas(gcf,[path,'\groupPlots\diff_heading.png']);

%% Look at bump magnitude evolution in 60 sec bouts in the different blocks

figure('Position',[100 100 800 1000]),

for fly = 1:length(data)
    
    subplot(5,1,1)
    plot(mean_bout_BM_pre_panels{fly},'-o','color',[.5 .5 .5])
    hold on
    title('Initial panels');
    ylabel('Bump magnitude');
    ylim([0 3]);
    
    subplot(5,1,2)
    plot(mean_bout_BM_pre_wind{fly},'-o','color',[.5 .5 .5])
    hold on
    title('Initial wind');
    ylabel('Bump magnitude');
    ylim([0 3]);
    
    subplot(5,1,3)
    plot(mean_bout_BM_cue_combination{fly},'-o','color',[.5 .5 .5])
    hold on
    title('Cue combination');
    ylabel('Bump magnitude');
    xlim([1 20]);
    ylim([0 3]);
    
    subplot(5,1,4)
    plot(mean_bout_BM_post_panels{fly},'-o','color',[.5 .5 .5])
    hold on
    title('Final panels');
    ylabel('Bump magnitude');
    ylim([0 3]);
    
    subplot(5,1,5)
    plot(mean_bout_BM_post_wind{fly},'-o','color',[.5 .5 .5])
    hold on
    title('Final wind');
    ylabel('Bump magnitude');
    ylim([0 3]);
    
end

%Add means
subplot(5,1,1)
max_length_of_ip_bouts = max(cellfun(@length,mean_bout_BM_pre_panels));
all_mean_bout_BM_ip = NaN(length(data),max_length_of_ip_bouts);
for fly = 1:length(data)
    all_mean_bout_BM_ip(fly,1:length(mean_bout_BM_pre_panels{fly})) = mean_bout_BM_pre_panels{fly};
end
plot(nanmean(all_mean_bout_BM_ip),'-ko','linewidth',2)

subplot(5,1,2)
max_length_of_iw_bouts = max(cellfun(@length,mean_bout_BM_pre_wind));
all_mean_bout_BM_iw = NaN(length(data),max_length_of_iw_bouts);
for fly = 1:length(data)
    all_mean_bout_BM_iw(fly,1:length(mean_bout_BM_pre_wind{fly})) = mean_bout_BM_pre_wind{fly};
end
plot(nanmean(all_mean_bout_BM_iw),'-ko','linewidth',2)

subplot(5,1,3)
max_length_of_cc_bouts = max(cellfun(@length,mean_bout_BM_cue_combination));
all_mean_bout_BM_cc = NaN(length(data),max_length_of_cc_bouts);
for fly = 1:length(data)
    all_mean_bout_BM_cc(fly,1:length(mean_bout_BM_cue_combination{fly})) = mean_bout_BM_cue_combination{fly};
end
plot(nanmean(all_mean_bout_BM_cc),'-ko','linewidth',2)

subplot(5,1,4)
max_length_of_fp_bouts = max(cellfun(@length,mean_bout_BM_post_panels));
all_mean_bout_BM_fp = NaN(length(data),max_length_of_fp_bouts);
for fly = 1:length(data)
    all_mean_bout_BM_fp(fly,1:length(mean_bout_BM_post_panels{fly})) = mean_bout_BM_post_panels{fly};
end
plot(nanmean(all_mean_bout_BM_fp),'-ko','linewidth',2)

subplot(5,1,5)
max_length_of_fw_bouts = max(cellfun(@length,mean_bout_BM_post_wind));
all_mean_bout_BM_fw = NaN(length(data),max_length_of_fw_bouts);
for fly = 1:length(data)
    all_mean_bout_BM_fw(fly,1:length(mean_bout_BM_post_wind{fly})) = mean_bout_BM_post_wind{fly};
end
plot(nanmean(all_mean_bout_BM_fw),'-ko','linewidth',2)

saveas(gcf,[path,'\groupPlots\bump_mag_across_session_type.png']);

%% Repeat dividing by block number instead of type

figure('Position',[100 100 800 1000]),
for fly = 1:length(data)
    
    %Plot initial block
    subplot(5,1,1)
    max_length_of_first_bout = max([cellfun(@length,mean_bout_BM_pre_panels),cellfun(@length,mean_bout_BM_pre_wind)]);
    all_mean_bout_BM_1 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_mean_bout_BM_1(fly,1:length(mean_bout_BM_pre_panels{fly})) = mean_bout_BM_pre_panels{fly};
        else
            all_mean_bout_BM_1(fly,1:length(mean_bout_BM_pre_wind{fly})) = mean_bout_BM_pre_wind{fly};            
        end
            
    end
    plot(all_mean_bout_BM_1','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BM_1),'-ko','linewidth',2)
    ylim([0 3]);
    title('First block'); ylabel('Bump magnitude');
    
    subplot(5,1,2)
    all_mean_bout_BM_2 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_mean_bout_BM_2(fly,1:length(mean_bout_BM_pre_panels{fly})) = mean_bout_BM_pre_panels{fly};
        else
            all_mean_bout_BM_2(fly,1:length(mean_bout_BM_pre_wind{fly})) = mean_bout_BM_pre_wind{fly};            
        end
    end
    plot(all_mean_bout_BM_2','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BM_2),'-ko','linewidth',2)
    ylim([0 3]);
    title('Second block'); ylabel('Bump magnitude');
    
    subplot(5,1,3)
    max_length_of_cc_bouts = max(cellfun(@length,mean_bout_BM_cue_combination));
    all_mean_bout_BM_cc = NaN(length(data),max_length_of_cc_bouts);
    for fly = 1:length(data)
        all_mean_bout_BM_cc(fly,1:length(mean_bout_BM_cue_combination{fly})) = mean_bout_BM_cue_combination{fly};
    end
    plot(all_mean_bout_BM_cc','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BM_cc),'-ko','linewidth',2)
    ylim([0 3]); xlim([1 13]);
    title('Third block'); ylabel('Bump magnitude');
    
    subplot(5,1,4)
    max_length_of_fourth_bout = max([cellfun(@length,mean_bout_BM_post_panels),cellfun(@length,mean_bout_BM_post_wind)]);
    all_mean_bout_BM_4 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_mean_bout_BM_4(fly,1:length(mean_bout_BM_post_panels{fly})) = mean_bout_BM_post_panels{fly};
        else
            all_mean_bout_BM_4(fly,1:length(mean_bout_BM_post_wind{fly})) = mean_bout_BM_post_wind{fly};            
        end
    end
    plot(all_mean_bout_BM_4','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BM_4),'-ko','linewidth',2)
    ylim([0 3]); xlim([1 5]);
    title('Fourth block'); ylabel('Bump magnitude');
    
    subplot(5,1,5)
    all_mean_bout_BM_5 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_mean_bout_BM_5(fly,1:length(mean_bout_BM_post_panels{fly})) = mean_bout_BM_post_panels{fly};
        else
            all_mean_bout_BM_5(fly,1:length(mean_bout_BM_post_wind{fly})) = mean_bout_BM_post_wind{fly};            
        end
    end
    plot(all_mean_bout_BM_5','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BM_5),'-ko','linewidth',2)
    ylim([0 3]);  xlim([1 5]);
    title('Fifth block'); ylabel('Bump magnitude');
    
end

saveas(gcf,[path,'\groupPlots\bump_mag_across_session_number.png']);

%% Repeat for bump width

figure('Position',[100 100 800 1000]),
for fly = 1:length(data)
    
    %Plot initial block
    subplot(5,1,1)
    all_mean_bout_BW_1 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_mean_bout_BW_1(fly,1:length(mean_bout_BW_pre_panels{fly})) = mean_bout_BW_pre_panels{fly};
        else
            all_mean_bout_BW_1(fly,1:length(mean_bout_BW_pre_wind{fly})) = mean_bout_BW_pre_wind{fly};            
        end
    end
    plot(all_mean_bout_BW_1','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BW_1),'-ko','linewidth',2)
    ylim([0 3]);
    title('First block'); ylabel('Bump width');
    
    subplot(5,1,2)
    all_mean_bout_BW_2 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_mean_bout_BW_2(fly,1:length(mean_bout_BW_pre_panels{fly})) = mean_bout_BW_pre_panels{fly};
        else
            all_mean_bout_BW_2(fly,1:length(mean_bout_BW_pre_wind{fly})) = mean_bout_BW_pre_wind{fly};            
        end
    end
    plot(all_mean_bout_BW_2','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BW_2),'-ko','linewidth',2)
    ylim([0 3]);
    title('Second block'); ylabel('Bump width');
    
    subplot(5,1,3)
    all_mean_bout_BW_cc = NaN(length(data),max_length_of_cc_bouts);
    for fly = 1:length(data)
        all_mean_bout_BW_cc(fly,1:length(mean_bout_BW_cue_combination{fly})) = mean_bout_BW_cue_combination{fly};
    end
    plot(all_mean_bout_BW_cc','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BW_cc),'-ko','linewidth',2)
    ylim([0 3]); xlim([1 13]);
    title('Third block'); ylabel('Bump width');
    
    subplot(5,1,4)
    all_mean_bout_BW_4 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_mean_bout_BW_4(fly,1:length(mean_bout_BW_post_panels{fly})) = mean_bout_BW_post_panels{fly};
        else
            all_mean_bout_BW_4(fly,1:length(mean_bout_BW_post_wind{fly})) = mean_bout_BW_post_wind{fly};            
        end
    end
    plot(all_mean_bout_BW_4','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BW_4),'-ko','linewidth',2)
    ylim([0 3]); xlim([1 5]);
    title('Fourth block'); ylabel('Bump width');
    
    subplot(5,1,5)
    all_mean_bout_BW_5 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_mean_bout_BW_5(fly,1:length(mean_bout_BW_post_panels{fly})) = mean_bout_BW_post_panels{fly};
        else
            all_mean_bout_BW_5(fly,1:length(mean_bout_BW_post_wind{fly})) = mean_bout_BW_post_wind{fly};            
        end
    end
    plot(all_mean_bout_BW_5','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_BW_5),'-ko','linewidth',2)
    ylim([0 3]);  xlim([1 5]);
    title('Fifth block'); ylabel('Bump width');
    
end

saveas(gcf,[path,'\groupPlots\bump_width_across_session_number.png']);

%% Single cue plasticity based on initial offset differences

figure('Position',[100 100 1000 1000]),
subplot(2,1,1)
plot(abs(initial_cue_diff(configuration == 1)),abs(bar_offset_diff(1,configuration == 1)),'ko','MarkerFaceColor','k', 'MarkerSize',8)
hold on
plot(abs(initial_cue_diff(configuration == 2)),abs(bar_offset_diff(1,configuration == 2)),'ro','MarkerFaceColor','r', 'MarkerSize',8)
ylabel('Bar offset difference (post vs pre cue combination)','fontweight','bold');
legend('Bar first','Wind first','location','best');

subplot(2,1,2)
plot(abs(initial_cue_diff(configuration == 1)),abs(wind_offset_diff(1,configuration == 1)),'ko','MarkerFaceColor','k', 'MarkerSize',8)
hold on
plot(abs(initial_cue_diff(configuration == 2)),abs(wind_offset_diff(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
xlabel('Abs of initial bar and wind offset difference','fontweight','bold');
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
legend('Bar was first stim','Wind was first stim','location','best');

saveas(gcf,[path,'\groupPlots\single_cue_diff_vs_initial_diff.png']);

%% Repeat using only the last portion of each block

figure('Position',[100 100 1000 1000]),
subplot(2,1,1)
plot(abs(initial_cue_diff_last_part(configuration == 1)),abs(bar_offset_diff_last_part(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(initial_cue_diff_last_part(configuration == 2)),abs(bar_offset_diff_last_part(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Bar offset difference (post vs pre cue combination)','fontweight','bold');
legend('Bar first','Wind first','location','best');

subplot(2,1,2)
plot(abs(initial_cue_diff_last_part(configuration == 1)),abs(wind_offset_diff_last_part(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(initial_cue_diff_last_part(configuration == 2)),abs(wind_offset_diff_last_part(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
xlabel('Abs of initial bar and wind offset difference','fontweight','bold');
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
legend('Bar first','Wind first','location','best');

saveas(gcf,[path,'\groupPlots\single_cue_diff_vs_initial_diff_last_part.png']);


%% Plot probability of flies to be stopped across sessions

figure('Position',[100 100 800 1000]),
for fly = 1:length(data)
    
    %Plot initial block
    subplot(5,1,1)
    max_length_of_first_bout = max([cellfun(@length,p_stop_pre_panels),cellfun(@length,p_stop_pre_wind)]);
    all_mean_bout_p_stop_1 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_mean_bout_p_stop_1(fly,1:length(p_stop_pre_panels{fly})) = p_stop_pre_panels{fly};
        else
            all_mean_bout_p_stop_1(fly,1:length(p_stop_pre_wind{fly})) = p_stop_pre_wind{fly};            
        end
    end
    plot(all_mean_bout_p_stop_1','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_p_stop_1),'-ko','linewidth',2)
    ylim([0 1]);
    title('First block'); ylabel('P(stopped)');
    
    subplot(5,1,2)
    all_mean_bout_p_stop_2 = NaN(length(data),max_length_of_first_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_mean_bout_p_stop_2(fly,1:length(p_stop_pre_panels{fly})) = p_stop_pre_panels{fly};
        else
            all_mean_bout_p_stop_2(fly,1:length(p_stop_pre_wind{fly})) = p_stop_pre_wind{fly};            
        end
    end
    plot(all_mean_bout_p_stop_2','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_p_stop_2),'-ko','linewidth',2)
    ylim([0 1]);
    title('First block'); ylabel('P(stopped)');
    
    subplot(5,1,3)
    max_length_of_cc_bouts = max(cellfun(@length,p_stop_cue_combination));
    all_mean_bout_p_stop_cc = NaN(length(data),max_length_of_cc_bouts);
    for fly = 1:length(data)
        all_mean_bout_p_stop_cc(fly,1:length(p_stop_cue_combination{fly})) = p_stop_cue_combination{fly};
    end
    plot(all_mean_bout_p_stop_cc','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_p_stop_cc),'-ko','linewidth',2)
    ylim([0 1]); xlim([1 13]);
    title('Third block'); ylabel('P(stopped)');
    
    subplot(5,1,4)
    max_length_of_fourth_bout = max([cellfun(@length,p_stop_post_panels),cellfun(@length,p_stop_post_wind)]);
    all_mean_bout_p_stop_4 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 1
            all_mean_bout_p_stop_4(fly,1:length(p_stop_post_panels{fly})) = p_stop_post_panels{fly};
        else
            all_mean_bout_p_stop_4(fly,1:length(p_stop_post_wind{fly})) = p_stop_post_wind{fly};            
        end
    end
    plot(all_mean_bout_p_stop_4','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_p_stop_4),'-ko','linewidth',2)
    ylim([0 1]); xlim([1 5]);
    title('Fourth block'); ylabel('P(stopped)');
    
    subplot(5,1,5)
    all_mean_bout_p_stop_5 = NaN(length(data),max_length_of_fourth_bout);
    for fly = 1:length(data)
        if configuration(fly) == 2
            all_mean_bout_p_stop_5(fly,1:length(p_stop_post_panels{fly})) = p_stop_post_panels{fly};
        else
            all_mean_bout_p_stop_5(fly,1:length(p_stop_post_wind{fly})) = p_stop_post_wind{fly};            
        end
    end
    plot(all_mean_bout_p_stop_5','-o','color',[.5 .5 .5])
    hold on
    plot(nanmean(all_mean_bout_p_stop_5),'-ko','linewidth',2)
    ylim([0 1]); xlim([1 5]);
    title('Fifth block'); ylabel('P(stopped)');
    
end
saveas(gcf,[path,'\groupPlots\p_stop_across_sessions.png']);

%% Mean probability of stopping per session

figure('Position',[100 100 1200 800]),
plot(p_stopped','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(p_stopped),'-ko','linewidth',2)
xlabel('Block type','fontweight','bold');
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('P(stopped)','fontweight','bold');

%% Repeat with zscored data

zscored_pstopped = zscore(p_stopped,[],2);
figure,
plot(zscored_pstopped','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_pstopped),'-ko','linewidth',3,'MarkerFaceColor','k')
xlabel('Block type','fontweight','bold','fontsize',14);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Zscored P(stopped)','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')

saveas(gcf,[path,'\groupPlots\mean_pstopped.png']);
%the flies appear to stop more in the block following the cue combination
%period.

%% Combine mvt and pstop plots for presentation display

figure,
subplot(1,2,1)
plot(zscored_mvt','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_mvt),'-ko','linewidth',3,'MarkerFaceColor','k')
xlabel('Block type','fontweight','bold','fontsize',14);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Zscored total movement','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')

subplot(1,2,2)
plot(zscored_pstopped','-o','color',[.5 .5 .5]);
xlim([0 6]); xticks([1 2 3 4 5]);
hold on
plot(nanmean(zscored_pstopped),'-ko','linewidth',3,'MarkerFaceColor','k')
xlabel('Block type','fontweight','bold','fontsize',14);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
ylabel('Zscored P(stopped)','fontweight','bold','fontsize',14);
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',14,'FontWeight','bold')

%% Changes in bump parameters as a function of the inital offset difference

change_BM = [];
change_BW = [];
for fly = 1:length(data)
    if configuration(fly) == 1
        change_BM(fly) = mean_bout_BM_cue_combination{fly}(1)-mean_bout_BM_pre_wind{fly}(end);
        change_BW(fly) = mean_bout_BW_cue_combination{fly}(1)-mean_bout_BW_pre_wind{fly}(end);
    else
        change_BM(fly) = mean_bout_BM_cue_combination{fly}(1)-mean_bout_BM_pre_panels{fly}(end);
        change_BW(fly) = mean_bout_BW_cue_combination{fly}(1)-mean_bout_BW_pre_panels{fly}(end);
    end   
end

figure('Position',[100 100 1000 800]),
subplot(2,1,1)
plot(abs(initial_cue_diff),change_BM,'o')
xlabel('Abs of initial cue difference','fontweight','bold');
ylabel({'Change in bump magnitude';'(beginning cue combination - end single cue)'},'fontweight','bold');

subplot(2,1,2)
plot(abs(initial_cue_diff),change_BW,'o')
xlabel('Abs of initial cue difference','fontweight','bold');
ylabel({'Change in bump width';'(beginning cue combination - end single cue)'},'fontweight','bold');

saveas(gcf,[path,'\groupPlots\change_bump_par_vs_initial_cue_diff.png']);

%there are no clear trends

%% Repeat using the last part of the session to determine the offset

figure('Position',[100 100 1000 800]),
subplot(2,1,1)
plot(abs(initial_cue_diff_last_part),change_BM,'o')
xlabel('Abs of initial cue difference','fontweight','bold');
ylabel({'Change in bump magnitude';'(beginning cue combination - end single cue)'},'fontweight','bold');

subplot(2,1,2)
plot(abs(initial_cue_diff_last_part),change_BW,'o')
xlabel('Abs of initial cue difference','fontweight','bold');
ylabel({'Change in bump width';'(beginning cue combination - end single cue)'},'fontweight','bold');

saveas(gcf,[path,'\groupPlots\change_bump_par_vs_initial_cue_diff_last_part.png']);

%again, there are no clear trends

%% Inverse analysis, looking at the change in post-pre offset for single cues as a function of the change in bump parameters

%For the bar
figure('Position',[100 100 1000 800]),
subplot(2,1,1)
plot(change_BM,abs(bar_offset_diff),'o')
xlabel({'Change in bump magnitude';'(beginning cue combination - end single cue)'},'fontweight','bold');
ylabel('Abs difference in bar offset (post - pre)');

subplot(2,1,2)
plot(change_BW,abs(bar_offset_diff),'o')
xlabel({'Change in bump width';'(beginning cue combination - end single cue)'},'fontweight','bold');
ylabel('Abs difference in bar offset (post - pre)');

saveas(gcf,[path,'\groupPlots\bar_plasticity_vs_change_bump_par.png']);


%For the wind
figure('Position',[100 100 1000 800]),
subplot(2,1,1)
plot(change_BM,abs(wind_offset_diff),'o')
xlabel({'Change in bump magnitude';'(beginning cue combination - end single cue)'},'fontweight','bold');
ylabel('Abs difference in wind offset (post - pre)');

subplot(2,1,2)
plot(change_BW,abs(wind_offset_diff),'o')
xlabel({'Change in bump width';'(beginning cue combination - end single cue)'},'fontweight','bold');
ylabel('Abs difference in wind offset (post - pre)');

saveas(gcf,[path,'\groupPlots\wind_plasticity_vs_change_bump_par.png']);


%% Repeat computing change in bump parameters as the difference in the mean bump parameters across blocks

change_mean_BM = [];
change_mean_BW = [];
for fly = 1:length(data)
        change_mean_BM(fly) = thresh_BM_mean(fly,3) - thresh_BM_mean(fly,2);
        change_mean_BW(fly) = thresh_BW_mean(fly,3) - thresh_BW_mean(fly,2);  
end

%For the bar
figure('Position',[100 100 1000 800]),
subplot(2,1,1)
plot(change_mean_BM,abs(bar_offset_diff),'o')
xlabel({'Change in mean bump magnitude'},'fontweight','bold');
ylabel('Abs difference in bar offset (post - pre)');

subplot(2,1,2)
plot(change_mean_BW,abs(bar_offset_diff),'o')
xlabel({'Change in mean bump width'},'fontweight','bold');
ylabel('Abs difference in bar offset (post - pre)');


%For the wind
figure('Position',[100 100 1000 800]),
subplot(2,1,1)
plot(change_mean_BM,abs(wind_offset_diff),'o')
xlabel({'Change in mean bump magnitude'},'fontweight','bold');
ylabel('Abs difference in wind offset (post - pre)');

subplot(2,1,2)
plot(change_mean_BW,abs(wind_offset_diff),'o')
xlabel({'Change in mean bump width'},'fontweight','bold');
ylabel('Abs difference in wind offset (post - pre)');


%% Cue combination and single cue difference vs initial cue difference

%Determine difference between cue combination and initial cue
offset_diff_cc_w = [];
offset_diff_cc_v = [];
for fly = 1:length(data)
    if configuration(fly) == 1
        offset_diff_cc_w(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,2)));
        offset_diff_cc_v(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,1)));
    else
        offset_diff_cc_w(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,1)));
        offset_diff_cc_v(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,2)));
    end
end

figure('Position',[100 100 1000 1000]),
subplot(2,1,1)
plot(abs(initial_cue_diff(configuration == 1)),abs(offset_diff_cc_v(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(initial_cue_diff(configuration == 2)),abs(offset_diff_cc_v(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Abs of cue combination and initial bar offset difference','fontweight','bold');
legend('Bar first','Wind first','location','best');

subplot(2,1,2)
plot(abs(initial_cue_diff(configuration == 1)),abs(offset_diff_cc_w(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(initial_cue_diff(configuration == 2)),abs(offset_diff_cc_w(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
xlabel('Abs of initial bar and wind offset difference','fontweight','bold');
ylabel('Abs of cue combination and initial wind offset difference','fontweight','bold');
legend('Bar first','Wind first','location','best');

saveas(gcf,[path,'\groupPlots\cue_comb_and_single_cue_diff_vs_initial_diff.png']);


%% Repeat analysis combining wind and bar data and color coding by last cue

for fly = 1:length(data)
    single_cue_cc_diff_1(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,1)));
    single_cue_cc_diff_2(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,2)));
end

figure('Position',[100 100 1000 600]),

plot(abs(initial_cue_diff),abs(single_cue_cc_diff_1),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(initial_cue_diff),abs(single_cue_cc_diff_2),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Abs of cue combination and initial single cue offset difference','fontweight','bold');
xlabel('Abs of initial bar and wind offset difference','fontweight','bold');
legend('Single cue was first cue','Single cue was second cue','location','best');

saveas(gcf,[path,'\groupPlots\cue_comb_and_single_cue_diff_vs_initial_diff_combined.png']);



%% Repeat focusing on the last part for the offsets

for fly = 1:length(data)
    single_cue_cc_diff_1_final(fly) = rad2deg(circ_dist(offset_mean_final(fly,3),offset_mean_final(fly,1)));
    single_cue_cc_diff_2_final(fly) = rad2deg(circ_dist(offset_mean_final(fly,3),offset_mean_final(fly,2)));
end

figure('Position',[100 100 1000 600]),

plot(abs(initial_cue_diff_last_part),abs(single_cue_cc_diff_1_final),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(initial_cue_diff_last_part),abs(single_cue_cc_diff_2_final),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Abs of cue combination and initial single cue offset difference','fontweight','bold');
xlabel('Abs of initial bar and wind offset difference','fontweight','bold');
legend('Single cue was first cue','Single cue was second cue','location','best');

saveas(gcf,[path,'\groupPlots\cue_comb_and_single_cue_diff_vs_initial_diff_last_part.png']);


%% Single cue plasticity for visual cues as a function of cue combination offset

figure('Position',[100 100 1500 600]),
subplot(1,3,1)
plot(abs(offset_diff_cc_v(1,configuration == 1)),abs(bar_offset_diff(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(offset_diff_cc_v(1,configuration == 2)),abs(bar_offset_diff(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Bar offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial bar offset difference','fontweight','bold');
legend('Bar first','Wind first','location','best');
correlation = corrcoef(abs(offset_diff_cc_v),abs(bar_offset_diff));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,2)
plot(abs(offset_diff_cc_v(1,configuration == 1)),abs(bar_offset_diff(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
xlabel('Cue combination and initial bar offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_v(1,configuration == 1)),abs(bar_offset_diff(1,configuration == 1)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,3)
plot(abs(offset_diff_cc_v(1,configuration == 2)),abs(bar_offset_diff(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
xlabel('Cue combination and initial bar offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_v(1,configuration == 2)),abs(bar_offset_diff(1,configuration == 2)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

saveas(gcf,[path,'\groupPlots\single_cue_v_diff_vs_cue_combination_plasticity.png']);

%% Repeat for wind cues

figure('Position',[100 100 1500 600]),
subplot(1,3,1)
plot(abs(offset_diff_cc_w(1,configuration == 1)),abs(wind_offset_diff(1,configuration == 1)),'ko', 'MarkerFaceColor','k', 'MarkerSize',8)
hold on
plot(abs(offset_diff_cc_w(1,configuration == 2)),abs(wind_offset_diff(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
legend('Bar first','Wind first','location','best');
correlation = corrcoef(abs(offset_diff_cc_w),abs(wind_offset_diff));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,2)
plot(abs(offset_diff_cc_w(1,configuration == 1)),abs(wind_offset_diff(1,configuration == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_w(1,configuration == 1)),abs(wind_offset_diff(1,configuration == 1)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,3)
plot(abs(offset_diff_cc_w(1,configuration == 2)),abs(wind_offset_diff(1,configuration == 2)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_w(1,configuration == 2)),abs(wind_offset_diff(1,configuration == 2)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

saveas(gcf,[path,'\groupPlots\single_cue_w_diff_vs_cue_combination_plasticity.png']);

%% Repeat combining both cues

for fly = 1:length(data)
    single_cue_diff_1(fly) = rad2deg(circ_dist(offset_mean(fly,4),offset_mean(fly,1)));
    single_cue_diff_2(fly) = rad2deg(circ_dist(offset_mean(fly,5),offset_mean(fly,2)));
end

figure('Position',[100 100 1500 600]),
subplot(1,3,1)
plot(abs(single_cue_cc_diff_1),abs(single_cue_diff_1),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(single_cue_cc_diff_2),abs(single_cue_diff_2),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Single cue offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
legend('Single cue went first','Single cue went second','location','best');
[correlation, pval] = corrcoef([abs(single_cue_cc_diff_1),abs(single_cue_cc_diff_2)],[abs(single_cue_diff_1),abs(single_cue_diff_2)]);
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,2)
plot(abs(single_cue_cc_diff_1),abs(single_cue_diff_1),'ko', 'MarkerFaceColor','k','MarkerSize',8)
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
[correlation, pval] = corrcoef(abs(single_cue_cc_diff_1),abs(single_cue_diff_1));
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,3)
plot(abs(single_cue_cc_diff_2),abs(single_cue_diff_2),'ro', 'MarkerFaceColor','r','MarkerSize',8)
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
[correlation, pval] = corrcoef(abs(single_cue_cc_diff_2),abs(single_cue_diff_2));
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);

saveas(gcf,[path,'\groupPlots\single_cue_diff_vs_cue_combination_plasticity.png']);

%% Color-code plasticity by offset variability in the single-cue-phase

figure('Position',[100 100 1200 800]),
subplot(1,2,1)
scatter(single_cue_cc_diff_1,single_cue_diff_1,58,offset_var_r(:,4),'filled','MarkerEdgeColor','k')
colormap(gray); colorbar;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
[correlation, pval] = corrcoef(single_cue_cc_diff_1,single_cue_diff_1);
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([-180 180]); ylim([-180 180]);

subplot(1,2,2)
scatter(single_cue_cc_diff_2,single_cue_diff_2,58,offset_var_r(:,5),'filled','MarkerEdgeColor','k')
colormap(gray); colorbar;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
[correlation, pval] = corrcoef(single_cue_cc_diff_2,single_cue_diff_2);
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([-180 180]); ylim([-180 180]);

saveas(gcf,[path,'\groupPlots\single_cue_diff_vs_cue_combination_plasticity_colored_by_off_var.png']);


%% Color-code plasticity by bump magnitude in the cue combination period

figure('Position',[100 100 1200 800]),
subplot(1,2,1)
scatter(single_cue_cc_diff_1,single_cue_diff_1,58,thresh_BM_mean(:,3),'filled','MarkerEdgeColor','k')
colormap(gray); colorbar;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
[correlation, pval] = corrcoef(single_cue_cc_diff_1,single_cue_diff_1);
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([-180 180]); ylim([-180 180]);

subplot(1,2,2)
scatter(single_cue_cc_diff_2,single_cue_diff_2,58,thresh_BM_mean(:,5),'filled','MarkerEdgeColor','k')
colormap(gray); colorbar;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold');
[correlation, pval] = corrcoef(single_cue_cc_diff_2,single_cue_diff_2);
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([-180 180]); ylim([-180 180]);

saveas(gcf,[path,'\groupPlots\single_cue_diff_vs_cue_combination_plasticity_colored_by_BM.png']);

%% Color-code by the change in bump mag between blocks 3 and 2

figure('Position',[100 100 1200 800]),
subplot(1,2,1)
scatter(abs(single_cue_cc_diff_1),abs(single_cue_diff_1),58,change_mean_BM,'filled','MarkerEdgeColor','k')
colormap(gray); h = colorbar;
ylabel(h,'bump magnitude change','FontSize',12,'Rotation',270);
h.Label.Position(1) = 3;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold','fontsize',12);
ylabel('Single cue plasticity (post - pre offset)','fontweight','bold','fontsize',12);
[correlation, pval] = corrcoef(abs(single_cue_cc_diff_1),abs(single_cue_diff_1));
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);

subplot(1,2,2)
scatter(abs(single_cue_cc_diff_2),abs(single_cue_diff_2),58,change_mean_BM,'filled','MarkerEdgeColor','k')
colormap(gray); h = colorbar;
ylabel(h,'bump magnitude change','FontSize',12,'Rotation',270);
h.Label.Position(1) = 3;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold','fontsize',12);
[correlation, pval] = corrcoef(abs(single_cue_cc_diff_2),abs(single_cue_diff_2));
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);


%% Repeat color-coding by bump width

figure('Position',[100 100 1200 800]),
subplot(1,2,1)
scatter(abs(single_cue_cc_diff_1),abs(single_cue_diff_1),58,change_mean_BW,'filled','MarkerEdgeColor','k')
colormap(gray); h = colorbar;
ylabel(h,'bump width change','FontSize',12,'Rotation',270);
h.Label.Position(1) = 3;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold','fontsize',12);
ylabel('Single cue plasticity (post - pre offset)','fontweight','bold','fontsize',12);
[correlation, pval] = corrcoef(abs(single_cue_cc_diff_1),abs(single_cue_diff_1));
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);

subplot(1,2,2)
scatter(abs(single_cue_cc_diff_2),abs(single_cue_diff_2),58,change_mean_BW,'filled','MarkerEdgeColor','k')
colormap(gray); h = colorbar;
ylabel(h,'bump width change','FontSize',12,'Rotation',270);
h.Label.Position(1) = 3;
xlabel('Cue combination and initial single cue offset difference','fontweight','bold','fontsize',12);
[correlation, pval] = corrcoef(abs(single_cue_cc_diff_2),abs(single_cue_diff_2));
title([['Corr = ',num2str(round(correlation(1,2),2))],['    Pval = ',num2str(round(pval(1,2),2))]]);
xlim([0 180]); ylim([0 180]);


%% Repeat previous analysis dividing into cases when flies were on average walking upwind (or -/+ 90 deg), vs downwind

%determine if flies are on average walking upwind or downwind during cue
%combination bout
upwind_walking = [];
for fly = 1:length(data)
        if (heading_mean(fly,3) < pi/2 & heading_mean(fly,3) > -pi/2)
            upwind_walking(fly) = 1;
        else
            upwind_walking(fly) = 0;
        end
end

%For visual data
figure('Position',[100 100 1500 600]),
subplot(1,3,1)
plot(abs(offset_diff_cc_v(1,upwind_walking == 1)),abs(bar_offset_diff(1,upwind_walking == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(offset_diff_cc_v(1,upwind_walking == 0)),abs(bar_offset_diff(1,upwind_walking == 0)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
legend('Upwind walking','non-upwind walking','location','best');
correlation = corrcoef(abs(offset_diff_cc_v),abs(bar_offset_diff));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,2)
plot(abs(offset_diff_cc_v(1,upwind_walking == 1)),abs(bar_offset_diff(1,upwind_walking == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_v(1,upwind_walking == 1)),abs(bar_offset_diff(1,upwind_walking == 1)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,3)
plot(abs(offset_diff_cc_v(1,upwind_walking == 0)),abs(bar_offset_diff(1,upwind_walking == 0)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_v(1,upwind_walking == 0)),abs(bar_offset_diff(1,upwind_walking == 0)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

saveas(gcf,[path,'\groupPlots\single_cue_v_diff_vs_cue_combination_plasticity_including_heading.png']);


%% For wind data

figure('Position',[100 100 1500 600]),
subplot(1,3,1)
plot(abs(offset_diff_cc_w(1,upwind_walking == 1)),abs(wind_offset_diff(1,upwind_walking == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(offset_diff_cc_w(1,upwind_walking == 0)),abs(wind_offset_diff(1,upwind_walking == 0)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
legend('Upwind walking','non-upwind walking','location','best');
correlation = corrcoef(abs(offset_diff_cc_w),abs(wind_offset_diff));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,2)
plot(abs(offset_diff_cc_w(1,upwind_walking == 1)),abs(wind_offset_diff(1,upwind_walking == 1)),'ko', 'MarkerFaceColor','k','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_w(1,upwind_walking == 1)),abs(wind_offset_diff(1,upwind_walking == 1)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

subplot(1,3,3)
plot(abs(offset_diff_cc_w(1,upwind_walking == 0)),abs(wind_offset_diff(1,upwind_walking == 0)),'ro', 'MarkerFaceColor','r','MarkerSize',8)
ylabel('Wind offset difference (post vs pre cue combination)','fontweight','bold');
xlabel('Cue combination and initial wind offset difference','fontweight','bold');
correlation = corrcoef(abs(offset_diff_cc_w(1,upwind_walking == 0)),abs(wind_offset_diff(1,upwind_walking == 0)));
title(['Corr = ',num2str(round(correlation(1,2),2))]);
xlim([0 180]); ylim([0 180]);

saveas(gcf,[path,'\groupPlots\single_cue_w_diff_vs_cue_combination_plasticity_including_heading.png']);

%% Dividing the data by heading during the test period

%Divide test period into 60 sec bouts

%Compute mean heading in that bout

%Determine if mean heading in that bout is upwind or downwind

%Determine offset plasticity in that bout

%Correlate offset plasticity to difference between cc and single cue, and
%divide by upwind or downwind walking


%% Offset variability in 60-sec bouts


%%
close all; clear all;