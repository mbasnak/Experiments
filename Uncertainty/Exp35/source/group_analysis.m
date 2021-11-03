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

offset_var = [];
offset_var_r = [];
offset_mean = [];
heading_mean = [];
BM_mean = [];
BW_mean = [];
thresh_BM_mean = [];
thresh_BW_mean = [];
total_mvt = [];

for fly = 1:length(data)
    offset_var = [offset_var;data(fly).offset_var];
    offset_var_r = [offset_var_r;data(fly).offset_var_r];
    offset_mean = [offset_mean;data(fly).offset_mean];
    heading_mean = [heading_mean;data(fly).heading_mean];
    BM_mean = [BM_mean;data(fly).allBM_thresh'];
    BW_mean = [BW_mean;data(fly).allBW_thresh'];  
    thresh_BM_mean = [thresh_BM_mean;data(fly).allBM_thresh_final'];
    thresh_BW_mean = [thresh_BW_mean;data(fly).allBW_thresh_final'];  
    total_mvt = [total_mvt;data(fly).all_total_mvt_thresh'];
end

%% Plot offset variation per block

set(0,'DefaultTextInterpreter','none')

figure,
plot(offset_var_r','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(offset_var_r),'-ko','linewidth',2)
xlabel('Block #');
ylabel('circ_std(offset)');

%note that the durations of the blocks are different, so this will affect
%circ_std.

saveas(gcf,[path,'\groupPlots\offset_var.png']);

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

figure,
plot(abs(diff_offset_mean_deg'),'-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(abs(diff_offset_mean_deg)),'-ko','linewidth',2)
xlabel('Block #');
ylabel([{'Absolute difference in offset mean'};{'with respect to previous block'}]);

saveas(gcf,[path,'\groupPlots\diff_offset.png']);

%% Comparison in offset between block 1 and 3

figure('Position',[100 100 1000 800]),
subplot(1,4,[1 3]),
plot(rad2deg([offset_mean(:,1),offset_mean(:,3)]'),'-o','color',[.5 .5 .5])
xlim([0 3]);
ylim([-180 180]);
xticks([1 2])
xticklabels({'block 1', 'block 3'});
title('Mean offset');

diff_offset = rad2deg(abs(circ_dist(offset_mean(:,1),offset_mean(:,3))));
subplot(1,4,4),
boxplot(diff_offset)
title('Abs difference');

saveas(gcf,[path,'\groupPlots\diff_offset_1_and_3.png']);

%% Comparison in offset betwwn block 2 and 4

figure('Position',[100 100 1000 800]),
subplot(1,4,[1 3]),
plot(rad2deg([offset_mean(:,2),offset_mean(:,4)]'),'-o','color',[.5 .5 .5])
xlim([0 3]);
ylim([-180 180]);
xticks([1 2])
xticklabels({'block 2', 'block 4'});
title('Mean offset');

diff_offset = rad2deg(abs(circ_dist(offset_mean(:,2),offset_mean(:,4))));
subplot(1,4,4),
boxplot(diff_offset)
title('Abs difference');

saveas(gcf,[path,'\groupPlots\diff_offset_2_and_4.png']);

%% Comparison between blocks 3 and 5

figure('Position',[100 100 1000 800]),
subplot(1,4,[1 3]),
plot(rad2deg([offset_mean(:,3),offset_mean(:,5)]'),'-o','color',[.5 .5 .5])
xlim([0 3]);
ylim([-180 180]);
xticks([1 2])
xticklabels({'block 3', 'block 5'});
title('Mean offset');

diff_offset = rad2deg(abs(circ_dist(offset_mean(:,3),offset_mean(:,5))));
subplot(1,4,4),
boxplot(diff_offset)
title('Abs difference');

saveas(gcf,[path,'\groupPlots\diff_offset_3_and_5.png']);

%% Comparison between blocks 1 and 4 (same cue)

figure('Position',[100 100 1000 800]),
subplot(1,4,[1 3]),
plot(rad2deg([offset_mean(:,1),offset_mean(:,4)]'),'-o','color',[.5 .5 .5])
xlim([0 3]);
ylim([-180 180]);
xticks([1 2])
xticklabels({'block 1', 'block 4'});
title('Mean offset');

diff_offset = rad2deg(abs(circ_dist(offset_mean(:,1),offset_mean(:,4))));
subplot(1,4,4),
boxplot(diff_offset)
title('Abs difference');

saveas(gcf,[path,'\groupPlots\diff_offset_1_and_4.png']);


%% Comparison between blocks 2 and 5 (same cue)

figure('Position',[100 100 1000 800]),
subplot(1,4,[1 3]),
plot(rad2deg([offset_mean(:,2),offset_mean(:,5)]'),'-o','color',[.5 .5 .5])
xlim([0 3]);
ylim([-180 180]);
xticks([1 2])
xticklabels({'block 2', 'block 5'});
title('Mean offset');

diff_offset = rad2deg(abs(circ_dist(offset_mean(:,2),offset_mean(:,5))));
subplot(1,4,4),
boxplot(diff_offset)
title('Abs difference');

saveas(gcf,[path,'\groupPlots\diff_offset_2_and_5.png']);

%% Plot bump parameters

figure,
subplot(1,2,1)
plot(BM_mean','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(BM_mean),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump magnitude');

subplot(1,2,2)
plot(BW_mean','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(BW_mean),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump width');

saveas(gcf,[path,'\groupPlots\bump_parameters_evolution.png']);

%% Plot thresholded bump parameters

figure,
subplot(1,2,1)
plot(thresh_BM_mean','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmean(thresh_BM_mean),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump magnitude');

subplot(1,2,2)
plot(thresh_BW_mean','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmean(thresh_BW_mean),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump width');

saveas(gcf,[path,'\groupPlots\thresh_bump_parameters_evolution.png']);

%% Plot zscored-bump parameters

zscored_BM = zscore(BM_mean,[],2);
zscored_BW = zscore(BW_mean,[],2);

figure,
subplot(1,2,1)
plot(zscored_BM','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(zscored_BM),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean zscored bump magnitude');

subplot(1,2,2)
plot(zscored_BW','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(zscored_BW),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean zscored bump width');

saveas(gcf,[path,'\groupPlots\zscored_bump_parameters_evolution.png']);

%% Plot mean total movement to check that bump parameter evolution is not driven by that


figure,
plot(total_mvt','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(total_mvt),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Total movement (deg/s)');

saveas(gcf,[path,'\groupPlots\total_mvt.png']);

%% Plot change in heading mean per block

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
xlim([0 6]);
hold on
plot(mean(abs(diff_heading_mean_deg)),'-ko','linewidth',2)
xlabel('Block #');
ylabel([{'Absolute difference in heading mean'};{'with respect to previous block'}]);

saveas(gcf,[path,'\groupPlots\diff_heading.png']);

%%
close all; clear all;