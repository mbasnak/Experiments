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
offset_mean = [];
heading_mean = [];
thresh_BM_mean = [];
thresh_BW_mean = [];
total_mvt = [];
initial_cue_diff = [];
bar_offset_diff = [];
wind_offset_diff = [];
configuration = [];
p_stop_pre_panels = {};
p_stop_pre_wind = {};
p_stop_cue_combination = {};
p_stop_post_panels = {};
p_stop_post_wind = {};
p_stopped = [];
allSummaryData = array2table(zeros(0,4),'VariableNames', {'BumpMag','BumpWidth','TotalMvt','BlockType'});
flyNumber = [];

for fly = 1:length(data)
    
    offset_var_r = [offset_var_r;data(fly).offset_var_r];
    offset_mean = [offset_mean;data(fly).offset_mean];
    heading_mean = [heading_mean;data(fly).heading_mean];
    thresh_BM_mean = [thresh_BM_mean;data(fly).allBM_thresh'];
    thresh_BW_mean = [thresh_BW_mean;data(fly).allBW_thresh'];  
    total_mvt = [total_mvt;data(fly).all_total_mvt_thresh'];
    initial_cue_diff(fly) = data(fly).initial_cue_diff;
    bar_offset_diff(fly) = data(fly).bar_offset_diff;
    wind_offset_diff(fly) = data(fly).wind_offset_diff;
    configuration(fly) = data(fly).configuration;
    p_stop_pre_panels{fly} = data(fly).p_stop_pre_panels;
    p_stop_pre_wind{fly} = data(fly).p_stop_pre_wind;
    p_stop_cue_combination{fly} = data(fly).p_stop_cue_combination;
    p_stop_post_panels{fly} = data(fly).p_stop_post_panels;
    p_stop_post_wind{fly} = data(fly).p_stop_post_wind;
    p_stopped = [p_stopped;data(fly).p_stopped];  
    allSummaryData = [allSummaryData;data(fly).summary_data];
    flyNumber = [flyNumber,repelem(fly,length(data(fly).summary_data.BumpMag))];
end
allSummaryData = addvars(allSummaryData,flyNumber');
allSummaryData.Properties.VariableNames{'Var5'} = 'Fly';


%% Plot offset variation per block

set(0,'DefaultTextInterpreter','none')

figure('Position',[100 100 800 800]),
plot(offset_var_r','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmedian(offset_var_r),'-ko','linewidth',3,'MarkerFaceColor','k')
ylabel('Offset variability (rad)','fontweight','bold','fontsize',16);
xticks([1 2 3 4 5]);
xticklabels({'single cue','single cue','cue combination','single cue','single cue'})
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',16,'FontWeight','bold')

%we downsampled the blocks to compute the offset variability so that we
%used the number of frames of the shorter block

saveas(gcf,[path,'\groupPlots\offset_var.png']);
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\offset_var.svg');


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
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\diff_offset_with_respect_to_previous_block.svg');

%% Plot thresholded bump parameters

figure('Position',[100 100 1200 1000]),
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
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\thresh_bump_parameters_evolution.svg');

%% Plot zscored-bump parameters

zscored_BM = zscore(thresh_BM_mean,[],2);
zscored_BW = zscore(thresh_BW_mean,[],2);

figure('Position',[100 100 1200 1000]),
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
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\zscored_bump_parameters_evolution.svg');

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
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\total_mvt.svg');

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

%% Plot the distribution of cue offset differences with respect to the cue combination offset 

for fly = 1:length(data)
    single_cue_cc_diff_1(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,1)));
    single_cue_cc_diff_2(fly) = rad2deg(circ_dist(offset_mean(fly,3),offset_mean(fly,2)));
end

%As box plot
figure,
cue_order = [repelem(1,length(single_cue_cc_diff_1)),repelem(2,length(single_cue_cc_diff_2))];
boxplot([abs(single_cue_cc_diff_1),abs(single_cue_cc_diff_2)],cue_order,'color','k');
hold on
scatter(cue_order,[abs(single_cue_cc_diff_1),abs(single_cue_cc_diff_2)],[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Cue combination offset - single cue offset','fontsize',12);
xticklabels({'First cue','Second cue'});
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12);
ylim([0 180]);


%As line plot
figure,
plot([abs(single_cue_cc_diff_1);abs(single_cue_cc_diff_2)],'-o','color',[.5 .5 .5])
hold on
plot([mean(abs(single_cue_cc_diff_1));mean(abs(single_cue_cc_diff_2))],'-ko','MarkerFaceColor','k','linewidth',2)
xlim([0 3]);
xticks([1 2]);
xticklabels({'First cue','Second cue'});
ylabel('Cue combination offset - single cue offset (deg)','fontsize',12);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12);

saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\single_cue_cc_offset_diff.svg');


%Color coding by the type of cue that came first
figure,
offset_diff = [abs(single_cue_cc_diff_1);abs(single_cue_cc_diff_2)]';
wind_color = [255 36 85]/256;
plot(offset_diff(configuration == 1,:)','-o','color','b')
hold on
plot(offset_diff(configuration == 2,:)','-o','color',wind_color)
plot([mean(abs(single_cue_cc_diff_1));mean(abs(single_cue_cc_diff_2))],'-ko','MarkerFaceColor','k','linewidth',2)
xlim([0 3]);
xticks([1 2]);
xticklabels({'First cue','Second cue'});
ylabel('Cue combination offset - single cue offset','fontsize',12);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12);


%% Cue plasticity vs cc pull (exp 35)

for fly = 1:length(data)
    single_cue_diff_1(fly) = rad2deg(circ_dist(offset_mean(fly,4),offset_mean(fly,1)));
    single_cue_diff_2(fly) = rad2deg(circ_dist(offset_mean(fly,5),offset_mean(fly,2)));
end

figure,
plot(abs(single_cue_cc_diff_1),abs(single_cue_diff_1),'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(abs(single_cue_cc_diff_2),abs(single_cue_diff_2),'ko', 'MarkerFaceColor','k','MarkerSize',8)
ylabel({'Single cue offset difference', '(post vs pre cue combination)'},'fontweight','bold');
xlabel({'Cue combination and initial single cue','offset difference (deg)'},'fontweight','bold');
xticks([0 30 60 90 120 150 180]);
xticklabels({'0','30','60','90','120','150','180'});
yticks([0 30 60 90 120 150 180]);
yticklabels({'0','30','60','90','120','150','180'});
xlim([0 180]); ylim([0 180]);
refline(1,0);

%signed, doubled
figure('Position',[100 100 1000 1000]),
plot(single_cue_cc_diff_1,single_cue_diff_1,'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(single_cue_cc_diff_1,single_cue_diff_1+360,'ko', 'MarkerFaceColor','k','MarkerSize',8)
plot(single_cue_cc_diff_2,single_cue_diff_2,'ko', 'MarkerFaceColor','k','MarkerSize',8)
plot(single_cue_cc_diff_2,single_cue_diff_2+360,'ko', 'MarkerFaceColor','k','MarkerSize',8)
ylabel({'Post - pre (deg)'},'fontweight','bold','fontsize',14);
xlabel({'Training - pre (deg)'},'fontweight','bold','fontsize',14);
xticks([-180 -90 0 90 180]);
xticklabels({'-180','-90','0','90','180'});
yticks([-180 0 180 360 540]);
yticklabels({'-180','0','180','360','540'});
xlim([-180 180]); ylim([-180 540]);
refline(1,0);
refline(1,360);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12);

saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Block-Experiment\plasticity_scales_with_conflict.svg');


%signed
figure,
plot(single_cue_cc_diff_1,single_cue_diff_1,'ko', 'MarkerFaceColor','k','MarkerSize',8)
hold on
plot(single_cue_cc_diff_2,single_cue_diff_2,'ko', 'MarkerFaceColor','k','MarkerSize',8)
ylabel({'Post - pre (deg)'},'fontweight','bold','fontsize',12);
xlabel({'Training - pre (deg)'},'fontweight','bold','fontsize',12);
xticks([-180 -90 0 90 180]);
xticklabels({'-180','-90','0','90','180'});
yticks([-180 -90 0 90 180]);
yticklabels({'-180','-90','0','90','180'});
xlim([-180 180]); ylim([-180 180]);
refline(1,0);
a = get(gca,'YTickLabel');
set(gca,'YTickLabel',a,'fontsize',12);

%% Combine all differences

all_single_cue_cc_diff = [single_cue_cc_diff_1,single_cue_cc_diff_2];
all_single_cue_diff = [single_cue_diff_1,single_cue_diff_2];

%%
close all; clear all;