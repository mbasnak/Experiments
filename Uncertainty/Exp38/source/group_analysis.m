%Code to analyze the group results

clear all; close all;

%% Load data

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp38\data\';

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\data.mat']);
   end
end

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%% Combine all data

thresh_BM_mean = [];
thresh_BW_mean = [];
total_mvt = [];
all_BM_data = [];
all_BW_data = [];
stim_ID = [];
fly_ID = [];
wind_offset_diff = [];
bar_offset_diff = [];
BM_pre_jump = [];
BW_pre_jump = [];

for fly = 1:length(data)
    thresh_BM_mean = [thresh_BM_mean;data(fly).meanBM_perblock_thresh];
    thresh_BW_mean = [thresh_BW_mean;data(fly).meanBW_perblock_thresh];  
    total_mvt = [total_mvt;data(fly).mean_total_mvt];
    all_BM_data = [all_BM_data,data(fly).alldata_BM_thresh.allData_thresh];
    all_BW_data = [all_BW_data,data(fly).alldata_BW_thresh.allDataBW_thresh];
    stim_ID = [stim_ID,data(fly).stim_ID_thresh];
    fly_ID = [fly_ID,repelem(fly,1,length(data(fly).stim_ID_thresh))];
    wind_offset_diff = [wind_offset_diff,data(fly).wind_offset_diff'];
    bar_offset_diff = [bar_offset_diff,data(fly).bar_offset_diff'];
    BM_pre_jump = [BM_pre_jump;data(fly).pre_BM_prejump];
    BW_pre_jump = [BW_pre_jump;data(fly).pre_BW_prejump];
end

%% Plot parameters distribution

figure('Position',[100 100 1200 800]),
subplot(1,2,1)
boxplot(all_BM_data,stim_ID,'color','k')
set(findobj(gca,'type','line'),'linew',2)
ylim([0 5]);
xticks([1 2 3 4])
xticklabels({'Bar only','Wind only','Both cues pre','Both cues post'})
ylabel('Bump magnitude');

subplot(1,2,2)
boxplot(all_BW_data,stim_ID,'color','k')
set(findobj(gca,'type','line'),'linew',2)
ylim([0 5]);
xticks([1 2 3 4])
xticklabels({'Bar only','Wind only','Both cues pre','Both cues post'})
ylabel('Bump width');

saveas(gcf,[path,'\groupPlots\bump_parameters_distribution.png']);


%% Plot parameters mean trends

figure('Position',[100 100 1200 800]),
for fly = 1:length(data)
    subplot(1,2,1)
    plot([mean(all_BM_data(stim_ID == 1 & fly_ID == fly )),mean(all_BM_data(stim_ID == 2 & fly_ID == fly )),mean(all_BM_data(stim_ID == 3 & fly_ID == fly )),mean(all_BM_data(stim_ID == 4 & fly_ID == fly ))],'-','color',[.5 .5 .5])
    hold on
    xlim([0 5]); ylim([0 3]);
    
    subplot(1,2,2)
    plot([mean(all_BW_data(stim_ID == 1 & fly_ID == fly )),mean(all_BW_data(stim_ID == 2 & fly_ID == fly )),mean(all_BW_data(stim_ID == 3 & fly_ID == fly )),mean(all_BW_data(stim_ID == 4 & fly_ID == fly ))],'-','color',[.5 .5 .5])
    hold on
    xlim([0 5]); ylim([0 3]);   
end
subplot(1,2,1)
plot([mean(all_BM_data(stim_ID == 1)),mean(all_BM_data(stim_ID == 2)),mean(all_BM_data(stim_ID == 3)),mean(all_BM_data(stim_ID == 4))],'-ko','linewidth',2)
ylabel('Bump magnitude');

subplot(1,2,2)
plot([mean(all_BW_data(stim_ID == 1)),mean(all_BW_data(stim_ID == 2)),mean(all_BW_data(stim_ID == 3)),mean(all_BW_data(stim_ID == 4))],'-ko','linewidth',2)
ylabel('Bump width');

saveas(gcf,[path,'\groupPlots\bump_parameters_mean.png']);

%% Plot thresholded bump parameters for first 5 bouts

figure,
subplot(1,2,1)
plot(thresh_BM_mean(:,1:5)','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmean(thresh_BM_mean(:,1:5)),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump magnitude');

subplot(1,2,2)
plot(thresh_BW_mean(:,1:5)','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(nanmean(thresh_BW_mean(:,1:5)),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump width');

saveas(gcf,[path,'\groupPlots\thresh_bump_parameters_evolution.png']);


%% Repeat for all bouts

figure,
subplot(1,2,1)
plot(thresh_BM_mean','-o','color',[.5 .5 .5]);
xlim([0 18]);
hold on
plot(nanmean(thresh_BM_mean),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump magnitude');

subplot(1,2,2)
plot(thresh_BW_mean','-o','color',[.5 .5 .5]);
xlim([0 18]);
hold on
plot(nanmean(thresh_BW_mean),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean bump width');

saveas(gcf,[path,'\groupPlots\thresh_bump_parameters_evolution_all_bouts.png']);

%% Zscoring the data for the first 5 bouts

zscored_BM = zscore(thresh_BM_mean,[],2);
zscored_BW = zscore(thresh_BW_mean,[],2);

figure,
subplot(1,2,1)
plot(zscored_BM(:,1:5)','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(zscored_BM(:,1:5)),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean zscored bump magnitude');

subplot(1,2,2)
plot(zscored_BW(:,1:5)','-o','color',[.5 .5 .5]);
xlim([0 6]);
hold on
plot(mean(zscored_BW(:,1:5)),'-ko','linewidth',2)
xlabel('Block #');
ylabel('Mean zscored bump width');

saveas(gcf,[path,'\groupPlots\zscored_bump_parameters_evolution.png']);

%% Repeat for all bouts color coding by bout type

single_cue_bouts = [1,2,4,5,8,9,12,13,16,17];
cue_combination_bouts = [3,6,7,10,11,14,15];

color_cue = ones(1,17);
color_cue(cue_combination_bouts) = 0;

figure('Position',[100 100 1600 500]),
subplot(1,2,1)
plot(zscored_BM','-o','color',[.5 .5 .5]);
xlim([0 18]);
hold on
plot(mean(zscored_BM),'-ko','linewidth',2)
scatter([1:17],mean(zscored_BM),45,color_cue,'filled');
colormap(gca,'cool')
xlabel('Block #');
ylabel('Mean zscored bump magnitude');

subplot(1,2,2)
plot(zscored_BW','-o','color',[.5 .5 .5]);
xlim([0 18]);
hold on
plot(mean(zscored_BW),'-ko','linewidth',2)
scatter([1:17],mean(zscored_BW),45,color_cue,'filled');
colormap(gca,'cool')
xlabel('Block #');
ylabel('Mean zscored bump width');

saveas(gcf,[path,'\groupPlots\zscored_bump_parameters_evolution_all_bouts.png']);

%% Look at wind plasticity vs bump parameters

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(BM_pre_jump',abs(wind_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'pre jump'}); ylabel('Wind offset difference');
ylim([0 180]);

subplot(1,2,2)
plot(BW_pre_jump',abs(wind_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'pre jump'});
ylim([0 180]);

%% Zscore the data

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(zscore(BM_pre_jump)',abs(wind_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'pre jump'}); ylabel('Wind offset difference');
ylim([0 180]);

subplot(1,2,2)
plot(zscore(BW_pre_jump)',abs(wind_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'pre jump'});
ylim([0 180]);

%% Look at bar plasticity vs bump parameters

figure('Position',[100 100 1000 800]),
subplot(1,2,1)
plot(BM_pre_jump',abs(bar_offset_diff),'o')
xlabel({'Bump magnitude in preceding cue combination';'pre jump'}); ylabel('Bar offset difference');
ylim([0 180]);

subplot(1,2,2)
plot(BW_pre_jump',abs(bar_offset_diff),'o')
xlabel({'Bump width in preceding cue combination';'pre jump'});
ylim([0 180]);