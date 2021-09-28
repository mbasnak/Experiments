clear all; close all

%% Load data

folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\pilots\high_contrast');

change_bm = [];
change_bw = [];

for folder = 1:length(folderNames)
   if contains(folderNames(folder).name,'60D05') == 1
       load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\bump_parameter_change.mat'))
       mean_bm_aj = mean(bm_aj,2);
       mean_bw_aj = mean(bw_aj,2);
       change_bm = [change_bm,mean_bm_aj];
       change_bw = [change_bw,mean_bw_aj];
   end
end

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\pilots\high_contrast\group_plots';

%% Plot BM

figure,
plot(change_bm,'-o','color',[.5 .5 .5])
hold on
plot(mean(change_bm,2),'-ko','linewidth',2)
xlim([0 3]);
ylim([min(min(change_bm))-0.5 max(max(change_bm))+0.5]);
xticks([1 2])
xticklabels({'pre-jump','post-jump'})
ylabel('Mean bump magnitude (max-min)');
saveas(gcf,[path,'\change_in_bm_AJ.png'])


%% Plot BW

figure,
plot(change_bw,'-o','color',[.5 .5 .5])
hold on
plot(mean(change_bw,2),'-ko','linewidth',2)
xlim([0 3]);
ylim([min(min(change_bw))-0.5 max(max(change_bw))+0.5]);
xticks([1 2])
xticklabels({'pre-jump','post-jump'})
ylabel('Mean bump width');
saveas(gcf,[path,'\change_in_bw_AJ.png'])

