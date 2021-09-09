clear all; close all

%% Load data

folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\pilots\high_contrast');

change_bm = [];

for folder = 1:length(folderNames)
   if contains(folderNames(folder).name,'60D05') == 1
       load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\bm_change.mat'))
       mean_bm_aj = mean(bm_aj,2);
       change_bm = [change_bm,mean_bm_aj];
   end
end

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp30\data\pilots\high_contrast\group_plots';

%% Plot

figure,
plot(change_bm,'-o','color',[.5 .5 .5])
hold on
plot(mean(change_bm,2),'-ko','linewidth',2)
xlim([0 3]);
ylim([0.5 1.3]);
xticks([1 2])
xticklabels({'pre-jump','post-jump'})
ylabel('Mean bump magnitude (max-min)');
saveas(gcf,[path,'\change_in_bm_AJ.png'])