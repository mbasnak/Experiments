%Code to analyze the pooled data for the offset control

clear all; close all;

%% Load data

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp40\data';

folderNames = dir(path);

for folder = 1:length(folderNames)
   if (contains(folderNames(folder).name,'60D05') ==1)
       data(folder) = load(fullfile(path,folderNames(folder).name,'analysis\offset_var_data.mat'));
   end
end

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));
data2 =  squeeze(cell2mat(struct2cell(data)));

%Reorder to make it bar - wind - empty
reordered_data = data2;
reordered_data(2,:) = data2(3,:);
reordered_data(3,:) = data2(2,:);

%% Plot

figure,
plot(reordered_data,'-o','color',[.5 .5 .5])
hold on
errorbar(1:3,mean(reordered_data,2),std(reordered_data,[],2)/sqrt(size(reordered_data,2)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',6)
xticks(1:3);
xlim([0 4]);
xticklabels({'Bar','Wind','Empty'});
ylim([0 3]);
ylabel('Offset variability (rad)');

saveas(gcf,[path,'\groupPlots\offset_var.png']);
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\Offset-Control\offset_var.svg');