%Code to see starvation time for each fly in this experiment:

clear all; close all;

%% Load metadata

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp35\data\low_reliability';

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       metadata(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\ball\fly_metadata.mat']);
   end
end

%delete empty cells
all_data = squeeze(struct2cell(metadata));
all_data = all_data(~cellfun('isempty',all_data));

%% Recover starvation times

for fly = 1:length(all_data)
    starvation_time(fly) = str2num(all_data{fly,1}{1,2}{1,1});
end

%% Plot

figure,
boxplot(starvation_time,'color','k')
hold on
yline(0);
scatter(repelem(1,length(starvation_time),1),starvation_time,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Starvation time');

max(starvation_time)
min(starvation_time)

