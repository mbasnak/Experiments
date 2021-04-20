%code to compare the inverted gain portion to the empty trial to see how
%the heading or bar offset is much more stable during the inverted gain
%portion

clear all; close all;

%% Load data

%list all of the folders in the directory
folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');

%initialize necessary empty arrays
gain_change_offset_var = [];
empty_offset_var = [];

%for the fly folders
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'60D05')==1 & contains(folderNames(folder).name,'20210127')==0 & contains(folderNames(folder).name,'20210218_60D05_7f_fly2')==0) %I'm excluding the first fly because I didn't run a darkness trial for it
        
        path = [folderNames(folder).folder,'\',folderNames(folder).name];
        %get the sessions info
        load([path,'\sessions_info.mat'])
        
        %load the empty trial and inverted gain trial data
        empty_trial = load([path,'\analysis\analysis_sid_',num2str(sessions_info.empty_trial),'_tid_0.mat']);
        gain_change = load([path,'\analysis\analysis_sid_',num2str(sessions_info.gain_change),'_tid_0.mat']);
        
        %load info about whether this is a type 1 or type 2 fly, to know which
        %offset to compare to
        load([path,'\analysis\gain_change_data.mat'])
        
        % Compute mean offset variability during the inverted gain period
        %get hdf5 files in directory
        file_names = dir(fullfile([path,'\ball\'],'*hdf5'));
        for file = 1:length(file_names)
            if contains(file_names(file).name,['sid_',num2str(sessions_info.gain_change),'_'])
                hdf5_file_to_read = fullfile(file_names(file).folder,file_names(file).name);
            end
        end       
        gain_yaw = double(h5read(hdf5_file_to_read,'/gain_yaw'));
        %downsample to match data length
        gain_yaw_ds = resample(gain_yaw,length(gain_change.data.time),length(gain_yaw));
        gain_yaw_ds(gain_yaw_ds<0) = -1;
        gain_yaw_ds(gain_yaw_ds>0) = 1;        
        %determine gain changes
        gain_changes = find(abs(diff(gain_yaw_ds))>0.5);
        gain_changes = gain_changes(1:2);
        
        % Set block limits      
        blockLimits{1} = [1,gain_changes(1)-1];
        blockLimits{2} = [gain_changes(1),gain_changes(2)];
        blockLimits{3} = [gain_changes(2)+1,length(gain_change.data.time)];
        
        if type_of_fly == 1
            [var empty_trial_offset_var] = circ_std(deg2rad(empty_trial.data.heading_offset));
            empty_offset_var = [empty_offset_var,empty_trial_offset_var];
            [var gain_change_trial_offset_var] = circ_std(deg2rad(gain_change.data.heading_offset(blockLimits{2}(1):blockLimits{2}(2))));
            gain_change_offset_var = [gain_change_offset_var,gain_change_trial_offset_var];
        else
            [var empty_trial_offset_var] = circ_std(deg2rad(empty_trial.data.bar_offset));
            empty_offset_var = [empty_offset_var,empty_trial_offset_var];
            [var gain_change_trial_offset_var] = circ_std(deg2rad(gain_change.data.bar_offset(blockLimits{2}(1):blockLimits{2}(2)))); 
            gain_change_offset_var = [gain_change_offset_var,gain_change_trial_offset_var];
        end
        
    end
end

%% Plot

allData = [empty_offset_var',gain_change_offset_var'];

%plot
figure('Position',[100 100 600 800]),
plot(allData','color',[.5 .5 .5])
hold on
errorbar(mean(allData),std(allData)/13,'-ko','markerfacecolor','k','linewidth',2)
xlim([0 3]);
ylim([0 3]);
xticks([1 2])
xticklabels({'Empty trial', 'Gain change trial'})
ylabel('Circular std of offset');

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\empty_trial_offset_comparison.png');


%% Now let's focus only on flies type I, in the last third of the inverted gain block

clear all; close all;
folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data');

%initialize necessary empty arrays
gain_change_offset_var = [];
empty_offset_var = [];

%for the fly folders
for folder = 1:length(folderNames)
    if (contains(folderNames(folder).name,'60D05')==1 & contains(folderNames(folder).name,'20210127')==0 & contains(folderNames(folder).name,'20210218_60D05_7f_fly2')==0) %I'm excluding the first fly because I didn't run a darkness trial for it
        
        path = [folderNames(folder).folder,'\',folderNames(folder).name];
        %get the sessions info
        load([path,'\sessions_info.mat'])
        
        %load the empty trial and inverted gain trial data
        empty_trial = load([path,'\analysis\analysis_sid_',num2str(sessions_info.empty_trial),'_tid_0.mat']);
        gain_change = load([path,'\analysis\analysis_sid_',num2str(sessions_info.gain_change),'_tid_0.mat']);
        
        %load info about whether this is a type 1 or type 2 fly, to know which
        %offset to compare to
        load([path,'\analysis\gain_change_data.mat'])
        
        % Compute mean offset variability during the inverted gain period
        %get hdf5 files in directory
        file_names = dir(fullfile([path,'\ball\'],'*hdf5'));
        for file = 1:length(file_names)
            if contains(file_names(file).name,['sid_',num2str(sessions_info.gain_change),'_'])
                hdf5_file_to_read = fullfile(file_names(file).folder,file_names(file).name);
            end
        end       
        gain_yaw = double(h5read(hdf5_file_to_read,'/gain_yaw'));
        %downsample to match data length
        gain_yaw_ds = resample(gain_yaw,length(gain_change.data.time),length(gain_yaw));
        gain_yaw_ds(gain_yaw_ds<0) = -1;
        gain_yaw_ds(gain_yaw_ds>0) = 1;        
        %determine gain changes
        gain_changes = find(abs(diff(gain_yaw_ds))>0.5);
        gain_changes = gain_changes(1:2);
        
        % Set block limits      
        blockLimits{1} = [1,gain_changes(1)-1];
        blockLimits{2} = [gain_changes(1),gain_changes(2)];
        blockLimits{3} = [gain_changes(2)+1,length(gain_change.data.time)];
        
        if type_of_fly == 1
            [var empty_trial_offset_var] = circ_std(deg2rad(empty_trial.data.heading_offset));
            empty_offset_var = [empty_offset_var,empty_trial_offset_var];
            [var gain_change_trial_offset_var] = circ_std(deg2rad(gain_change.data.heading_offset(blockLimits{2}(2)-2500:blockLimits{2}(2))));
            gain_change_offset_var = [gain_change_offset_var,gain_change_trial_offset_var];
        else
        end
        
    end
end

allData = [empty_offset_var',gain_change_offset_var'];

%plot
figure('Position',[100 100 600 800]),
plot(allData','color',[.5 .5 .5])
hold on
errorbar(mean(allData),std(allData)/13,'-ko','markerfacecolor','k','linewidth',2)
xlim([0 3]);
ylim([0 3]);
xticks([1 2])
xticklabels({'Empty trial', 'Gain change trial'})
ylabel('Circular std of offset');
title('Focusing on type I flies, last third of inverted gain')

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp28\data\groupPlots\empty_trial_offset_comparison_typeI_lastthird.png');
