%Code to analyze the results across flies, using the continuous dff method

%Clean workspace
clear all; close all;

%% Import data

%Define main directory
exp_dir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\GCaMP_control\data';

%List all the folders
folderNames = dir(exp_dir);

for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')
        load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\open_loop_data.mat'));
        summary_data{folder} = summarydata;
        flyID = repelem(folder,1,size(summarydata,1));
        summary_data{folder} = addvars(summary_data{folder},flyID');
    end
end

%Remove empty cells
summary_data = summary_data(~cellfun('isempty',summary_data));
all_data = table;
for fly = 1:length(summary_data)
    all_data = [all_data;summary_data{fly}];
end
all_data.Properties.VariableNames{'Var6'} = 'Fly';

%Sort by stimulus velocity
all_data = sortrows(all_data,{'Fly','stim_vel'},{'ascend','ascend'});
%Average data with same stimulus speed
all_data = varfun(@nanmean,all_data,'InputVariables',{'offset_var','bump_mag','bump_width','total_mvt'},...
       'GroupingVariables',{'stim_vel','Fly'});

mean_offset_data = varfun(@mean,all_data,'InputVariables','nanmean_offset_var',...
       'GroupingVariables',{'stim_vel'});
mean_bm_data = varfun(@mean,all_data,'InputVariables','nanmean_bump_mag',...
       'GroupingVariables',{'stim_vel'});
mean_bw_data = varfun(@mean,all_data,'InputVariables','nanmean_bump_width',...
       'GroupingVariables',{'stim_vel'});
mean_mvt_data = varfun(@nanmedian,all_data,'InputVariables','nanmean_total_mvt',...
       'GroupingVariables',{'stim_vel'});
   
%% Plot

figure('Position',[100 100 1400 800]),
    
for fly = 3:6
    subplot(1,4,1)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.nanmean_offset_var(all_data.Fly == fly),'-o','color',[.5 .5 .5])
    hold on
    
    subplot(1,4,2)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.nanmean_bump_mag(all_data.Fly == fly),'-o','color',[.5 .5 .5])
    hold on
    
    subplot(1,4,3)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.nanmean_bump_width(all_data.Fly == fly),'-o','color',[.5 .5 .5])
    hold on
    
    subplot(1,4,4)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.nanmean_total_mvt(all_data.Fly == fly),'-o','color',[.5 .5 .5])
    hold on
    
end
%Add mean trends
subplot(1,4,1)
plot(mean_offset_data.stim_vel,mean_offset_data.mean_nanmean_offset_var,'-ko','linewidth',2)
ylim([0 1.5]);
title('Offset variability');
xlabel('Stimulus velocity (deg/s)');
ylabel('Circular std of offset');

subplot(1,4,2)
plot(mean_bm_data.stim_vel,mean_bm_data.mean_nanmean_bump_mag,'-ko','linewidth',2)
ylim([0 2.5]);
title('Bump magnitude');
xlabel('Stimulus velocity (deg/s)');
ylabel('Bump magnitude');

subplot(1,4,3)
plot(mean_bw_data.stim_vel,mean_bw_data.mean_nanmean_bump_width,'-ko','linewidth',2)
ylim([0 4.5]);
title('Bump width');
xlabel('Stimulus velocity (deg/s)');
ylabel('Bump width');

subplot(1,4,4)
plot(mean_mvt_data.stim_vel,mean_mvt_data.nanmedian_nanmean_total_mvt,'-ko','linewidth',2)
title('Total fly movement');
ylim([0 250]);
xlabel('Stimulus velocity (deg/s)');
ylabel('Total fly movement (deg/s)');

saveas(gcf,[exp_dir,'\groupPlots\parameters_vs_stim_speed.png'])
saveas('C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\GCaMP-control\parameters_vs_stim_speed.svg')
