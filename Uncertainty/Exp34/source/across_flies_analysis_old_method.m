%Code to analyze the results across flies, using the continuous dff method

%Clean workspace
clear all; close all;

%% Import data

%Define main directory
exp_dir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp34\data';

%List all the folders
folderNames = dir(exp_dir);

for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')
        load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\open_loop_data_old_method.mat'));
        summary_data{folder} = summarydata;
        flyID = repelem(folder,1,size(summarydata,1));
        summary_data{folder} = addvars(summary_data{folder},flyID');
        
        fly = repelem(folder,1,size(bump_mag_data,1));
        bump_mag{folder} = bump_mag_data;
        bump_mag{folder} = addvars(bump_mag_data,fly');
        
        bump_width{folder} = bump_width_data;
        bump_width{folder} = addvars(bump_width_data,fly');
    end
end

%Remove empty cells
summary_data = summary_data(~cellfun('isempty',summary_data));
bump_mag = bump_mag(~cellfun('isempty',bump_mag));
bump_width = bump_width(~cellfun('isempty',bump_width));
all_data = table;
model_data = table;
model_data_bw = table;
for fly = 1:length(summary_data)
    all_data = [all_data;summary_data{fly}];
    model_data = [model_data;bump_mag{fly}];
    model_data_bw = [model_data_bw;bump_width{fly}];    
end
all_data.Properties.VariableNames{'Var6'} = 'Fly';
model_data.Properties.VariableNames{'Var6'} = 'Fly';
model_data_bw.Properties.VariableNames{'Var6'} = 'Fly';

mean_offset_data = varfun(@mean,all_data,'InputVariables','offset_var',...
       'GroupingVariables',{'stim_vel'});
mean_bm_data = varfun(@mean,all_data,'InputVariables','bump_mag',...
       'GroupingVariables',{'stim_vel'});
mean_bw_data = varfun(@mean,all_data,'InputVariables','bump_width',...
       'GroupingVariables',{'stim_vel'});
mean_mvt_data = varfun(@nanmedian,all_data,'InputVariables','total_mvt',...
       'GroupingVariables',{'stim_vel'});
   
%% Plot

figure('Position',[100 100 1400 800]),
    
for fly = 1:length(summary_data)
    subplot(1,4,1)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.offset_var(all_data.Fly == fly),'o','color',[.5 .5 .5])
    hold on
    
    subplot(1,4,2)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.bump_mag(all_data.Fly == fly),'o','color',[.5 .5 .5])
    hold on
    
    subplot(1,4,3)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.bump_width(all_data.Fly == fly),'o','color',[.5 .5 .5])
    hold on
    
    subplot(1,4,4)
    plot(all_data.stim_vel(all_data.Fly == fly),all_data.total_mvt(all_data.Fly == fly),'o','color',[.5 .5 .5])
    hold on
    
end
%Add mean trends
subplot(1,4,1)
plot(mean_offset_data.stim_vel,mean_offset_data.mean_offset_var,'-ko','linewidth',2)
ylim([0 1.5]);
title('Offset variability');
xlabel('Stimulus velocity (deg/s)');
ylabel('Circular std of offset');

subplot(1,4,2)
plot(mean_bm_data.stim_vel,mean_bm_data.mean_bump_mag,'-ko','linewidth',2)
ylim([0 1.5]);
title('Bump magnitude');
xlabel('Stimulus velocity (deg/s)');
ylabel('Bump magnitude (from von Mises fit)');

subplot(1,4,3)
plot(mean_bw_data.stim_vel,mean_bw_data.mean_bump_width,'-ko','linewidth',2)
ylim([0 4.5]);
title('Bump width');
xlabel('Stimulus velocity (deg/s)');
ylabel('Bump width');

subplot(1,4,4)
plot(mean_mvt_data.stim_vel,mean_mvt_data.nanmedian_total_mvt,'-ko','linewidth',2)
title('Total movement');
ylim([0 250]);
xlabel('Stimulus velocity (deg/s)');
ylabel('Total movement (deg/s)');


%% Run model

model_data.Fly = categorical(model_data.Fly);
model_data_bw.Fly = categorical(model_data_bw.Fly);

mdl_BM = fitlme(model_data,'BumpMagnitude ~ stimVel+TotalMovement+(1|Fly)')
mdl_BW = fitlme(model_data_bw,'BumpWidth ~ stimVel+TotalMovement+(1|Fly)')
