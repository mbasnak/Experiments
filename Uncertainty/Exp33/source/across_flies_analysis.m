%Code to analyze the results of the closed-loop wind experiment at
%different intensities across flies

%Clean workspace
clear all; close all;

%% Import data

%Define main directory
exp_dir = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp33\data';

%List all the folders
folderNames = dir(exp_dir);

%For the experimental folders, import the summary data
for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')
        data{folder} = load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\all_summary_data.mat'));
    end
end

%Remove empty cells
data = data(~cellfun('isempty',data));

%% Create arrays combining all the important variables


for fly = 1:length(data)
    airflow{fly} = data{1,fly}.summary_data.airflow;
    heading_var{fly} = data{1,fly}.summary_data.heading_var;
    offset_var{fly} = data{1,fly}.summary_data.offset_var;
    mean_bm{fly} = data{1,fly}.summary_data.mean_bm;
    mean_bw{fly} = data{1,fly}.summary_data.mean_bw;
end
 
%% Replicate the individual plot with all the data points

%Create a table with the variables
for folder = 1:length(folderNames)
    if contains(folderNames(folder).name,'60D05')
        load(fullfile(folderNames(folder).folder,folderNames(folder).name,'analysis\all_summary_data.mat'));
        alldata{folder} = array2table([summary_data.airflow',summary_data.heading_var',summary_data.offset_var',summary_data.mean_bm',summary_data.mean_bw'],'VariableNames', {'Airflow','HeadingVar','OffsetVar','MeanBM','MeanBW'});
        flyID = repelem(folder,1,length(summary_data.airflow));
        alldata{folder} = addvars(alldata{folder},flyID');
    end
end

%Remove empty cells
alldata = alldata(~cellfun('isempty',alldata));
all_data = table;

for fly = 1:length(data)
    all_data = [all_data;alldata{fly}];
end

all_data.Properties.VariableNames{'Var6'} = 'Fly';

%% Get mean variables by airflow

mean_offset_data = varfun(@mean,all_data,'InputVariables','OffsetVar',...
       'GroupingVariables',{'Airflow'});
   
mean_heading_data = varfun(@mean,all_data,'InputVariables','HeadingVar',...
       'GroupingVariables',{'Airflow'});
   
mean_bm_data = varfun(@mean,all_data,'InputVariables','MeanBM',...
       'GroupingVariables',{'Airflow'});
   
mean_bw_data = varfun(@mean,all_data,'InputVariables','MeanBW',...
       'GroupingVariables',{'Airflow'});
   

%%

figure('Position',[100 100 800 800]),
    
for fly = 1:length(data)
    
    subplot(4,1,1)
    plot(airflow{fly},heading_var{fly}, 'o','color',[.5 .5 .5])
    hold on
    title('Heading variability');
    xlim([-0.1 0.8]);
    ylim([0 2.5]);
    
    subplot(4,1,2)
    plot(airflow{fly},offset_var{fly},'o','color',[.5 .5 .5])
    hold on
    title('Offset variability');
    xlim([-0.1 0.8]);
    ylim([0 2.5]);
    
    subplot(4,1,3)
    plot(airflow{fly},mean_bm{fly},'o','color',[.5 .5 .5])
    hold on
    title('Mean bump magnitude')
    xlim([-0.1 0.8]);
    ylim([0 2]);
    
    subplot(4,1,4)
    plot(airflow{fly},mean_bw{fly},'o','color',[.5 .5 .5])
    hold on
    title('Mean bump width');
    xlabel('Airflow (L/min)');
    xlim([-0.1 0.8]);
    ylim([1 5]);

end
%Add mean trends
subplot(4,1,1)
plot(mean_offset_data.Airflow,mean_offset_data.mean_OffsetVar,'-ko','linewidth',2)

subplot(4,1,2)
plot(mean_heading_data.Airflow,mean_heading_data.mean_HeadingVar,'-ko','linewidth',2)

subplot(4,1,3)
plot(mean_bm_data.Airflow,mean_bm_data.mean_MeanBM,'-ko','linewidth',2)

subplot(4,1,4)
plot(mean_bw_data.Airflow,mean_bw_data.mean_MeanBW,'-ko','linewidth',2)


save('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp33\data\groupPlots\across_flies_analysis.png');

