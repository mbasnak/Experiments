%%Group analysis for the closed-loop wind data

clear all; close all;

%% Load data

%Get directory you're interested in
path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp33\data');

%I will create a .txt file for each fly with the session numbers
%Read the .txt file containing the ID of the open loop sessions.
fileID = fopen([path,'\2p\wind_sessions.txt'],'r');
formatSpec = '%d';
wind_sessions = fscanf(fileID,formatSpec);

for session = 1:length(wind_sessions)
   data{session} = load([path,'\analysis\summary_data_old_method_sid_',num2str(wind_sessions(session)),'_.mat']);    
end

%% Combine the variables

heading_var = [];
offset_var = [];
meanBM = [];
meanBW = [];
airflow = [];

for session = 1:length(data)
    heading_var = [heading_var,data{1,session}.heading_var];
    offset_var = [offset_var,data{1,session}.offset_var];
    meanBM = [meanBM,data{1,session}.meanBM];
    meanBW = [meanBW,data{1,session}.meanBW];
    airflow = [airflow,data{1,session}.airflow];
end

%% Plot each variable as a function of the airflow

figure('Position',[100 100 800 800]),
subplot(4,1,1)
plot(airflow,heading_var, 'o')
title('Heading variability');
xlim([-0.1 0.8]);
ylim([0 2.5]);

subplot(4,1,2)
plot(airflow,offset_var,'o')
title('Offset variability');
xlim([-0.1 0.8]);
ylim([0 2.5]);

subplot(4,1,3)
plot(airflow,meanBM,'o')
title('Mean bump magnitude')
xlim([-0.1 0.8]);
ylim([0 2]);

subplot(4,1,4)
plot(airflow,meanBW,'o')
title('Mean bump width');
xlabel('Airflow (L/min)');
xlim([-0.1 0.8]);
ylim([1 4]);


save([path,'\group_closed_loop_wind_old_method.png']);


%% Adding movement

for session = 1:length(wind_sessions)
   all_data{session} = load([path,'\analysis\analysis_sid_',num2str(wind_sessions(session)),'_tid_0.mat']);    
end

%%
%ignore the gof, but add a movement threshold
mvt_thresh = 25;
meanBM_thresh = [];
meanBW_thresh = [];
offset_var_all = [];
heading_var_all = [];
for session = 1:length(wind_sessions)
    meanBM_thresh = [meanBM_thresh,mean(all_data{1,session}.data.bump_magnitude(all_data{1,session}.data.total_mvt_ds>mvt_thresh))];
    bump_width = compute_bump_width(all_data{1,session}.data.mean_dff_EB);
    meanBW_thresh = [meanBW_thresh,mean(bump_width(all_data{1,session}.data.total_mvt_ds>mvt_thresh))];    
    [~,offset_var_] = circ_std(deg2rad(all_data{1,session}.data.heading_offset));
    offset_var_all = [offset_var_all,offset_var_];
    [~,heading_var_] = circ_std(all_data{1,session}.data.motor_pos,[],[],2);
    heading_var_all = [heading_var_all,heading_var_];    
end


figure('Position',[100 100 800 800]),
subplot(4,1,1)
plot(airflow,heading_var_all, 'o')
title('Heading variability');
xlim([-0.1 0.8]);
ylim([0 2.5]);

subplot(4,1,2)
plot(airflow,offset_var_all,'o')
title('Offset variability');
xlim([-0.1 0.8]);
ylim([0 2.5]);

subplot(4,1,3)
plot(airflow,meanBM_thresh,'o')
title('Mean bump magnitude')
xlim([-0.1 0.8]);
ylim([0 2]);

subplot(4,1,4)
plot(airflow,meanBW_thresh,'o')
title('Mean bump width');
xlabel('Airflow (L/min)');
xlim([-0.1 0.8]);
ylim([1 4]);

save([path,'\group_closed_loop_wind_mvt_thresh_old_method.jpg']);


%% Save data for across fly comparison

summary_data.airflow = airflow;
summary_data.offset_var = offset_var_all;
summary_data.heading_var = heading_var_all;
summary_data.mean_bm = meanBM_thresh;
summary_data.mean_bw = meanBW_thresh;

save([path,'\analysis\all_summary_data_old_method.mat'],'summary_data')

%% Clean up

clear all; close all