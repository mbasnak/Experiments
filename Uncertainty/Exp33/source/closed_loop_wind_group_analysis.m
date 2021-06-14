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
   data{session} = load([path,'\analysis\summary_data_sid_',num2str(wind_sessions(session)),'_.mat']);    
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

figure,
subplot(4,1,1)
plot(airflow,heading_var, 'o')
title('Heading variability');

subplot(4,1,2)
plot(airflow,offset_var,'o')
title('Offset variability');

subplot(4,1,3)
plot(airflow,meanBM,'o')
title('Mean bump magnitude');

subplot(4,1,4)
plot(airflow,meanBW,'o')
title('Mean bump width');
xlabel('Airflow (L/min)');

