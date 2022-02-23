%Code to analyze the imaging and behavior data in EB-DAN experiments

function EBandbulb_data_analysis(path,sid)

tid = 0; %I only ever run 1 trial per session

global slash;
if isunix() == 1 %if running in Linux or mac
    slash = '/'; %define the slash this way
else %if running in windows
    slash = '\'; %define the other way
end
set(0,'DefaultTextInterpreter','none');

%% Create the data analysis folder

data_analysis_dir = [path slash 'analysis' slash];

% if(~exist(data_analysis_dir, 'dir'))
%     mkdir(data_analysis_dir);
% end

%% Load the imaging data

%Move to the folder of interest
cd(path)

%Load the roi data 
load(['2p' slash 'ROI' slash 'ROI_midline_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']);
%Load the registered imaging stack
load(['2p' slash 'sid_' num2str(sid) '_tid_' num2str(tid) slash 'rigid_sid_' num2str(sid) '_tid_' num2str(tid) '_Chan_1_sessionFile.mat']);


%% Sort the ROIs alphabetically

T = struct2table(roi); % convert the struct array to a table
sortedT = sortrows(T, 'name'); % sort the table by 'DOB'
sorted_roi = table2struct(sortedT); % change it back to struct array if necessary


%% Get the summed GCaMP7f data

summedData = squeeze(sum(regProduct,3)); %add the data across the z layers making up each volume, to obtain 1 PB image per timepoint


%% Get baseline and compute DF/F

%For each ROI
figure,
for i = 1:length(sorted_roi)
    for timepoint = 1:size(summedData,3)
        roi_data = summedData(:,:,timepoint);
        trace(timepoint) = sum(roi_data(sorted_roi(i).BW));
    end
    sorted = sort(trace);
    %figure,plot(sorted);
    baseline = prctile(sorted,10);
    dff{i} = (trace - baseline)./baseline;
    
    subplot(length(sorted_roi),1,i)
    plot(dff{i})
    if i == 1
        title('EB');
    elseif i == 2
        title('Left bulb');
    else
        title('Right bulb');
    end
end
saveas(gcf,[path,'\analysis\activity_in_all_neuropiles.png'])

%% Load the ball data

ball_dir = [path slash 'ball' slash];
expression = ['bdata' '*sid_' num2str(sid) '_tid_' num2str(tid) '*.mat'];
ball_file = dir(fullfile(ball_dir, expression));
ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
runobjFile = dir(fullfile([ball_dir 'runobj' slash], ['*_sid_' num2str(sid) '_*']));
load(fullfile([ball_dir 'runobj' slash], runobjFile.name)); %load run_obj


%% Convert the ball data panel position into an angle, bar position; fictrac position/velocities

% 1)Import behavior file
bdata_raw = ballData.trial_bdata; %get the ball data
bdata_time = ballData.trial_time; %get the trial time

% 2)Use an auxiliary function to get the different components of the behavior data
number_x = 96;
[smoothed, bdata_time_out, panel_angle, flyPosRad, motor_pos] = get_data_360(bdata_time, bdata_raw, number_x);

% 3)Recover relevant movement parameters
vel_for = smoothed.xVel';
vel_yaw = smoothed.angularVel;
vel_for_deg = smoothed.xVelDeg';
vel_side_deg = smoothed.yVelDeg';
total_mvt = smoothed.total_mvt;

%If you want to check what the movement data looks like, uncomment the
%following line
%plot_mvt_parameters(vel_for_deg,vel_side_deg,vel_yaw,total_mvt)
panel_y = downsample(bdata_raw(:,7), floor(4000/50));
wind_valve = downsample(bdata_raw(:,6), floor(4000/50));

% 4)Subsample all the variables to have the length of the number of
%volumes scanned:
volumes = length(dff{1}); %get the volume number
time_ds = bdata_time_out(round(linspace(1, length(bdata_time_out), volumes))); 
vel_for_ds = vel_for(round(linspace(1, length(vel_for), volumes)));
vel_yaw_ds = vel_yaw(round(linspace(1, length(vel_yaw), volumes)));
vel_for_deg_ds = vel_for_deg(round(linspace(1, length(vel_for_deg), volumes)));
vel_side_deg_ds = vel_side_deg(round(linspace(1, length(vel_side_deg), volumes)));
total_mvt_ds = total_mvt(round(linspace(1, length(total_mvt), volumes)));
panel_y_ds = panel_y(round(linspace(1, length(panel_y), volumes)));
wind_valve_ds = wind_valve(round(linspace(1, length(wind_valve), volumes)));

%With this convention, a positive change in panel_angle_ds implies a
%clockwise change in bar position.
panel_angle_ds = panel_angle(round(linspace(1, length(panel_angle), volumes)));
motor_pos_ds = motor_pos(round(linspace(1, length(motor_pos), volumes)));

%With this convention, a positive change in flyPosRad_ds implies a
%clockwise change in heading.
%flyPosRad_ds = resample(wrapToPi(flyPosRad),volumes,length(flyPosRad));
flyPosRad_ds = wrapToPi(flyPosRad(round(linspace(1, length(flyPosRad), volumes))));


%%

figure,
yyaxis left
plot(dff{1})
ylabel('EB-DAN activity (dff)');
yyaxis right
plot(total_mvt_ds)
ylabel('Total movement (deg/s)');
saveas(gcf,[path,'\analysis\activity_and_movement.png'])


%% Save the data into the analysis folder

% General experiment info
data.parentDir = path;
data.sid = sid;
data.tid = tid;
data.time = time_ds;
data.run_obj = run_obj;
data.trial_dur = run_obj.trial_t;
data.fr_y_ds = panel_y_ds;

% Movement parameters
data.vel_for = vel_for;
data.vel_yaw = vel_yaw;
data.vel_for_ds = vel_for_ds;
data.vel_yaw_ds = vel_yaw_ds;
data.vel_side_deg_ds = vel_side_deg_ds;
data.vel_for_deg_ds = vel_for_deg_ds;
data.total_mvt_ds = total_mvt_ds;
data.heading = flyPosRad_ds;
data.heading_deg = rad2deg(flyPosRad_ds);

%Devices
data.panel_angle = panel_angle_ds;
data.motor_pos = motor_pos_ds;
data.wind_valve = wind_valve_ds;

% Imaging data
data.volumes = volumes;
data.dff = dff;


% Write file
filename = [path slash 'analysis' slash 'analysis_sid_' num2str(sid) '_tid_' num2str(tid) '.mat'];
save(filename, 'data');

close all;
end