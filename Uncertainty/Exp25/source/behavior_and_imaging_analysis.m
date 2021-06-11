function [data,smoothed] = behavior_and_imaging_analysis(parentDir, sid, tid)

%Code to combine the imaging and behavior data



%% Create the data analysis folder
data_analysis_dir = [parentDir,'\analysis\'];
if(~exist(data_analysis_dir, 'dir'))
    mkdir(data_analysis_dir);
end



%% Import behavior file

expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*.mat'];
ball_dir = [parentDir,'\ball\'];
ball_file = dir(fullfile(ball_dir, expression));
ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
runobjFile = dir(fullfile([ball_dir,'\runobj\'], ['*_sid_' num2str(sid) '_*']));
load(fullfile([ball_dir,'\runobj\'], runobjFile.name)); %load run_obj


%% Convert the ball data into appropriate units, and compute velocities

%1)Import behavior file
bdata_raw = ballData.trial_bdata; %get the ball data
bdata_time = ballData.trial_time; %get the trial time

%2)Use an auxiliary function to get the different components of the behavior data
number_x_dim = 96;
number_y_dim = 3;
[smoothed, bdata_time_out, stim_pos, flyPosRad] = get_data_360(bdata_time, bdata_raw, number_x_dim, number_y_dim);

%3)Recover relevant movement parameters
vel_for = smoothed.xVel';
vel_yaw = smoothed.angularVel;
vel_for_deg = smoothed.xVelDeg';
vel_side_deg = smoothed.yVelDeg';
total_mvt = smoothed.total_mvt;

%If you want to check what the movement data looks like, uncomment the
%following line
%plot_mvt_parameters(vel_for_deg,vel_side_deg,vel_yaw,total_mvt)

panel_y = downsample(bdata_raw(:,6), floor(4000/50));

% 4)Subsample all the variables to have the length of the number of
%volumes scanned:
volumes = 
time_ds = bdata_time_out(round(linspace(1, length(bdata_time_out), volumes))); 
vel_for_ds = vel_for(round(linspace(1, length(vel_for), volumes)));
vel_yaw_ds = vel_yaw(round(linspace(1, length(vel_yaw), volumes)));
vel_for_deg_ds = vel_for_deg(round(linspace(1, length(vel_for_deg), volumes)));
vel_side_deg_ds = vel_side_deg(round(linspace(1, length(vel_side_deg), volumes)));
total_mvt_ds = total_mvt(round(linspace(1, length(total_mvt), volumes)));
panel_y_ds = panel_y(round(linspace(1, length(panel_y), volumes)));

%With this convention, a positive change in panel_angle_ds implies a
%clockwise change in bar position.
stim_pos_ds = resample(stim_pos,volumes,length(stim_pos));

%With this convention, a positive change in flyPosRad_ds implies a
%clockwise change in heading.
flyPosRad_ds = resample(flyPosRad,volumes,length(flyPosRad));

%% Compute the offset
%we will compute and define the offset as the
%circular distance between the negative of dff_pva and the heading
%(since a clockwise movement in the EB should elicit a counterclockwise movement of the fly)
heading_offset = wrapTo360(rad2deg(circ_dist(-dff_pva,flyPosRad_ds)));
bar_offset = wrapTo360(rad2deg(circ_dist(dff_pva,deg2rad(stim_pos_ds))));


%% Save the data into the analysis folder

% general experiment info
data.parentDir = parentDir;
data.sid = sid;
data.tid = tid;
data.time = time_ds;
data.run_obj = run_obj;
data.trial_dur = run_obj.trial_t;
data.fr_y_ds = panel_y_ds;

% movement parameters
data.vel_for = vel_for;
data.vel_yaw = vel_yaw;
data.vel_for_ds = vel_for_ds;
data.vel_yaw_ds = vel_yaw_ds;
data.vel_side_deg_ds = vel_side_deg_ds;
data.vel_for_deg_ds = vel_for_deg_ds;
data.total_mvt_ds = total_mvt_ds;
data.heading = flyPosRad_ds;
data.heading_deg = rad2deg(flyPosRad_ds);
data.panel_angle = stim_pos_ds;
data.bar_offset = bar_offset;
data.heading_offset = heading_offset;
data.volumes = volumes;

% imaging data
continuous_data.dff_matrix = dff;
continuous_data.bump_magnitude = bump_mag;
continuous_data.bump_width = bump_width;
continuous_data.adj_rs = adj_rs;
continuous_data.bump_pos = u;

% write file
filename = [data_analysis_dir, ['\continuous_analysis_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']];
save(filename, 'data', 'smoothed');

disp('data analysis done!');

end