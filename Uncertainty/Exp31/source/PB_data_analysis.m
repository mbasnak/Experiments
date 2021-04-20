function [data,smoothed] = PB_data_analysis(parentDir, sid, tid)

%% PB imaging and behavior data processing and analysis
% Performs the following actions:
% Takes behavior data
% Compiles basic experiment data for given trial
% Saves data to analysis folder

%% Look for session Files post registration. Find if there is one or two
% channels.

% find the ROI_analysis folder
roi_analysis_dir = [parentDir,'\2p\ROI_analysis\'];

% find the behavior ball data folder
ball_dir = [parentDir,'\ball\'];

% create the data analysis folder
data_analysis_dir = [parentDir,'\analysis\'];
if(~exist(data_analysis_dir, 'dir'))
    mkdir(data_analysis_dir);
end


%% Get ROI_analysis file, load

expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*.mat'];
roi_analysis_file = dir(fullfile(roi_analysis_dir, expression));
load(fullfile(roi_analysis_dir, roi_analysis_file.name)); %load ROI data (dff, zscore,...)

%Sort by proper ROI order using PB glomeruli order (this way the 1st
%will be the glomerulus on the far right, which will get plotted on the
%top of the imagesc)
order = [8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16];
for i = 1:length(order)
    roi_data(i).order = order(i);%add field with the ROI order
end
T = struct2table(roi_data); % convert the struct array to a table
sortedT = sortrows(T, 'order'); % sort the table by ROI order
roi_data = table2struct(sortedT);


%% Create dF/F matrix (1-16)

volumes = length(roi_data(1).dff); %get the number of imaging volumes
dff_matrix = zeros(16, volumes);

for i = 1:16
    dff_matrix(i, :) = roi_data(i).dff; 
end

%% Create zscored matrix (1-16)

%Next we obtain the zscored values, where the intensity for a given
%glomerulus for each timepoint is subtracted the mean intensity across the
%full session for that glomerulus, and that is divided by the std for that
%same glomerulus across the full session
z_matrix = zeros(16, volumes);
for i = 1:16
    z_matrix(i, :) = roi_data(i).zscore; 
end

%% Calculate power spectrum and phase using Fourier transform

%1) Shift the data to match EB coordinates, i.e., such that the two halves of the glomeruli appear
%as two consecutive equal halves (move the G1 on the le to the end)
shifted_dff_matrix = dff_matrix([9,16,15,14,13,12,11,10,8,7,6,5,4,3,2,1],:);

%2) Apply a fast fourier transform to the dff data
transform = fft(shifted_dff_matrix, 16, 1);
transform(1,:) = []; %Remove the zeroth order component of the transform, which does not contribute to the phase

%3) Compute the phase and amplitude of relevant FFT components
phase = angle(transform(1:8,:));  %take the phase of the first 8 components
phase_value = -phase(2,:); %recover the phase from the second component, which is the one we care about when imaging in the PB
%In this convention, positive changes in phase imply counterclockwise
%movement of the bump around the PB and EB

amplitude =  2*abs(transform/16); %obtain the amplitude of the different FFT components
amplitude = amplitude(2,:); %save the amplitude of the second component, which is the one we care about and corresponds to a bump magnitude estimate

%4) Calculate the power of the 8 glomerulus signal associated with the bump
%(follows the bump magnitude trend)
s = .5.^(1:.025:5); % create a geometric series for the periodogram -- this is the FREQUENCY
pxx = periodogram(shifted_dff_matrix, [], s, 1, 'power'); % POWER SPECTRA
position_eight = find(s == .125); % Look for Period of 8 glomeruli
power_value = pxx(position_eight, :); %find the power of the period=8 component


%5) If you want to check what the past data looks like, uncomment the
%following line
%plot_fft_results(dff_matrix,phase_value,amplitude,power_value)


%% Combine dF/F from L&R glomeruli and calculate PVA

% 1)Separate the PB into the complimentary halves
left_dff = dff_matrix([9,16,15,14,13,12,11,10],:);
right_dff = dff_matrix([8,7,6,5,4,3,2,1],:);

% 2)Average both halves
mean_dff = (left_dff + right_dff)./2;

% 3)Shift the mean to obtain the top of the EB (i.e., glomeruli 5 of the PB))
%at the top
mean_dff_EB = circshift(mean_dff,4);
dff_pva = circ_mean(repmat([pi/8:pi/4:15*pi/8], size(mean_dff_EB,2),1), mean_dff_EB', 2);
%In this convention, positive changes in dff_pva imply clockwise
%movement of the bump around the PB and EB

% 4)Get the mean data for the zscored version too
left_z = z_matrix([9,16,15,14,13,12,11,10],:);
right_z = z_matrix([8,7,6,5,4,3,2,1],:);
mean_z = (left_z + right_z)./2;
mean_z_EB = circshift(mean_z,4);

%% Convert the ball data panel position into an angle, bar position; fictrac position/velocities

% 1)Import behavior file
ball_file = dir(fullfile(ball_dir, expression));
ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
runobjFile = dir(fullfile([ball_dir,'\runobj\'], ['*_sid_' num2str(sid) '_*']));
load(fullfile([ball_dir,'\runobj\'], runobjFile.name)); %load run_obj
bdata_raw = ballData.trial_bdata; %get the ball data
bdata_time = ballData.trial_time; %get the trial time

% 2)Use an auxiliary function to get the different components of the behavior data
[smoothed, bdata_time_out, stim_pos, flyPosRad] = get_data_360(bdata_time, bdata_raw);

% 3)Recover relevant movement parameters
vel_for = smoothed.xVel';
vel_yaw = smoothed.angularVel;
vel_for_deg = smoothed.xVelDeg';
vel_side_deg = smoothed.yVelDeg';
total_mvt = smoothed.total_mvt;

%If you want to check what the movement data looks like, uncomment the
%following line
%plot_mvt_parameters(vel_for_deg,vel_side_deg,vel_yaw,total_mvt)

% 4)Subsample all the variables to have the length of the number of
%volumes scanned:
time_ds = bdata_time_out(round(linspace(1, length(bdata_time_out), volumes))); 
vel_for_ds = vel_for(round(linspace(1, length(vel_for), volumes)));
vel_yaw_ds = vel_yaw(round(linspace(1, length(vel_yaw), volumes)));
vel_for_deg_ds = vel_for_deg(round(linspace(1, length(vel_for_deg), volumes)));
vel_side_deg_ds = vel_side_deg(round(linspace(1, length(vel_side_deg), volumes)));
total_mvt_ds = total_mvt(round(linspace(1, length(total_mvt), volumes)));

stim_pos_ds = stim_pos(round(linspace(1, length(stim_pos), volumes)));
%With this convention, a positive change in stim_pos_ds implies a
%clockwise change in wind position.

flyPosRad_ds = flyPosRad(round(linspace(1, length(flyPosRad), volumes)));
%With this convention, a positive change in flyPosRad_ds implies a
%clockwise change in heading.


%% Offset calculation
%we will compute and define the offset as the
%circular distance between the negative of dff_pva and the heading
%(since a clockwise movement in the EB should elicit a counterclockwise movement of the fly)

heading_offset = wrapTo360(rad2deg(circ_dist(-dff_pva,flyPosRad_ds)));
wind_offset = wrapTo360(rad2deg(circ_dist(dff_pva,deg2rad(stim_pos_ds))));


%% Save the data into the analysis folder

% general experiment info
data.parentDir = parentDir;
data.sid = sid;
data.tid = tid;
data.time = time_ds;
data.run_obj = run_obj;
data.trial_dur = run_obj.trial_t;

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
data.motor_angle = stim_pos_ds;
data.wind_offset = wind_offset;
data.heading_offset = heading_offset;
data.volumes = volumes;

% imaging data
data.z_matrix = z_matrix;
data.dff_matrix = dff_matrix;
data.shifted_dff_matrix = shifted_dff_matrix;
data.ftpower = power_value;
data.phase = phase_value;
data.bump_magnitude = amplitude;
data.mean_dff_EB = mean_dff_EB;
data.mean_z_EB = mean_z_EB;
data.dff_pva = dff_pva;

filename = [data_analysis_dir, ['\analysis_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']];
save(filename, 'data', 'smoothed');

