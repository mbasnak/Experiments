function [data,smoothed] = PB_data_analysis(parentDir, sid, tid_list)

%% PB imaging and behavior data processing and analysis
% Performs the following actions:
% Takes behavior data
% Cleans fictrac dropout problems and imaging dropout problems (parameter to include or exclude datapoint)
% Removes first 3 seconds of trial for pattern start
% Compiles basic experiment data for given trial
% Saves data to analysis folder: 
%     data.parentDir = parentDir;
%     data.sid = sid;
%     data.tid = tid;
%     data.time = time_ds;
%     data.vel_for = vel_for;
%     data.vel_side = vel_side;
%     data.vel_yaw = vel_yaw;
%     data.vel_for_ds = vel_for_ds;
%     data.vel_side_ds = vel_side_ds;
%     data.vel_yaw_ds = vel_yaw_ds;
%     data.fr = fr;
%     data.panel_angle = panel_angle_rad;
%     data.volumes = volumes;
%     data.z_matrix = z_matrix;
%     data.dff_matrix = dff_matrix;
%     data.ftpower = power_value;
%     data.phase = phase_value;
%     data.mean_dff = mean_dff;
%     data.dff_pva_rad = dff_pva_rad;
%     data.dff_pva_glomerulus = dff_pva;
%     data.offset = offset;
%     data.mean_offset = mean_offset;
%     data.power_period = 1./s;
%     data.power_specta = pxx;
%     data.pattern = run_obj.pattern_number;
%     data.run_obj = run_obj;
%     data.exclude_timepoint_ds = exclude_timepoint_ds;
%     data.include_timepoint_ds = ones(1, length(exclude_timepoint_ds)) - exclude_timepoint_ds;

%% Look for session Files post registration. Find if there is one or two
% channels.

global slash;
if isunix() == 1
    slash = '/';
else
    slash = '\';
end

% find the ROI_analysis folder
roi_analysis_dir = [parentDir slash '2p' slash 'ROI_analysis' slash];

% find the behavior ball data folder
ball_dir = [parentDir slash 'ball' slash];

% create the data analysis folder
data_analysis_dir = [parentDir slash 'analysis' slash];
if(~exist(data_analysis_dir, 'dir'))
    mkdir(data_analysis_dir);
end

disp('starting analysis..');

num_tids = length(tid_list);

for i = 1:num_tids %for each trial
    
    tid = tid_list(i); %get the trial ID
    
    %% Get ROI_analysis file, load
    expression = ['*sid_' num2str(sid) '_tid_' num2str(tid) '*.mat'];
    roi_analysis_file = dir(fullfile(roi_analysis_dir, expression));
    load(fullfile(roi_analysis_dir, roi_analysis_file.name)); %load ROI data (dff, zscore,...)
    
    %Sort by proper ROI order (using PB glomeruli order)-I have added this
    %MB
    order = [8,9,7,10,6,11,5,12,4,13,3,14,2,15,1,16];
    for i = 1:length(order)
        roi_data(i).order = order(i);%add field with the ROI order
    end
    T = struct2table(roi_data); % convert the struct array to a table
    sortedT = sortrows(T, 'order'); % sort the table by ROI order
    roi_data = table2struct(sortedT);   
    
    %% Import behavior file
    ball_file = dir(fullfile(ball_dir, expression));
    ballData = load(fullfile(ball_dir, ball_file.name)); %load ballData
    runobjFile = dir(fullfile([ball_dir slash 'runobj' slash], ['*_sid_' num2str(sid) '_*']));
    load(fullfile([ball_dir slash 'runobj' slash], runobjFile.name)); %load run_obj
    
    %% Create dF/F matrix (1-16)    
    volumes = length(roi_data(1).dff); %get the volume number
    
    dff_matrix = zeros(16, volumes);
    
    for i = 1:16
        dff_matrix(i, :) = roi_data(i).dff; %I changed this because the field 'num' does not exist - MB 20190917
    end
    
    %% Create z matrix (1-16)
    
    z_matrix = zeros(16, volumes);
    for i = 1:16
        z_matrix(i, :) = roi_data(i).zscore; %same
    end
     
    %% Calculate power spectrum and phase using Fourier transform
      
    transform = fft(dff_matrix, 16, 1); %apply a fast fourier transform to the dff data
    %when one or more of the glomeruli give you nans, all the values in the
    %transform will, and then the phase will just be a nan.
    
    transform(1,:) = []; %The first row is taken away. This doesn't mean that the first ROI is being taken,
    %but the first row of the transform. 
    %I think the first row might correspond to the zero-frequency component
    %of the transform, so it's not useful to us.    
   
    n = 16;
    phase = angle(transform(1:8,:));  %take the phase of the FFT of period 8 glomeruli
    phase_value = -squeeze(phase(2,:)); %% This is actually the negative of the phase. However, in order to maintain consistency with Green et al., I realized I had to make the FFT phase negative.
    
    s = .5.^(1:.025:5); % create a geometric series for the periodogram -- this is the FREQUENCY
    pxx = periodogram(dff_matrix, [], s, 1, 'power'); % POWER SPECTRA
    %pxx = periodogram(z_matrix, [], s, 1, 'power'); % I'm changing the power to be calculated with the zscored data to avoid nans.
    position_eight = find(s == .125); % Look for Period of 8 glomeruli
    power_value = pxx(position_eight, :);
    
    %% Combine dF/F from L&R glomeruli and calculate PVA
    
    left_dff = dff_matrix(1:8,:);
    right_dff = dff_matrix(9:16,:);
       
    mean_dff = (left_dff + right_dff)./2;
    mu = circ_mean(repmat([-7*pi/8:pi/4:7*pi/8], size(mean_dff,2),1), mean_dff', 2);
    dff_pva = ((mu)/pi*4+4.5);
    dff_pva = dff_pva';
    dff_pva_rad = mu';
    %dff_pva_rad = circ_mean(mean_dff,[],2);
    
    %% Convert the ball data panel position into an angle, bar position; fictrac position/velocities
    
    bdata_raw = ballData.trial_bdata; %get the ball data
    bdata_time = ballData.trial_time; %get the trial time
    number_y = 2; %I've changed this to 2 because my stimuli have 2 y dimension - MB 20191004
    
    %use a function to get the different components of the behavior data
    %change channel settings!
    [smoothed, bdata_time_out, ~, d_yaw, ~, fr, panel_angle, vel_for, vel_side, vel_yaw, ~,flyPosRad] = get_data_360(bdata_time, bdata_raw, run_obj.number_frames, number_y);    
   
    %I'm adding lines to have the fwd and side velocity use my vel
    %calculation method instead of Jenny's - MB 202002210
    vel_for = smoothed.xVel';
    vel_yaw = smoothed.angularVel;
    
    panel_y = downsample(bdata_raw(:,6), floor(4000/50));

    n = length(bdata_time_out);
    time_ds = bdata_time_out(round(linspace(1, n, volumes))); %get a vector of length the total trial time and as many points as the volumes scanned
    %subsample all the variables to have the length of the number of
    %volumes scanned:
    %vel_for_ds = vel_for(round(linspace(1, n, volumes)));
    vel_for_ds = resample(vel_for,volumes,length(vel_for));
    vel_side_ds = vel_side(round(linspace(1, n, volumes)));
    %vel_yaw_ds = vel_yaw(round(linspace(1, n, volumes)));
    vel_yaw_ds = resample(vel_yaw,volumes,length(vel_yaw));
    panel_angle_ds = panel_angle(round(linspace(1, n, volumes)));
    fr_ds = fr(round(linspace(1, n, volumes)));
    panel_y_ds = panel_y(round(linspace(1, n, volumes)));
    %flyPosRad_ds = wrapToPi(flyPosRad(round(linspace(1, n, volumes))));
    flyPosRad_ds = resample(wrapToPi(flyPosRad),volumes,length(flyPosRad));
    fr = fr';
    
    %% Calculate continuous phase offset throughout the trial
    
    panel_angle_rad = panel_angle_ds'.*2*pi./360;
    %obtain offset between neural data and behavioral data
    offset = wrapToPi(phase_value - panel_angle_rad);
    offset2 = offset(~isnan(offset));
    %I'm adding the next line with the new fly heading
    offset3 = wrapToPi(phase_value - flyPosRad_ds');
    mean_offset = circ_mean(offset2'); %this can be used to shift the phase for plotting.
    
    
    %% DATA ANALYSIS CLEAN: EVALUATING WHETHER TO INCLUDE OR EXCLUDE TIMEPOINT
    % Criterion 1: Fictrac dropout filtering
    % These thresholds were empirically determined.
    % Delete frames with too abrupt vel changes.
    ft_dropout = zeros(length(vel_yaw_ds)-1,1);
    ft_dropout(diff(vel_yaw_ds)>10) = 1;
    ft_dropout(diff(vel_yaw_ds)<-10) = 1;
    ft_dropout(diff(vel_for_ds)>10) = 1;
    ft_dropout(diff(vel_for_ds)<-10) = 1;
    ft_dropout = [1; ft_dropout]; 
    ft_dropout(movvar(vel_yaw_ds, 10) < 0.001) = 1;
    ft_dropout(movvar(vel_for_ds, 10) < 0.001) = 1;
    ft_dropout_smooth_ds = (movmean(ft_dropout,50)>.2)';
%     if size(ft_dropout_smooth_ds, 2)>size(vel_yaw_ds, 1)
%         ft_dropout_smooth_ds = ft_dropout_smooth_ds(1:end-1);
%     end
    
    ft_dropout_smooth_nonds = logical(interp1(time_ds, double(ft_dropout_smooth_ds), bdata_time_out));
    
    % Criterion 2: Imaging dropout or Power of fourier insufficient (using
    % empirical value of .1)
    
    imaging_dropout_ds = power_value < 0.1;
    
    % Criterion 3: Exclude the first 3 seconds of every trial for pattern
    % startup
    exclude_timepoint_ds = ft_dropout_smooth_ds + imaging_dropout_ds;
    exclude_timepoint_ds(exclude_timepoint_ds >= 2) = 1;
    exclude_timepoint_ds(time_ds<3) = 1;

    
    %% Save the data into the analysis folder
    data.parentDir = parentDir;
    data.sid = sid;
    data.tid = tid;
    data.time = time_ds;
    data.vel_for = vel_for;
    data.vel_side = vel_side;
    data.vel_yaw = vel_yaw;
    data.vel_for_ds = vel_for_ds;
    data.vel_side_ds = vel_side_ds;
    data.vel_yaw_ds = vel_yaw_ds;
    data.fr = fr;
    data.fr_y_ds = panel_y_ds;
    data.fr_ds = fr_ds;
    data.panel_angle = panel_angle_rad;
    data.flyPosRad = flyPosRad_ds;
    data.volumes = volumes;
    data.z_matrix = z_matrix;
    data.dff_matrix = dff_matrix;
    data.ftpower = power_value;
    data.phase = phase_value;
    data.mean_dff = mean_dff;
    data.dff_pva_rad = dff_pva_rad;
    data.dff_pva_glomerulus = dff_pva;
    data.offset = offset;
    data.offset3 = offset3;
    data.mean_offset = mean_offset;
    data.power_period = 1./s;
    data.power_spectra = pxx;
    data.pattern = run_obj.pattern_number;
    data.run_obj = run_obj;
    data.exclude_timepoint_ds = exclude_timepoint_ds;
    data.include_timepoint_ds = ones(1, length(exclude_timepoint_ds)) - exclude_timepoint_ds;
    data.ft_dropout = ft_dropout_smooth_nonds;
    data.jump_time = run_obj.bar_jump_time;
    data.trial_dur = run_obj.trial_t;
    
    filename = [data_analysis_dir slash ['analysis_sid_' num2str(sid) '_tid_' num2str(tid) '.mat']];
    
    save(filename, 'data', 'smoothed');

    
disp('data analysis done!');
end