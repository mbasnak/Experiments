function [smoothed, t, motor_angle, flyPosRad] = get_data_360(trial_time, trial_data)

%% Processes data from ball and panels for the wind experiments

% DAQ channels and sensor settings, this will change if using a different
% setup
settings = sensor_settings;
settings.fictrac_yaw_gain_DAQ_AI = 1;
settings.fictrac_x_DAQ_AI = 5;
settings.fictrac_yaw_DAQ_AI = 3;
settings.fictrac_y_DAQ_AI = 2;
settings.motor = 8;

% Assign fictrac data and compute position and velocities in proper units
data.ficTracIntx = trial_data(:,settings.fictrac_x_DAQ_AI);
data.ficTracAngularPosition = trial_data(:,settings.fictrac_yaw_gain_DAQ_AI);
data.ficTracAngularPosition_uncorrected = trial_data(:,settings.fictrac_yaw_DAQ_AI);
data.ficTracInty = trial_data(:,settings.fictrac_y_DAQ_AI);

[smoothed] = singleTrialVelocityAnalysis9mm(data,settings.sampRate);
flyPosRad = smoothed.angularPosition;

% Get stimulus position
motor = trial_data( :, settings.motor);
[motor_angle] = process_motor_360(motor);

% Get time in proper units
[ t ] = process_time(trial_time);

end

