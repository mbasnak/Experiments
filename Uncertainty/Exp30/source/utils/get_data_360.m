function [smoothed, t, angle, flyPosRad] = get_data_360(trial_time, trial_data)

%% Processes data from ball and panels for a 360 degree panel arena.

% DAQ channels and sensor settings, this will change if using a different
% setup
settings = sensor_settings;
settings.fictrac_yaw_gain_DAQ_AI = 1;
settings.fictrac_x_DAQ_AI = 4;
settings.fictrac_yaw_DAQ_AI = 3;
settings.fictrac_y_DAQ_AI = 2;
settings.panels_DAQ_AI = 5;
settings.panels_DAQ_AI_Y = 6;

% Assign fictrac data and compute position and velocities in proper units
data.ficTracIntx = trial_data(:,settings.fictrac_x_DAQ_AI);
data.ficTracAngularPosition = trial_data(:,settings.fictrac_yaw_gain_DAQ_AI);
data.ficTracAngularPosition_uncorrected = trial_data(:,settings.fictrac_yaw_DAQ_AI);
data.ficTracInty = trial_data(:,settings.fictrac_y_DAQ_AI);

[smoothed] = singleTrialVelocityAnalysis9mm(data,settings.sampRate);
flyPosRad = smoothed.angularPosition;

% Get stimulus position
panels = trial_data( :, settings.panels_DAQ_AI);
[angle] = process_panel_360(panels);

% Get time in proper units
[ t ] = process_time(trial_time);

end

