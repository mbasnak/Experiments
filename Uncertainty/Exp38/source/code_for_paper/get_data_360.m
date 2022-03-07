function [smoothed, t, visual_stim, fly_pos_rad, motor_pos] = get_data_360(trial_time, trial_data, num_x_px)

%% Processes data from ball and panels for a 360 degree panel arena.

% DAQ channels and sensor settings
settings = sensor_settings;
settings.fictrac_yaw_gain_DAQ_AI = 4;
settings.fictrac_x_DAQ_AI = 1;
settings.fictrac_y_DAQ_AI = 2;
settings.panels_DAQ_AI = 5;
settings.panels_DAQ_AI_Y = 7;
settings.wind_valve = 6;
settings.motor = 9;

% Assign fictrac positions on the three directions, and stimulus position,
% based on the channels' nature
panels = trial_data( :, settings.panels_DAQ_AI );
motor = trial_data( :, settings.motor );
data.ficTracIntx = trial_data( :, settings.fictrac_x_DAQ_AI );
data.ficTracInty = trial_data( :, settings.fictrac_y_DAQ_AI );
data.ficTracAngularPosition = trial_data( :, settings.fictrac_yaw_gain_DAQ_AI );

% Run auxiliary function to compute the velocity
smoothed = singleTrialVelocityAnalysis9mm(data,settings.sampRate);
fly_pos_rad = smoothed.angularPosition;

% Run auxiliary function to process the stimulus and time data
visual_stim = process_panel_360(panels, num_x_px);
motor_pos = process_motor_360(motor);
t = process_time(trial_time);

end

