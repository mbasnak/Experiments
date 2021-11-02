function [smoothed, t, angle, flyPosRad, motor_pos] = get_data_360( trial_time, trial_data, num_frames_x)

%% Processes data from ball, panels and wind device

% DAQ channels and sensor settings 
settings = sensor_settings;
settings.fictrac_yaw_gain_DAQ_AI = 4;
settings.fictrac_x_DAQ_AI = 1;
settings.fictrac_y_DAQ_AI = 2;
settings.panels_DAQ_AI = 5;
settings.panels_DAQ_AI_Y = 7;
settings.wind_valve = 6;
settings.motor = 9;

% When assigning fictrac directions, x will be forward, y will be side
ft_for = trial_data( :, settings.fictrac_x_DAQ_AI );
ft_yaw = trial_data( :, settings.fictrac_yaw_gain_DAQ_AI );
ft_side = trial_data( :, settings.fictrac_y_DAQ_AI );
panels = trial_data( :, settings.panels_DAQ_AI );
panels_y_raw = trial_data( :, settings.panels_DAQ_AI_Y );
motor = trial_data( :, settings.motor );
wind_valve = trial_data( :, settings.wind_valve);

%I'm adding the following lines to make my velocity analysis - MB 20191005
data.ficTracIntx = ft_for;
data.ficTracInty = ft_side;
%I'm inverting the sign in the angular position, for the fly's movement to
%be the opposite of the ball's movement
data.ficTracAngularPosition = ft_yaw;
[smoothed] = singleTrialVelocityAnalysis9mm(data,settings.sampRate);
%%%%%%%%%%%%%%%%%%%%

[angle] = process_panel_360(panels, num_frames_x);
flyPosRad = smoothed.angularPosition;
motor_pos = process_motor_360(motor);

%
[t] = process_time(trial_time);
end

