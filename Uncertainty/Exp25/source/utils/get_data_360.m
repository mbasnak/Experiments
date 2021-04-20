function [smoothed, t, angle, flyPosRad] = get_data_360( trial_time, trial_data, num_frames_x, num_frames_y)

%% Processes data from ball and panels for a 360 degree panel arena.

% DAQ channels and sensor settings % check the channel settings, especially
% the panels - MB 20191005
settings = sensor_settings;
settings.fictrac_x_DAQ_AI = 3;
settings.fictrac_yaw_DAQ_AI = 1;
settings.fictrac_y_DAQ_AI = 2;
settings.panels_DAQ_AI = 5;
settings.panels_DAQ_AI_Y = 6;

% When assigning fictrac directions, x will be forward, y will be side
ft_for = trial_data( :, settings.fictrac_x_DAQ_AI );
ft_yaw = trial_data( :, settings.fictrac_yaw_DAQ_AI );
ft_side = trial_data( :, settings.fictrac_y_DAQ_AI );
panels = trial_data( :, settings.panels_DAQ_AI );
panels_y_raw = trial_data( :, settings.panels_DAQ_AI_Y );

%I'm adding the following lines to make my velocity analysis - MB 20191005
data.ficTracIntx = ft_for;
data.ficTracInty = ft_side;
%I'm inverting the sign in the angular position, for the fly's movement to
%be the opposite of the ball's movement
data.ficTracAngularPosition = ft_yaw;
[smoothed] = singleTrialVelocityAnalysis9mm(data,settings.sampRate);
%I'm adding the next line to get the heading the way I do
flyPosRad = smoothed.angularPosition;
%%%%%%%%%%%%%%%%%%%%

[angle] = process_panel_360(panels, num_frames_x);
[panel_y] = process_panel_y(panels_y_raw, num_frames_y);
[t] = process_time(trial_time);
end

