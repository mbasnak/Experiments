%I'm commenting this code and changing it to make sense with my set-up - MB
%20190808

function [stim_pos] = process_motor_360(rawData)

% Initial variables
settings = sensor_settings;
maxVal = 10;
minVal = 0;

%invert data
inv_data = 10-rawData;

%compute in deg
deg_data = inv_data*360/10;

%unwrap
unwrapped_data = rad2deg(unwrap(deg2rad(deg_data)));

%downsample data
ds_data = downsample(unwrapped_data, settings.sampRate/25);

%wrap
stim_pos = wrapTo360(ds_data);