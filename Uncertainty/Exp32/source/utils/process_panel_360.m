%I'm commenting this code and changing it to make sense with my set-up - MB
%20190808

function [stim_pos] = process_panel_360(rawData)

% Initial variables
settings = sensor_settings;
maxVal = 10;
minVal = 0;
initialPx = 3; %my x=1 positon is 3 pixels to the right of the animal
initialAngle = initialPx*360/96;

%invert the direction of the data to match my convention
inv_data = 10 - rawData;

%compute in deg
deg_data = inv_data*360/10;

%unwrap
unwrapped_data = rad2deg(unwrap(deg2rad(deg_data)));

%shift according to the panels pos
shifted_data = unwrapped_data + initialAngle;

%downsample data
ds_data = downsample(shifted_data, settings.sampRate/25);

%wrap
stim_pos = wrapTo360(ds_data);