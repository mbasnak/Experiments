%I'm commenting this code and changing it to make sense with my set-up - MB
%20190808

function [stim_pos] = process_panel_360(rawData, frames)

% Initial variables
settings = sensor_settings;
maxVal = 10;
minVal = 0;
initialPx = 3; %my x=1 positon is 3 pixels to the right of the animal
initialAngle = initialPx*360/frames;

%invert the direction of the data to match my convention
inv_data = 10 - rawData;

%compute in deg
deg_data = inv_data*360/10;

%shift according to the panels pos
shifted_data = deg_data + initialAngle;

%wrap
stim_pos = wrapTo360(shifted_data);

%downsample data
stim_pos = resample(stim_pos, 25, settings.sampRate);