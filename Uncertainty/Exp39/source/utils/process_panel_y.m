function [fr] = process_panel_y(rawData, frames)
%% 360 arena July 2018

% Initial variables
settings = sensor_settings;
maxVal = 10;
minVal = 0;
voltsPerStep = (maxVal-minVal)./(frames-1);

rate = 2*(50/settings.sampRate);
[kb, ka] = butter(2,rate);

% Set limits on voltage; then filter
rawData(rawData < minVal) = minVal;
rawData(rawData > maxVal) = maxVal;

smoothedData = filtfilt(kb, ka, rawData);

% Calculate the frame number (round to nearest integer), calculate the
% pixel angle of the bar given the grame number.
fr = round((smoothedData - minVal)./voltsPerStep);

%% 4-5-19 downsampling
settings = sensor_settings;

n = floor(settings.sampRate/settings.sensorPollFreq);
fr_downsampled = downsample(fr, n);
fr = fr_downsampled;