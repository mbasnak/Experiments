%I'm commenting this code and changing it to make sense with my set-up - MB
%20190808

function [angle] = process_panel_360(rawData, frames)

% Initial variables
settings = sensor_settings;
maxVal = 10;
minVal = 0;
%initialAngle = -3; %I think this means when the x voltage for the panels reads 0 her bar is 15 pixels to one side. In mine it is 3 px to the other side
initialAngle = 0;
barWidth = 2;
voltsPerStep = (maxVal-minVal)./(frames-1); %how much voltage corresponds to a pixel movement

rate = 2*(50/settings.sampRate); %I don't know why this would be the rate
[kb, ka] = butter(2,rate);

% Set limits on voltage; then filter
rawData(rawData < minVal) = minVal; %remove voltages under 0 if there was any
rawData(rawData > maxVal) = maxVal; %remove voltages over 10 if there was any

smoothedData = filtfilt(kb, ka, rawData); %filter the data.

% Calculate the frame number (round to nearest integer), calculate the
% pixel angle of the bar given the grame number.
fr = round((smoothedData - minVal)./voltsPerStep);
pixelAngle = 360./96;
arenaAngle = frames*pixelAngle; %this is just 360
angle = (initialAngle-((fr-1)+barWidth/2).*pixelAngle); % accounts for the bar width
%angle = (fr+1+initialAngle).*pixelAngle;

% Wrap to 180 
angle = wrapTo180(angle);
if arenaAngle < 360
    halfArena = arenaAngle./2;
    indexOver = angle < -halfArena;
    angle = angle + indexOver.*arenaAngle;
end

%% 4-5-19 downsampling
settings = sensor_settings;

n = floor(settings.sampRate/settings.sensorPollFreq);

angle_downsampled = downsample(angle, n);
angle = angle_downsampled;
