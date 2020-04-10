function [smoothed] = posDataDecoding(data, sampleRate)

%OUTPUT
%smoothed = has the smoothed position velocity information from FicTrac's

%INPUT
%data = the raw data from FicTrac
%sampleRate = the acquisition rate of the NiDaq

%% (1) Downsample (I downsample to 25 Hz because my FicTrac usually runs at 50 Hz)
% I'm using the function 'resample' instead of 'downsample' because it
% applies anti-aliasing filter before sampling

 downs.angularPosition = resample(data.ficTracAngularPosition,25,sampleRate);
 downs.Intx = resample(data.ficTracIntx,25,sampleRate);
 downs.Inty = resample(data.ficTracInty,25,sampleRate);


%% (2) Transform units to rad

 downsRad.angularPosition = downs.angularPosition .* 2 .* pi ./ 10;
 downsRad.Intx =  downs.Intx.* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
 downsRad.Inty = downs.Inty .* 2 .* pi ./ 10;
 
 %% 2) Unwrap
 
 unwrapped.angularPosition = unwrap(downsRad.angularPosition);
 unwrapped.Intx = unwrap(downsRad.Intx);
 unwrapped.Inty = unwrap(downsRad.Inty);
 
    
 %% (3) Low pass filter: remove any potentially artifactual high frequency data
    
 smoothed.angularPosition = lowPassFilter(unwrapped.angularPosition, 25, sampleRate);
 smoothed.Intx = lowPassFilter(unwrapped.Intx, 25, sampleRate);
 smoothed.Inty = lowPassFilter(unwrapped.Inty, 25, sampleRate);

%% (5) Transform to degrees or mm

 deg.angularPosition = (smoothed.angularPosition / (2*pi)) * 360;
 deg.xPosition = smoothed.Intx * 4.5; %4.5 is the radius of the ball I use
 deg.yPosition = smoothed.Inty * 4.5;

%% (6) Take the derivative to get the velocity

 smoothed.angularVel = gradient(deg.angularPosition).* 25; %I multiply the derivative by 25 to get the velocity is relevant units (mm/s)
 smoothed.xVel = gradient(deg.xPosition).* 25;
 smoothed.yVel = gradient(deg.yPosition).* 25;


end
