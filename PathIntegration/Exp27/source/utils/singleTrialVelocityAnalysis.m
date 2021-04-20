function [smoothed] = singleTrialVelocityAnalysis(data, sampleRate)
% This function processes the raw position data obtained from FicTrac to
% give back the velocity in degrees.

% It takes the rawData (data) and the sample rate (sampleRate) as inputs
% and gives back a struct names "smoothed" with the 3 final velocities for
% x, y and angular velocity.


%% Tranform signal from voltage to radians for unwrapping

    dataRad.Intx = data.ficTracIntx .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
    dataRad.angularPosition = data.ficTracAngularPosition .* 2 .* pi ./ 10;
    
% Now the position is going between 0 and 2 pi.

%% Unwrapping 

    unwrapped.Intx = unwrap(dataRad.Intx);
    unwrapped.angularPosition = unwrap(dataRad.angularPosition);

% Now the position is unwrapped, so it doesn't jump when moving from 0 to
% 2pi and vice versa
%% Downsample the position data to match FicTrac's output


    % Downsample to match FicTrac's output
    downsampled.Intx = resample(unwrapped.Intx,25,sampleRate); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames
    downsampled.angularPosition = resample(unwrapped.angularPosition,25,sampleRate);
          
% The output is downsampled. It isn't noticeable when plotting solid lines, and
% it is barely noticeable when plotting dotted lines.


%% Smooth the data

    smoothed.Intx = smoothdata(downsampled.Intx,'rlowess',25); %to 25 Hz
    %smoothed.Inty = smoothdata(unwrapped.Inty,'rlowess',25); 
    smoothed.angularPosition = smoothdata(downsampled.angularPosition,'rlowess',25);
    
     
%% Transform to useful systems 
    
    deg.Intx = smoothed.Intx * 4.5; % wer tranform the pos to mm by scaling the value by the sphere's radius
    %deg.Inty = smoothed.Inty * 4.5;
    deg.angularPosition = (smoothed.angularPosition / (2*pi)) * 360; % we transform the angular position to degrees
    
%% Take the derivative

    diffe.Intx = gradient(deg.Intx).* 25; %we multiply by 25 because we have downsampled to 25 Hz
    %diff.Inty = gradient(deg.Inty).* 25; 
    diffe.angularPosition = gradient(deg.angularPosition).* 25;

%% Calculate the distribution and take away values that are below 2.5% and above 97.5%
    
    percentile25AV = prctile(diffe.angularPosition,2.5);
    percentile975AV = prctile(diffe.angularPosition,97.5);
    boundedDiffAngularPos = diffe.angularPosition;
    boundedDiffAngularPos(diffe.angularPosition<percentile25AV | diffe.angularPosition>percentile975AV) = NaN;
    
    percentile25FV = prctile(diffe.Intx,2.5);
    percentile975FV = prctile(diffe.Intx,97.5);
    boundedDiffIntx = diffe.Intx;
    boundedDiffIntx(boundedDiffIntx<percentile25FV | boundedDiffIntx>percentile975FV) = NaN;
    
%     percentile25SV = prctile(diff.Inty,2.5);
%     percentile975SV = prctile(diff.Inty,97.5);
%     boundedDiffInty = diff.Inty;
%     boundedDiffInty(boundedDiffInty<percentile25SV | boundedDiffInty>percentile975SV) = NaN;
%     
    

 %% Linearly interpolate to replace the NaNs with values.
 
    [pointsVectorAV] = find(~isnan(boundedDiffAngularPos));
    valuesVectorAV = boundedDiffAngularPos(pointsVectorAV);
    xiAV = 1:length(boundedDiffAngularPos);
    interpAngVel = interp1(pointsVectorAV,valuesVectorAV,xiAV);
    
    [pointsVectorFV] = find(~isnan(boundedDiffIntx));
    valuesVectorFV = boundedDiffIntx(pointsVectorFV);
    xiFV = 1:length(boundedDiffIntx);
    interpxVel = interp1(pointsVectorFV,valuesVectorFV,xiFV);
    
%     [pointsVectorSV] = find(~isnan(boundedDiffInty));
%     valuesVectorSV = boundedDiffInty(pointsVectorSV);
%     xiSV = 1:length(boundedDiffInty);
%     interpyVel = interp1(pointsVectorSV,valuesVectorSV,xiSV);
       
 %%  Smooth again
 
    %smoothed.xVel = smoothdata(diff.Intx,'rlowess',15); 
    smoothed.xVel = smoothdata(interpxVel,'rlowess',15);
    %smoothed.yVel = smoothdata(interpyVel,'rlowess',15);
    smoothed.angularVel = smoothdata(interpAngVel,'rlowess',15);


end