function [smoothed] = singleTrialVelocityAnalysis(data, sampleRate)
% This function processes the raw position data obtained from FicTrac to
% give back the velocity in degrees.

% It takes the rawData (data) and the sample rate (sampleRate) as inputs
% and gives back a struct names "smoothed" with the 3 final velocities for
% x, y and angular velocity.


%% Tranform signal from voltage to radians for unwrapping

    dataRad_Intx = data.ficTracIntx .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq
    dataRad_Inty = data.ficTracInty .* 2 .* pi ./ 10;
    dataRad_angularPosition = data.ficTracAngularPosition .* 2 .* pi ./ 10;
    dataRad_angularPosition_uncorrected = data.ficTracAngularPosition_uncorrected .* 2 .* pi ./ 10;
    

% Now the position is going between 0 and 2 pi.

%% Unwrapping 

    unwrapped_Intx = unwrap(dataRad_Intx);
    unwrapped_Inty = unwrap(dataRad_Inty);
    unwrapped_angularPosition = unwrap(dataRad_angularPosition);
    unwrapped_angularPosition_uncorrected = unwrap(dataRad_angularPosition_uncorrected);
% Now the position is unwrapped, so it doesn't jump when moving from 0 to
% 2pi and vice versa


%% Downsample the position data to match FicTrac's output


    % Downsample to match FicTrac's output
    downsampled_Intx = downsample(unwrapped_Intx,160); 
    downsampled_Inty = downsample(unwrapped_Inty,160); %For a 4000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 160 frames  
    downsampled_angularPosition = downsample(unwrapped_angularPosition,160);
    downsampled_angularPosition_uncorrected = downsample(unwrapped_angularPosition_uncorrected,160);
    
        
% The output is downsampled. It isn't noticeable when plotting solid lines, and
% it is barely noticeable when plotting dotted lines.
%% Smooth the data

    smoothed.Intx = smoothdata(downsampled_Intx,'rlowess',25); 
    smoothed.Inty = smoothdata(downsampled_Inty,'rlowess',25); 
    smoothed.angularPosition = smoothdata(downsampled_angularPosition,'rlowess',25);
    smoothed.angularPosition_uncorrected = smoothdata(downsampled_angularPosition_uncorrected,'rlowess',25);    
    
     
%% Transform to useful systems 
    
    Intx_mm = smoothed.Intx * 4.5; % wer tranform the pos to mm by scaling the value by the sphere's radius
    Inty_mm = smoothed.Inty * 4.5;
    angularPosition_deg = rad2deg(smoothed.angularPosition); % we transform the angular position to degrees
    %subtract the initial angle from the set up (the positon 1 is actually
    %3 pixels to the right of the fly)
    angularPosition_deg_uncorrected = rad2deg(smoothed.angularPosition_uncorrected);
    initialPx = 3; %my x=1 positon is 3 pixels to the right of the animal
    initialAngle = initialPx*360/96;
    %shift according to the panels pos
    shifted_data = angularPosition_deg - initialAngle;
    %wrap
    angularPosition_deg = wrapTo360(shifted_data);
    smoothed.angularPosition = deg2rad(angularPosition_deg);
    
    %add position for x and y in degrees to later compute total movement
    Intx_deg = rad2deg(smoothed.Intx);
    Inty_deg = rad2deg(smoothed.Inty);    
    
%% Take the derivative

    diff_Intx = gradient(Intx_mm).* 25; %we multiply by 25 because we have downsampled to 25 Hz
    diff_Inty = gradient(Inty_mm).* 25; 
    diff_angularPosition = gradient(angularPosition_deg_uncorrected).* 25;
    
    %%add the x and y in deg/s
    diff_Intx_deg = gradient(Intx_deg).* 25; %we multiply by 25 because we have downsampled to 25 Hz
    diff_Inty_deg = gradient(Inty_deg).* 25;

%% Calculate the distribution and take away values that are below 2.5% and above 97.5%
    
    percentile25AV = prctile(diff_angularPosition,2.5);
    percentile975AV = prctile(diff_angularPosition,97.5);
    boundedDiffAngularPos = diff_angularPosition;
    boundedDiffAngularPos(diff_angularPosition<percentile25AV | diff_angularPosition>percentile975AV) = NaN;
    
    percentile25FV = prctile(diff_Intx,2.5);
    percentile975FV = prctile(diff_Intx,97.5);
    boundedDiffIntx = diff_Intx;
    boundedDiffIntx(boundedDiffIntx<percentile25FV | boundedDiffIntx>percentile975FV) = NaN;
    
    percentile25SV = prctile(diff_Inty,2.5);
    percentile975SV = prctile(diff_Inty,97.5);
    boundedDiffInty = diff_Inty;
    boundedDiffInty(boundedDiffInty<percentile25SV | boundedDiffInty>percentile975SV) = NaN;
    
    %add the deg/s components
    percentile25FVDeg = prctile(diff_Intx_deg,2.5);
    percentile975FVDeg = prctile(diff_Intx_deg,97.5);
    boundedDiffIntxDeg = diff_Intx_deg;
    boundedDiffIntxDeg(boundedDiffIntxDeg<percentile25FVDeg | boundedDiffIntxDeg>percentile975FVDeg) = NaN;
    
    percentile25SVDeg = prctile(diff_Inty_deg,2.5);
    percentile975SVDeg = prctile(diff_Inty_deg,97.5);
    boundedDiffIntyDeg = diff_Inty_deg;
    boundedDiffIntyDeg(boundedDiffIntyDeg<percentile25SVDeg | boundedDiffIntyDeg>percentile975SVDeg) = NaN;

 %% Linearly interpolate to replace the NaNs with values.
 
    [pointsVectorAV] = find(~isnan(boundedDiffAngularPos));
    valuesVectorAV = boundedDiffAngularPos(pointsVectorAV);
    xiAV = 1:length(boundedDiffAngularPos);
    interpAngVel = interp1(pointsVectorAV,valuesVectorAV,xiAV);
    
    [pointsVectorFV] = find(~isnan(boundedDiffIntx));
    valuesVectorFV = boundedDiffIntx(pointsVectorFV);
    xiFV = 1:length(boundedDiffIntx);
    interpxVel = interp1(pointsVectorFV,valuesVectorFV,xiFV);
    
    [pointsVectorSV] = find(~isnan(boundedDiffInty));
    valuesVectorSV = boundedDiffInty(pointsVectorSV);
    xiSV = 1:length(boundedDiffInty);
    interpyVel = interp1(pointsVectorSV,valuesVectorSV,xiSV);
    
    %for the deg/s velocities
    [pointsVectorFVDeg] = find(~isnan(boundedDiffIntxDeg));
    valuesVectorFVDeg = boundedDiffIntxDeg(pointsVectorFVDeg);
    xiFVDeg = 1:length(boundedDiffIntxDeg);
    interpxVelDeg = interp1(pointsVectorFVDeg,valuesVectorFVDeg,xiFVDeg);
    
    [pointsVectorSVDeg] = find(~isnan(boundedDiffIntyDeg));
    valuesVectorSVDeg = boundedDiffIntyDeg(pointsVectorSVDeg);
    xiSVDeg = 1:length(boundedDiffIntyDeg);
    interpyVelDeg = interp1(pointsVectorSVDeg,valuesVectorSVDeg,xiSVDeg);
       
 %%  Smooth again
 
    smoothed.xVel = smoothdata(interpxVel,'rlowess',15);
    smoothed.yVel = smoothdata(interpyVel,'rlowess',15);
    smoothed.angularVel = smoothdata(interpAngVel,'rlowess',15);
    % add the deg/s components
    smoothed.xVelDeg = smoothdata(interpxVelDeg,'rlowess',15);
    smoothed.yVelDeg = smoothdata(interpyVelDeg,'rlowess',15);
    %add total movement
    smoothed.total_mvt = abs(smoothed.angularVel)+abs(smoothed.xVelDeg)+abs(smoothed.yVelDeg);


end