%looking at the velocity for the full experiment


%raw data
figure,
plot(abs(rawData.trial_bdata(:,xPanels)-5));
hold on
plot(rawData.trial_bdata(:,OptoTrigger));
    allExpVel = rawData.trial_bdata(:,xFly);
    downsRad = allExpVel .* 2 .* pi ./ 10;
    unwrapped = unwrap(downsRad);
    smoothedVel = smoothdata(unwrapped,'rlowess',25);
    degx = smoothedVel * 4.5; 
    diffx = gradient(degx).* 25;
    smoothedxVel = smoothdata(diffx,'rlowess',15);
plot(smoothedxVel)
legend('x panels', 'opto trigger', 'fwd velocity');



%downsampling it
sampleRate = 4000;

figure,
plot(downsample(abs(rawData.trial_bdata(:,xPanels)-5),sampleRate/25));
hold on
plot(downsample(rawData.trial_bdata(:,OptoTrigger),sampleRate/25));
    allExpVel = rawData.trial_bdata(:,xFly);
    downs = downsample(allExpVel,sampleRate/25);
    downsRad = downs .* 2 .* pi ./ 10;
    unwrapped = unwrap(downsRad);
    smoothedVel = smoothdata(unwrapped,'rlowess',25);
    degx = smoothedVel * 4.5; 
    diffx = gradient(degx).* 25;
    smoothedxVel = smoothdata(diffx,'rlowess',15);
plot(smoothedxVel)
legend('x panels', 'opto trigger', 'fwd velocity');


%looking at the individual trips
for i = 1:length(tripData)
    figure,
    plot(downsample(tripData{1,i}(:,xPanels),sampleRate/25))
    hold on
    plot(downsample(tripData{1,i}(:,OptoTrigger),sampleRate/25))
    plot(smoothed{1,i}.xVel)
end

%% Trying different velocity analysis to see how they turn out.

% Steps that need to be taken for sure:
%1) change from V to rad.
%2) unwrap
%3) either filter or downsample
%4) change the position to mm using scaling factor
%5) calculate the velocity using the position data

% Some decisions to be made
%1) Will we use filtering or downsampling, and if filtering, what should be
%the cutoff frequency.
%2) Are we also going to smooth the data, if so, how many times and with
%what method?
%3) Are we using diff or gradient for the velocity calculation?
%4) Do we want to exclude extreme velocities? Which I currently do.


%Plotting output of current method
allXPos = rawData.trial_bdata(2.8E5:3.18E5,xFly);
totalTime = linspace(1, length(allXPos)/sampleRate,length(allXPos));

figure('Position',[100 100 1800 1000]),
subplot(3,4,1)
plot(totalTime, abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime, rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,allXPos)
ylabel('Voltage'); xlabel('Time (sec)');
legend('x panels', 'opto trigger', 'x data');
title('Raw x position voltage data');
ylim([-1 10]);

%1) First step: downsampling
downsampledIntx = downsample(allXPos,sampleRate/25); %For a 1000 rate acquisition frame rate from the NiDaq, downsampling to 25 Hz equals taking 1 every 40 frames
downsTime = linspace(1, length(allXPos)/sampleRate, length(downsampledIntx));

subplot(3,4,2),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,downsampledIntx)
ylabel('Voltage'); xlabel('Time (sec)');
title('Downsampled x position voltage data');
ylim([-1 10]);

%2) Second step: transform signal to radians
downsRadIntx = downsampledIntx .* 2 .* pi ./ 10; %10 is for the max voltage outputed by the daq

subplot(3,4,3),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,downsRadIntx)
ylabel('Rad'); xlabel('Time (sec)');
title('Downsampled x position in rad');
ylim([-1 10]);

%3) Third step: unwrapping
unwrappedIntx = unwrap(downsRadIntx);

subplot(3,4,5),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,unwrappedIntx)
ylabel('Rad'); xlabel('Time (sec)');
title('Unrapped x position');
ylim([-1 10]);

%4) Fourth step: smoothing
smoothedIntx = smoothdata(unwrappedIntx,'rlowess',25);

subplot(3,4,6),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,smoothedIntx)
ylabel('Rad'); xlabel('Time (sec)');
title('Smoothed x position');
ylim([-1 10]);

%5) Fifth step: transform to mm
degIntx = smoothedIntx * 4.5;

subplot(3,4,7),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,degIntx)
ylabel('mm'); xlabel('Time (sec)');
title('Smoothed x position in mm');
ylim([-1 40]);

%6) Sixth step: take the derivative
diffIntx = gradient(degIntx).* 25; %we multiply by 25 because we have downsampled to 25 Hz

subplot(3,4,9),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,diffIntx)
ylabel('mm/s'); xlabel('Time (sec)');
title('x Vel in mm');
ylim([-1 10]);

%7) Seventh step: remove values that are below 2.5% and above 97.5%
percentile25FV = prctile(diffIntx,2.5);
percentile975FV = prctile(diffIntx,97.5);
boundedDiffIntx = diffIntx;
boundedDiffIntx(boundedDiffIntx<percentile25FV | boundedDiffIntx>percentile975FV) = NaN;

subplot(3,4,10),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,boundedDiffIntx)
ylabel('mm/s'); xlabel('Time (sec)');
title('x Vel in mm without extreme values');
ylim([-1 10]);

%8) Eigth step: linearly interpolate to replace the missing values
[pointsVectorFV] = find(~isnan(boundedDiffIntx));
valuesVectorFV = boundedDiffIntx(pointsVectorFV);
xiFV = 1:length(boundedDiffIntx);
interpxVel = interp1(pointsVectorFV,valuesVectorFV,xiFV);

subplot(3,4,11),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,interpxVel)
ylabel('mm/s'); xlabel('Time (sec)');
title('x Vel in mm with interpolated values');
ylim([-1 10]);

%9) Last step: smoothing again
smoothedxVel = smoothdata(interpxVel,'rlowess',15);

subplot(3,4,12),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(downsTime,smoothedxVel)
ylabel('mm/s'); xlabel('Time (sec)');
title('Smoothed x vel');
ylim([-1 10]);


%% Overlaying the different velocities

figure,
plot(diffIntx)
hold on
plot(boundedDiffIntx)
plot(interpxVel)
plot(smoothedxVel)
legend('x Vel','no extremes','interp values','smoothed')
xlabel('Time'); ylabel('Vel (mm/s)');
title('Different steps in vel calculation');


%% Diff vs gradient

diffx = diff(degIntx).* 25;

figure,
plot(diffx)
hold on
plot(diffIntx)
legend('diff function', 'gradient');
xlabel('Time'); ylabel('Vel (mm/s)');
title('Different ways to get the velocity');

%% Using Yvette's protocol

figure('Position',[100 100 1800 1000]),
subplot(3,3,1)
plot(totalTime, abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime, rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,allXPos)
ylabel('Voltage'); xlabel('Time (sec)');
legend('x panels', 'opto trigger', 'x data');
title('Raw x position voltage data');
ylim([-1 10]);

%1) First step: transform signal to radians
RadIntx = allXPos .* 2 .* pi ./ 10;

subplot(3,3,2),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,RadIntx)
ylabel('Rad'); xlabel('Time (sec)');
title('X position in rad');
ylim([-1 10]);

%2) Second step: unwrapping
unwrappedIntx = unwrap(RadIntx);

subplot(3,3,3),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,unwrappedIntx)
ylabel('Rad'); xlabel('Time (sec)');
title('Unrapped x position');
ylim([-1 10]);

% Third step: find indexes where the unwrapping happened and put NaNs
upwrappedIndexes = find (abs(diff(RadIntx)) > pi); 
NUM_SAMPLES_FROM_WRAP_TO_REPLACE = 2;
upwrappedIndexes = upwrappedIndexes( upwrappedIndexes > NUM_SAMPLES_FROM_WRAP_TO_REPLACE & upwrappedIndexes < (length (unwrappedIntx) - NUM_SAMPLES_FROM_WRAP_TO_REPLACE) ); 
cleanedPos = unwrappedIntx;
for i = 1: length ( upwrappedIndexes )
    index_start = upwrappedIndexes(i) -  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ; 
    index_end = upwrappedIndexes(i) +  NUM_SAMPLES_FROM_WRAP_TO_REPLACE ;    
    cleanedPos( index_start : index_end ) = NaN;
end

subplot(3,3,4),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,cleanedPos)
ylabel('Rad'); xlabel('Time (sec)');
title('Unrapped x position');
ylim([-1 10]);

% Fourth step: replace NaN values with the last preceding value that was a real number
nanIDX = find(isnan(cleanedPos)); % find NaN indexes
% replace with preceeding value
while(~isempty(nanIDX))
    cleanedPos(nanIDX) = cleanedPos(nanIDX - 1); 
    nanIDX  = find(isnan(cleanedPos));
end

subplot(3,3,5),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,cleanedPos)
ylabel('Rad'); xlabel('Time (sec)');
title('Unrapped x position');
ylim([-1 10]);

% Fifth step: low pass filter the position array (she uses 25 hz)
filteredPosition = lowPassFilter(cleanedPos, 25, sampleRate);

subplot(3,3,6),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,filteredPosition)
ylabel('Rad'); xlabel('Time (sec)');
title('Filtered x position');
ylim([-1 10]);
 
% Sixth step: low pass filter the position array again to be more aggressize
filteredPosition = lowPassFilter(filteredPosition, 25, sampleRate);

subplot(3,3,7),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,filteredPosition)
ylabel('Rad'); xlabel('Time (sec)');
title('Filtered x position');
ylim([-1 10]);

% Seventh step: transform from radians into mm
accumulatedPositionOut = filteredPosition * 4.5;

subplot(3,3,8),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,accumulatedPositionOut)
ylabel('mm'); xlabel('Time (sec)');
title('Filtered x position');
ylim([-1 40]);

% take derivative and ajust for sample rate to solve for deg/s
velocityOut = gradient(accumulatedPositionOut) .* sampleRate ; 

subplot(3,3,9),
plot(totalTime,abs(rawData.trial_bdata(2.8E5:3.18E5,xPanels)-5));
hold on
plot(totalTime,rawData.trial_bdata(2.8E5:3.18E5,OptoTrigger));
plot(totalTime,velocityOut)
ylabel('Vel (mm/s)'); xlabel('Time (sec)');
title('Velocity');
ylim([-1 10]);

%% Overlaying Yvette's velocity and mine

figure,
plot(totalTime,velocityOut)
hold on
plot(downsTime,diffIntx)

%adding a smoothed version of the data with Yvette's protocol
smoothingData = smoothdata(velocityOut,'rlowess',50);
plot(totalTime,smoothingData)
legend('Yvette protocol', 'my protocol', 'smoothing Yvette protocol')
xlabel('Time'); ylabel('Velocity (mm/s');
