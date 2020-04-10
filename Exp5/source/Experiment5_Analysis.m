% Experiment 5 analysis

%%% This code analyses the output of the panels and the FicTrac data
%for experiment 5, in which the fly gets thre blocks of bar jumps of
%different error level (low-high-low)
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\Experiments\Exp5\data'
[file,path] = uigetfile();
load([path,file]);

rawData = daq_data'; %transpose matrix for future use

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;
PanelStatus = 6; %this signal tells whether the panels are on or off.

%Define the error block we're in according to the file name
if findstr('Low',file) == 1
    error = 1; %(low)
else
    error = 2; %(high)
end
    
%% Subset acquisition of x and y pos, as well as FicTrac data

data.xpanelVolts =  rawData (:,xPanels); 
VOLTAGE_RANGE_x = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)

data.yPanelVolts =  rawData (:, yPanels);
VOLTAGE_RANGE_y = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
maxValY = 96;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2

%FicTrac data
data.fictracAngularPosition = rawData ( : , headingFly); 
%I'm negating the x data because my axes were inverted in fictrac
data.ficTracIntx = -rawData ( : , xFly); 
data.ficTracInty = rawData ( : , yFly); 

%% Determine bar jumps
%This is going to be done in two ways:
%1) using the signal from channel 6 to determine when the panels were on
%and off, and adding jumps every 200 sec in between
%2) taking the derivative of the voltage signal from the ypanels channel to identify jumps

% Plot Panel acquisition time
figure,
subplot(2,1,1)
plot(rawData(:,6))
ylabel('Voltage signal (V)');
title('Panels ON and OFF');
xlim([0 size(daq_data,2)]);
ylim([-0.1 10.1]);

% Define when the panels turn on and off
%take the derivative of the panels' signal
panelON = diff(rawData(:,6));
subplot(2,1,2)
plot(panelON)
xlim([0 size(daq_data,2)]);
xlabel('Time (frames)');
ylabel('Diff of voltage signal (V)');
%Find the frame when the panels turn on
[M, I] = max(panelON);
%Find the frame when they turn off.
[M2, I2] = min(panelON);
startFrame = I(1);
endFrame = I2(1);
%add text indicating start and end frames
text(startFrame+35000,5,strcat('startFrame',num2str(startFrame)))
text(endFrame-100000,5,strcat('endFrame',num2str(endFrame)))


% get the jumps using the signal from the yPanels channel
Jumps = diff(data.yPanelVolts);
Jumps(abs(Jumps)>0.4 & abs(Jumps)<1)=1;
Jumps = round(Jumps);

figure,
suptitle('Bar jumps');
subplot(3,1,1)
plot(data.yPanelVolts)
ylabel('Voltage (V)');xlabel('Time');
subplot(3,1,2)
plot(Jumps);
ylabel('Voltage difference (V)');xlabel('Time');

j = find(Jumps); %indices of the actual bar jumps, taken from the y signal
if error == 2
 j = j(2:end);
end
jsec = j/1000;

%plot the data from the yPanels and add lines of the previously determined
%bar jumps
%plot the panels y dimension signal
subplot(3,1,3), plot(data.yPanelVolts)
title('Bar jumps');
xlabel('Time (frames)'); ylabel('Voltage (V)');
hold on
%add the bar jumps
for i = 1:length(j)
    plot([j(i) j(i)],[0 10],'r');
end

%% Fixing the data relative to the bar jumps

%The x position of the panels and the FicTrac output are "ignorant" of the
%bar jumps. The x position of the bar will move with the angular position
%of the fly, but the coordinate system changes every time the bar jumps.
%We need to make sure the xpos and heading of the fly are corrected to take
%this coordinate change into account.

yVoltsBJ = data.yPanelVolts(j-1);
yVoltsAJ = data.yPanelVolts(j+1);
ydiff = yVoltsAJ-yVoltsBJ; %this is the offset that I need to adjust by

%Correct the data after every jump except for the last
data.xPanelVoltsUW = data.xpanelVolts;
data.ficTracAngularPositionUW = data.fictracAngularPosition;
for i = 1:size(j)-1
   data.xPanelVoltsUW(j(i)+1:j(i+1)) = (data.xpanelVolts(j(i)+1:j(i+1)))+sum(ydiff(1:i));  
   data.ficTracAngularPositionUW(j(i)+1:j(i+1)) = (data.fictracAngularPosition(j(i)+1:j(i+1)))+sum(ydiff(1:i));
end

%Correct the data after the last jump
data.xPanelVoltsUW(j(end)+1:end) = data.xpanelVolts(j(end)+1:end)+sum(ydiff(1:end));
data.ficTracAngularPositionUW(j(end)+1:end) = data.fictracAngularPosition(j(end)+1:end)+sum(ydiff(1:end));

%I now have to wrap this data to get it to be between 0 and 10 V. 
data.xPanelVolts = data.xPanelVoltsUW;
data.ficTracAngularPosition = data.ficTracAngularPositionUW;

for i = 1:size(data.xPanelVolts)
    
    if data.xPanelVoltsUW(i) > 10
    data.xPanelVolts(i) = data.xPanelVoltsUW(i)-10;
    else
    data.xPanelVolts(i) = data.xPanelVoltsUW(i);
    end
    
    if data.ficTracAngularPositionUW(i) > 10
    data.ficTracAngularPosition(i) = data.ficTracAngularPositionUW(i)-10;
    else
    data.ficTracAngularPosition(i) = data.ficTracAngularPositionUW(i);
    end

end

% Getting the data in x and y position from the voltage info.
data.xPanelPos = round ((data.xPanelVolts  * maxValX ) /VOLTAGE_RANGE_x); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE_y);

%check the wrapping graphically
figure('Position', [100 100 1200 900]),
subplot(2,1,1),
plot(data.yPanelVolts,'r')
hold on
plot(data.xPanelVoltsUW,'b')
plot(data.xPanelVolts,'g')
ylabel('Voltage (V)');
title('x panels voltage signals before and after wrapping');
ylim([0 max(data.xPanelVoltsUW)]);
legend({'yPanels','Unwrapped xPanels', 'Wrapped xPanels'});

subplot(2,1,2),
plot(data.yPanelVolts,'r')
hold on
plot(data.ficTracAngularPositionUW,'b')
plot(data.ficTracAngularPosition,'g')
xlabel('Time (frame)');ylabel('Voltage (V)');
title('angular position voltage signals before and after wrapping');
ylim([0 max(data.ficTracAngularPositionUW)]);
legend({'yPanels','Unwrapped angular position', 'Wrapped angular position'});


% Getting the degrees that each jump represents and checking they are
% correct
figure('Position', [100 100 1200 900]),

jumpPos = data.yPanelPos(j+1)-data.yPanelPos(j-1);
degJumps = wrapTo180(jumpPos*(360/96));
subplot(1,3,1), plot(degJumps,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from yPanels');
ylabel('deg');xlabel('Trial #');xlim([1 TrialNum]);
% compare that to the jump function we have stored
hold on
%Define the jumps according to the error level
if error == 1
    jumps = jumps;
elseif error == 2
    jumps = [-90,-90,-90,90,90,90,90,90,-90,90,90,90,90,-90,-90,-90,90,90,90,90,90,-90,-90,90];
end
plot(jumps(1:TrialNum),'b')
legend({'Jumps from Y data','Jump function used'});

%Check if the jump magnitude appears ok in the x panel data
jumpMag = data.xPanelPos(j+1)-data.xPanelPos(j-1); 
degMag = wrapTo180(jumpMag*(360/96));
subplot(1,3,2), plot(degMag,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from xPanels');
ylabel('deg');xlabel('Trial #'); xlim([1 TrialNum]);
% compare that to the jump function we have stored
hold on
plot(jumps(1:TrialNum),'b')
legend({'Jumps from X data','Jump function used'});

%check if the jump magnitude appears ok in the angular position data
jumpMag2 = data.ficTracAngularPosition(j+1)-data.ficTracAngularPosition(j-1); 
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180
subplot(1,3,3), plot(degMag2,'ro')
ylim([-180 180]);
title('Jump magnitude, taken from angular position');
ylabel('deg');xlabel('Trial #');xlim([1 TrialNum]);
% compare that to the jump function we have stored
hold on
plot(jumps(1:TrialNum),'b')
legend({'Jumps from angular position','Jump function used'});

%% Downsample, unwrap and smooth position data, then get velocity and smooth

sizeBall = 9;
sampleRate = 1000;
[smoothed] = singleTrialVelocityAnalysis9mm(data,sampleRate);

%% 2D trajectories

%Uncorrected trajectory (without correcting the heading). This is to
%visualize the fly's response to the jumps, when they do very well.
[posx,posy]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.AngularPosition);

time = linspace(1,2400,length(posx));

figure,
subplot(1,2,1)
scatter(posx,posy,0.5,time)
h = colorbar
title(h,'Time (s)');
axis equal
axis tight
title('Uncorrected 2D trajectory');
set(gca,'xtick',[])
set(gca,'ytick',[])


%Corrected trajectory: this is the actual trajectory the fly took during
%the block
[posx2,posy2]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.angularPosition);

subplot(1,2,2), scatter(posx2,posy2,0.5,time);
h2 = colorbar
title(h2, 'Time (s)')
axis equal
axis tight
title('Corrected 2D trajectory');
set(gca,'xtick',[])
set(gca,'ytick',[])

saveas(gcf,strcat(path,'plots\Trajectories',file(1:end-4),'.png'));

%% Forward velocity analysis

% The forward velocity is a good indicative of whether the fly is walking
% well in that trial or not

%To look at the velocity, I should do so without adding the offsets of the
%jumps, because the offsets will give me weird jumps in velocity that the
%fly might not actually have made. For the forward velocity, I can use the
%one calculates by the singleTrialVelocityAnalysis function, because it
%uses the variable Intx, that I haven't corrected for the offset

forwardVelocity = smoothed.xVel;
meanVelocity = mean(forwardVelocity);
time = linspace(0,(length(rawData)/1000),length(forwardVelocity));

%Global forward velocity analysis
figure ('Position', [100 100 1200 900]),
subplot(2,1,1)
plot(time,forwardVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
ylim([min(forwardVelocity)-10 max(forwardVelocity)+10]);
hold on
%add the jumps
for i = 1:length(jsec)
     plot([jsec(i) jsec(i)],[min(forwardVelocity)-10 max(forwardVelocity)+10],'g');
end
hline = refline([0 meanVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title('Forward velocity of the fly', 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (mm/s)')

subplot(2,1,2)
histogram(forwardVelocity,'FaceColor',[.4 .2 .6])
title('Distribution of forward velocities');
xlabel('Forward velocity (mm/s)');
ylabel('Frequency');
xlim([min(forwardVelocity) max(forwardVelocity)]);

saveas(gcf,strcat(path,'plots\ForwardVelocity',file(1:end-4),'.png'));

%% Angular velocity

%To look at the velocity, I should do so without adding the offsets of the
%jumps, because the offsets will give me weird jumps in velocity that the
%fly might not actually have made. For the angular velocity, I'm using a new 
%function that I made to downsample, unwrap, smooth and get the velocity
%with the angular position uncorrected for the offset

AngularPosition = rawData ( : , headingFly);
angularVelocity = getAngVel(AngularPosition);

meanAngVelocity = mean(angularVelocity);
time = linspace(0,(length(rawData)/1000),length(angularVelocity));

figure('Position', [100 100 1200 900]),
subplot(2,1,1)
plot(time,angularVelocity,'k','HandleVisibility','off')
xlim([0 time(end)]);
ylim([min(angularVelocity)-50 max(angularVelocity)+50]);
hold on
for i = 1:length(j) %add the bar jumps
     plot([jsec(i) jsec(i)],[min(angularVelocity)-50 max(angularVelocity)+50],'g');
end
hline = refline([0 meanAngVelocity]);
hline.Color = 'r'; hline.LineStyle = '--';
rline = refline([0 0]);
rline.Color = [.5 .5 .5]; rline.LineWidth = 1.5;
title('Angular velocity of the fly', 'Interpreter', 'none');
xlabel('Time (s)')
ylabel('Velocity (deg/s)')

subplot(2,1,2)
histogram(angularVelocity,'FaceColor',[.2 .8 .6])
title('Distribution of angular velocities');
xlabel('Angular velocity (deg/s)');
ylabel('Frequency');
xlim([min(angularVelocity) max(angularVelocity)]);

saveas(gcf,strcat(path,'plots\SmoothedAngularVelocity',file(1:end-4),'.png'))


%%  Activity levels

% We are going to decide whether a fly is moving or not based on the
% forward velocity. If it's above 1 mm/s we will consider it is moving
% We will work with the downsampled data

downsampled.xPanelPos = downsample(data.xPanelPos,1000/25); %downsample the panels position
dataMoving.xPanelPos = downsampled.xPanelPos(forwardVelocity>1); %keep the position frames during which the fly moved
moving = smoothed.angularPosition(forwardVelocity>1); %keep the angular position frames during which the fly moved
% I need to think more carefully about whether there is a problem in using
% the forward velocity obtained from smoothed data to choose the frames
% during which the fly is or not moving in unsmoothed data...

percentageActivity = 100*size(moving)/size(smoothed.angularPosition);
activity = zeros(length(forwardVelocity),1);

for i = 1:length(forwardVelocity)
    if forwardVelocity(1,i) > 1
        activity(i,1) = 1;
    else
        activity(i,1) = 0;
    end
end

time = linspace(0,(length(rawData)/1000),length(activity));

figure,
set(gcf, 'Position', [500, 500, 1000, 100])
newMap = flipud(gray);
xaxis = time;
trials = [0:1];
imagesc(xaxis,trials,activity')
colormap(newMap)
title(strcat('Activity raster plot, percentage activity:', num2str(percentageActivity), '%'));
ylabel('Activity');
xlabel('Time (s)');

save(strcat(path,'dataFromAnalysis\Activity',file(1:end-4),'.mat'),'time','activity','percentageActivity');
saveas(gcf,strcat(path,'plots\ActivityRP',file(1:end-4),'.png'))

%% Output in degrees of the Panels position

pxToDeg = 360/96; % There are 96 possible positions and this represents 360 deg
posToDeg = zeros(1,length(downsampled.xPanelPos));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
for i=1:length(downsampled.xPanelPos)
    if downsampled.xPanelPos(i) == 70
        posToDeg(i) = 0;
    elseif downsampled.xPanelPos(i) >70 
        posToDeg(i) = (downsampled.xPanelPos(i)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        posToDeg(i) = (downsampled.xPanelPos(i)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

%% Probability density function of the stimulus position

% Remapping the positions to span -180 to 180 deg
remapPosToDeg = wrapTo180(posToDeg);

figure('Position', [100 100 1600 900]),

% Stimulus position in time, with every frame
time = linspace(0,(length(rawData)/1000),length(remapPosToDeg));
subplot(2,4,[1,5])
scatter(remapPosToDeg, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(remapPosToDeg, time,'k','HandleVisibility','off')
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [jsec(i) jsec(i)],'r--');
end
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
subplot(2,4,[4,8])
scatter(remapPosToDeg(forwardVelocity>1), time(forwardVelocity>1));
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [jsec(i) jsec(i)],'r--');
end
ax = gca;
ax.YDir = 'reverse'; 


% Plot the histogram and probability density
%With every frame
edges = [-180:20:180];
[counts] = histcounts(remapPosToDeg,edges);
probabilities = counts./sum(counts);
degs = linspace(-180,180,length(counts));

subplot(2,4,2)
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
xlim([-180 180]); ylim([0 max(probabilities)+0.05]);
title('Stimulus position, all frames');
ylabel('Probability density'); xlabel('Stimulus position (deg)');

%With only frames with velocity>1
[countsMoving] = histcounts(remapPosToDeg(forwardVelocity>1),edges);
probabilitiesMoving = countsMoving./sum(countsMoving);
degsMoving = linspace(-180,180,length(countsMoving));

subplot(2,4,3)
plot(degsMoving,probabilitiesMoving,'k')
set(0, 'DefaulttextInterpreter', 'none')
xlim([-180 180]); ylim([0 max(probabilitiesMoving)+0.05]);
title('Stimulus position, only moving frames');
ylabel('Probability density'); xlabel('Stimulus position (deg)');


% Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);
% some statistics...
CircularStats = circ_stats(posToRad);
[pval,z] = circ_rtest(posToRad);
circLength = circ_r(posToRad,[],2);

%Plot the histogram in polar coordinates
circedges = [0:20:360];
circedges = deg2rad(circedges);
subplot(2,4,6)
polarhistogram(posToRad,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

subplot(2,4,7)
polarhistogram(posToRad(forwardVelocity>1),circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top';

saveas(gcf,strcat(path,'plots\BarPosition',file(1:end-4),'.png'))

%% Fly's heading thoughout the experiment

%I think for the fly's heading I don't need to remap anything, cause it
%should be based on fictrac's calibration and have the front be 0?

flyPosPx = rad2deg(smoothed.angularPosition)/pxToDeg;
FlyPosDeg = zeros(1,length(smoothed.angularPosition));

% Convert from xpos to degrees, knowing that xpos 70 = 0 deg
for i=1:length(smoothed.angularPosition)
    if flyPosPx(i) == 70
        FlyPosDeg(i) = 0;
    elseif flyPosPx(i) >70 
        FlyPosDeg(i) = (flyPosPx(i)-70)*pxToDeg; % Correct the offset and multiply by factor to get deg
    else
        FlyPosDeg(i) = (flyPosPx(i)+27)*pxToDeg; % Correct the offset and multiply by factor to get deg
    end
end

FlyPos360 = wrapTo360(FlyPosDeg);
flyPos180 = wrapTo180(FlyPos360);

figure('Position', [100 100 1600 900]),

% Heading in time, with every frame
time = linspace(0,(length(rawData)/1000),length(flyPos180));
subplot(2,4,[1,5])
scatter(flyPos180, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(flyPos180, time,'k','HandleVisibility','off')
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [jsec(i) jsec(i)],'r--');
end
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
%time = linspace(0,(length(rawData)/1000),length(flyPos180(forwardVelocity>1)));
subplot(2,4,[4,8])
scatter(flyPos180(forwardVelocity>1), time(forwardVelocity>1));
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [jsec(i) jsec(i)],'r--');
end
ax = gca;
ax.YDir = 'reverse'; 
 

% Plot the histogram and probability density
edges = [-180:20:180];
[countsFly] = histcounts(flyPos180,edges);
[countsFlyMoving] = histcounts(flyPos180(forwardVelocity>1),edges);
probabilitiesFly = countsFly./sum(countsFly);
probabilitiesFlyMoving = countsFlyMoving./sum(countsFlyMoving);
degsFly = linspace(-180,180,length(countsFly));
degsFlyMoving = linspace(-180,180,length(countsFlyMoving));

subplot(2,4,2) %with every frame
plot(degsFly,probabilitiesFly,'k')
xlim([-180 180]); ylim([0 max(probabilitiesFly)+0.05]);
title('Fly heading, every frame');
ylabel('Probability density'); xlabel('Fly heading (deg)');
subplot(2,4,3) %with moving frames only
plot(degsFlyMoving,probabilitiesFlyMoving,'k')
xlim([-180 180]); ylim([0 max(probabilitiesFlyMoving)+0.05]);
title('Fly heading. only moving frames');
ylabel('Probability density'); xlabel('Fly heading (deg)');

% In polar coordinates...
%Taking every frame into account
posToRadFly = deg2rad(FlyPos360);
CircularStatsFly = circ_stats(posToRadFly);
circedges = [0:20:360];
circedges = deg2rad(circedges);
subplot(2,4,6)
polarhistogram(posToRadFly,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%with only moving frames
CircularStatsFlyMoving = circ_stats(posToRadFly(forwardVelocity>1));
subplot(2,4,7)
posToRadFlyMoving = posToRadFly(forwardVelocity>1);
polarhistogram(posToRadFlyMoving,circedges,'Normalization','probability','FaceColor',[1,0.2,0.7],'HandleVisibility','off');
ax = gca;
ax.ThetaDir='clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

save(strcat(path,'dataFromAnalysis\FlyPosition',file(1:end-4),'.mat'),'posToRadFlyMoving','circedges');
saveas(gcf,strcat(path,'plots\FlyPosition',file(1:end-4),'.png'))

%% Look at the goal and calculate the distance to it...

%Taking all the experiment
figure,
%with every frame
goal = circ_mean(posToRadFly,[],2);
dist2goal2 = circ_dist(posToRadFly,goal);
dist2goal = wrapTo180(rad2deg(dist2goal2));
[countsDist] = histcounts(dist2goal,edges);
probabilitiesDist = countsDist./sum(countsDist);
degsFlyDist = linspace(-180,180,length(countsDist));
subplot(1,2,1), plot(degsFlyDist,probabilitiesDist,'r')
title('Distance to the goal with every frame');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%taking only 'moving frames'
goalMoving = circ_mean(posToRadFly(forwardVelocity>1),[],2);
devGoalMoving = circ_std(posToRadFly(forwardVelocity>1),[],[],2);
dist2goalMoving2 = circ_dist(posToRadFly(forwardVelocity>1),goalMoving);
dist2goalMoving = wrapTo180(rad2deg(dist2goalMoving2));
[countsDistMoving] = histcounts(dist2goalMoving,edges);
probabilitiesDistMoving = countsDistMoving./sum(countsDistMoving);
degsFlyDistMoving = linspace(-180,180,length(countsDistMoving));
subplot(1,2,2), plot(degsFlyDistMoving,probabilitiesDistMoving,'r')
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

saveas(gcf,strcat(path,'plots\Dist2Goal',file(1:end-4),'.png'))
save(strcat(path,'dataFromAnalysis\goals',file(1:end-4),'.mat'),'goal','goalMoving','dist2goal','dist2goalMoving','degsFlyDistMoving','probabilitiesDistMoving','devGoalMoving');


%% Looking at the goal in 10 segments throughout the experiment

posToRadFlyMoving = posToRadFly(forwardVelocity>1);

stepSize = round(length(posToRadFlyMoving)/10)-1;

count = 1;
for i = 1:10
    
    goalMovingTrial(i) = circ_mean(posToRadFlyMoving(count:count+stepSize),[],2);
    goalMovingTrial(i) = wrapTo180(rad2deg(goalMovingTrial(i)));
    count = count + stepSize;

end

figure, plot(goalMovingTrial,'ko')

save(strcat(path,'dataFromAnalysis\goalTrial',file(1:end-4),'.mat'),'goalMovingTrial');

%% Using the 100 sec before and after each jump to compute several goal distributions
%per experiment...if it's a low error block

if error == 1
%I get the 100 sec before and after the jump for every jump except for the
%last
sec = 100;
sizeBall = 9;
j2 = j(1:end-1);
shortData = getDataAroundJump(rawData,j2,sec,sizeBall,sampleRate);
%apend jump data and save for using in the pooled flies analysis
shortData.jumpMag = jumps(1:end-1);
shortData.angPos = getPosAroundJump(data.ficTracAngularPosition,j(1:end-1),sec);


%Calculate and plot distribution of dist2goal around every jump
%with every frame
figure('Position', [100 100 1600 900]),
for i = 1:length(j)-1
    %For every frame
    goal100sec = circ_mean(deg2rad(shortData.angPos));
    dist2goal2100sec{i} = circ_dist(deg2rad(shortData.angPos(:,i)),goal100sec(i));
    dist2goal100sec{i} = wrapTo180(rad2deg(dist2goal2100sec{i}));
    [countsDist100sec{i}] = histcounts(dist2goal100sec{i},edges);
    probabilitiesDist100sec{i} = countsDist100sec{i}./sum(countsDist100sec{i});
    degsFlyDist100sec{i} = linspace(-180,180,length(countsDist100sec{i}));
    
    %For moving frames only
    datamoving{i} = shortData.angPos(:,i);
    datamoving{i}(shortData.forwardVel(:,i)<=1) = NaN;
    datamoving{i} = datamoving{i}(~isnan(datamoving{i}))';
    goal100secMoving{i} = circ_mean(deg2rad(datamoving{i}),[],2);
    dist2goal2100secMoving{i} = circ_dist(deg2rad(datamoving{i}),goal100secMoving{i});
    dist2goal100secMoving{i} = wrapTo180(rad2deg(dist2goal2100secMoving{i}));
    [countsDist100secMoving{i}] = histcounts(dist2goal100secMoving{i},edges);
    probabilitiesDist100secMoving{i} = countsDist100secMoving{i}./sum(countsDist100secMoving{i});
    degsFlyDist100secMoving{i} = linspace(-180,180,length(countsDist100secMoving{i}));
    
    %plot
    subplot(3,4,i),
    plot(degsFlyDist100sec{i},probabilitiesDist100sec{i},'r')
    hold on
    plot(degsFlyDist100secMoving{i},probabilitiesDist100secMoving{i},'k')
    xlim([-180, 180]); xlabel('Distance to the goal (deg)');
    ylabel('Probability');
end
    suptitle('Distance to the goal 100 sec around the jumps');
    
saveas(gcf,strcat(path,'plots\Dist2goal100sec',file(1:end-4),'.png'))

    
%Plot them using colorscales
probaDist100sec = cell2mat(probabilitiesDist100sec);
probaDist100sec = reshape(probaDist100sec,length(probaDist100sec)/length(probabilitiesDist100sec),length(probabilitiesDist100sec));
probaDist100secMoving = cell2mat(probabilitiesDist100secMoving);
probaDist100secMoving = reshape(probaDist100secMoving,length(probaDist100secMoving)/length(probabilitiesDist100secMoving),length(probabilitiesDist100secMoving));

%create new colormap
newMap = flipud(gray);
xaxis = [-180:360/17:180];
trials = [1:length(jumps)];
figure, imagesc(xaxis,trials,probaDist100secMoving')
colormap(newMap)
colorbar
saveas(gcf,strcat(path,'plots\HeatMapDist2goal100sec',file(1:end-4),'.png'))

%save data
save(strcat(path,'dataFromAnalysis\shortData',file(1:end-4),'.mat'),'shortData','probaDist100secMoving');

end

%% Using 10 sec around jump

%1) Use the 10 sec preceding each jump to compute the goal by taking
%the mean heading
sec = 10;
shortData10sec = getDataAroundJump(rawData,j,sec,sizeBall,sampleRate);
shortData10sec.jumpMag = jumps;
shortData10sec.angPos = getPosAroundJump(data.ficTracAngularPosition,j,sec);

for i = 1:length(j)
    goal10sec(i) = circ_mean(deg2rad(shortData10sec.angPos(1:250,i)));
end

%2) Look at the polar distribution of headings whithin those 10 sec
%in every case, and maybe discard those where the SD is too big?
figure('Position', [100 100 1600 900]),
for i = 1:length(j)
    distribution10sec(:,i) = circ_dist(deg2rad(shortData10sec.angPos(1:250,i)),goal10sec(i));
    [s(i) s0(i)] = circ_std(distribution10sec(:,i));
    subplot(4,ceil(length(jumps)/4),i)
    polarhistogram(distribution10sec(:,i),circedges,'Normalization','probability','FaceColor',[1,0,0],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
    title(strcat('SD = ',num2str(rad2deg(s0(i)))));
end
suptitle('Distribution of the distance to the goal 10 sec prior to jumps')

%For moving frames only
figure('Position', [100 100 1600 900]),
for i = 1:length(j)
    datamoving10sec{i} = shortData10sec.angPos(1:250,i);
    datamoving10sec{i}(shortData10sec.forwardVel(1:250,i)<=1) = NaN;
    datamoving10sec{i} = datamoving10sec{i}(~isnan(datamoving10sec{i}))';
    goal10secMoving{i} = circ_mean(deg2rad(datamoving10sec{i}),[],2);
    distribution10secMoving{i} = circ_dist(deg2rad(datamoving10sec{i}),goal10secMoving{i});
    [s(i) s010sec(i)] = circ_std(distribution10secMoving{i},[],[],2);
    subplot(4,length(jumps)/4,i)
    polarhistogram(distribution10secMoving{i},circedges,'Normalization','probability','FaceColor',[0,0,0],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
    title(strcat('SD = ',num2str(rad2deg(s010sec(i)))));
end
suptitle('Distribution of the distance to the goal 10 sec prior to jumps, moving frames')     
saveas(gcf,strcat(path,'plots\Distribution10secBJ',file(1:end-4),'.png'))
    
%3) Calculate the distance to the goal in the seconds after the
%jump
figure('Position', [100 100 1600 900]),
for i = 1:length(j)
    datamoving10secAJ{i} = shortData10sec.angPos(:,i);
    datamoving10secAJ{i}(shortData10sec.forwardVel(:,i)<=1) = NaN;
    datamoving10secAJ{i} = datamoving10secAJ{i}(~isnan(datamoving10secAJ{i}))';
    dist2goal10secMoving{i} = circ_dist(deg2rad(datamoving10secAJ{i}),goal10secMoving{i});
    dist2goal210secMoving{i} = wrapTo180(rad2deg(dist2goal10secMoving{i}));
    [countsDist10secMoving{i}] = histcounts(dist2goal210secMoving{i},edges);
    probabilitiesDist10secMoving{i} = countsDist10secMoving{i}./sum(countsDist10secMoving{i});
    degsFlyDist10secMoving{i} = linspace(-180,180,length(countsDist10secMoving{i}));
    
    %plot
    subplot(4,length(jumps)/4,i),
    plot(degsFlyDist10secMoving{i},probabilitiesDist10secMoving{i},'k')
    xlim([-180, 180]); xlabel('Distance to the goal (deg)');
    ylabel('Probability');
end
suptitle('Distance to the goal, 10 sec around jumps')     
saveas(gcf,strcat(path,'plots\Dist2goal10sec',file(1:end-4),'.png'))

%Plot them using colorscales
probaDist10secMoving = cell2mat(probabilitiesDist10secMoving);
probaDist10secMoving = reshape(probaDist10secMoving,length(probaDist10secMoving)/length(probabilitiesDist10secMoving),length(probabilitiesDist10secMoving));

%create new colormap
newMap = flipud(gray);
xaxis = [-180:360/17:180];
trials = [1:length(jumps)];
figure, imagesc(xaxis,trials,probaDist10secMoving')
colormap(newMap)
colorbar
saveas(gcf,strcat(path,'plots\HeatmapDist2goal10sec',file(1:end-4),'.png'))

%save data
save(strcat(path,'dataFromAnalysis\shortData10sec',file(1:end-4),'.mat'),'shortData10sec','probaDist10secMoving');
%% Per 'trial'

%I'm using a function to get and smooth data around the jumps
sec = 10; %how many sec before and after the jump I want to look at
perTrialData = getDataAroundJump(rawData,j,sec,sizeBall,sampleRate);

%apend jump data and save for using in the pooled flies analysis
perTrialData.jumpMag = jumps;
perTrialData.angPos = getPosAroundJump(data.ficTracAngularPosition,j,sec);

save(strcat(path,'dataFromAnalysis\perTrialData',file(1:end-4),'.mat'),'perTrialData');

%% Velocity around the jumps

time = linspace(-sec,sec,length(perTrialData.forwardVel));

%Individually
perTrialData.forwardVel(:,:) = -perTrialData.forwardVel(:,:); %I need to negate this for it to make sense given my fictrac axes at the moment this exp was run

for i = 1:length(j)
    figure,
    subplot(1,3,1)
    plot(time,perTrialData.forwardVel(:,i),'.')
    hold on
    plot(time,perTrialData.forwardVel(:,i))
    line([0,0],[min(perTrialData.forwardVel(:,i)), max(perTrialData.forwardVel(:,i))],'Color','black');
    ylim([min(perTrialData.forwardVel(:,i)), max(perTrialData.forwardVel(:,i))]);
    title('Forward velocity around the bar jumps');
    xlabel('Time(s)');
    ylabel('Velocity (mm/s)');
    
    subplot(1,3,2)
    plot(time,perTrialData.angVel(:,i),'.')
    hold on
    plot(time,perTrialData.angVel(:,i))
    line([0,0],[min(perTrialData.angVel(:,i)), max(perTrialData.angVel(:,i))],'Color','black');
    ylim([min(perTrialData.angVel(:,i)), max(perTrialData.angVel(:,i))]);
    title('Angular velocity around the bar jumps');
    xlabel('Time(s)');
    ylabel('Velocity (deg/s)');

    % Position around the jumps
%for ths, I need to use the corrected angular position of the fly in time

    time2 = linspace(-sec,sec,length(perTrialData.angPos));
    subplot(1,3,3)
    plot(time2,wrapTo360(perTrialData.angPos(:,i)),'.')
    hold on
    plot(time2,wrapTo360(perTrialData.angPos(:,i)))
    line([0,0],[0, 360],'Color','black');
    title('Angular Position around the bar jumps');
    xlabel('Time(s)');
    ylabel('Angular Position');
    ylim([0 360]);      

end


%Pooling results from similar magnitude jumps
%1) make a jump vector using the preloaded jump vector from the experiment
%and taking only as many elements as trials there were
trials = jumps(1:size(j,1));
%2) identify elements in that vector belonging to the 4 different groups,
%and put the perTrialData into those groups
Data90.forwardVel = perTrialData.forwardVel(:,trials == 90);
Data90.angVel = perTrialData.angVel(:,trials == 90);
DataNeg90.forwardVel = perTrialData.forwardVel(:,trials == -90);
DataNeg90.angVel = perTrialData.angVel(:,trials == -90);

%plot forward and angular velocity for every group

figure,
subplot(1,2,1)
plot(time,Data90.forwardVel,'.')
hold on
plot(time,Data90.forwardVel)
title('Forward velocity for 90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');    
subplot(1,2,2)
plot(time,Data90.angVel,'.')
hold on
plot(time,Data90.angVel)
title('Angular velocity for 90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (deg/s)');


figure,
subplot(1,2,1)
plot(time,DataNeg90.forwardVel,'.')
hold on
plot(time,DataNeg90.forwardVel)
title('Forward velocity for -90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (mm/s)');    
subplot(1,2,2)
plot(time,DataNeg90.angVel,'.')
hold on
plot(time,DataNeg90.angVel)
title('Angular velocity for -90 deg jumps');
xlabel('Time(s)');
ylabel('Velocity (deg/s)');


%Use heatmaps to see the evolution
figure,
trials = size(DataNeg90.angVel,2);
imagesc(trials,time,DataNeg90.angVel)
colorbar


%plot mean forward and angular velocity per group
meanForwardVel90 = mean(Data90.forwardVel,2);
meanForwardVelNeg90 = mean(DataNeg90.forwardVel,2);

meanAngVel90 = mean(Data90.angVel,2);
meanAngVelNeg90 = mean(DataNeg90.angVel,2);

figure('Position', [100 100 1600 900]),
subplot(1,2,1)
plot(time,meanForwardVel90,'r')
hold on
plot(time,meanForwardVelNeg90,'k')
title('Mean forward velocity');
legend({'90','-90'});
ylabel('Forward velocity (mm/s)'); xlabel('Time from bar jump (s)');
xlim([-10,10]);
plot([-10,10],[0,0],'-.k','HandleVisibility','off');
subplot(1,2,2)
plot(time,meanAngVel90,'r')
hold on
plot(time,meanAngVelNeg90,'k')
title('Mean angular velocity');
legend({'90','-90'});
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');
xlim([-5, 5]); ylim([min(meanAngVel90)-10, max(meanAngVelNeg90)+10]);
plot([0,0],[min(meanAngVel90)-10, max(meanAngVelNeg90)+10],'k','HandleVisibility','off');
plot([-5,5],[0,0],'-.k','HandleVisibility','off');

saveas(gcf,strcat(path,'plots\MeanAJvelocities',file(1:end-4),'.png'))