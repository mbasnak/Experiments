function singleTrialAnalysisExp14(path,file)

load([path,file]);

rawData = daq_data'; %transpose matrix for future use

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;
PanelStatus = 6; %this signal tells whether the panels are on or off.

    
%% Subset acquisition of x and y pos, as well as FicTrac data

data.xpanelVolts =  rawData (:,xPanels); 
VOLTAGE_RANGE_x = 10;
maxValX =  96 ;

data.yPanelVolts =  rawData (:, yPanels);
VOLTAGE_RANGE_y = 10;
maxValY = 96;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2

%FicTrac data
data.fictracAngularPosition = rawData ( : , headingFly); 
data.ficTracIntx = rawData ( : , xFly); 
data.ficTracInty = rawData ( : , yFly); 


%% Determine bar jumps
%This is going to be done in two ways:
%1) using the signal from channel 6 to determine when the panels were on
%and off, and adding jumps every 200 sec in between
%2) taking the derivative of the voltage signal from the ypanels channel to identify jumps


% Define when the panels turn on and off
%take the derivative of the panels' signal
panelON = diff(rawData(:,6));
%Find the frame when the panels turn on
[M, I] = max(panelON);
%Find the frame when they turn off.
[M2, I2] = min(panelON);
startFrame = I(1);
endFrame = I2(1);

% get the jumps using the signal from the yPanels channel
Jumps = diff(data.yPanelVolts);
Jumps(abs(Jumps)>0.4 & abs(Jumps)<1)=1;
Jumps = round(Jumps);

%Identify the changes in block using the y panel voltage, since in block 2,
%it goes above 5
block2 = find(data.yPanelVolts>5);
block2Start = block2(1);
block2End = block2(end);

j = find(Jumps); %indices of the actual bar jumps, taken from the y signal
%remove spurious consecutive jumps
j(diff(j)<4000)=[];
if size(j,1)>48
    j = j(2:end);
end
jsec = j/1000;


%Find the jumps that correspond to the changes in block
[~,change1] = min(abs(j-block2Start));
changeBlock1 = j(change1);

[~,change2] = min(abs(j-block2End));
changeBlock2 = j(change2);
%% Correct voltage

voltageCorr = data.yPanelVolts(j(change1)+1) - data.yPanelVolts(j(change1)-1);

data.yPanelVolts(j(change1)+1:j(change2)) = data.yPanelVolts(j(change1)+1:j(change2)) - voltageCorr;
data.yPanelVolts = 2*data.yPanelVolts;

%remove values from jumps vector
indices = [change1,change2];
j(indices) = [];
jsec = j/1000;

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
data.xPanelPos = round ((data.xPanelVolts  * maxValX) /VOLTAGE_RANGE_x); % Convert from what it reads in volts from the Ni-Daq to an X position in pixels in the panels
data.yPanelPos = round ((data.yPanelVolts  * maxValY) /VOLTAGE_RANGE_y);

jumpPos = data.yPanelPos(j+1)-data.yPanelPos(j-1);
degJumps = wrapTo180(jumpPos*(360/96));
%Check if the jump magnitude appears ok in the x panel data
jumpMag = data.xPanelPos(j+1)-data.xPanelPos(j-1); 
degMag = wrapTo180(jumpMag*(360/96));
%check if the jump magnitude appears ok in the angular position data
jumpMag2 = data.ficTracAngularPosition(j+1)-data.ficTracAngularPosition(j-1); 
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180

% If the jumps look different because the function got cut off because of
% the frame number, we create a newJumps vector with the real jump values:

 newJumps = zeros(size(j));
 
 for i = 1:length(newJumps)
    if degJumps(i)>-100 & degJumps(i)<-80
        newJumps(i) = -90;
    else
        newJumps(i) = 90;
    end
 end
 
 
%Check if the jump magnitude appears ok in the x panel data
jumpMag = data.xPanelPos(j+1)-data.xPanelPos(j-1); 
degMag = wrapTo180(jumpMag*(360/96));

%check if the jump magnitude appears ok in the angular position data
jumpMag2 = data.ficTracAngularPosition(j+1)-data.ficTracAngularPosition(j-1); 
radMag2 = jumpMag2.* 2 .* pi ./ 10; %convert from voltage to radians
degMag2 = wrapTo180(rad2deg(radMag2)); %convert from radians to degrees and wrap 180

%% Downsample, unwrap and smooth position data, then get velocity and smooth

sizeBall = 9;
sampleRate = 1000;
[smoothed] = singleTrialVelocityAnalysis9mm(data,sampleRate);
%[smoothed] = posDataDecoding(data,sampleRate);

%% 2D trajectories

%Uncorrected trajectory (without correcting the heading). This is to
%visualize the fly's response to the jumps, when they do very well.
[posx,posy]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.AngularPosition);

time = linspace(1,1000,length(posx));

figure, 
subplot(1,2,1)
scatter(posx,posy,0.5,time)
colorbar
axis equal
axis tight
title('Uncorrected 2D trajectory');


%Corrected trajectory: this is the actual trajectory the fly took during
%the block
[posx2,posy2]=FlyTrajectory(smoothed.Intx,smoothed.Inty,smoothed.angularPosition);

subplot(1,2,2), scatter(posx2,posy2,0.5,time);
colorbar
axis equal
axis tight
title('Corrected 2D trajectory');

saveas(gcf,strcat(path,'plots\Trajectories',file(1:end-4),'.png'));
close;

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
close;

%% Angular velocity

%To look at the velocity, I should do so without adding the offsets of the
%jumps, because the offsets will give me weird jumps in velocity that the
%fly might not actually have made. For the angular velocity, I'm using a new 
%function that I made to downsample, unwrap, smooth and get the velocity
%with the angular position uncorrected for the offset

AngularPosition = rawData ( : , headingFly);
%angularVelocity = getAngVel(AngularPosition);
angularVelocity = smoothed.AngularVel;

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
close;

%% Plot forward and angular velocities throughout the experiment using a raster plot

newMap = flipud(gray);
figure
set(gcf, 'Position', [300, 500, 1600, 500]),
subplot(2,1,1)
imagesc(time,[],forwardVelocity)
colormap(hot)
colorbar
xlabel('Time (s)');
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Forward Velocity (mm/s)')

subplot(2,1,2)
imagesc(time,[],angularVelocity')
colormap(hot)
colorbar
xlabel('Time (s)');
set(gca,'ytick',[])
set(gca,'yticklabel',[])
title('Angular Velocity (deg/s)')

saveas(gcf,strcat(path,'plots\VelocityRP',file(1:end-4),'.png'))
close;

%%  Activity levels

% We are going to decide whether a fly is moving or not based on the
% forward velocity. If it's above 1 mm/s we will consider it is moving
% We will work with the downsampled data

downsampled.xPanelPos = resample(data.xPanelPos,25,sampleRate); %downsample the panels position
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
close;

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
subplot(3,4,[1,5,9])
scatter(remapPosToDeg, time, [], forwardVelocity); %add the forward velocity as a color
hold on
plot(remapPosToDeg, time,'k','HandleVisibility','off')
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [changeBlock1/1000 changeBlock1/1000],'r','LineWidth',2);
     plot([-180 180], [changeBlock2/1000 changeBlock2/1000],'r','LineWidth',2);
end
ax = gca;
ax.YDir = 'reverse';

 % Heading in time, with only moving frames.
subplot(3,4,[4,8,12])
scatter(remapPosToDeg(forwardVelocity>1), time(forwardVelocity>1));
xlabel('Heading angle (deg)'); ylabel('Time (s)');
xlim([-180 180]); ylim([0 max(time)]);
hold on
for i = 1:length(j)
     plot([-180 180], [changeBlock1/1000 changeBlock1/1000],'r','LineWidth',2);
     plot([-180 180], [changeBlock2/1000 changeBlock2/1000],'r','LineWidth',2);
end
ax = gca;
ax.YDir = 'reverse'; 


% Plot the polar histograms

% Polar coordinates analysis of the stimulus position
posToRad = deg2rad(posToDeg);

%Plot the histogram in polar coordinates
circedges = [0:20:360];
circedges = deg2rad(circedges);


FB = posToRad(1:(changeBlock1/40));
SB = posToRad((changeBlock1/40)+1:(changeBlock2/40));
TB = posToRad((changeBlock2/40)+1:end);

forwardVelFB = forwardVelocity(1:(changeBlock1/40));
forwardVelSB = forwardVelocity((changeBlock1/40)+1:(changeBlock2/40));
forwardVelTB = forwardVelocity((changeBlock2/40)+1:end);


%Block1
subplot(3,4,2)
polarhistogram(FB,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block2
subplot(3,4,6)
polarhistogram(SB,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block3
subplot(3,4,10)
polarhistogram(TB,circedges,'Normalization','probability','FaceColor',[0.2,0.5,1],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top


%with moving frames
%Block1
subplot(3,4,3)
polarhistogram(FB(forwardVelFB>1),circedges,'Normalization','probability','FaceColor',[0,0.5,0.3],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block2
subplot(3,4,7)
polarhistogram(SB(forwardVelSB>1),circedges,'Normalization','probability','FaceColor',[0,0.5,0.3],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

%Block3
subplot(3,4,11)
polarhistogram(TB(forwardVelTB>1),circedges,'Normalization','probability','FaceColor',[0,0.5,0.3],'HandleVisibility','off');
ax = gca;
ax.ThetaDir = 'clockwise';
ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top

saveas(gcf,strcat(path,'plots\BarPosition',file(1:end-4),'.png'))
close;

%% Look at the goal and calculate the distance to it...

movingPos = posToRad(forwardVelocity>1);

edges = [-180:20:180];
%Taking all the experiment
figure,
%with every frame
goal = circ_mean(posToRad,[],2);
dist2goal2 = circ_dist(posToRad,goal);
dist2goal = wrapTo180(rad2deg(dist2goal2));
[countsDist] = histcounts(dist2goal,edges);
probabilitiesDist = countsDist./sum(countsDist);
degsFlyDist = linspace(-180,180,length(countsDist));
subplot(1,2,1), plot(degsFlyDist,probabilitiesDist,'r')
title('Distance to the goal with every frame');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

%taking only 'moving frames'
goalMoving = circ_mean(movingPos,[],2);
devGoalMoving = circ_std(movingPos,[],[],2);
dist2goalMoving2 = circ_dist(movingPos,goalMoving);
dist2goalMoving = wrapTo180(rad2deg(dist2goalMoving2));
[countsDistMoving] = histcounts(dist2goalMoving,edges);
probabilitiesDistMoving = countsDistMoving./sum(countsDistMoving);
degsFlyDistMoving = linspace(-180,180,length(countsDistMoving));
subplot(1,2,2), plot(degsFlyDistMoving,probabilitiesDistMoving,'r')
title('Distance to the goal with moving frames');
xlim([-180, 180]); xlabel('Distance to the goal (deg)');
ylabel('Probability');

saveas(gcf,strcat(path,'plots\Dist2Goal',file(1:end-4),'.png'))
close;
save(strcat(path,'dataFromAnalysis\goals',file(1:end-4),'.mat'),'goal','goalMoving','dist2goal','dist2goalMoving','degsFlyDistMoving','probabilitiesDistMoving','devGoalMoving');

%% Using 3 sec around jump

%1) Use the 3 sec preceding each jump to compute the goal by taking
%the mean heading
sec = 3;
shortData3sec = getDataAroundJump(rawData,j,sec,sizeBall,sampleRate);
shortData3sec.jumpMag = newJumps;
shortData3sec.angPos = getPosAroundJump(data.ficTracAngularPosition,j,sec);

for i = 1:length(j)
    goal3sec(i) = circ_mean(deg2rad(shortData3sec.angPos(1:75,i)));
end

%2) Look at the polar distribution of headings whithin those 10 sec
%in every case, and maybe discard those where the SD is too big?
figure('Position', [100 100 1600 900]),
for i = 1:length(j)
    distribution3sec(:,i) = circ_dist(deg2rad(shortData3sec.angPos(1:75,i)),goal3sec(i));
    [s(i) s0(i)] = circ_std(distribution3sec(:,i));
    subplot(4,ceil(length(newJumps)/4),i)
    polarhistogram(distribution3sec(:,i),circedges,'Normalization','probability','FaceColor',[1,0,0],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
    title(strcat('SD = ',num2str(rad2deg(s0(i)))));
end
suptitle('Distribution of the distance to the goal 3 sec prior to jumps')


%For moving frames only
figure('Position', [100 100 1600 900]),
for i = 1:length(j)
    datamoving3sec{i} = shortData3sec.angPos(1:75,i);
    datamoving3sec{i}(shortData3sec.forwardVel(1:75,i)<=1) = NaN;
    datamoving3sec{i} = datamoving3sec{i}(~isnan(datamoving3sec{i}))';
    goal3secMoving{i} = circ_mean(deg2rad(datamoving3sec{i}),[],2);
    distribution3secMoving{i} = circ_dist(deg2rad(datamoving3sec{i}),goal3secMoving{i});
    [s(i) s03sec(i)] = circ_std(distribution3secMoving{i},[],[],2);
    subplot(4,ceil(length(newJumps)/4),i)
    polarhistogram(distribution3secMoving{i},circedges,'Normalization','probability','FaceColor',[0,0,0],'HandleVisibility','off');
    ax = gca;
    ax.ThetaDir='clockwise';
    ax.ThetaZeroLocation = 'top'; %rotate the plot so that 0 is on the top
    title(strcat('SD = ',num2str(rad2deg(s03sec(i)))));
end
suptitle('Distribution of the distance to the goal 3 sec prior to jumps, moving frames')     
saveas(gcf,strcat(path,'plots\Distribution3secBJ',file(1:end-4),'.png'))
close;

%3) Calculate the distance to the goal in the seconds after the
%jump
figure('Position', [100 100 1600 900]),
for i = 1:length(j)
    datamoving3secAJ{i} = shortData3sec.angPos(:,i);
    datamoving3secAJ{i}(shortData3sec.forwardVel(:,i)<=1) = NaN;
    datamoving3secAJ{i} = datamoving3secAJ{i}(~isnan(datamoving3secAJ{i}))';
    dist2goal3secMoving{i} = circ_dist(deg2rad(datamoving3secAJ{i}),goal3secMoving{i});
    dist2goal23secMoving{i} = wrapTo180(rad2deg(dist2goal3secMoving{i}));
    [countsDist3secMoving{i}] = histcounts(dist2goal23secMoving{i},edges);
    probabilitiesDist3secMoving{i} = countsDist3secMoving{i}./sum(countsDist3secMoving{i});
    degsFlyDist3secMoving{i} = linspace(-180,180,length(countsDist3secMoving{i}));
    
    %plot
    subplot(4,ceil(length(newJumps)/4),i),
    plot(degsFlyDist3secMoving{i},probabilitiesDist3secMoving{i},'k')
    xlim([-180, 180]); xlabel('Distance to the goal (deg)');
    ylabel('Probability');
end
suptitle('Distance to the goal, 3 sec around jumps')     
saveas(gcf,strcat(path,'plots\Dist2goal3sec',file(1:end-4),'.png'))
close;

%Plot them using colorscales
probaDist3secMoving = cell2mat(probabilitiesDist3secMoving);
probaDist3secMoving = reshape(probaDist3secMoving,length(probaDist3secMoving)/length(probabilitiesDist3secMoving),length(probabilitiesDist3secMoving));
%create new colormap
newMap = flipud(gray);
xaxis = [-180:360/17:180];
trials = [1:length(newJumps)];
figure, imagesc(xaxis,trials,probaDist3secMoving')
colormap(newMap)
colorbar
saveas(gcf,strcat(path,'plots\HeatmapDist2goal3sec',file(1:end-4),'.png'))
close;

%save data
save(strcat(path,'dataFromAnalysis\shortData3sec',file(1:end-4),'.mat'),'shortData3sec','probaDist3secMoving');


%% Per 'trial'

%I'm using a function to get and smooth data around the jumps
sec = 3; %how many sec before and after the jump I want to look at
%perTrialData = getDataAroundJump(rawData,j,sec,sizeBall);
perTrialData = getDataAroundJumpFiltered(rawData,j,sec,sizeBall);

%apend jump data and save for using in the pooled flies analysis
perTrialData.jumpMag = newJumps;
%perTrialData.angPos = getPosAroundJump(data.ficTracAngularPosition,j,sec);
perTrialData.angPos = getPosAroundJumpFiltered(data.ficTracAngularPosition,j,sec);
save(strcat(path,'dataFromAnalysis\perTrialData',file(1:end-4),'.mat'),'perTrialData');


%% Velocity and around the jumps

time = linspace(-sec,sec,length(perTrialData.forwardVel));

%Pooling results from similar magnitude jumps
%1) make a jump vector using the preloaded jump vector from the experiment
%and taking only as many elements as trials there were
trials = newJumps;
%2) identify elements in that vector belonging to the 4 different groups,
%and put the perTrialData into those groups
Data90.forwardVel = perTrialData.forwardVel(:,trials == 90);
Data90.angVel = perTrialData.angVel(:,trials == 90);
DataNeg90.forwardVel = perTrialData.forwardVel(:,trials == -90);
DataNeg90.angVel = perTrialData.angVel(:,trials == -90);


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
xlim([-3,3]);
plot([-3,3],[0,0],'-.k','HandleVisibility','off');
subplot(1,2,2)
plot(time,meanAngVel90,'r')
hold on
plot(time,meanAngVelNeg90,'k')
title('Mean angular velocity');
legend({'90','-90'});
ylabel('Angular velocity (deg/s)'); xlabel('Time from bar jump (s)');
xlim([-3, 3]); %ylim([min(meanAngVel90)-10, max(meanAngVelNeg90)+10]);
plot([0,0],[min(meanAngVel90)-10, max(meanAngVelNeg90)+10],'k','HandleVisibility','off');
plot([-3,3],[0,0],'-.k','HandleVisibility','off');

saveas(gcf,strcat(path,'plots\MeanAJvelocities',file(1:end-4),'.png'))
close;

end