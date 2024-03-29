
function singleTrialAnalysisExp9(path,file)

load([path,file]);

rawData = daq_data'; %transpose matrix for future use

% Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xPanels = 4;
yPanels = 5;
PanelStatus = 6; %this signal tells whether the panels are on or off.


%% Determine trial change

% Define when the panels turn on and off
%take the derivative of the panels' signal
panelON = diff(rawData(:,6));

% get the trial changes using the signal from the 6th channel
panelsON = find(panelON>2);
panelsOFF = find(panelON<-2);

%% Separating the data per trial type

%All of the odd trials are stripe fixation trials (i.e., between the first
%panelON and the first panelOFF, between the 3rd panelON and the third
%panelOFF, etc,...)
oddJumps = [1:2:length(panelsOFF)];

for i = 1:length(oddJumps)
    stripeFixationTrials{i} = daq_data(:,panelsON(oddJumps(i)):panelsOFF(oddJumps(i)));
end

%First take all of the optomotor response trials
evenJumps = [2:2:length(panelsOFF)];

for i = 1:length(evenJumps)
    optomotorResponseTrials{i} = daq_data(:,panelsON(evenJumps(i)):panelsOFF(evenJumps(i)));
end

%Use the jumpFunctions stored to determine the direction of rotation of the
%stimulus
clockwiseTrials = optomotorResponseTrials(jumpFunction(1:length(optomotorResponseTrials)) == 50);
counterclockwiseTrials = optomotorResponseTrials(jumpFunction(1:length(optomotorResponseTrials)) == 51);

    
%% Subset acquisition of x and y pos, as well as FicTrac data

    VOLTAGE_RANGE_x = 9.77; % This should be 10 V, but empirically I measure 0.1 V for pos x=1 and 9.87 V for pos x=96
    maxValX =  96 ;% pattern.x_num (I am using 96 for every pattern now, but if it wasn't the case I would need to adjust it)
    VOLTAGE_RANGE_y = 9.86; %likewise, empirically this should be 10V, but I am getting 9.86
    maxValY = 96;% I think I am using 1 for my Y dimension for every pattern except the 4px grating, which uses 2

    
for i = 1:length(stripeFixationTrials)
    dataSF{i}.xpanelVolts =  stripeFixationTrials{i}(xPanels,:); 
    dataSF{i}.yPanelVolts =  stripeFixationTrials{i}(yPanels,:);
    dataSF{i}.ficTracAngularPosition = stripeFixationTrials{i}(headingFly,:); 
    dataSF{i}.ficTracIntx = -stripeFixationTrials{i}(xFly,:); 
    dataSF{i}.ficTracInty = stripeFixationTrials{i}(yFly,:); 
end

for i = 1:length(clockwiseTrials)
    dataClockwise{i}.xpanelVolts =  clockwiseTrials{i}(xPanels,:); 
    dataClockwise{i}.yPanelVolts =  clockwiseTrials{i}(yPanels,:);
    dataClockwise{i}.ficTracAngularPosition = clockwiseTrials{i}(headingFly,:); 
    dataClockwise{i}.ficTracIntx = -clockwiseTrials{i}(xFly,:); 
    dataClockwise{i}.ficTracInty = clockwiseTrials{i}(yFly,:); 
end
    
for i = 1:length(counterclockwiseTrials)
    dataCounterclockwise{i}.xpanelVolts =  counterclockwiseTrials{i}(xPanels,:); 
    dataCounterclockwise{i}.yPanelVolts =  counterclockwiseTrials{i}(yPanels,:);
    dataCounterclockwise{i}.ficTracAngularPosition = counterclockwiseTrials{i}(headingFly,:); 
    dataCounterclockwise{i}.ficTracIntx = -counterclockwiseTrials{i}(xFly,:); 
    dataCounterclockwise{i}.ficTracInty = counterclockwiseTrials{i}(yFly,:); 
end
    
    
%% Downsample, unwrap and smooth position data, then get velocity and smooth

sizeBall = 9;
sampleRate = 1000;

for i = 1:length(stripeFixationTrials)
    [smoothedSF{i}] = posDataDecoding(dataSF{1,i},sampleRate);
end

for i = 1:length(clockwiseTrials)
    [smoothedClockwise{i}] = posDataDecoding(dataClockwise{1,i},sampleRate);
end

for i = 1:length(counterclockwiseTrials)
    [smoothedCounterclockwise{i}] = posDataDecoding(dataCounterclockwise{1,i},sampleRate);
end


%% Forward velocity analysis

%(1) Stripe fixation trials
for i = 1:length(smoothedSF)
    forwardVelocity{i} = smoothedSF{1,i}.xVel(1:201);
end

forwardVelocity = cell2mat(forwardVelocity);
forwardVelocity = reshape(forwardVelocity,[length(smoothedSF{1,i}.xVel(1:201)),length(forwardVelocity)/length(smoothedSF{1,i}.xVel(1:201))]);

meanVelocity = mean(forwardVelocity,2);
time = linspace(0,8,length(forwardVelocity));

figure('Position', [100 100 1400 900]),
subplot(1,3,1)
plot(time,forwardVelocity)
hold on
plot(time,meanVelocity,'k','lineWidth',2)
title('Forward velocity for stripe fixation trials')
xlabel('Time (sec)'); ylabel('Velocity (mm/s)');
xlim([0 8]); ylim([min(min(forwardVelocity))-2, max(max(forwardVelocity))+2]);


%(2) Clockwise optomotor trials
for i = 1:length(smoothedClockwise)
    forwardVelocityClockwise{i} = smoothedClockwise{i}.xVel(1:24);
end

forwardVelocityClockwise = cell2mat(forwardVelocityClockwise);
forwardVelocityClockwise = reshape(forwardVelocityClockwise,[length(smoothedClockwise{1}.xVel),length(forwardVelocityClockwise)/length(smoothedClockwise{1}.xVel)]);

meanVelocityClockwise = mean(forwardVelocityClockwise,2);
time = linspace(0,3,24);

subplot(1,3,2)
plot(time,forwardVelocityClockwise)
hold on
plot(time,meanVelocityClockwise,'k','lineWidth',2)
title('Forward velocity for clockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(forwardVelocity))-2, max(max(forwardVelocity))+2]);


%(3) Counterclockwise optomotor trials
for i = 1:length(smoothedCounterclockwise)
    forwardVelocityCounterclockwise{i} = smoothedCounterclockwise{i}.xVel(1:24);
end

forwardVelocityCounterclockwise = cell2mat(forwardVelocityCounterclockwise);
forwardVelocityCounterclockwise = reshape(forwardVelocityCounterclockwise,[length(smoothedCounterclockwise{1}.xVel),length(forwardVelocityCounterclockwise)/length(smoothedCounterclockwise{1}.xVel)]);

meanVelocityCounterclockwise = mean(forwardVelocityCounterclockwise,2);
time = linspace(0,3,24);

subplot(1,3,3)
plot(time,forwardVelocityCounterclockwise)
hold on
plot(time,meanVelocityCounterclockwise,'k','lineWidth',2)
title('Forward velocity for counterclockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(forwardVelocity))-2, max(max(forwardVelocity))+2]);

saveas(gcf,strcat(path,'plots\ForwardVelocity',file(1:end-4),'.png'));
close;

%% Angular velocity

%(1) Stripe fixation trials
for i = 1:length(smoothedSF)
    angularVelocity{i} = smoothedSF{i}.angularVel(1:201);
end

angularVelocity = cell2mat(angularVelocity);
angularVelocity = reshape(angularVelocity,[length(smoothedSF{1}.angularVel(1:201)),length(angularVelocity)/length(smoothedSF{1}.angularVel(1:201))]);

meanVelocity = mean(angularVelocity,2);
time = linspace(0,8,length(smoothedSF{1}.angularVel(1:201)));

figure('Position', [100 100 1400 900]),
subplot(1,3,1)
plot(time,angularVelocity)
hold on
plot(time,meanVelocity,'k','lineWidth',2)
title('Angular velocity for stripe fixation trials')
xlabel('Time (sec)'); ylabel('Velocity (mm/s)');
xlim([0 8]); ylim([min(min(angularVelocity))-2, max(max(angularVelocity))+2]);


%(2) Clockwise optomotor trials
for i = 1:length(smoothedClockwise)
    angularVelocityClockwise{i} = smoothedClockwise{i}.angularVel(1:24);
end

angularVelocityClockwise = cell2mat(angularVelocityClockwise);
angularVelocityClockwise = reshape(angularVelocityClockwise,[length(smoothedClockwise{1}.angularVel),length(angularVelocityClockwise)/length(smoothedClockwise{1}.angularVel)]);

meanVelocityClockwise = mean(angularVelocityClockwise,2);

%(3) Counterclockwise optomotor trials
for i = 1:length(smoothedCounterclockwise)
    angularVelocityCounterclockwise{i} = smoothedCounterclockwise{i}.angularVel(1:24);
end

angularVelocityCounterclockwise = cell2mat(angularVelocityCounterclockwise);
angularVelocityCounterclockwise = reshape(angularVelocityCounterclockwise,[length(smoothedCounterclockwise{1}.angularVel),length(angularVelocityCounterclockwise)/length(smoothedCounterclockwise{1}.angularVel)]);

meanVelocityCounterclockwise = mean(angularVelocityCounterclockwise,2);
time = linspace(0,3,24);

subplot(1,3,2)
plot(time,angularVelocityClockwise)
hold on
plot(time,meanVelocityClockwise,'k','lineWidth',2)
title('Angular velocity for clockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(angularVelocityClockwise))-2, max(max(angularVelocityCounterclockwise))+2]);

subplot(1,3,3)
plot(time,angularVelocityCounterclockwise)
hold on
plot(time,meanVelocityCounterclockwise,'k','lineWidth',2)
title('Angular velocity for counterclockwise optomotor trials')
xlabel('Time (sec)');
xlim([0 3]); ylim([min(min(angularVelocityClockwise))-2, max(max(angularVelocityCounterclockwise))+2]);

saveas(gcf,strcat(path,'plots\AngularVelocity',file(1:end-4),'.png'));
close;

%Look at their evolution with a heatmap.
figure('Position', [100 100 1000 900]),
colormap(hot)
subplot(1,2,1),
imagesc(time,[1:length(clockwiseTrials)],angularVelocityClockwise')
colorbar
title('Angular Velocity in clockwise trials');

subplot(1,2,2),
imagesc(time,[1:length(counterclockwiseTrials)],angularVelocityCounterclockwise')
colorbar
title('Angular Velocity in counterclockwise trials');

saveas(gcf,strcat(path,'plots\AngularVelocityEvolution',file(1:end-4),'.png'));
close;

end