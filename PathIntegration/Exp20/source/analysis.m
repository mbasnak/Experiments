%Code to analyze exp 20: the virtual hallway that spins 180 when the fly reaches
%the end, gives the reward in the same location for both trial types and
%the gain is 1 (shorter hallway)
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp20\data\experimental flies'

[fileName,pathName] = uigetfile();   
rawData = load([pathName,fileName]);

sid = fileName(end-10);

%load the run object
run_obj_dir = [pathName,'runobj'];
cd (run_obj_dir)
run_obj_name = dir(strcat('*',sid,'_runobj.mat'));
run_obj = load(strcat(run_obj_dir,'\',run_obj_name.name));

%%  Define Ni-Daq channels ID
  
headingFly = 1;
%yFly = 2; %right now I'm not collecting the y signal
%through the NiDaq
yawFly =  2;
xFly = 3;
xFlyGain = 7;
xPanels = 4;
yPanels = 5;
PanelStatus = 6;
OptoTrigger = 8;

samplingRate = 4000;

%% Use the signal from the turns in the y panels voltage to determine the changes in trip

%convert from frames to seconds
time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(rawData.trial_bdata));

%smooth a little the x panel data to better detect the voltage thresholds
testData = movmedian(rawData.trial_bdata(:,xPanels),500);

%set the conditions and get the turns: turns are hapenning when the x panel voltage
%is between 9.87 and 5, or between 0 and 5
turns = (testData > 5 & testData < 9.84) | (testData > 0 & testData < 4.84);

%get the turn points taking the derivative and look for peaks bigger than
%0.5
test = abs(diff(turns));
test = [test;test(end)];
turnPoints = find(test > 0.5);

%Plot to see if they're accurate
figure,
plot(time,testData)
hold on
plot(time(turnPoints),testData(turnPoints),'r.')

%get the trips as trip 1 = turnpoint1 (start):turnpoint 2, trip 2 =
%turnpoint3:turnpoint4, trip 3 = turnpoint5:turnpoint6
tripData = {};
for i = 1:floor(length(turnPoints)/2)
   tripData{i} = rawData.trial_bdata(turnPoints(i+i-1):turnPoints(i+i),:);
end

%We are going to take a few steps to verify that the trips and turns have
%been properly identified
%(1) check the median variance in the voltage to see if it actually
%identified the trips, or it these were turns
for i = 1:length(tripData)
    Varian(i) = var(tripData{i}(:,4));
end
% For the trips, the median voltage variance should be around 0. For the
% turns, it should be around 2.
if median(Varian) > 1
    turnData = tripData; %if the median voltage is too high, these trips were actually turns
    clear tripData;
    for i = 1:floor(length(turnPoints)/2)-1
        tripData{i} = rawData.trial_bdata(turnPoints(i+i):turnPoints(i+i+1),:); %get the trips as you would usually get the turns
    end
else
    turnData = {};
    for i = 1:floor(length(turnPoints)/2)-1
       turnData{i} = rawData.trial_bdata(turnPoints(i+i):turnPoints(i+i+1),:); 
    end
end

%Plot the x panel bouts of the trip bouts: they should all look flat
figure('Position',[100 200 1400 800]),
for i = 1:length(tripData)
    subplot(5,ceil(length(tripData)/5),i)
    plot(tripData{i}(:,4))
    title(['Trip # ',num2str(i)]);
end
%suptitle('Changes in panel x dimension during the trips');

%Plot the x panel bouts of the turn bouts: they should all look oblique
figure%('Position',[100 200 1400 800]),
for i = 1:length(turnData)
    subplot(5,ceil(length(turnData)/5),i)
    plot(turnData{i}(:,4))
    title(['Turn # ',num2str(i)]);
end
%suptitle('Changes in panel x dimension during the turns');

%Remove single spurious turns or trials
for i = 1:length(tripData)
    Varian(i) = var(tripData{i}(:,4));
    if Varian(i) > 1
        tripData{i} = [];
    end
end
tripData = tripData(~cellfun('isempty',tripData));

for i = 1:length(turnData)
    Varian2(i) = var(turnData{i}(:,4));
    if Varian2(i) < 1
        turnData{i} = [];
    end
end
turnData = turnData(~cellfun('isempty',turnData));


%Plot the x panel bouts of the trip bouts: they should all look flat
figure%('Position',[100 200 1400 800]),
for i = 1:length(tripData)
    subplot(5,ceil(length(tripData)/5),i)
    plot(tripData{i}(:,4))
    title(['Trip # ',num2str(i)]);
end
%suptitle('Changes in panel x dimension during the trips');

%Plot the x panel bouts of the trip bouts: they should all look oblique
figure%('Position',[100 200 1400 800]),
for i = 1:length(turnData)
    subplot(5,ceil(length(turnData)/5),i)
    plot(turnData{i}(:,4))
    title(['Turn # ',num2str(i)]);
end
%suptitle('Changes in panel x dimension during the turns');
%% Separate the full data into trips and smooth it

sizeBall = 9;

for i = 1:length(tripData)
    %Panels data
    data{i}.xPanelVolts =  tripData{1,i}(:,xPanels); 
    VOLTAGE_RANGE_x = 10;
    maxValX =  96;

    data{i}.yPanelVolts =  tripData{1,i}(:, yPanels);
    VOLTAGE_RANGE_y = 10;
    maxValY = 92; 
           
    %FicTrac data
    data{i}.ficTracAngularPosition = tripData{1,i}( : , yawFly); 
    data{i}.ficTracIntx = tripData{1,i}( : , xFly); 
    %data{i}.ficTracInty = tripData{1,i}( : , yFly); 

    %smooth the data and get the velocity
    smoothed{i} = singleTrialVelocityAnalysis(data{i},samplingRate);
    %smoothed{i} = newVelocityAnalysis(data{i},samplingRate);

    ypanels{i} = resample(tripData{1,i}(:,yPanels),25,samplingRate);
    optoPulse{i} = resample(tripData{1,i}(:,OptoTrigger),25,samplingRate);
    trialDuration(i) = length(data{i}.xPanelVolts)/samplingRate;
end


%% Trial dur vs trial number

probeTrials = sort(str2num(rawData.probe_trials));
trial_num = length(tripData);
trialNum = 1:trial_num;
optoTrials = sort(setdiff(trialNum,probeTrials));
probeTrials = probeTrials(probeTrials<=length(trialDuration)); %remove trials that excede the total number of trials the fly got
probeTrials = probeTrials(probeTrials > 5); %remove the first 5 trials since she didn't get an opto experience up to that point
optoTrials = optoTrials(optoTrials<=length(trialDuration)); %remove trials that excede the total number of trials the fly got

PT = ones(1,length(probeTrials));
OT = ones(1,length(optoTrials));
group = [PT, 2*OT];
durPT = trialDuration(probeTrials);
durOT = trialDuration(optoTrials);

% Plot the trial duration in time
figure('Position', [400 200 1000 800])
subplot(1,2,1)
plot(trialNum(probeTrials),trialDuration(probeTrials),'ko','MarkerFaceColor', 'k','DisplayName','Probe trials')
hold on
plot(trialNum(optoTrials),trialDuration(optoTrials),'ro','MarkerFaceColor', 'r','DisplayName','Opto trials')
xlabel('Trial number'); ylabel('Trial duration (s)'); 
legend('Probe trials', 'Opto trials');

subplot(1,2,2)
boxplot([durPT';durOT'],group)
ax = gca;
ax.XTick = [1 2];
ax.XTickLabel={'Probe Trials','Opto Trials'};
title('Distribution of trial duration');ylabel('Trial duration (s)');

%saveas(gcf,[pathName(1:end-5),'plots\TrialDurVsTrialNum_sid_',sid,'.png'])


%% Remove outlier that have a really big variance in the xVel, because it probably indicates an artifact

for i = 1:size(smoothed,2)
    varian(i) = nanvar(smoothed{1,i}.xVel);
end

%exclude trials with var > 8 in the xVel.
excludetrials=[];
for i = 1:size(smoothed,2)
    if varian(i) > 8
        excludetrials=[excludetrials i];
    end
end

%remove outlier trials 
optoTrials=setdiff(optoTrials,excludetrials);
probeTrials=setdiff(probeTrials,excludetrials);

%% Look at something analogous to spike triggered average

%looking at the opto pulse, there are trials that should be opto
%trials that are not showing any pulse, in which the pulse seems to have failed.
%I'm getting rid of those for the following analyses.

for i = 1:length(optoTrials)    
    if (sum(diff(tripData{1,optoTrials(i)}(:,OptoTrigger))<-1)>0 && sum(diff(tripData{1,optoTrials(i)}(:,OptoTrigger))>1)>0)
        WorkingOptoPulse(i) = optoTrials(i); %keep the trial
        startPulse(i) = floor((find(diff(tripData{1,optoTrials(i)}(:,OptoTrigger))>1,1))*25/4000); %find the start point of the pulse
        EndPulse{i} = find(diff(tripData{1,optoTrials(i)}(:,OptoTrigger))<-1); %find the end point of the pulse
        endPulse(i) = floor(EndPulse{i}(end)*25/4000);
    else
        WorkingOptoPulse(i) = NaN;
        startPulse(i) = NaN;
        endPulse(i) = NaN;
    end
end
WorkingOptoPulse = WorkingOptoPulse(~isnan(WorkingOptoPulse));
startPulse = startPulse(~isnan(startPulse));
endPulse = endPulse(~isnan(endPulse));

%calculate the mean velocity before and after the pulse start and pulse end
velAroundPulse = {};
for i = 1:length(WorkingOptoPulse)  
    if startPulse(i)>14
        velAroundPulse{i} = smoothed{1,WorkingOptoPulse(i)}.xVel(startPulse(i)-13:startPulse(i)+25);
    else
        velAroundPulse{i} = smoothed{1,WorkingOptoPulse(i)}.xVel(1:startPulse(i)+25);
    end
end

trialL = cellfun(@length,velAroundPulse);
shortT = find(trialL<39);
for i = 1:length(shortT)
    velAroundPulse{shortT(i)} = [repelem(nan,1,14-startPulse(shortT(i))),velAroundPulse{shortT(i)}];
end    
velAroundPulse = cell2mat(velAroundPulse);
velAroundPulse = reshape(velAroundPulse,39,length(trialL));

meanVelBeforePulse = nanmean(velAroundPulse(1:13,:));
meanVelDuringPulse = nanmean(velAroundPulse(14:26,:));
meanVelAfterPulse = nanmean(velAroundPulse(27:39,:));
medianVelAroundPulse = [meanVelBeforePulse',meanVelDuringPulse',meanVelAfterPulse'];

figure('Position', [200 200 1200 800])
subplot(1,2,1)
time = linspace(-0.5,1,39);
boundedline(time,nanmean(velAroundPulse,2),nanstd(velAroundPulse,[],2)/sqrt(length(WorkingOptoPulse)))
hold on
line([0 0],[min(min(nanmean(velAroundPulse,2)))-0.3 max(max(nanmean(velAroundPulse,2)))+0.3],'Color','r','LineWidth',2)
line([0.5 0.5],[min(min(nanmean(velAroundPulse,2)))-0.3 max(max(nanmean(velAroundPulse,2)))+0.3],'Color','r','LineWidth',2)
ylim([min(min(nanmean(velAroundPulse,2)))-0.3 max(max(nanmean(velAroundPulse,2)))+0.3]);
ylabel('Forward velocity (mm/s)');
xlabel('Time around the pulse (s)');
xlim([-0.5 1]);
text(-0.25,max(max(nanmean(velAroundPulse,2)))+0.2,'1','FontSize',14)
text(0.25,max(max(nanmean(velAroundPulse,2)))+0.2,'2','FontSize',14)
text(0.75,max(max(nanmean(velAroundPulse,2)))+0.2,'3','FontSize',14)

%Profile plots of the velocity in different time points
subplot(1,2,2)
plot(medianVelAroundPulse','color',[0.5 0.5 0.5])
hold on
err = iqr(medianVelAroundPulse);
errorbar([1,2,3],mean(medianVelAroundPulse),err,'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
xticks([1 2 3])
xticklabels({'1','2','3'})
xlim([0 4]);

%suptitle('Forward velocity around the pulse');

%saveas(gcf,[pathName(1:end-5),'plots\VelAroundPulseOpto_sid_',sid,'.png'])


%% Plotting velocity as a function of distance travelled

WorkingProbeTrials = probeTrials; %in this configuration of the hallway, all the probetrials are useful trials
%load the pattern voltages to then define the voltage steps for each
%distance point
load('Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp20\source\utils\stimulusVoltages.mat');


%bin the forward velocity data into voltages using the ypanel data
for i = 1:length(smoothed)
    for j = 1:48
        step{i,j} = find(ypanels{1,i} > medianVoltages(j)-0.04 & ypanels{1,i} < medianVoltages(j)+0.04);
        stepMeans(i,j) = mean(smoothed{1,i}.xVel(step{i,j}));
        angularStepMeans(i,j) = mean(abs(smoothed{1,i}.angularVel(step{i,j})));
    end
end

%Linearly interpolate the data to fill in the gaps
for i = 1:size(stepMeans,1)
    stepMeans(i,:) = fillmissing(stepMeans(i,:),'linear');
    angularStepMeans(i,:) = fillmissing(angularStepMeans(i,:),'linear');
end

%Divide trips in left and right
%define the first one as left (=1)
%alternate left:right (1:0) for the following
tripType = zeros(1,length(tripData));
for i = 1:2:length(tripData)
    tripType(i) = 1;
end

WorkingOptoLeft = WorkingOptoPulse(tripType(WorkingOptoPulse) == 1);
WorkingProbeLeft = WorkingProbeTrials(tripType(WorkingProbeTrials) == 1);
WorkingOptoRight = WorkingOptoPulse(tripType(WorkingOptoPulse) == 0);
WorkingProbeRight = WorkingProbeTrials(tripType(WorkingProbeTrials) == 0);

%Get medians and error for each trip type
optoLeftMeans = nanmedian(stepMeans(WorkingOptoLeft,:));
probeLeftMeans = nanmedian(stepMeans(WorkingProbeLeft,:));
optoLeftError = (nanstd(stepMeans(WorkingOptoLeft,:)))./sqrt(length(WorkingOptoLeft));
probeLeftError = (nanstd(stepMeans(WorkingProbeLeft,:)))./sqrt(length(WorkingProbeLeft));

optoRightMeans = nanmedian(stepMeans(WorkingOptoRight,:));
probeRightMeans = nanmedian(stepMeans(WorkingProbeRight,:));
optoRightError = (nanstd(stepMeans(WorkingOptoRight,:)))./sqrt(length(WorkingOptoRight));
probeRightError = (nanstd(stepMeans(WorkingProbeRight,:)))./sqrt(length(WorkingProbeRight));

optoLeftAngMeans = nanmedian(angularStepMeans(WorkingOptoLeft,:));
probeLeftAngMeans = nanmedian(angularStepMeans(WorkingProbeLeft,:));
optoLeftAngError = (nanstd(angularStepMeans(WorkingOptoLeft,:)))./sqrt(length(WorkingOptoLeft));
probeLeftAngError = (nanstd(angularStepMeans(WorkingProbeLeft,:)))./sqrt(length(WorkingProbeLeft));

optoRightAngMeans = nanmedian(angularStepMeans(WorkingOptoRight,:));
probeRightAngMeans = nanmedian(angularStepMeans(WorkingProbeRight,:));
optoRightAngError = (nanstd(angularStepMeans(WorkingOptoRight,:)))./sqrt(length(WorkingOptoRight));
probeRightAngError = (nanstd(angularStepMeans(WorkingProbeRight,:)))./sqrt(length(WorkingProbeRight));

%Get the distance data
ydimension = [1:48]; %number of stimulus dimensions we're using
%convert the ydimensions to distance using the xgain
distance = 9*pi*(ydimension/92)/run_obj.run_obj.gain_x;
%get reward location
rewardDimension = round(run_obj.run_obj.reward_distance*48/100);

figure('Position', [100 100 1600 1000]),
ax1 = subplot(2,2,1);
[opT] = boundedline(distance,optoLeftMeans,optoLeftError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeLeftMeans,probeLeftError,'-ko','alpha');
%xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Forward velocity (mm/s)');
title('Left trips')

ax2 = subplot(2,2,2);
[opT] = boundedline(distance,flip(optoRightMeans),flip(optoRightError),'-ro','alpha');
hold on
[prT] = boundedline(distance,flip(probeRightMeans),flip(probeRightError),'-ko','alpha');
%xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
title('Right trips')

ax3 = subplot(2,2,3);
[opT] = boundedline(distance,optoLeftAngMeans,optoLeftAngError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeLeftAngMeans,probeLeftAngError,'-ko','alpha');
%xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Angular speed (deg/s)');
 xlabel('Distance travelled (mm)');

ax4 = subplot(2,2,4);
[opT] = boundedline(distance,flip(optoRightAngMeans),flip(optoRightAngError),'-ro','alpha');
hold on
[prT] = boundedline(distance,flip(probeRightAngMeans),flip(probeRightAngError),'-ko','alpha');
%xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance travelled (mm)');

%get the min and max velocities to set the axes limits
minOFwdVel = min(min([optoLeftMeans,optoRightMeans]));
minPFwdVel = min(min([probeLeftMeans,probeRightMeans]));
minFwdVels = [minOFwdVel;minPFwdVel];
minFwdVel = min(minFwdVels);
maxOFwdVel = max(max([optoLeftMeans,optoRightMeans]));
maxPFwdVel = max(max([probeLeftMeans,probeRightMeans]));
maxFwdVels = [maxOFwdVel;maxPFwdVel];
maxFwdVel = max(maxFwdVels);
ylimits = ([minFwdVel-2,maxFwdVel+2]);
ylim(ax1,ylimits); ylim(ax2,ylimits);

maxOAngVel = max(max([optoLeftAngMeans,optoRightAngMeans]));
maxPAngVel = max(max([probeLeftAngMeans,probeRightAngMeans]));
maxAngVels = [maxOAngVel;maxPAngVel];
maxAngVel = max(maxAngVels);
ylimits2 = ([0,maxAngVel+10]);
ylim(ax3,ylimits2); ylim(ax4,ylimits2);

%saveas(gcf,[pathName(1:end-5),'plots\MeanVelvsPanels_sid_',sid,'.png'])

%% Compare to shifted control

figure('Position',[100 100 1400 1000]),
subplot(2,2,1)
shiftedControl(distance,stepMeans,WorkingOptoLeft,rewardDimension,1,1)
subplot(2,2,2)
shiftedControl(distance,stepMeans,WorkingOptoRight,rewardDimension,1,0)
subplot(2,2,3)
shiftedControl(distance,stepMeans,WorkingProbeLeft,rewardDimension,0,1)
subplot(2,2,4)
shiftedControl(distance,stepMeans,WorkingProbeRight,rewardDimension,0,0)

%saveas(gcf,[pathName(1:end-5),'plots\ShiftedControl_sid_',sid,'.png'])

%% Forward velocity in individual trials

figure('Position', [100 100 1600 1000]),
subplot(2,2,1)
for i = 1:length(WorkingOptoLeft)
    plot(distance,stepMeans(WorkingOptoLeft,:))
    hold on
end
plot(distance,optoLeftMeans,'r','LineWidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Left trips');
ylabel('Forward velocity (mm/s)');
ylim([0,7]);

subplot(2,2,2)
for i = 1:length(WorkingOptoRight)
    plot(distance,flip(stepMeans(WorkingOptoRight,:),2))
    hold on
end
plot(distance,flip(optoRightMeans),'r','LineWidth',2)
xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
title('Right trips');
ylim([0,7]);

subplot(2,2,3)
for i = 1:length(WorkingProbeLeft)
    plot(distance,stepMeans(WorkingProbeLeft,:))
    hold on
end
plot(distance,probeLeftMeans,'k','LineWidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');
ylabel('Forward velocity (mm/s)');
ylim([0,7]);

subplot(2,2,4)
for i = 1:length(WorkingProbeRight)
    plot(distance,flip(stepMeans(WorkingProbeRight,:),2))
    hold on
end
plot(distance,flip(probeRightMeans),'k','LineWidth',2)
xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');
ylim([0,7]);

saveas(gcf,[pathName(1:end-5),'plots\VelvsPanelsSingleTrials_sid',sid,'.png'])


%% Look at individual opto and probe trials to detect anticipation


%For left opto trips only
figure('Position', [100 100 300 1000]),
for i = 1:length(WorkingOptoLeft)
    subplot(length(WorkingOptoLeft),1,i),
    plot(distance,stepMeans(WorkingOptoLeft(i),:),'-ro')
    hold on
    xline(distance(rewardDimension),'lineWidth',2,'color','b');
end
suptitle('Left opto trips');
saveas(gcf,[pathName(1:end-5),'plots\SingleTrialsOptoLeft_sid',sid,'.png'])


%For left probe trips only
figure('Position', [100 100 300 1000]),
for i = 1:length(WorkingProbeLeft)
    subplot(length(WorkingProbeLeft),1,i),
    plot(distance,stepMeans(WorkingProbeLeft(i),:),'-ko')
    hold on
    xline(distance(rewardDimension),'lineWidth',2,'color','b');
end
suptitle('Left probe trips');
saveas(gcf,[pathName(1:end-5),'plots\SingleTrialsProbeLeft_sid',sid,'.png'])


%For right opto trips only
figure('Position', [100 100 300 1000]),
for i = 1:length(WorkingOptoRight)
    subplot(length(WorkingOptoRight),1,i),
    plot(distance,flip(stepMeans(WorkingOptoRight(i),:)),'-ro')
    hold on
    xline(distance(rewardDimension),'lineWidth',2,'color','b');
end
suptitle('Right opto trips');
saveas(gcf,[pathName(1:end-5),'plots\SingleTrialsOptoRight_sid',sid,'.png'])


%For right probe trips only
figure('Position', [100 100 300 1000]),
for i = 1:length(WorkingProbeRight)
    subplot(length(WorkingProbeRight),1,i),
    plot(distance,flip(stepMeans(WorkingProbeRight(i),:)),'-ko')
    hold on
    xline(distance(rewardDimension),'lineWidth',2,'color','b');
end
suptitle('Right probe trips');
saveas(gcf,[pathName(1:end-5),'plots\SingleTrialsProbeRight_sid',sid,'.png'])


%% Angular speed in individual trials

figure('Position', [100 100 1600 1000]),
subplot(2,2,1)
for i = 1:length(WorkingOptoLeft)
    plot(distance,angularStepMeans(WorkingOptoLeft,:))
    hold on
end
plot(distance,optoLeftAngMeans,'r','LineWidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Left trips');
ylabel('Angular speed (deg/s)');
ylim([0,100]);

subplot(2,2,2)
for i = 1:length(WorkingOptoRight)
    plot(distance,flip(angularStepMeans(WorkingOptoRight,:),2))
    hold on
end
plot(distance,flip(optoRightAngMeans),'r','LineWidth',2)
xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
title('Right trips');
ylim([0,100]);

subplot(2,2,3)
for i = 1:length(WorkingProbeLeft)
    plot(distance,angularStepMeans(WorkingProbeLeft,:))
    hold on
end
plot(distance,probeLeftAngMeans,'k','LineWidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');
ylabel('Angular speed (deg/s)');
ylim([0,100]);

subplot(2,2,4)
for i = 1:length(WorkingProbeRight)
    plot(distance,flip(angularStepMeans(WorkingProbeRight,:),2))
    hold on
end
plot(distance,flip(probeRightAngMeans),'k','LineWidth',2)
xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');
ylim([0,100]);

saveas(gcf,[pathName(1:end-5),'plots\AngSpeedvsPanelsSingleTrials_sid',sid,'.png'])


%% Get the turns and check their duration

for i = 1:length(turnData)
    %Panels data
    turndata{i}.xPanelVolts = turnData{1,i}(:,xPanels); 
    VOLTAGE_RANGE_x = 10;
    maxValX =  96;

    turndata{i}.yPanelVolts = turnData{1,i}(:, yPanels);
    VOLTAGE_RANGE_y = 10;
    maxValY = 92; 
           
    %FicTrac data
    turndata{i}.ficTracAngularPosition = turnData{1,i}( : , yawFly); 
    turndata{i}.ficTracIntx = turnData{1,i}( : , xFly); 
    %turndata{i}.ficTracInty = turnData{1,i}( : , yFly); 

    %smooth the data and get the velocity
    turnsmoothed{i} = singleTrialVelocityAnalysis(turndata{i},samplingRate);
    turnDuration(i) = length(turndata{i}.xPanelVolts)/samplingRate;
end

%verify the turn duration
figure, plot(turnDuration,'o')
ylim([1 5]), xlim([1 length(turnData)]);
title('Turn duration');
ylabel('Duration (sec)'); xlabel('Turn number');

%% Get clockwise and counterclockwise turns and plot fly's behavior

clockwise = true(1,length(turnData));
counterclockwise = true(1,length(turnData));

%In clockwise turns, the voltage at the end of the turn is bigger than it
%was at the beginning, and the opposite is true for counterclockwise turns.
for i = 1:length(turnData)
    if turndata{1,i}.xPanelVolts(end) > turndata{1,i}.xPanelVolts(1)
        clockwise(i) = true;
        counterclockwise(i) = false;
    else
        clockwise(i) = false;
        counterclockwise(i) = true;
    end
end

turnNum = [1:1:length(turnData)];
clockwiseTurns = turnNum(clockwise);
counterclockwiseTurns = turnNum(counterclockwise);

%get minimum size of the turns, and have them all be that length
for i = 1:length(turnsmoothed)
    sizeTurns(i) = length(turnsmoothed{1,i}.angularVel);
end
[minSizeTurn,indMinSizeTurn] = min(sizeTurns);
TurnTime = linspace(0,3,length(turnsmoothed{1,indMinSizeTurn}.angularVel));

figure,
ax1 = subplot(1,2,1);
for i = 1:sum(clockwise)
    clockAngVel(:,i) = turnsmoothed{1,clockwiseTurns(i)}.angularVel(1:length(turnsmoothed{1,indMinSizeTurn}.angularVel));
    plot(TurnTime,clockAngVel(:,i))
    hold on
end
plot(TurnTime,mean(clockAngVel,2),'k','lineWidth',2)
ylabel('Angular velocity (deg/s)');
xlabel('Time from turn start (sec)');
title('Clockwise turns');

ax2 = subplot(1,2,2);
for i = 1:sum(counterclockwise)
    counterclockAngVel(:,i) = turnsmoothed{1,counterclockwiseTurns(i)}.angularVel(1:length(turnsmoothed{1,indMinSizeTurn}.angularVel));
    plot(TurnTime,counterclockAngVel(:,i))
    hold on
end
plot(TurnTime,mean(counterclockAngVel,2),'k','lineWidth',2)
xlabel('Time from turn start (sec)');
title('Counterclockwise turns');
%get min and max angular velocities to set plot limits
minCAngVel = min(min(clockAngVel));
minCCAngVel = min(min(counterclockAngVel));
minAngVels = [minCAngVel;minCCAngVel];
minAngVel = min(minAngVels);
maxCAngVel = max(max(clockAngVel));
maxCCAngVel = max(max(counterclockAngVel));
maxAngVels = [maxCAngVel;maxCCAngVel];
maxAngVel = max(maxAngVels);
limitsAngVel = [abs(minAngVel),maxAngVel];
limitAngVel = max(limitsAngVel);
ylimits = ([-limitAngVel-10,limitAngVel+10]);
ylim(ax1,ylimits); ylim(ax2,ylimits);

saveas(gcf,[pathName(1:end-5),'plots\TurnBehavior_sid',sid,'.png'])

%% Correlation between optomotor response in previous turn and anticipatory behavior

%get angular velocity magnitude for each turn
for i = 1:length(turnsmoothed)
    turnMag(i) = abs(max(turnsmoothed{1,i}.angularVel)-min(turnsmoothed{1,i}.angularVel));
end

%calculate error in estimation for each trial
%(I will do this by taking the local min and assigning that as the
%estimate, and then taking the difference between that estimation and the
%actual reward location)

for i = 1:length(WorkingProbeTrials)
    [localmin{i},Prom{i}] = islocalmin(stepMeans(WorkingProbeTrials(i),:)); %find local minima  
    if sum(localmin{i}) > 0
        [M{i},I{i}] = max(Prom{1,i}); %save the index of the biggest valley
    else
        I{i} = [];
    end
%     figure,
%     plot(ydimension,stepMeans(probeTrials(i),:))
%     hold on
%     plot(ydimension(I{i}),stepMeans(probeTrials(i),I{i}),'ro')
%     xlabel('y dimension'); ylabel('Forward vel (mm/s)');
    
    if sum(localmin{i}) > 0
        estimate(i) = distance(I{i});
        estimateError(i) = abs(estimate(i)-distance(rewardDimension));
%         if mod(i,2) == 1
%             estimateError(i) = estimate(i)-distance(rewardDimension);
%         else
%             estimateError(i) = estimate(i)-distance(end)-distance(rewardDimension);
%         end
        
    else
        estimate(i) = nan;
        estimateError(i) = nan;
    end
end

%Plot the correlation
figure,
plot(turnMag(WorkingProbeTrials-1),abs(estimateError),'o')
xlabel('Optomotor response magnitude in preceding open loop bout');
ylabel('Error in distance estimation');
hold on
lsline

saveas(gcf,[pathName(1:end-5),'plots\OptoResponseVsDistEstimation_sid',sid,'.png'])

%% Plotting the animal's trajectories

%load the HDF5 file data
files = dir(fullfile(pathName, '*.hdf5'));
pause(10)
fileName = files(1).name;

headingData = h5read([pathName,'\',fileName],'/heading');
timeData = h5read([pathName,'\',fileName],'/time');
posx = h5read([pathName,'\',fileName],'/posx');
posy = h5read([pathName,'\',fileName],'/posy');
open_loop_value = h5read([pathName,'\',fileName],'/open loop x');

%use the open loop value to determine the turning points
test = abs(diff(open_loop_value));
turningPoints = find(test > 0.4);
timing = timeData(2:end);

figure, plot(timing, abs(diff(open_loop_value)))
hold on
plot(timing(turningPoints),test(turningPoints),'ro')

%get the turns
for i = 1:floor(length(turningPoints)/2)
   TurnData(i).posx = posx(turningPoints(i+i-1):turningPoints(i+i),:);
   TurnData(i).posy = posy(turningPoints(i+i-1):turningPoints(i+i),:);
   TurnData(i).time = timeData(turningPoints(i+i-1):turningPoints(i+i),:);
end

%get the trips
TripData(1).posx =  posx(1:turningPoints(1),:);
TripData(1).posy =  posy(1:turningPoints(1),:);
TripData(1).time =  timeData(1:turningPoints(1),:);
for i = 1:floor(length(turningPoints)/2)-1
   MostTripData(i).posx = posx(turningPoints(i+i):turningPoints(i+i+1),:); 
   MostTripData(i).posy = posy(turningPoints(i+i):turningPoints(i+i+1),:); 
   MostTripData(i).time = timeData(turningPoints(i+i):turningPoints(i+i+1),:); 
end
for i = 2:length(MostTripData)+1
   TripData(i).posx = MostTripData(i-1).posx;
   TripData(i).posy = MostTripData(i-1).posy;
   TripData(i).time = MostTripData(i-1).time;
end

scaleFactor = 4.5; %ball's radius in mm

%Calculate the 'curviness' of the trajectory, as the sum of the total
%angular speed of a given trial
for i = 1:length(smoothed)
    curviness(i) = abs(nansum(smoothed{1,i}.angularVel));
end

%trips arranged by straightness
[sortedCurviness,I] = sort(curviness);

%2nd criterion for straightness: D/L with D being the distance from the
%start to the end of the trajectory, and L being the length of the
%trajectory.
D = distance(end);
%getting the length of the trajectories as the sum of the distances between
%consecutive points (using pythagoras).
for i = 1:length(smoothed)
    TripData(i).distx = diff(TripData(i).posx.*scaleFactor);
    TripData(i).disty = diff(TripData(i).posy.*scaleFactor);
    TripData(i).dist = sqrt((TripData(i).distx).^2 + (TripData(i).disty).^2);
    L(i) = sum(TripData(i).dist);
    tortuosity(i) = D/L(i);
end

%trips arranged by tortuosity
[sortedTortuosity,I2] = sort([0-tortuosity]);

%% Aligning the opto data and plotting the trajectories

angle = rawData.trial_bdata( : , yawFly);
downsampledAngle = downsample(angle,4000/25);
angle = downsampledAngle.* 2 .* pi ./ 10;
angle = wrapTo180(rad2deg(angle));
time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(angle));
interp_panel_heading = wrapTo180(rad2deg(interp1(timeData, unwrap(deg2rad(headingData)), time)));
interp_panel_heading(isnan(interp_panel_heading)) = 0;

[c, lags] = xcorr(angle, -interp_panel_heading, 100);
[value1, idx1] = max(c);
[value2, idx2] = min(c);
if idx1 < 100 && idx1 > 5
    idx = idx1;
else
    idx = idx2;
end
offset_time = -lags(idx); % offset_time is in the time spacing of processed bdata


%Plot them from 'straighter' to less straight
figure('Position',[0 0 1800 1100]),
for i = 1:length(smoothed)
    subplot(10,ceil(length(smoothed)/10),i),
    time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(TripData(I2(i)).posx));
    timePulse = linspace(1,length(rawData.trial_bdata)/samplingRate,length(optoPulse{1,I2(i)}));
    optoPulse{1,I2(i)}(optoPulse{1,I2(i)} < 3) = 0;
    optoPulse{1,I2(i)}(optoPulse{1,I2(i)} > 3) = 5;
    interp_opto_pulse = interp1(timePulse,optoPulse{1,I2(i)},time);
    interp_opto_pulse2 = [zeros(abs(offset_time),1)',interp_opto_pulse];
    opto_pulse = interp_opto_pulse2(1:length(interp_opto_pulse));
    scatter((TripData(I2(i)).posx-TripData(I2(i)).posx(1))*scaleFactor,(TripData(I2(i)).posy-TripData(I2(i)).posy(1))*scaleFactor,[],opto_pulse);
    title(num2str(tortuosity(I2(i))));
end
suptitle('Trajectories sorted by tortuosity');
saveas(gcf,[pathName(1:end-5),'plots\TrajectoriesTortuosity_sid',sid,'.png'])

%Plot them from 'straighter' to less straight
figure('Position',[0 0 1800 1100]),
for i = 1:length(smoothed)
    subplot(10,ceil(length(smoothed)/10),i),
    time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(TripData(I(i)).posx));
    timePulse = linspace(1,length(rawData.trial_bdata)/samplingRate,length(optoPulse{1,I(i)}));
    optoPulse{1,I(i)}(optoPulse{1,I(i)} < 3) = 0;
    optoPulse{1,I(i)}(optoPulse{1,I(i)} > 3) = 5;
    interp_opto_pulse = interp1(timePulse,optoPulse{1,I(i)},time);
    interp_opto_pulse2 = [zeros(abs(offset_time),1)',interp_opto_pulse];
    opto_pulse = interp_opto_pulse2(1:length(interp_opto_pulse));
    scatter((TripData(I(i)).posx-TripData(I(i)).posx(1))*scaleFactor,(TripData(I(i)).posy-TripData(I(i)).posy(1))*scaleFactor,[],opto_pulse);
    title(num2str(curviness(I(i))));
end
suptitle('Trajectories sorted by curviness');
saveas(gcf,[pathName(1:end-5),'plots\TrajectoriesCurviness_sid',sid,'.png'])


%% Correlate straightness of trips to slowing down behavior

% 1) Criterion 1 for straightness = summed angular speed for the trip. The
% bigger the value, the less straight the trip.

figure('Position',[300 500 1000 600]),
subplot(1,2,1)
plot(-curviness(WorkingProbeTrials),abs(estimateError),'o')
xlabel('Straightness in walking behavior');
ylabel('Error in distance estimation');
hold on
lsline

% 2) Criterion 2 for straightness: D/L

subplot(1,2,2)
plot(tortuosity(WorkingProbeTrials),abs(estimateError),'ro')
xlabel('Straightness in walking behavior (using trajectory length)');
ylabel('Error in distance estimation');
hold on
lsline

saveas(gcf,[pathName(1:end-5),'plots\StraightnessVsError_sid',sid,'.png'])

%% Animal's trajectories during the open loop turns

figure('Position',[0 0 1800 1100]),
for i = 1:length(TurnData)
    subplot(10,ceil(length(TripData)/10),i),
    scatter((TurnData(i).posx-TurnData(i).posx(1))*scaleFactor,(TurnData(i).posy-TurnData(i).posy(1))*scaleFactor,[],(TurnData(i).time-TurnData(i).time(1)));
    hold on
    plot((TurnData(i).posx-TurnData(i).posx(1))*scaleFactor,(TurnData(i).posy-TurnData(i).posy(1))*scaleFactor,'k');
end
%% Analysis of initial probe trials

%get initial probe trials
LeftInitialProbeTrials = [1,3,5];
RightInitialProbeTrials = [2,4,6];

figure('Position', [100 100 1600 1000]),
subplot(2,1,1)
for i = 1:length(LeftInitialProbeTrials)
    plot(distance,stepMeans(LeftInitialProbeTrials,:))
    hold on
end
plot(distance,median(stepMeans(LeftInitialProbeTrials,:)),'k','LineWidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Left trips');
ylabel('Forward velocity (mm/s)');
ylim([0,7]);

subplot(2,1,2)
for i = 1:length(RightInitialProbeTrials)
    plot(distance,flip(stepMeans(RightInitialProbeTrials,:),2))
    hold on
end
plot(distance,median(flip(stepMeans(RightInitialProbeTrials,:),2)),'k','LineWidth',2)
xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
title('Right trips');
ylabel('Forward velocity (mm/s)');
xlabel('Distance walked (mm)');
ylim([0,7]);

suptitle('Pre-training probe block');

saveas(gcf,[pathName(1:end-5),'plots\PreTrainingBlock_sid',sid,'.png'])
