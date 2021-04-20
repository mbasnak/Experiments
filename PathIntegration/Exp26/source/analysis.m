%Code to analyze exp 26: the inbound/outbound hallway with 360 loops

close all; clear all;

%% Prompt the user to select the file to open and load it.

cd 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp26\data'

%Select the file name that starts with 'bdata'
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
yawFly =  2;
xFly = 3;
xFlyGain = 7;
xPanels = 4;
yPanels = 5;
PanelStatus = 6;
OptoTrigger = 8;
cameraTrigger = 9;

samplingRate = 4000;

%% Check if the change in trip type and trip counter make sense

figure,
%plot the matlab data
subplot(2,1,1)
%Plot the panel's y dimension
plot(rawData.trial_time,rawData.trial_bdata(:,yPanels))
hold on
%Add the change in trip
for changeInTrip = 1:length(rawData.trialData)
    [M I] = min(abs(rawData.trial_time-rawData.trialData(1,changeInTrip)));
    if rawData.trialData(2,changeInTrip) == 2
        plot(rawData.trialData(1,changeInTrip),rawData.trial_bdata(I,yPanels),'ro')
    else
        plot(rawData.trialData(1,changeInTrip),rawData.trial_bdata(I,yPanels),'bo')
    end    
    
end
%Add text for the trip counter
for changeInTrip = 1:length(rawData.trialData)
    text(rawData.trialData(1,changeInTrip),9,num2str(rawData.trialData(3,changeInTrip)))
end
xlabel('Time (sec)');
ylabel('Y panels voltage (V)');
ylim([0 10]);
title('NiDaq data')

% Compare to the hdf5 file to see if python and matlab are matching

%Get file name
fileList = dir(fullfile(pathName, '*.hdf5'));
file_name = [pathName,fileList.name];

%Import relevant variables
timestamp = h5read(file_name,'/time');
voltage_out = h5read(file_name,'/output_voltage_x_gain');
trip_type = h5read(file_name,'/trip_type');
trial_counter = h5read(file_name,'/trial_counter');

%Get the times of changes in trip_type
change_trip = diff(trip_type);
change_trip_frames = find(abs(change_trip)>0.5);
change_to_right = [];
change_to_left = [];
for change = 1:length(change_trip_frames)
    if trip_type(change_trip_frames(change)+1) == 2
        change_to_right = [change_to_right,change_trip_frames(change)];
    elseif (trip_type(change_trip_frames(change)+1) == 1 | trip_type(change_trip_frames(change)+1) == 3)
        change_to_left = [change_to_left,change_trip_frames(change)];
    end
end

%Get the times of changes in trial_counter
change_trial = diff(trial_counter);
change_trial_frames = find(change_trial>0.5);

subplot(2,1,2)
%Plot the panel's y dimension
plot(timestamp,voltage_out)
hold on
%Add the change in trip
for changeInTrip = 1:length(change_to_right)
    plot(timestamp(change_to_right(changeInTrip)),voltage_out(change_to_right(changeInTrip)),'ro') 
end
for changeInTrip = 1:length(change_to_left)
    plot(timestamp(change_to_left(changeInTrip)),voltage_out(change_to_left(changeInTrip)),'bo') 
end
%Add text for the trip counter
for changeInTrip = 1:length(change_trial_frames)
    text(timestamp(change_trial_frames(changeInTrip)),9,num2str(trial_counter(change_trial_frames(changeInTrip)+1)))
end
xlabel('Time (sec)');
ylabel('Y panels voltage (V)');
ylim([0 10]);
title('Python data');

%% Use the signal from the turns in the y panels voltage to determine the changes in trip

%convert from frames to seconds
time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(rawData.trial_bdata));

%smooth a little the x panel data to better detect the voltage thresholds
panelxData = movmedian(rawData.trial_bdata(:,xPanels),500);

%using the panel data, get the turns
Test = diff(panelxData);
Test(abs(Test) >1) = 0;
turning_times = abs(Test)>0.04;

%Find the start and end of each turn by looking at the number of zeroes
%before and after
turn_points = zeros(length(turning_times),1);
for i = 8000:length(turning_times)-7999
    if ((turning_times(i) == 1) & ((sum(turning_times(i-7999:i-1))==0) | (sum(turning_times(i+1:i+7999))==0)))
        turn_points(i) = 1;
    else
        turn_points(i) = 0;
    end
end
turnPoints = find(turn_points > 0.5);

%Plot to see if they're accurate
figure,
plot(time,panelxData)
hold on
plot(time(turnPoints),panelxData(turnPoints),'r.')

%get the trips as trip 1 = start:turnpoint1,
%trip n = (nx2)-2:(nx2)
tripData = {};
tripData{1} = rawData.trial_bdata(1:turnPoints(1),:);
for i = 2:floor(length(turnPoints)/2)
   tripData{i} = rawData.trial_bdata(turnPoints(2*i-2):turnPoints(2*i-1),:);
end

%Plot the x panel bouts of the trip bouts: they should all look flat
figure('Position',[100 200 1400 800]),
for i = 1:length(tripData)
    subplot(5,ceil(length(tripData)/5),i)
    plot(tripData{i}(:,4))
    title(['Trip # ',num2str(i)]);
end
suptitle('Changes in panel x dimension during the trips');

%get the turns as tn = n*2-1:n*2
turnData = {};
for i = 1:floor(length(turnPoints)/2)
   turnData{i} = rawData.trial_bdata(turnPoints(2*i-1):turnPoints(2*i),:);
end


%Plot the x panel bouts of the turn bouts: they should all look oblique
figure('Position',[100 200 1400 800]),
for i = 1:length(turnData)
    subplot(5,ceil(length(turnData)/5),i)
    plot(turnData{i}(:,4))
    title(['Turn # ',num2str(i)]);
end
suptitle('Changes in panel x dimension during the turns');

%% Add possibility to remove or add turning points by clicking




%% Merge the left360_firsthalf and left360_secondhalf using the panels data.

%Get the turning range as the voltage range between the start and end of
%each turn
turning_range = [];
for i = 1:length(turnData)
    turning_range(i) = range(turnData{1,i}(500:end-500,4));
end

%find 360 turns by finding those turns where the voltage had a bigger range
%(close to 10 instead of 5 V)
full_turns = find(turning_range > 6);

%tag the trips preceding those turns as left360
for trip = 1:length(full_turns)
   left360(trip) = full_turns(trip)-trip+1;   
end

%combine the first and second half of left360 trips
TripData = {};
TripData{1} = tripData{1};

for trip = 2:length(tripData)
   if ~isempty(find(full_turns == trip))
       TripData{trip} = [tripData{trip};tripData{trip+1}];
   elseif ~isempty(find(full_turns == trip-1))
       TripData{trip} = [];
   else
       TripData{trip} = tripData{trip};
   end
end

TripData = TripData(~cellfun('isempty',TripData));



%% Separate the full data into trips and smooth it

sizeBall = 9;

for i = 1:length(TripData)
    
    %Panels data
    data{i}.xPanelVolts =  TripData{1,i}(:,xPanels); 
    VOLTAGE_RANGE_x = 10;
    maxValX =  96; %our pattern has 96 x dimensions

    data{i}.yPanelVolts =  TripData{1,i}(:, yPanels);
    VOLTAGE_RANGE_y = 10;
    maxValY = 92; %our pattern has 92 y dimensions
           
    %FicTrac data
    data{i}.ficTracAngularPosition = TripData{1,i}( : , yawFly); 
    data{i}.ficTracIntx = TripData{1,i}( : , xFly); 

    %smooth the data and get the velocity
    smoothed{i} = singleTrialVelocityAnalysis(data{i},samplingRate);

    %downsample and split into trips other useful data
    ypanels{i} = resample(TripData{1,i}(:,yPanels),25,samplingRate);
    xpanels{i} = resample(TripData{1,i}(:,xPanels),25,samplingRate);
    optoPulse{i} = resample(TripData{1,i}(:,OptoTrigger),25,samplingRate);
    
    %calculate the trial duration
    trialDuration(i) = length(data{i}.xPanelVolts)/samplingRate;
    
end


%% Trial dur vs trial number

probeTrials = sort(str2num(rawData.probe_trials)); %get probe trials from run object
trial_num = length(TripData);
trialNum = 1:trial_num;
optoTrials = sort(setdiff(trialNum,probeTrials)); %get opto trials as the non-probe trials

probeTrials = probeTrials(probeTrials<=trial_num); %remove trials that excede the total number of trials the fly got
probeTrials = probeTrials(probeTrials > 5); %remove the first 5 trials since she didn't get an opto experience up to that point
optoTrials = optoTrials(optoTrials<=trial_num); %remove trials that excede the total number of trials the fly got

PT = ones(1,length(probeTrials));
OT = ones(1,length(optoTrials));
group = [PT, 2*OT];
durPT = trialDuration(probeTrials);
durOT = trialDuration(optoTrials);

% Plot the trial duration in time
figure%('Position', [400 200 1000 800])
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


%% Determine 'working opto pulses'

%Define working opto trials as though in which the LED pulse was actually
%triggered, using the 'optoPulse' signal

for i = 1:length(optoTrials)    
    if (sum(diff(TripData{1,optoTrials(i)}(:,OptoTrigger))<-1)>0 && sum(diff(TripData{1,optoTrials(i)}(:,OptoTrigger))>1)>0)
        WorkingOptoPulse(i) = optoTrials(i); %keep the trial
        startPulse(i) = floor((find(diff(TripData{1,optoTrials(i)}(:,OptoTrigger))>1,1))*25/4000); %find the start point of the pulse
        EndPulse{i} = find(diff(TripData{1,optoTrials(i)}(:,OptoTrigger))<-1); %find the end point of the pulse
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


%% Plot forward velocity around the time of the pulse in the opto trials

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

%% Identify long left and right trips

right360 = left360 + 1;

%% Plot forward velocity vs distance for each trial, along with trial type, reward location, and whether it was opto or probe

%load the pattern voltages to then define the voltage steps for each
%distance point
load('Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp20\source\utils\stimulusVoltages.mat');

%I'm using 92 dimensions in y, so there are 92 possible voltage steps the
%panels can take
%bin the forward velocity data into voltages using the ypanel data
for i = 1:length(smoothed)
    for j = 1:length(medianVoltages)
        step{i,j} = find(ypanels{1,i} > medianVoltages(j)-0.04 & ypanels{1,i} < medianVoltages(j)+0.04);
        stepMeans(i,j) = nanmedian(smoothed{1,i}.xVel(step{i,j}));
    end
end


%Linearly interpolate the data to fill in the gaps
for i = 1:size(stepMeans,1)
    if (any(left360 == i) | any(right360 == i)) %for long left and right trips
        StepMeans{i} = stepMeans(i,:);
    else
        StepMeans{i} = stepMeans(i,1:46); %If this is not one of the longer trials, restrict the voltage between 0 and 5 V
    end
    StepMeans{i} = fillmissing(StepMeans{i},'linear');
end

%Plot individual trips
ydimension = [1:92];
distance = 9*pi*(ydimension/92)/run_obj.run_obj.gain_x;
%get reward location
rewardDimension = round(run_obj.run_obj.reward_distance*46/100);

for trip = 1:length(smoothed)
    figure,
    if ~isempty(find(probeTrials == trip))
        plot(distance(1:length(StepMeans{trip})),StepMeans{trip},'-ko')
        hold on
        %xline(distance(rewardDimension),'lineWidth',2,'color','b');
        line([distance(rewardDimension) distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
        ylabel('Forward velocity (mm/s)');
        xlabel('Distance (mm)');
        title(['Trip #',num2str(trip)]);
    else
        plot(distance(1:length(StepMeans{trip})),StepMeans{trip},'-ro')
        ylabel('Forward velocity (mm/s)');
        hold on
        if any(right360 == trip)
            %xline(distance(rewardDimension)+15,'lineWidth',2,'color','b');
            line([distance(rewardDimension)+15 distance(rewardDimension)+15], [0 5],'lineWidth',2,'color','b');
        else
            %xline(distance(rewardDimension),'lineWidth',2,'color','b');  
            line([distance(rewardDimension) distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
        end
        xlabel('Distance (mm)');
        title(['Trip #',num2str(trip)]);
    end
end


%% Plotting velocity as a function of distance travelled

WorkingProbeTrials = probeTrials; %in this configuration of the hallway, all the probetrials are useful trials


%bin the forward velocity data into voltages using the ypanel data
for i = 1:length(smoothed)
    if ~(any(left360 == i) | any(right360 == i)) %for regular left and right trips
        for j = 1:46
            step{i,j} = find(ypanels{1,i} > medianVoltages(j)-0.04 & ypanels{1,i} < medianVoltages(j)+0.04);
            stepMeansRegularTrials(i,j) = mean(smoothed{1,i}.xVel(step{i,j}));
            angularStepMeansRegularTrials(i,j) = mean(abs(smoothed{1,i}.angularVel(step{i,j})));
            stepMeansSpecialTrials(i,j) = NaN;
            angularStepMeansSpecialTrials(i,j) = NaN;
        end                
    else
        for j = 1:92 %for long left and right trips
            step{i,j} = find(ypanels{1,i} > medianVoltages(j)-0.04 & ypanels{1,i} < medianVoltages(j)+0.04);
            stepMeansSpecialTrials(i,j) = mean(smoothed{1,i}.xVel(step{i,j}));
            angularStepMeansSpecialTrials(i,j) = mean(abs(smoothed{1,i}.angularVel(step{i,j})));
            stepMeansRegularTrials(i,j) = NaN;
            angularStepMeansRegularTrials(i,j) = NaN;
        end
    end
end

%Linearly interpolate the data to fill in the gaps
for i = 1:size(stepMeansRegularTrials,1)
    stepMeansRegularTrials(i,1:46) = fillmissing(stepMeansRegularTrials(i,1:46),'linear');
    angularStepMeansRegularTrials(i,1:46) = fillmissing(angularStepMeansRegularTrials(i,1:46),'linear');
end

%Remove everything past dimension 46
stepMeansRegularTrials = stepMeansRegularTrials(:,1:46);
angularStepMeansRegularTrials = angularStepMeansRegularTrials(:,1:46);



%Divide trips in left and right
%define the first one as left (=1)
%alternate left:right (1:0) for the following
tripType = zeros(1,length(TripData));
for i = 1:length(TripData)
    if ((mod(i,2) ~= 0) & (~any(left360 == i)))
        tripType(i) = 1;
    elseif ((mod(i,2) ~= 0) & (any(left360 == i)))
        tripType(i) = 3;
    elseif ((mod(i,2) == 0) & (~any(right360 == i)))
        tripType(i) = 0;
    else
        tripType(i) = 2;
    end
end

WorkingOptoLeft = WorkingOptoPulse(tripType(WorkingOptoPulse) == 1);
WorkingProbeLeft = WorkingProbeTrials(tripType(WorkingProbeTrials) == 1);
WorkingOptoRight = WorkingOptoPulse(tripType(WorkingOptoPulse) == 0);
WorkingProbeRight = WorkingProbeTrials(tripType(WorkingProbeTrials) == 0);

%Get medians and error for each trip type
optoLeftMeans = nanmedian(stepMeansRegularTrials(WorkingOptoLeft,:));
probeLeftMeans = nanmedian(stepMeansRegularTrials(WorkingProbeLeft,:));
optoLeftError = (nanstd(stepMeansRegularTrials(WorkingOptoLeft,:)))./sqrt(length(WorkingOptoLeft));
probeLeftError = (nanstd(stepMeansRegularTrials(WorkingProbeLeft,:)))./sqrt(length(WorkingProbeLeft));

optoRightMeans = nanmedian(stepMeansRegularTrials(WorkingOptoRight,:));
probeRightMeans = nanmedian(stepMeansRegularTrials(WorkingProbeRight,:));
optoRightError = (nanstd(stepMeansRegularTrials(WorkingOptoRight,:)))./sqrt(length(WorkingOptoRight));
probeRightError = (nanstd(stepMeansRegularTrials(WorkingProbeRight,:)))./sqrt(length(WorkingProbeRight));

optoLeftAngMeans = nanmedian(angularStepMeansRegularTrials(WorkingOptoLeft,:));
probeLeftAngMeans = nanmedian(angularStepMeansRegularTrials(WorkingProbeLeft,:));
optoLeftAngError = (nanstd(angularStepMeansRegularTrials(WorkingOptoLeft,:)))./sqrt(length(WorkingOptoLeft));
probeLeftAngError = (nanstd(angularStepMeansRegularTrials(WorkingProbeLeft,:)))./sqrt(length(WorkingProbeLeft));

optoRightAngMeans = nanmedian(angularStepMeansRegularTrials(WorkingOptoRight,:));
probeRightAngMeans = nanmedian(angularStepMeansRegularTrials(WorkingProbeRight,:));
optoRightAngError = (nanstd(angularStepMeansRegularTrials(WorkingOptoRight,:)))./sqrt(length(WorkingOptoRight));
probeRightAngError = (nanstd(angularStepMeansRegularTrials(WorkingProbeRight,:)))./sqrt(length(WorkingProbeRight));

%Get the distance data
ydimension = [1:46]; %number of stimulus dimensions we're using
%convert the ydimensions to distance using the xgain
distance = 9*pi*(ydimension/92)/run_obj.run_obj.gain_x;
%get reward location
rewardDimension = round(run_obj.run_obj.reward_distance*46/100);

figure('Position', [100 100 1600 1000]),
ax1 = subplot(2,2,1);
[opT] = boundedline(distance,optoLeftMeans,optoLeftError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeLeftMeans,probeLeftError,'-ko','alpha');
%xline(distance(rewardDimension),'lineWidth',2,'color','b');
line([distance(rewardDimension) distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
ylabel('Forward velocity (mm/s)');
title('Left trips')

ax2 = subplot(2,2,2);
[opT] = boundedline(distance,flip(optoRightMeans),flip(optoRightError),'-ro','alpha');
hold on
[prT] = boundedline(distance,flip(probeRightMeans),flip(probeRightError),'-ko','alpha');
%xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
line([distance(end)-distance(rewardDimension) distance(end)-distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');


title('Right trips')

ax3 = subplot(2,2,3);
[opT] = boundedline(distance,optoLeftAngMeans,optoLeftAngError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeLeftAngMeans,probeLeftAngError,'-ko','alpha');
%xline(distance(rewardDimension),'lineWidth',2,'color','b');
line([distance(rewardDimension) distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
ylabel('Angular speed (deg/s)');
xlabel('Distance travelled (mm)');

ax4 = subplot(2,2,4);
[opT] = boundedline(distance,flip(optoRightAngMeans),flip(optoRightAngError),'-ro','alpha');
hold on
[prT] = boundedline(distance,flip(probeRightAngMeans),flip(probeRightAngError),'-ko','alpha');
%xline(distance(end)-distance(rewardDimension),'lineWidth',2,'color','b');
line([distance(end)-distance(rewardDimension) distance(end)-distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
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

%% Combine all the trips, since the reward is in the middle

%combine trips
optoTrips = [stepMeansRegularTrials(WorkingOptoLeft,:);flip(stepMeansRegularTrials(WorkingOptoRight,:))];
probeTrips = [stepMeansRegularTrials(WorkingProbeLeft,:);flip(stepMeansRegularTrials(WorkingProbeRight,:))];
optoAngTrips = [angularStepMeansRegularTrials(WorkingOptoLeft,:);flip(angularStepMeansRegularTrials(WorkingOptoRight,:))];
probeAngTrips = [angularStepMeansRegularTrials(WorkingProbeLeft,:);flip(angularStepMeansRegularTrials(WorkingProbeRight,:))];

%get means and error
optoMeans = nanmean(optoTrips);
optoError = (nanstd(optoTrips./sqrt(length(optoTrips))));
probeMeans = nanmean(probeTrips);
probeError = (nanstd(probeTrips./sqrt(length(probeTrips))));

optoAngMeans = nanmean(optoAngTrips);
optoAngError = (nanstd(optoAngTrips./sqrt(length(optoAngTrips))));
probeAngMeans = nanmean(probeAngTrips);
probeAngError = (nanstd(probeAngTrips./sqrt(length(probeAngTrips))));

figure,
ax1 = subplot(2,1,1);
[opT] = boundedline(distance,optoMeans,optoError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeMeans,probeError,'-ko','alpha');
%xline(distance(rewardDimension),'lineWidth',2,'color','b');
line([distance(rewardDimension) distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
ylabel('Forward velocity (mm/s)');

ax1 = subplot(2,1,2);
[opT] = boundedline(distance,optoAngMeans,optoAngError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeAngMeans,probeAngError,'-ko','alpha');
%xline(distance(rewardDimension),'lineWidth',2,'color','b');
line([distance(rewardDimension) distance(rewardDimension)], [0 100],'lineWidth',2,'color','b');
ylabel('Angular speed (deg/s)');


%% Plot for each trip what happens temporally and what happens distance wise

%CHECK FLIPPING!


for trip = 1:length(WorkingOptoPulse)
    figure,
    %Time wise
    subplot(2,1,1)
    plot(smoothed{1,WorkingOptoPulse(trip)}.xVel)
    hold on
    line([startPulse(trip) startPulse(trip)],[0 5], 'color', 'r', 'LineWidth',2)
    ylim([-5 5]);
    
    %Distance wise
    ydimension = [1:92];
    distance = 9*pi*(ydimension/92)/run_obj.run_obj.gain_x;
    subplot(2,1,2)
    %if this a right trip, flip the data
    if (tripType(WorkingOptoPulse(trip) == 0) | tripType(WorkingOptoPulse(trip) == 2))
        plot(distance(1:length(StepMeans{WorkingOptoPulse(trip)})),flip(StepMeans{WorkingOptoPulse(trip)}),'-ko')        
    %if not, don't
    else
        plot(distance(1:length(StepMeans{WorkingOptoPulse(trip)})),StepMeans{WorkingOptoPulse(trip)},'-ko')
    end
    ylim([-5 5]);
    if any(right360 == WorkingOptoPulse(trip))
        line([distance(rewardDimension)+15 distance(rewardDimension)+15], [0 5],'lineWidth',2,'color','b');
    else
        line([distance(rewardDimension) distance(rewardDimension)], [0 5],'lineWidth',2,'color','b');
    end
    
end