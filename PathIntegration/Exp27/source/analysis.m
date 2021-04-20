%Code to analyze exp 27: the virtual hallway that spins 180 when the fly reaches
%the end, gives the reward in the same location for both trial types and
%the gain is 1 (shorter hallway)
close all; clear all;

% prompt the user to select the file to open and load it.
cd 'Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp27\data'

[fileName,pathName] = uigetfile();   
rawData = load([pathName,fileName]);

sid = fileName(end-10);

%load the run object
run_obj_dir = [pathName,'runobj'];
cd (run_obj_dir)
run_obj_name = dir(strcat('*',sid,'_runobj.mat'));
run_obj = load(strcat(run_obj_dir,'\',run_obj_name.name));

%% Make directory to save plots if it doesn't already exist
  
%Move to the analysis folder
cd(pathName)
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(pathName,'plots'); 
end
%List the contents of the 'plots' folder
cd([pathName,'plots\'])

%%  Define Ni-Daq channels ID
  
headingFly = 1;
%yFly = 2; %right now we're not collecting the y signal
%through the NiDaq
yawFly =  2;
xFly = 3;
xFlyGain = 7;
xPanels = 4;
yPanels = 5;
PanelStatus = 6;
OptoTrigger = 8;
CameraTrigger = 9;

samplingRate = 4000;
% 
% 
% %% Check if the change in trip type and trip counter make sense and if matlab and python data match
% 
% figure('Position',[100 100 1600 800]),
% %plot the matlab data
% subplot(2,1,1)
% %Plot the panel's y dimension
% plot(rawData.trial_time,rawData.trial_bdata(:,yPanels))
% hold on
% %Add the change in trip
% for changeInTrip = 1:length(rawData.trialData)
%     [M I] = min(abs(rawData.trial_time-rawData.trialData(1,changeInTrip)));
%     if rawData.trialData(2,changeInTrip) == 2
%         plot(rawData.trialData(1,changeInTrip),rawData.trial_bdata(I,yPanels),'ro')
%     else
%         plot(rawData.trialData(1,changeInTrip),rawData.trial_bdata(I,yPanels),'bo')
%     end    
%     
% end
% %Add text for the trip counter, every 10 trips
% for changeInTrip = 1:10:length(rawData.trialData)
%     text(rawData.trialData(1,changeInTrip),9,num2str(rawData.trialData(3,changeInTrip)))
% end
% xlabel('Time (sec)');
% ylabel('Y panels voltage (V)');
% ylim([0 10]);
% xlim([0 rawData.trial_time(end)]);
% title('NiDaq data');
% %add custom legend
% h = zeros(2, 1);
% h(1) = plot(NaN,NaN,'ro');
% h(2) = plot(NaN,NaN,'bo');
% legend(h, 'Right trip start','Left trip start');
% 
% % Compare to the hdf5 file to see if python and matlab are matching
% 
%Get file name
if exist('fileList','var') == 0
    fileList = dir(fullfile(pathName, '*.hdf5'));
    for file = 1:length(fileList)
        if contains(fileList(file).name,fileName(end-14:end-6))
            file_name{file} = [pathName,fileList(file).name];
        else
            file_name{file} = [];
        end
    end
    file_name = char(file_name(~cellfun(@isempty, file_name)));
end

%Import relevant variables
timestamp = h5read(file_name,'/time');
voltage_out = h5read(file_name,'/output_voltage_x_gain');
trip_type = h5read(file_name,'/trip_type');
trial_counter = h5read(file_name,'/trial_counter');
% 
% %Get the times of changes in trip_type
% change_trip = diff(trip_type);
% change_trip_frames = find(abs(change_trip)>0.5);
% change_to_right = [];
% change_to_left = [];
% for change = 1:length(change_trip_frames)
%     if trip_type(change_trip_frames(change)+1) == 2
%         change_to_right = [change_to_right,change_trip_frames(change)];
%     elseif (trip_type(change_trip_frames(change)+1) == 1 | trip_type(change_trip_frames(change)+1) == 3)
%         change_to_left = [change_to_left,change_trip_frames(change)];
%     end
% end
% 
% %Get the times of changes in trial_counter
% change_trial = diff(trial_counter);
% change_trial_frames = find(change_trial>0.5);
% 
% subplot(2,1,2)
% %Plot the panel's y dimension
% plot(timestamp,voltage_out)
% hold on
% %Add the change in trip
% for changeInTrip = 1:length(change_to_right)
%     plot(timestamp(change_to_right(changeInTrip)),voltage_out(change_to_right(changeInTrip)),'ro') 
% end
% for changeInTrip = 1:length(change_to_left)
%     plot(timestamp(change_to_left(changeInTrip)),voltage_out(change_to_left(changeInTrip)),'bo') 
% end
% %Add text for the trip counter
% for changeInTrip = 1:10:length(change_trial_frames)
%     text(timestamp(change_trial_frames(changeInTrip)),9,num2str(trial_counter(change_trial_frames(changeInTrip)+1)))
% end
% xlim([0 timestamp(end)]);
% xlabel('Time (sec)');
% ylabel('Y panels voltage (V)');
% ylim([0 10]);
% title('Python data');
% %add custom legend
% h = zeros(2, 1);
% h(1) = plot(NaN,NaN,'ro');
% h(2) = plot(NaN,NaN,'bo');
% legend(h, 'Right trip start','Left trip start');
% 
% saveas(gcf,[pathName,'plots\checkingTripData_sid_',sid,'.png'])

%% Use the signal from the turns in the y panels voltage to determine the changes in trip

%convert from frames to seconds
time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(rawData.trial_bdata));

%smooth a little the x panel data to better detect the voltage thresholds
smooth_xpanel_data = movmedian(rawData.trial_bdata(:,xPanels),500);

%set the conditions and get the turns: turns are hapenning when the x panel voltage
%is between 9.84 and 5, or between 0 and 4.84
turns = (smooth_xpanel_data > 5 & smooth_xpanel_data < 9.84) | (smooth_xpanel_data > 0 & smooth_xpanel_data < 4.84);

%get the turn points taking the derivative and look for peaks bigger than
%0.5
diff_turns = abs(diff(turns));
diff_turns = [diff_turns;diff_turns(end)];
turnPoints = find(diff_turns > 0.5);

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


% %Plot the x panel bouts of the trip bouts: they should all look flat
% figure('Position',[100 100 1600 800]),
% for i = 1:length(tripData)
%     subplot(5,ceil(length(tripData)/5),i)
%     plot(tripData{i}(:,4),'linewidth',2)
%     title(['#',num2str(i)]);
%     set(gca,'xticklabel',{[]});
%     set(gca,'yticklabel',{[]});    
%     ylim([0 11]);
% end
% suptitle('Changes in panel x dimension during the trips');
% saveas(gcf,[pathName,'plots\tripVoltageData_sid_',sid,'.png'])
% 
% 
% %Plot the x panel bouts of the trip bouts: they should all look oblique
% figure('Position',[100 100 1600 800]),
% for i = 1:length(turnData)
%     subplot(5,ceil(length(turnData)/5),i)
%     plot(turnData{i}(:,4),'linewidth',2)
%     title(['#',num2str(i)]);
%     set(gca,'xticklabel',{[]});
%     set(gca,'yticklabel',{[]});    
%     ylim([0 11]);
% end
% suptitle('Changes in panel x dimension during the turns');
% 
% saveas(gcf,[pathName,'plots\turnVoltageData_sid_',sid,'.png'])

%% Separate the full data into trips, smooth it and get velocity

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
    data{i}.ficTracIntx = tripData{1,i}( : , xFly); %I use x and not xgain to get the actual velocity of the fly

    %smooth the data and get the velocity
    %smoothed{i} = singleTrialVelocityAnalysis(data{i},samplingRate);
    nonsmoothed{i} = singleTrialVelocityAnalysis2(data{i},samplingRate);

    ypanels{i} = resample(tripData{1,i}(:,yPanels),25,samplingRate);
    optoPulse{i} = resample(tripData{1,i}(:,OptoTrigger),25,samplingRate);
    trialDuration(i) = length(data{i}.xPanelVolts)/samplingRate;
end


%% Trial dur vs trial number

%get probe and opto trials
probeTrials = sort(str2num(rawData.probe_trials));
trial_num = length(tripData);
trialNum = 1:trial_num;
%compute the opto trials as all the trials that were not probe trials
optoTrials = sort(setdiff(trialNum,probeTrials));
probeTrials = probeTrials(probeTrials<=length(trialDuration)); %remove trials that excede the total number of trials the fly got
probeTrials = probeTrials(probeTrials > 6); %remove the first 6 trials since she didn't get an opto experience up to that point
optoTrials = optoTrials(optoTrials<=length(trialDuration)); %remove trials that excede the total number of trials the fly got

%Get trial length
PT = ones(1,length(probeTrials));
OT = ones(1,length(optoTrials));
group = [PT, 2*OT];
durPT = trialDuration(probeTrials);
durOT = trialDuration(optoTrials);

% Plot the trial duration in time
figure('Position', [100 100 1000 800])
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

saveas(gcf,[pathName,'plots\TrialDurVsTrialNum_sid_',sid,'.png'])

%% Look at the velocity temporally in individual trials

%make a figure every 10 trials
for group = 1:ceil((length(optoTrials)+length(probeTrials))/10)
    figure('Position',[50 50 300 900]),
    suptitle(['Trials ',num2str((group-1)*10+1),'to ',num2str(group*10)]);
    for trial = 1:10
        subplot(10,1,trial)
        if any(optoTrials==(group-1)*10+trial)
            plot(nonsmoothed{1,(group-1)*10+trial}.xVel,'k')
            hold on
            plot(optoPulse{1,(group-1)*10+trial}*2,'r')
        else
            plot(nonsmoothed{1,(group-1)*10+trial}.xVel,'k')
        end    
        xlim([0 length(nonsmoothed{1,(group-1)*10+trial}.xVel)]);
    end
    saveas(gcf,[pathName,'plots\ind_trials_temporally_group',num2str(group),'.png'])
end
close all;

%% Divide trips in left and right
%define the first one as left (=1)
%alternate left:right (1:0) for the following
tripType = zeros(1,length(tripData));
for i = 1:length(tripData)
    if (mod(i,2) ~= 0) 
        tripType(i) = 1;
    else
        tripType(i) = 2;
    end
end

% %% Remove outlier that have a really big variance in the xVel, because it probably indicates an artifact
% 
% for i = 1:size(nonsmoothed,2)
%     varian(i) = nanvar(nonsmoothed{1,i}.xVel);
% end
% 
% %exclude trials with var > 18 in the xVel.
% excludetrials=[];
% for i = 1:size(nonsmoothed,2)
%     if varian(i) > 18
%         excludetrials=[excludetrials i];
%     end
% end
% 
% %remove outlier trials 
% optoTrials=setdiff(optoTrials,excludetrials);
% probeTrials=setdiff(probeTrials,excludetrials);
% 
% %maybe instead, or on top of this, I should remove trials that have very
% %high velocity values (for example above 20). The smoothing might be
% %artificially increasing these.

%% Look at something analogous to spike triggered average

%looking at the opto pulse, there are sometimes trials that should be opto
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
    if (startPulse(i)>13 & length(nonsmoothed{1,WorkingOptoPulse(i)}.xVel)>startPulse(i)+25) 
        velAroundPulse{i} = nonsmoothed{1,WorkingOptoPulse(i)}.xVel(startPulse(i)-13:startPulse(i)+25);
    elseif (startPulse(i)<=13 & length(nonsmoothed{1,WorkingOptoPulse(i)}.xVel)>startPulse(i)+25) 
        velAroundPulse{i} = nonsmoothed{1,WorkingOptoPulse(i)}.xVel(1:startPulse(i)+25);
        velAroundPulse{i} = [repelem(nan,1,14-startPulse(i)),velAroundPulse{i}];
    elseif (startPulse(i)>13 & length(nonsmoothed{1,WorkingOptoPulse(i)}.xVel)<=startPulse(i)+25) 
        velAroundPulse{i} = nonsmoothed{1,WorkingOptoPulse(i)}.xVel(startPulse(i)-13:end);
        velAroundPulse{i} = [velAroundPulse{i},repelem(nan,1,25-(length(nonsmoothed{1,WorkingOptoPulse(i)}.xVel)-startPulse(i)))];       
    end
end

%remove empty cells (fix later!)
velAroundPulse = velAroundPulse(~cellfun('isempty',velAroundPulse));


%find trials in which the fly received the reward less than 0.5 sec into
%the trial and add nans to compensate for those timepoints and align to the
%rest of the trials
trialL = cellfun(@length,velAroundPulse);
% shortT = find(trialL<39);
% for i = 1:length(shortT)
%     velAroundPulse{shortT(i)} = [repelem(nan,1,14-startPulse(shortT(i))),velAroundPulse{shortT(i)}];
% end    
velAroundPulse = cell2mat(velAroundPulse);
velAroundPulse = reshape(velAroundPulse,39,length(trialL));

%compute the mean in 3 different sections around the pulse
meanVelBeforePulse = nanmean(velAroundPulse(1:13,:));
meanVelDuringPulse = nanmean(velAroundPulse(14:26,:));
meanVelAfterPulse = nanmean(velAroundPulse(27:39,:));
meanVelAroundPulse = [meanVelBeforePulse',meanVelDuringPulse',meanVelAfterPulse'];

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
plot(meanVelAroundPulse','color',[0.5 0.5 0.5])
hold on
err = std(meanVelAroundPulse);
errorbar([1,2,3],mean(meanVelAroundPulse),err,'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
xticks([1 2 3])
xticklabels({'1','2','3'})
xlim([0 4]);
suptitle('Forward velocity around the pulse');

saveas(gcf,[pathName,'plots\VelAroundPulseOpto_sid_',sid,'.png'])

%% Repeat the vel around opto stim overlaying individual trials to see their trend

figure
time = linspace(-0.5,1,39);
for trial = 1:size(velAroundPulse,2)
    plot(time,velAroundPulse(:,trial),'color',[.5,.5,.5])
    hold on
end
plot(time,nanmean(velAroundPulse,2),'k','linewidth',2)
line([0 0],[min(min(velAroundPulse))-0.3 max(max(velAroundPulse))+0.3],'Color','r','LineWidth',2)
line([0.5 0.5],[min(min(velAroundPulse))-0.3 max(max(velAroundPulse))+0.3],'Color','r','LineWidth',2)
ylim([min(min(velAroundPulse))-0.3 max(max(velAroundPulse))+0.3]);
ylabel('Forward velocity (mm/s)');
xlabel('Time around the pulse (s)');
xlim([-0.5 1]);
text(-0.25,max(max(nanmean(velAroundPulse,2)))+0.2,'1','FontSize',14)
text(0.25,max(max(nanmean(velAroundPulse,2)))+0.2,'2','FontSize',14)
text(0.75,max(max(nanmean(velAroundPulse,2)))+0.2,'3','FontSize',14)

saveas(gcf,[pathName,'plots\VelAroundPulseOpto_indTrials_sid_',sid,'.png'])

%% Binning data as a function of distance travelled

WorkingProbeTrials = probeTrials; 
%load the stimulus voltages to obtain the voltage steps for each dimension
load('Z:\Wilson Lab\Mel\Experiments\PathIntegration\Exp20\source\utils\stimulusVoltages.mat');

%bin the forward velocity data into voltages using the ypanel data
for i = 1:length(nonsmoothed)
    for j = 1:48
        step{i,j} = find(ypanels{1,i} > medianVoltages(j)-0.04 & ypanels{1,i} < medianVoltages(j)+0.04);
        stepMeans(i,j) = mean(nonsmoothed{1,i}.xVel(step{i,j}));
        angularStepMeans(i,j) = mean(abs(nonsmoothed{1,i}.angularVel(step{i,j})));
    end
end

%Linearly interpolate the data to fill in the gaps
for i = 1:size(stepMeans,1)
    stepMeans(i,:) = fillmissing(stepMeans(i,:),'linear');
    angularStepMeans(i,:) = fillmissing(angularStepMeans(i,:),'linear');
end

%flip the data for right trips
for trip = 1:length(tripData)
    if tripType(trip) == 2
        stepMeans(trip,:) = flip(stepMeans(trip,:));
        angularStepMeans(trip,:) = flip(angularStepMeans(trip,:));
    end
end

%Convert the data to distance in relevant units
ydimension = [1:48]; %number of stimulus dimensions we're using  (check if this is 48 or 46!)
%convert the ydimensions to distance using the xgain
distance = 9*pi*(ydimension/92)/run_obj.run_obj.gain_x;
%get reward location
rewardDimension = round(run_obj.run_obj.reward_distance*48/100);

%% Plot the velocity vs the distance for individual trips

%make a figure every 10 trials
for group = 1:ceil((length(WorkingOptoPulse)+length(WorkingProbeTrials))/10)
    figure('Position',[50 50 800 1000]),
    suptitle(['Trials ',num2str((group-1)*10+1),'to ',num2str(group*10)]);
    for trial = 1:10
        subplot(10,2,trial*2-1)
        if any(WorkingOptoPulse==(group-1)*10+trial)
            plot(distance,stepMeans((group-1)*10+trial,:),'-ro')
            hold on
            xline(distance(rewardDimension),'lineWidth',2,'color','b');
        else
            plot(distance,stepMeans((group-1)*10+trial,:),'-ko')
            xline(distance(rewardDimension),'lineWidth',2,'color','b','LineStyle','--');
        end
        
        if trial == 1
            title('Forward velocity (mm/s)');
        elseif trial == 10
            xlabel('Distance (mm)');
        end
        
        subplot(10,2,trial*2)
        if any(WorkingOptoPulse==(group-1)*10+trial)
            plot(distance,angularStepMeans((group-1)*10+trial,:),'-ro')
            hold on
            xline(distance(rewardDimension),'lineWidth',2,'color','b');
        else
            plot(distance,angularStepMeans((group-1)*10+trial,:),'-ko')
            xline(distance(rewardDimension),'lineWidth',2,'color','b','LineStyle','--');
        end
        if trial == 1
            title('Angular speed (deg/s)');
        elseif trial == 10
            xlabel('Distance (mm)');
        end
        
    end
    
    saveas(gcf,[pathName,'plots\ind_trials_distance_group',num2str(group),'.png'])
end

close all;
%% Plot the mean trends of probe and opto trials, pooling left and right trips since the reward is at 50%

optoMeans = nanmedian(stepMeans(WorkingOptoPulse,:));
probeMeans = nanmedian(stepMeans(WorkingProbeTrials,:));
optoError = (nanstd(stepMeans(WorkingOptoPulse,:)))./sqrt(length(WorkingOptoPulse));
probeError = (nanstd(stepMeans(WorkingProbeTrials,:)))./sqrt(length(WorkingProbeTrials));

optoAngMeans = nanmedian(angularStepMeans(WorkingOptoPulse,:));
probeAngMeans = nanmedian(angularStepMeans(WorkingProbeTrials,:));
optoAngError = (nanstd(angularStepMeans(WorkingOptoPulse,:)))./sqrt(length(WorkingOptoPulse));
probeAngError = (nanstd(angularStepMeans(WorkingProbeTrials,:)))./sqrt(length(WorkingProbeTrials));

figure('Position', [100 100 1400 800]),
ax1 = subplot(1,2,1);
[opT] = boundedline(distance,optoMeans,optoError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeMeans,probeError,'-ko','alpha');
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Forward velocity (mm/s)');
title('All trips')
xlabel('Distance (mm)');

ax1 = subplot(1,2,2);
[opT] = boundedline(distance,optoAngMeans,optoAngError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeAngMeans,probeAngError,'-ko','alpha');
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Angular speed (deg/s)');
title('All trips')
xlabel('Distance (mm)');

saveas(gcf,[pathName,'plots\meanVelVsDistance_sid_',sid,'.png'])

%save data
save([pathName,'distanceData.mat'],'optoMeans','probeMeans','optoAngMeans','probeAngMeans')

%% Compare to shifted control

figure('Position',[100 100 1400 1000]),
subplot(1,2,1)
shiftedControl(distance,stepMeans,WorkingOptoPulse,rewardDimension,1,1)
subplot(1,2,2)
shiftedControl(distance,stepMeans,WorkingProbeTrials,rewardDimension,0,1)

saveas(gcf,[pathName,'plots\meanVelVsDistanceShiftedControl_sid_',sid,'.png'])


%% Plot all individual trials and the mean

figure('Position',[100 100 1000 1000]),
subplot(2,1,1)
for optoTrial = 1:length(WorkingOptoPulse)
   lh = plot(distance,stepMeans(WorkingOptoPulse(optoTrial),:),'r');
   hold on
   lh.Color=[1,0,0,0.1];
end
plot(distance,optoMeans,'r','linewidth',3)
ylim([0 10]);
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Opto trials');
ylabel('Forward velocity (mm/s)');

subplot(2,1,2)
for probeTrial = 1:length(WorkingProbeTrials)
   lh = plot(distance,stepMeans(WorkingProbeTrials(probeTrial),:),'k');
   hold on
   lh.Color=[0,0,0,0.3];
end
plot(distance,probeMeans,'k','linewidth',3)
ylim([0 10]);
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Probe trials');
ylabel('Forward velocity (mm/s)');
xlabel('Distance (mm)');

saveas(gcf,[pathName,'plots\VelVsDistanceAllTrials_sid_',sid,'.png'])

%% Plotting the result as a heatmap

%zscore the velocity
zscored_vel = zscore(stepMeans,[],2);

figure,
subplot(2,1,1)
imagesc(zscored_vel(WorkingOptoPulse,:))

subplot(2,1,2)
imagesc(zscored_vel(WorkingProbeTrials,:))


saveas(gcf,[pathName,'plots\VelVsDistanceHeatmap_sid_',sid,'.png'])

%% Line plots for the zscored data

figure,
plot(distance,mean(zscored_vel(WorkingOptoPulse,:)),'r','linewidth',2)
hold on
plot(distance,mean(zscored_vel(WorkingProbeTrials,:)),'k','linewidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
xlabel('Distance (mm)');

saveas(gcf,[pathName,'plots\VelVsDistanceZscored_sid_',sid,'.png'])

%% Split the trials into left and right to see if the response is symmetric on both sides

%split trials
leftTrials = find(tripType==1);
rightTrials = find(tripType==2);
optoLeft = intersect(leftTrials,WorkingOptoPulse);
probeLeft = intersect(leftTrials,WorkingProbeTrials);
optoRight = intersect(rightTrials,WorkingOptoPulse);
probeRight = intersect(rightTrials,WorkingProbeTrials);

figure('Position',[100 100 1000 800]),
%left trials
subplot(2,2,1)
plot(distance,nanmean(stepMeans(optoLeft,:)),'r','linewidth',2)
hold on
plot(distance,nanmean(stepMeans(probeLeft,:)),'k','linewidth',2)
ylim([0 10]);
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Left trips');
ylabel('Forward velocity (mm/s)');

subplot(2,2,3)
plot(distance,nanmean(angularStepMeans(optoLeft,:)),'r','linewidth',2)
hold on
plot(distance,nanmean(angularStepMeans(probeLeft,:)),'k','linewidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylim([0 600]);
ylabel('Angular speed (deg/s)');
xlabel('Distance (mm)');

%right trials
subplot(2,2,2)
plot(distance,nanmean(stepMeans(optoRight,:)),'r','linewidth',2)
hold on
plot(distance,nanmean(stepMeans(probeRight,:)),'k','linewidth',2)
ylim([0 10]);
xline(distance(rewardDimension),'lineWidth',2,'color','b');
title('Right trips');

subplot(2,2,4)
plot(distance,nanmean(angularStepMeans(optoRight,:)),'r','linewidth',2)
hold on
plot(distance,nanmean(angularStepMeans(probeRight,:)),'k','linewidth',2)
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylim([0 600]);
xlabel('Distance (mm)');

saveas(gcf,[pathName,'plots\LeftVsRightTripsVel_sid_',sid,'.png'])

%% Trajectories

%load the HDF5 file data
headingData = h5read(file_name,'/heading');
posx = h5read(file_name,'/posx');
posy = h5read(file_name,'/posy');
open_loop_value = h5read(file_name,'/open loop x');

%use the open loop value to determine the turning points
test = abs(diff(open_loop_value));
turningPoints = find(test > 0.4);
timing = timestamp(2:end);

figure, plot(timing, abs(diff(open_loop_value)))
hold on
plot(timing(turningPoints),test(turningPoints),'ro')

%get the turns
for i = 1:floor(length(turningPoints)/2)
   TurnData(i).posx = posx(turningPoints(i+i-1):turningPoints(i+i),:);
   TurnData(i).posy = posy(turningPoints(i+i-1):turningPoints(i+i),:);
   TurnData(i).time = timestamp(turningPoints(i+i-1):turningPoints(i+i),:);
end

%get the trips
TripData(1).posx =  posx(1:turningPoints(1),:);
TripData(1).posy =  posy(1:turningPoints(1),:);
TripData(1).time =  timestamp(1:turningPoints(1),:);
for i = 1:floor(length(turningPoints)/2)-1
   MostTripData(i).posx = posx(turningPoints(i+i):turningPoints(i+i+1),:); 
   MostTripData(i).posy = posy(turningPoints(i+i):turningPoints(i+i+1),:); 
   MostTripData(i).time = timestamp(turningPoints(i+i):turningPoints(i+i+1),:); 
end
for i = 2:length(MostTripData)+1
   TripData(i).posx = MostTripData(i-1).posx;
   TripData(i).posy = MostTripData(i-1).posy;
   TripData(i).time = MostTripData(i-1).time;
end

scaleFactor = 4.5; %ball's radius in mm

%Calculate the 'curviness' of the trajectory, as the sum of the total
%angular speed of a given trial
for i = 1:length(TripData)
    curviness(i) = round(abs(nansum(nonsmoothed{1,i}.angularVel)));
end

%trips arranged by straightness
[sortedCurviness,I] = sort(curviness);

%2nd criterion for straightness: D/L with D being the distance from the
%start to the end of the trajectory, and L being the length of the
%trajectory.
D = distance(end);
%getting the length of the trajectories as the sum of the distances between
%consecutive points (using pythagoras).
for i = 1:length(TripData)
    TripData(i).distx = diff(TripData(i).posx.*scaleFactor);
    TripData(i).disty = diff(TripData(i).posy.*scaleFactor);
    TripData(i).dist = sqrt((TripData(i).distx).^2 + (TripData(i).disty).^2);
    L(i) = sum(TripData(i).dist);
    tortuosity(i) = round(D/L(i),2);
end

%trips arranged by tortuosity
[sortedTortuosity,I2] = sort([0-tortuosity]);

%% Aligning the opto data and plotting the trajectories

angle = rawData.trial_bdata( : , yawFly);
downsampledAngle = downsample(angle,4000/25);
angle = downsampledAngle.* 2 .* pi ./ 10;
angle = wrapTo180(rad2deg(angle));
time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(angle));
interp_panel_heading = wrapTo180(rad2deg(interp1(timestamp, unwrap(deg2rad(headingData)), time)));
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
figure('Units','normalized','Position',[0 0 1 1])
for i = 1:length(TripData)
    subplot(10,ceil(length(TripData)/10),i),
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
saveas(gcf,[pathName,'plots\TrajectoriesTortuosity_sid_',sid,'.png'])

%Plot them from 'straighter' to less straight
figure('Units','normalized','Position',[0 0 1 1])
for i = 1:length(TripData)
    subplot(10,ceil(length(TripData)/10),i),
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
saveas(gcf,[pathName,'plots\TrajectoriesCurviness_sid_',sid,'.png'])

%% Plot for each individual trial the velocities and the trajectory

for trip = 1:length(TripData)
    figure,
    subplot(2,2,1)
    if any(WorkingOptoPulse==trip)
        plot(distance,stepMeans(trip,:),'r')
    else
        plot(distance,stepMeans(trip,:),'k')
    end
    xline(distance(rewardDimension),'lineWidth',2,'color','b');
    ylabel('Forward velocity (mm/s)');
    
    subplot(2,2,3)
    if any(WorkingOptoPulse==trip)
        plot(distance,angularStepMeans(trip,:),'r')
    else
        plot(distance,angularStepMeans(trip,:),'k')
    end
    xline(distance(rewardDimension),'lineWidth',2,'color','b');
    ylabel('Angular speed (deg/s)');
    
    subplot(2,2,[2,4])
    time = linspace(1,length(rawData.trial_bdata)/samplingRate,length(TripData(trip).posx));
    timePulse = linspace(1,length(rawData.trial_bdata)/samplingRate,length(optoPulse{1,trip}));
    optoPulse{1,trip}(optoPulse{1,trip} < 3) = 0;
    optoPulse{1,trip}(optoPulse{1,trip} > 3) = 5;
    interp_opto_pulse = interp1(timePulse,optoPulse{1,trip},time);
    interp_opto_pulse2 = [zeros(abs(offset_time),1)',interp_opto_pulse];
    opto_pulse = interp_opto_pulse2(1:length(interp_opto_pulse));
    scatter((TripData(trip).posx-TripData(trip).posx(1))*scaleFactor,(TripData(trip).posy-TripData(trip).posy(1))*scaleFactor,[],opto_pulse);
    hold on
    plot([0 0],[0 0],'r+','linewidth',2)
    
    suptitle(['Trip #',num2str(trip)]);
    saveas(gcf,[pathName,'plots\ind_trajectory_trip',num2str(trip),'.png'])
end

close all;
%% Adding a fwdVel threshold for the angular speed analysis

%I will set a fwd vel threshold at 0.8 mm/s
mvt_thresh = 0.8;

for trial = 1:length(nonsmoothed)
    for timepoint = 1:length(nonsmoothed{1,trial}.angularVel)
        if nonsmoothed{1,trial}.xVel(timepoint)>mvt_thresh
            nonsmoothed{1,trial}.threshAngSpeed(timepoint) = abs(nonsmoothed{1,trial}.angularVel(timepoint));
        else
            nonsmoothed{1,trial}.threshAngSpeed(timepoint) = nan;
        end
    end
end


% get the mean for each distance step
%bin the forward velocity data into voltages using the ypanel data
for i = 1:length(nonsmoothed)
    for j = 1:48
        angularStepMeansThresh(i,j) = nanmean(nonsmoothed{1,i}.threshAngSpeed(step{i,j}));
    end
end


% %Linearly interpolate the data to fill in the gaps
% for i = 1:size(stepMeans,1)
%     angularStepMeansTresh(i,:) = fillmissing(angularStepMeansTresh(i,:),'linear');
% end

%flip the data for right trips
for trip = 1:length(tripData)
    if tripType(trip) == 2
        angularStepMeansThresh(trip,:) = flip(angularStepMeansThresh(trip,:));
    end
end

%Plot
optoAngTreshMeans = nanmean(angularStepMeansThresh(WorkingOptoPulse,:));
probeAngThreshMeans = nanmean(angularStepMeansThresh(WorkingProbeTrials,:));
optoAngThreshError = (nanstd(angularStepMeansThresh(WorkingOptoPulse,:)))./sqrt(length(WorkingOptoPulse));
probeAngThreshError = (nanstd(angularStepMeansThresh(WorkingProbeTrials,:)))./sqrt(length(WorkingProbeTrials));

figure('Position', [100 100 1400 800]),
ax1 = subplot(1,2,1);
[opT] = boundedline(distance,optoAngMeans,optoAngError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeAngMeans,probeAngError,'-ko','alpha');
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Angular speed (deg/s)');
title('Ang speed')

ax1 = subplot(1,2,2);
[opT] = boundedline(distance,optoAngTreshMeans,optoAngThreshError,'-ro','alpha');
hold on
[prT] = boundedline(distance,probeAngThreshMeans,probeAngThreshError,'-ko','alpha');
xline(distance(rewardDimension),'lineWidth',2,'color','b');
ylabel('Angular speed (deg/s)');
title('Ang speed thresholding fwd vel')

saveas(gcf,[pathName,'plots\angSpeedThresh_sid_',sid,'.png'])