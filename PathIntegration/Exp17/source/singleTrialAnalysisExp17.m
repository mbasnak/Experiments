
function singleTrialAnalysisExp17(dirName,sid)

cd(dirName)
fileNames = dir('*.mat');


%load data from all trials of a specific session.
for i = 1:length(fileNames)
    if regexp(fileNames(i).name,strcat('sid_',num2str(sid))) > 1        
        rawData{i} = load(strcat(fileNames(i).folder,'\',fileNames(i).name));
    end
end

%load the run object, which contains additional information about the
%session
run_obj_dir = strcat(dirName,'\runobj');
cd (run_obj_dir)
run_obj_name = dir(strcat('*',num2str(sid),'_runobj.mat'));
run_obj = load(strcat(run_obj_dir,'\',run_obj_name(1).name));

%remove empty cells
rawData = rawData(~cellfun('isempty',rawData));

%Define Ni-Daq channels ID
headingFly = 1;
yFly = 2;
xFly = 3;
xFlyGain = 7;
xPanels = 4;
yPanels = 5;
PanelStatus = 6;
OptoTrigger = 8;

samplingRate = 4000;

%determine the valence (sweet sensor rewarding, heat sensor aversive)
if contains(dirName, 'Gr28bd')
    valence = 'punishment';
else
    valence = 'reward';
end

%% Getting different trial types

trialNum = [1:length(rawData)];

probeTrials = sort(rawData{1,1}.probeTrials(rawData{1,1}.probeTrials < length(rawData)));
emptyTrials = rawData{1,1}.emptyTrials(rawData{1,1}.emptyTrials < length(rawData));
optoTrials = sort(setdiff(trialNum,[probeTrials,emptyTrials]));

%% Subset acquisition of x and y pos, as well as FicTrac data

%For each trial
for i = 1:length(rawData)
    
    %get the raw data
    Data = rawData{1,i}.trial_bdata;
    
    data.xPanelVolts =  Data(:,xPanels); 
    VOLTAGE_RANGE_x = 10;
    maxValX =  96 ;
      
    %In the non-empty trials, the python code lags behing the matlab code by a few seconds, so we
    %need to add a fix to cut the trials and show the traces from the
    %moment they started receiving a python signal
    
    %Find the data point when the acquisition starts and cut it so that it
    %starts showing us the info from there
        if sum(emptyTrials == i) == 0 %if this is not an empty trial
            StartPoint = find(diff(data.xPanelVolts)<-1 | diff(data.xPanelVolts)>1); %find when the panels are on
            if ~isempty(StartPoint)
                startPoints{i} = StartPoint;
                startPoint(i) = startPoints{i}(1); %take the point when they turn on
                Data = Data(startPoint(i):end,:); %cut the data before the panels were on
                trialsToRemove(i) = NaN; %this trial should be kept
            else %if there doesn't seem to be a starting signal
                trialsToRemove(i) = i; %label the trial to be removed
            end
        else
            trialsToRemove(i) = NaN; %this trial shouldn't be removed

        end   
        
        %FicTrac data
        data.ficTracAngularPosition = Data( : , headingFly); 
        data.ficTracIntx = Data( : , xFly); 
        data.ficTracInty = Data( : , yFly); 
      
        data.yPanelVolts =  Data(:, yPanels);
        VOLTAGE_RANGE_y = 10;
        maxValY = 92; %both virtual hallways have 92 y dimensions
        data.xPanelVolts = Data( : , xPanels);
    
        sizeBall = 9;
        %get the smoothed position and velocity with 2 different methods
        [smoothed{i}] = singleTrialVelocityAnalysis9mm(data,samplingRate);
        [forwardVel{i}, accumulatedx{i}] = ficTracSignalDecoding(data.ficTracIntx, samplingRate , 50, 10);
                
        ypanels{i} = resample(Data(:,yPanels),25,samplingRate); %downsampled panels data
        optoPulse{i} = resample(Data(:,OptoTrigger),25,samplingRate); %downsampled opto trigger
         
        trialDuration(i) = length(data.xPanelVolts)/samplingRate;

end

%getting rid of the NaNs.
trialsToRemove = trialsToRemove(~isnan(trialsToRemove));

%remove the trials defined in 'trialsToRemove' from the different trial
%categories
for i = 1:length(trialsToRemove)
    optoPulse{1,trialsToRemove(i)} = [];
    O = find(optoTrials == trialsToRemove(i));
    P = find(probeTrials == trialsToRemove(i));
    if (~isempty(O))
        optoTrials(O) = [];
    elseif (~isempty(P))
        probeTrials(P) = [];
    end   
end

%% Trial dur vs trial number

PT = ones(1,length(probeTrials));
ET = ones(1,length(emptyTrials));
OT = ones(1,length(optoTrials));
group = [PT, 2*ET, 3*OT];

durPT = trialDuration(probeTrials);
durET = trialDuration(emptyTrials);
durOT = trialDuration(optoTrials);


figure('Position', [300, 300, 1200, 800]),
subplot(1,2,1)
plot(trialNum(probeTrials),trialDuration(probeTrials),'ko','MarkerFaceColor', 'k','DisplayName','Probe trials')
hold on
plot(trialNum(emptyTrials),trialDuration(emptyTrials),'ko','DisplayName','Empty trials')
plot(trialNum(optoTrials),trialDuration(optoTrials),'ro','MarkerFaceColor', 'r','DisplayName','Opto trials')
line([0 trialNum(end)], [run_obj.run_obj.trial_t run_obj.run_obj.trial_t], 'Color', 'k', 'LineStyle', '--','DisplayName','Max trial duration')
ylim([0 run_obj.run_obj.trial_t+2]);
xlabel('Trial number'); ylabel('Trial duration (s)'); 
legend('Probe trials', 'Empty trials', 'Opto trials', 'Max trial duration');

subplot(1,2,2)
boxplot([durPT';durET';durOT'],group)
ax = gca;
ax.XTick = [1 2 3];
ax.XTickLabel={'Probe Trials','Empty Trials','Opto Trials'};
title('Distribution of trial duration');ylabel('Trial duration (s)');

saveas(gcf,[dirName(1:end-4),'\plots\TrialDurVsTrialNum.png'])
close;

%% remove trials with outlier velocity values 

%Remove trials that have particularly
%high or low velocity compared to the rest

% trialVelMeans=[]; %get mean forward velocity for each trial 
% for i=1:size(smoothed,2)
%     trialVelMeans=[trialVelMeans mean(smoothed{1,i}.xVel)];
% end
% %exclude trials with velocities in bottom 20% and top 5%
% percentile20 = prctile(trialVelMeans,20);
% percentile95 = prctile(trialVelMeans,95);
% excludetrials=[];
% for i=1:size(smoothed,2)
%     if (mean(smoothed{1,i}.xVel)<percentile20 | mean(smoothed{1,i}.xVel)>percentile95 | mean(smoothed{1,i}.xVel)>10 | mean(smoothed{1,i}.xVel)<-10)
%         excludetrials=[excludetrials i];
%     end
% end
% figure
% scatter(1:size(smoothed,2),trialVelMeans);
% hold on; scatter(excludetrials,trialVelMeans(excludetrials),'r');
% xlabel('Trial number');ylabel('mean forward velocity');
% yline(percentile20);yline(percentile95);
% title('Exclude velocity outliers');legend('trials to keep', 'trials to exclude');
% 
% 
% %remove outlier trials 
% optoTrials=setdiff(optoTrials,excludetrials);
% probeTrials=setdiff(probeTrials,excludetrials);
% emptyTrials=setdiff(emptyTrials,excludetrials);

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
emptyTrials=setdiff(emptyTrials,excludetrials);

%% Look at something analogous to spike triggered average

%looking at the opto pulse, there are many trials that should be opto
%trials that are not showing any pulse, (either those in which the fly
%never made it to the end of the hallway before the trial reset or a couple
%in which the pulse seems to have failed). I'm getting rid of those for the
%following analyses.

%use the opto trigger signal to figure out if the opto trials were actually
%opto and define the start and end of the opto pulse
for i = 1:length(optoTrials)
    if (~isempty(find(abs(diff(rawData{1,optoTrials(i)}.trial_bdata(:,8)))>1)))
        WorkingOptoPulse(i) = optoTrials(i);
        pulseFrames = find(abs(diff(rawData{1,optoTrials(i)}.trial_bdata(:,8)))>1); %we use the original data
        startPulse(i) = floor((pulseFrames(1)-startPoint(optoTrials(i)))*25/samplingRate); %we remove the chunk before the panels started, and downsample
        endPulse(i) = floor((pulseFrames(end)-startPoint(optoTrials(i)))*25/samplingRate);
    else
        WorkingOptoPulse(i) = NaN;
        startPulse(i) = NaN;
        endPulse(i) = NaN;
    end
end

WorkingOptoPulse = WorkingOptoPulse(~isnan(WorkingOptoPulse));
startPulse = startPulse(~isnan(startPulse));
endPulse = endPulse(~isnan(endPulse));

%get the velocity around the opto pulse
for i = 1:length(WorkingOptoPulse)
    if ((length(smoothed{1,WorkingOptoPulse(i)}.xVel) > startPulse(i)+25 )& startPulse(i) > 14)%if the trial goes on for an extra second after the max expansion
        velAroundPulse{i} = smoothed{1,WorkingOptoPulse(i)}.xVel(startPulse(i)-13:startPulse(i)+25); %collect the data until that timepoint
    elseif ((length(smoothed{1,WorkingOptoPulse(i)}.xVel) > startPulse(i)+25 )& startPulse(i) <= 14)      
        velAroundPulse{i} = smoothed{1,WorkingOptoPulse(i)}.xVel(1:startPulse(i)+25);
        velAroundPulse{i} = [repelem(nan,1,14-startPulse(i)),velAroundPulse{i}];
    else
        velAroundPulse{i} = smoothed{1,WorkingOptoPulse(i)}.xVel(startPulse(i)-13:end); %otherwise collect the data until the end of the trial
        velAroundPulse{i}(end+1:39) = nan;
    end        
end

trialLengths = cellfun(@length,velAroundPulse);  
velAroundPulse = cell2mat(velAroundPulse);
velAroundPulse = reshape(velAroundPulse,39,length(trialLengths));

%Divide time in 3 sections: before, during and after the pulse, and take
%the mean vel por each trial for each section
meanVelBeforePulse = nanmean(velAroundPulse(1:13,:));
meanVelDuringPulse = nanmean(velAroundPulse(14:26,:));
meanVelAfterPulse = nanmean(velAroundPulse(27:39,:));
meanVelAroundPulse = [meanVelBeforePulse',meanVelDuringPulse',meanVelAfterPulse'];

figure('Position', [100 100 1600 1000])
subplot(1,2,1)
time = linspace(-0.5,1,39);
boundedline(time,nanmean(velAroundPulse,2),nanstd(velAroundPulse,[],2)/sqrt(length(WorkingOptoPulse)))
hold on
line([0 0],[0 max(max(meanVelAroundPulse'))+0.3],'Color','r','LineWidth',2)
line([0.5 0.5],[0 max(max(meanVelAroundPulse'))+0.3],'Color','r','LineWidth',2)
ylim([0 max(max(meanVelAroundPulse'))+0.3]);title('Forward velocity around the pulse');
ylabel('Forward velocity (mm/s)');
xlabel('Time around the pulse (s)');
xlim([-0.5 1]);
text(-0.25,max(max(nanmean(velAroundPulse,2)))+0.2,'1','FontSize',14)
text(0.25,max(max(nanmean(velAroundPulse,2)))+0.2,'2','FontSize',14)
text(0.75,max(max(nanmean(velAroundPulse,2)))+0.2,'3','FontSize',14)

subplot(1,2,2)
plot(meanVelAroundPulse','color',[0.5 0.5 0.5])
hold on
err = iqr(meanVelAroundPulse);
errorbar([1,2,3],nanmean(meanVelAroundPulse),err,'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
xticks([1 2 3])
xticklabels({'1','2','3'})
xlim([0 4]);
ylim([0 max(max(meanVelAroundPulse'))+0.3]);

saveas(gcf,[dirName(1:end-4),'plots\VelAroundPulseOpto.png'])
close;

%% Look at that for the probe trials

%For the probe trials, I use the voltage of the y dimension of the panels to identify when
%the expanding object reached the size where the reward would have been
%given if it was an opto trial
  
for i = 1:length(probeTrials)
    expansionPoints{i} = find(ypanels{1,probeTrials(i)}>0.75 & ypanels{1,probeTrials(i)}<7.5);
    realExpansionPoints{i} = expansionPoints{i}(expansionPoints{i}>25);
    if (~isempty(realExpansionPoints{i}))
        WorkingProbeTrials(i) = probeTrials(i);
        startExpansion(i) = realExpansionPoints{i}(1);
    else
        WorkingProbeTrials(i) = NaN;
        startExpansion(i) = NaN; 
    end
end
WorkingProbeTrials = WorkingProbeTrials(~isnan(WorkingProbeTrials));
startExpansion = startExpansion(~isnan(startExpansion));

for i = 1:length(WorkingProbeTrials)
    if length(smoothed{1,WorkingProbeTrials(i)}.xVel) > startExpansion(i)+25 %if the trial goes on for an extra second after the max expansion
        velAroundExpansionPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startExpansion(i)-13:startExpansion(i)+25); %collect the data until that timepoint
    else
        velAroundExpansionPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startExpansion(i)-13:end); %otherwise collect the data until the end of the trial
    end        
end
    
%Find the trials for which the acquisition finished earlier and add
%NaNs at the end
trialLengths = cellfun(@length,velAroundExpansionPoint);
shortTrials = find(trialLengths<39);
if shortTrials
    for i = 1:length(shortTrials)
       velAroundExpansionPoint{shortTrials(i)}(end+1:39) = nan;
    end
end    
velAroundExpansionPoint = cell2mat(velAroundExpansionPoint);
velAroundExpansionPoint = reshape(velAroundExpansionPoint,39,length(trialLengths));
    
figure('Position', [100 100 1600 1000])
subplot(1,2,1)
time = linspace(-0.5,1,39);
boundedline(time,nanmean(velAroundExpansionPoint'),nanstd(velAroundExpansionPoint')/sqrt(length(WorkingProbeTrials)))
xlim([-0.5 1]);
xlabel('Time around expected reward (s)');
hold on
line([0 0],[0 5],'Color','k','LineWidth',2)
ylim([0 5]);
title(['Forward velocity around the expected ', valence]);
ylabel('Forward velocity (mm/s)');
    
subplot(1,2,2)
for i = 1:length(WorkingProbeTrials)
            
    if startExpansion(i)>20          
        velPrePoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startExpansion(i)-20:startExpansion(i)); %collect the data until that timepoint
    else
        velPrePoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(1:startExpansion(i));
        velPrePoint{i} = [repelem(nan,1,21-startExpansion(i)),velPrePoint{i}];
    end
        
    if length(smoothed{1,WorkingProbeTrials(i)}.xVel) > startExpansion(i)+20 %if the trial goes on for an extra second after the max expansion
        velPostPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startExpansion(i):startExpansion(i)+20); %collect the data until that timepoint
    else
        velPostPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startExpansion(i):end); %otherwise collect the data until the end of the trial
    end        
end

trialLength = cellfun(@length,velPostPoint);
shortTrial = find(trialLength<21);
if shortTrial
    for i = 1:length(shortTrial)
        velPostPoint{shortTrial(i)}(end+1:21) = nan;
    end
end 
velPrePoint = cell2mat(velPrePoint);
velPostPoint = cell2mat(velPostPoint);
velPostPoint = reshape(velPostPoint,21,length(trialLength));
velPrePoint = reshape(velPrePoint,21,length(trialLength));
    
profiles = [nanmedian(velPrePoint);nanmedian(velPostPoint)];   
plot(profiles,'color',[0.5 0.5 0.5])
hold on
err = iqr(profiles');
errorbar([1,2],median(profiles'),err,'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
xlim([0 3]);
ylabel('Forward velocity (mm/s)');
set(gca,'xTickLabel',[]);

saveas(gcf,[dirName(1:end-4),'plots\VelAroundExpected', valence, '.png'])
close;

%% add a random point to compare
%    
% for i = 1:length(WorkingProbeTrials)
%     RandomPoints{i} = find(ypanels{1,WorkingProbeTrials(i)}>0.45 & ypanels{1,WorkingProbeTrials(i)}<7.5);
%     realRandomPoints{i} = RandomPoints{i}(RandomPoints{i}>15);
%     if (~isempty(realRandomPoints{i}))
%         startRandom(i) = realRandomPoints{i}(1);
%     else
%         startRandom(i) = NaN; 
%     end
% end
% 
% startRandom = startRandom(~isnan(startRandom));
% 
% for i = 1:length(WorkingProbeTrials)
%     if length(smoothed{1,WorkingProbeTrials(i)}.xVel) > startRandom(i)+25 %if the trial goes on for an extra second after the max Random
%         velAroundRandomPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startRandom(i)-10:startRandom(i)+25); %collect the data until that timepoint
%     else
%         velAroundRandomPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startRandom(i)-10:end); %otherwise collect the data until the end of the trial
%     end        
% end
%     
%     %Find the trials for which the acquisition finished earlier and add
%     %NaNs at the end
% trialLengths = cellfun(@length,velAroundRandomPoint);
% shortTrials = find(trialLengths<36);
% if shortTrials
%     for i = 1:length(shortTrials)
%         velAroundRandomPoint{shortTrials(i)}(end+1:36) = nan;
%     end
% end    
% velAroundRandomPoint = cell2mat(velAroundRandomPoint);
% velAroundRandomPoint = reshape(velAroundRandomPoint,36,length(trialLengths));
%     
% figure('Position', [100 100 1200 1000])
% subplot(1,2,1)
% time = linspace(-0.5,1,36);
% boundedline(time,nanmedian(velAroundRandomPoint'),nanstd(velAroundRandomPoint')/sqrt(length(WorkingProbeTrials)))
% xlim([-0.5 1]); ylim([0 5]);
% xlabel('Time around random point (s)');
% hold on
% line([0 0],[0 5],'Color','k','LineWidth',2)
% ylim([0 5]);
% title('Forward velocity around random point');
% ylabel('Forward velocity (mm/s)');
%     
% subplot(1,2,2)
% for i = 1:length(WorkingProbeTrials)
%     velPreRandomPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startRandom(i)-15:startRandom(i)); %collect the data until that timepoint
%         
%     if length(smoothed{1,WorkingProbeTrials(i)}.xVel) > startRandom(i)+20 %if the trial goes on for an extra second after the max expansion
%         velPostRandomPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startRandom(i):startRandom(i)+20); %collect the data until that timepoint
%     else
%         velPostRandomPoint{i} = smoothed{1,WorkingProbeTrials(i)}.xVel(startRandom(i):end); %otherwise collect the data until the end of the trial
%     end        
% end
% 
% trialLength = cellfun(@length,velPostRandomPoint);
%     shortTrial = find(trialLength<21);
%     if shortTrial
%         for i = 1:length(shortTrial)
%             velPostRandomPoint{shortTrial(i)}(end+1:21) = nan;
%         end
%     end 
% velPreRandomPoint = cell2mat(velPreRandomPoint);
% velPostRandomPoint = cell2mat(velPostRandomPoint);
% velPostRandomPoint = reshape(velPostRandomPoint,21,length(trialLength));
% velPreRandomPoint = reshape(velPreRandomPoint,16,length(trialLength));
%     
% profilesRandom = [median(velPreRandomPoint);median(velPostRandomPoint)];
% 
% plot(profilesRandom,'color',[0.5 0.5 0.5])
% hold on
% errRandom = iqr(profiles');
% errorbar([1,2],median(profilesRandom'),errRandom,'-o','MarkerSize',5,'color','k','CapSize',0,'LineWidth',2,'MarkerFaceColor','k')
% xlim([-0 3]); ylim([0 5]);
% title('Velocity before and after random point');
% ylabel('Forward velocity (mm/s)');
% set(gca,'XTickLabel',[]);
% 
% saveas(gcf,[dirName(1:end-4),'plots\VelAroundRandomPoint', valence, '.png'])
% close;

%% Time estimate with fwd velocity
%We are mainly focusing on whether flies can estimate distances, but an
%alternative is that they might estimate the timing of a reward. The
%following analysis is to detect consistency in their responses with
%respect to time.

figure('Position', [100 100 1600 1000]),
subplot(1,2,1)
for i = 1:length(smoothed)
    if any(WorkingProbeTrials == i)
        Time{i} = linspace(1,length(smoothed{1,i}.xVel)/25,length(smoothed{1,i}.xVel));
        plot(Time{i},smoothed{1,i}.xVel,'k') %plot the fwd vel in time for a given trial
       [localmin{i},Prom{i}] = islocalmin(smoothed{1,i}.xVel); %find local minima of the fwd vel 
        if sum(localmin{i}) > 0
            [M{i},I{i}] = max(Prom{1,i}); %save the index of the biggest valley
            hold on
            plot(Time{i}(I{1,i}),smoothed{1,i}.xVel(I{i}),'ro') %add the biggest local min to the plot
        else
            I{i} = [];
        end
    end

    hold on
end
xlabel('Time (s)'); ylabel('Forward velocity (mm/s)');
title('Forward velocity local minima');

subplot(1,2,2)
minVel = cell2mat(I); %save the indexes of when in time each fly displayed a fwd vel min
boxplot(minVel/25,'orientation','horizontal');
title('Time estimate');
xlabel('Time (s)');
set(gca,'ytick',[]);

saveas(gcf,[dirName(1:end-4),'plots\TimeEstimateFwdVel.png'])
close;

%% Time estimate with ang speed
%They might be displaying behavioral changes in their angular speed instead
%of their fwd vel do this is to see if there are consistent changes in
%their ang speed with respect to time

figure('Position', [100 100 1600 1000]),
subplot(1,2,1)
for i = 1:length(smoothed)
    if any(WorkingProbeTrials == i)
        AngTime{i} = linspace(1,length(smoothed{1,i}.angularVel)/25,length(smoothed{1,i}.angularVel));
        plot(AngTime{i},abs(smoothed{1,i}.angularVel),'k')
       [localmax{i},Prom{i}] = islocalmax(abs(smoothed{1,i}.angularVel)); %find local max 
        if sum(localmax{i}) > 0
            [M{i},I{i}] = max(Prom{1,i}); %save the index of the biggest peak in ang speed
            hold on
            plot(AngTime{i}(I{1,i}),abs(smoothed{1,i}.angularVel(I{i})),'ro')
        else
            I{i} = [];
        end
    end

    hold on
end
xlabel('Time (s)'); ylabel('Angular speed (deg/s)');
title('Angular speed local maxima');

subplot(1,2,2)
maxAngVel = cell2mat(I);
boxplot(maxAngVel/25,'orientation','horizontal');
title('Time estimate');
xlabel('Time (s)');
set(gca,'ytick',[]);

saveas(gcf,[dirName(1:end-4),'plots\TimeEstimateAngSpeed.png'])
close;

%% Plotting velocity as a function of distance travelled
% we get the distance that the animals travelled looking at the y panels,
% since they move with the animal's forward walking

figure('Position', [100 100 1600 1000]),
subplot(1,2,1)
for i = 1:length(smoothed)
    if any(WorkingProbeTrials == i)
        plot(ypanels{1,i}, smoothed{1,i}.xVel,'r.')
    elseif any(emptyTrials == i)
        plot(ypanels{1,i}, smoothed{1,i}.xVel,'b.')
    elseif any (WorkingOptoPulse == i)
        plot(ypanels{1,i}, smoothed{1,i}.xVel,'k.')
    end
    hold on
end
xlabel('y panel voltage'); ylabel('Forward velocity (mm/s)');
xlim([0 1.50]);

subplot(1,2,2)
for i = 1:length(smoothed)
    if any(WorkingProbeTrials == i)
        plot(ypanels{1,i}, smoothed{1,i}.angularVel,'r.')
    elseif any(emptyTrials == i)
        plot(ypanels{1,i}, smoothed{1,i}.angularVel,'b.')
    elseif any (WorkingOptoPulse == i)
        plot(ypanels{1,i}, smoothed{1,i}.angularVel,'k.')
    end
    hold on
end
xlabel('y panel voltage'); ylabel('Angular velocity (deg/s)');
xlim([0 1.50]);

saveas(gcf,[dirName(1:end-4),'plots\VelVsDistSingleTrials.png'])
close;

%% Binning it

%The plots from the previous block are really hard to interpret because of
%the dot cloud. We will bin the data in different distance steps for more clarity

for i = 1:length(smoothed)
    step{i,1} = find(ypanels{1,i}>0.07 & ypanels{1,i}<0.15);
    step{i,2} = find(ypanels{1,i}>0.15 & ypanels{1,i}<0.26);
    step{i,3} = find(ypanels{1,i}>0.26 & ypanels{1,i}<0.36);
    step{i,4} = find(ypanels{1,i}>0.36 & ypanels{1,i}<0.48);
    step{i,5} = find(ypanels{1,i}>0.48 & ypanels{1,i}<0.58);
    step{i,6} = find(ypanels{1,i}>0.58 & ypanels{1,i}<0.68);
    step{i,7} = find(ypanels{1,i}>0.68 & ypanels{1,i}<0.80);
    step{i,8} = find(ypanels{1,i}>0.80 & ypanels{1,i}<0.90);
    step{i,9} = find(ypanels{1,i}>0.90 & ypanels{1,i}<1.01);
    step{i,10} = find(ypanels{1,i}>1.01 & ypanels{1,i}<1.13);
    step{i,11} = find(ypanels{1,i}>1.13 & ypanels{1,i}<1.23);
    step{i,12} = find(ypanels{1,i}>1.23 & ypanels{1,i}<1.35);
    step{i,13} = find(ypanels{1,i}>1.35 & ypanels{1,i}<1.41);
end

%get the mean for each step
for i = 1:length(smoothed)
    for j = 1:13
        stepMeans(i,j) = mean(smoothed{1,i}.xVel(step{i,j}));
        angularStepMeans(i,j) = mean(abs(smoothed{1,i}.angularVel(step{i,j})));
    end
end

%I linearly interpolate to fill in the gaps
for i = 1:length(stepMeans)
    stepMeans(i,:) = fillmissing(stepMeans(i,:),'linear');
    angularStepMeans(i,:) = fillmissing(angularStepMeans(i,:),'linear');
end

optoMeans = nanmedian(stepMeans(WorkingOptoPulse,:));
probeMeans = nanmedian(stepMeans(WorkingProbeTrials,:));
optoError = (nanstd(stepMeans(WorkingOptoPulse,:)))./sqrt(length(WorkingOptoPulse));
probeError = (nanstd(stepMeans(WorkingProbeTrials,:)))./sqrt(length(WorkingProbeTrials));

optoAngMeans = nanmedian(angularStepMeans(WorkingOptoPulse,:));
probeAngMeans = nanmedian(angularStepMeans(WorkingProbeTrials,:));
optoAngError = (nanstd(angularStepMeans(WorkingOptoPulse,:)))./sqrt(length(WorkingOptoPulse));
probeAngError = (nanstd(angularStepMeans(WorkingProbeTrials,:)))./sqrt(length(WorkingProbeTrials));

ydimension = [1:13];

figure('Position', [100 100 1600 1000]),
subplot(1,2,1)
[opT] = boundedline(ydimension,optoMeans,optoError,'-ro','alpha');
hold on
[prT] = boundedline(ydimension,probeMeans,probeError,'-ko','alpha');
rT = line([7 7],[0 max(max([probeMeans;optoMeans]))+1],'Color' ,'b','LineStyle', '-.');
ylabel('Forward velocity (mm/s)');
xlabel('Panels y dimension');
TrialTypes = [opT;prT;rT];
legend(TrialTypes,'opto trials', 'probe trials',[valence, ' location']);
ylim([0 max(max([probeMeans;optoMeans]))+1]);
xlim([1 13]);

subplot(1,2,2)
[opT] = boundedline(ydimension,optoAngMeans,optoAngError,'-ro','alpha');
hold on
[prT] = boundedline(ydimension,probeAngMeans,probeAngError,'-ko','alpha');
rT = line([7 7],[min(min([probeAngMeans;optoAngMeans]))-5 max(max([probeAngMeans;optoAngMeans]))+5],'Color' ,'b','LineStyle', '-.')
ylabel('Angular speed (deg/s)');
xlabel('Panels y dimension');
TrialTypes = [opT;prT;rT];
legend(TrialTypes,'opto trials', 'probe trials',[valence, ' location']);
ylim([min(min([probeAngMeans;optoAngMeans]))-5 max(max([probeAngMeans;optoAngMeans]))+5]);
xlim([1 13]);

saveas(gcf,[dirName(1:end-4),'plots\MeanVelvsPanels.png'])
close;

%save useful data for group analysis.
save([dirName(1:end-4),'dataFromAnalysis\data'],'smoothed','WorkingOptoPulse','WorkingProbeTrials','step','ypanels','optoPulse','velAroundExpansionPoint','velAroundPulse','velPostPoint','velPrePoint');

%% Getting that for early trials vs late trials to see if there is some evolution

earlyOptoTrials = WorkingOptoPulse(1:floor(length(WorkingOptoPulse)/2));
earlyProbeTrials = WorkingProbeTrials(1:floor(length(WorkingProbeTrials)/2));
lateOptoTrials = WorkingOptoPulse(floor(length(WorkingOptoPulse)/2)+1:end);
lateProbeTrials = WorkingProbeTrials(floor(length(WorkingProbeTrials)/2)+1:end);

earlyOptoMeans = nanmedian(stepMeans(earlyOptoTrials,:));
earlyProbeMeans = nanmedian(stepMeans(earlyProbeTrials,:));
earlyOptoError = (nanstd(stepMeans(earlyOptoTrials,:)))./sqrt(length(earlyOptoTrials));
earlyProbeError = (nanstd(stepMeans(earlyProbeTrials,:)))./sqrt(length(earlyProbeTrials));

lateOptoMeans = nanmedian(stepMeans(lateOptoTrials,:));
lateProbeMeans = nanmedian(stepMeans(lateProbeTrials,:));
lateOptoError = (nanstd(stepMeans(lateOptoTrials,:)))./sqrt(length(lateOptoTrials));
lateProbeError = (nanstd(stepMeans(lateProbeTrials,:)))./sqrt(length(lateProbeTrials));

figure('Position', [100 300 1000 600]),
subplot(1,2,1)
[opT] = boundedline(ydimension,earlyOptoMeans,earlyOptoError,'-ro','alpha');
hold on
[prT] = boundedline(ydimension,earlyProbeMeans,earlyProbeError,'-ko','alpha');
line([7 7],[0 max([earlyOptoMeans,earlyProbeMeans,lateOptoMeans,lateProbeMeans])+1],'Color' ,'b','LineStyle', '-.')
ylabel('Forward velocity (mm/s)');
xlabel('panels y dimension');
TrialTypes = [opT;prT];
legend(TrialTypes,'opto trials', 'probe trials');
title('Early trials')
ylim([0 max([earlyOptoMeans,earlyProbeMeans,lateOptoMeans,lateProbeMeans])+1]);
xlim([1 13]);

subplot(1,2,2)
[opT] = boundedline(ydimension,lateOptoMeans,lateOptoError,'-ro','alpha');
hold on
[prT] = boundedline(ydimension,lateProbeMeans,lateProbeError,'-ko','alpha');
line([7 7],[0 max([earlyOptoMeans,earlyProbeMeans,lateOptoMeans,lateProbeMeans])+1],'Color' ,'b','LineStyle', '-.')
xlabel('panels y dimension');
TrialTypes = [opT;prT];
legend(TrialTypes,'opto trials', 'probe trials');
title('Late trials');
ylim([0 max([earlyOptoMeans,earlyProbeMeans,lateOptoMeans,lateProbeMeans])+1]);
xlim([1 13]);

saveas(gcf,[dirName(1:end-4),'plots\MeanVelvsPanelsEarlyAndLateTrials.png'])
close;

%% Profile plots coloring looking at early and late probe trials

figure
subplot(1,2,1)
[prT] = boundedline(ydimension,earlyProbeMeans,earlyProbeError,'-ko','alpha');
hold on
plot(ydimension,stepMeans(earlyProbeTrials,:))%,'color',[0.5 0.5 0.5])
xline(7,'-.r','lineWidth',2)
ylabel('Forward velocity (mm/s)');
xlabel('panels y dimension');
title('Early trials')
ylim([0 max(max(stepMeans(probeTrials,:)))+1]);

subplot(1,2,2)
[prT] = boundedline(ydimension,lateProbeMeans,lateProbeError,'-ko','alpha');
hold on
plot(ydimension,stepMeans(lateProbeTrials,:))%,'color',[0.5 0.5 0.5])
xline(7,'-.r','lineWidth',2)
xlabel('panels y dimension');
title('Late trials');
ylim([0 max(max(stepMeans(probeTrials,:)))+1]);

saveas(gcf,[dirName(1:end-4),'plots\ProfilePlotsEarlyAndLateTrials.png'])
close;

%% Adding angular velocity

earlyAngProbeMeans = nanmedian(angularStepMeans(earlyProbeTrials,:));
lateAngProbeMeans = nanmedian(angularStepMeans(lateProbeTrials,:));
earlyAngProbeError = (nanstd(angularStepMeans(earlyProbeTrials,:)))./sqrt(length(earlyProbeTrials));
lateAngProbeError = (nanstd(angularStepMeans(lateProbeTrials,:)))./sqrt(length(lateProbeTrials));

figure
subplot(1,2,1)
[prT] = boundedline(ydimension,earlyAngProbeMeans,earlyAngProbeError,'-ko','alpha');
hold on
plot(ydimension,angularStepMeans(earlyProbeTrials,:))%,'color',[0.5 0.5 0.5])
xline(7,'-.r','lineWidth',2)
ylabel('Angular speed (deg/s)');
xlabel('panels y dimension');
title('Early trials')
ylim([0 max(max(angularStepMeans(probeTrials,:)))+1]);

subplot(1,2,2)
[prT] = boundedline(ydimension,lateAngProbeMeans,lateAngProbeError,'-ko','alpha');
hold on
plot(ydimension,angularStepMeans(lateProbeTrials,:))%,'color',[0.5 0.5 0.5])
xline(7,'-.r','lineWidth',2)
xlabel('panels y dimension');
title('Late trials');
ylim([0 max(max(angularStepMeans(probeTrials,:)))+1]);

saveas(gcf,[dirName(1:end-4),'plots\ProfilePlotsAVEarlyAndLateTrials.png'])
close;

%% Forward velocity and angular speed correlation

% figure,
% subplot(1,2,1)
% map = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 0 0; 0 0 0];
% 
% for i = 1:length(optoTrials)
% scatter(stepMeans(optoTrials(i),:),abs(angularStepMeans(optoTrials(i),:)),[],ydimension)
% xlabel('Forward velocity (deg/s)');ylabel('Angular speed (deg/s)');
% colormap(map)
% c = colorbar;
% c.Label.String = 'ydimension';
% hold on
% end
% xlim([0 5]);
% 
% subplot(1,2,2)
% for i = 1:length(probeTrials)
% scatter(stepMeans(probeTrials(i),:),abs(angularStepMeans(probeTrials(i),:)),[],ydimension)
% xlabel('Forward velocity (deg/s)');ylabel('Angular speed (deg/s)');
% colormap(map)
% c = colorbar;
% c.Label.String = 'ydimension';
% hold on
% end
% R2 = corr2(stepMeans(probeTrials,:),angularStepMeans(probeTrials,:));
% xlim([0 5]);
% 
% saveas(gcf,[dirName,'\VelocityCorrelation.png'])

%% Plotting the behavior in probe trials preceded by probe vs preceded by opto trials


for i = 1:length(WorkingProbeTrials)
    if any(WorkingOptoPulse == WorkingProbeTrials(i)-1)
        probeAfterOptoTrials(i) = WorkingProbeTrials(i);
        probeAfterProbeTrials(i) = nan;
    else
        probeAfterProbeTrials(i) = WorkingProbeTrials(i);
        probeAfterOptoTrials(i) = nan;
    end
end

probeAfterOptoTrials = probeAfterOptoTrials(~isnan(probeAfterOptoTrials));
probeAfterProbeTrials = probeAfterProbeTrials(~isnan(probeAfterProbeTrials));


probeAfterOptoMeans = nanmedian(stepMeans(probeAfterOptoTrials,:));
probeAfterProbeMeans = nanmedian(stepMeans(probeAfterProbeTrials,:));
probeAfterOptoError = (nanstd(stepMeans(probeAfterOptoTrials,:)))./sqrt(length(probeAfterOptoTrials));
probeAfterProbeError = (nanstd(stepMeans(probeAfterProbeTrials,:)))./sqrt(length(probeAfterProbeTrials));


figure,
subplot(1,2,1)
[opT] = boundedline(ydimension,probeAfterOptoMeans,probeAfterOptoError,'-ko','alpha');
ylabel('Forward velocity (mm/s)');
xlabel('panels y dimension');
title('Trials preceded by opto trials')
hold on
line([7 7],[0 max(probeAfterOptoMeans)+1],'Color' ,'b','LineStyle', '-.')
ylim([0 max(probeAfterOptoMeans)+1])
xlim([0 13])

subplot(1,2,2)
[prT] = boundedline(ydimension,probeAfterProbeMeans,probeAfterProbeError,'-ko','alpha');
xlabel('panels y dimension');
title('Trials preceded by probe trials')
line([7 7],[0 max(probeAfterOptoMeans)+1],'Color', 'b', 'LineStyle', '-.')
ylim([0 max(probeAfterOptoMeans)+1])
xlim([0 13])

saveas(gcf,[dirName(1:end-4),'plots\MeanVelvsPanelsPrecedingTrialType.png'])
close;

%% Evaluate response time in probe trials that followed opto trials in relationship with the time of reward in those trials
%If the animals are using time instead of distance, it might be hard to see
%when looking at every trial combined, since they walk at different speed
%in each trial. I decided here to look at the 'response' time (as the max
%slowing down) in probe trials with respect to the time of reward in the
%preceding opto trial, so see if I find a correlation

for i = 1:length(smoothed)
    
    if any(probeAfterOptoTrials == i) %for the probe trials following an opto trial
          [localminProbe{i},PromProbe{i}] = islocalmin(smoothed{1,i}.xVel); %find local minima  
    
    elseif any(probeAfterOptoTrials == i+1) %for the opto trials preceding a probe trial
          [localminOpto{i},PromOpto{i}] = islocalmin(smoothed{1,i}.xVel); %find local minima  

    end
       
end

localminProbe = localminProbe(~cellfun('isempty',localminProbe));
localminOpto = localminOpto(~cellfun('isempty',localminOpto));
PromProbe = PromProbe(~cellfun('isempty',PromProbe));
PromOpto = PromOpto(~cellfun('isempty',PromOpto));

for i = 1:length(localminProbe)
    if sum(localminProbe{i}) > 0
        [MProbe{i},IProbe{i}] = max(PromProbe{1,i}); %save the index of the biggest valley
    else
        MProbe{i} = NaN; 
        IProbe{i} = NaN; 
    end
    
    if sum(localminOpto{i}) > 0
        [MOpto{i},IOpto{i}] = max(PromOpto{1,i}); %save the index of the biggest valley
    else
        MOpto{i} = NaN; 
        IOpto{i} = NaN;
    end
end

IProbe = cell2mat(IProbe);
IOpto = cell2mat(IOpto);

%Calculate the Pearson correlation
R = corrcoef(IOpto,IProbe,'Rows','complete');

figure,
scatter(IOpto,IProbe)
ylabel('Time estimate in probe trial');
xlabel('Time of reward in preceding opto trial');
% add the correlation line
hold on
plot(IOpto,IOpto*R(2),'-.k')

saveas(gcf,[dirName(1:end-4),'plots\CorrTimeEstimate.png'])
close;

%% Local minimum analysis

%This analysis uses the local mins in fwd velocity as the animal's estimate
%of the reward location and we visualize them in a boxplot
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
        estimate(i) = ydimension(I{i});
    else
        estimate(i) = nan;
    end
end

%Proportion of probe trials with a local minimum
PropMin = sum(~isnan(estimate))*100/length(WorkingProbeTrials);


figure,
boxplot(estimate,'orientation','horizontal')
hold on
scatter(estimate,ones(1,length(estimate)).*(1+(rand(1,length(estimate))-0.5)/10),'k','filled')
xlim([1 13]); ylim([0 2]);
xlabel('y dimension');
xline(7,'--r','lineWidth',2)
text(1,1,[num2str(PropMin),' %']);
set(gca,'ytick',[])
title('Estimated reward distance (local minimum)');

saveas(gcf,[dirName(1:end-4),'plots\AllLocalMin.png'])
close;

%Dividing early and late trials
EarlyEstimate = estimate(1:floor(length(estimate)/2));
LateEstimate = estimate(floor(length(estimate)/2)+1:end);

AllEstimates = [EarlyEstimate,LateEstimate];
Groups = [ones(1,length(EarlyEstimate)),2*(ones(1,length(LateEstimate)))];

PropMinEarly = sum(~isnan(EarlyEstimate))*100/length(earlyProbeTrials);
PropMinLate = sum(~isnan(LateEstimate))*100/length(lateProbeTrials);

figure,
boxplot(AllEstimates(~isnan(AllEstimates)),Groups(~isnan(AllEstimates)),'orientation','horizontal')
hold on
jitter = [ones(1,length(EarlyEstimate)).*(1+(rand(1,length(EarlyEstimate))-0.5)/10),2*ones(1,length(LateEstimate)).*(1+(rand(1,length(LateEstimate))-0.5)/10)];
scatter(AllEstimates(~isnan(AllEstimates)),jitter(~isnan(AllEstimates)),'k','filled')
xlim([1 13]); ylim([0 3]);
xlabel('y dimension');
xline(7,'--r','lineWidth',2)
text(1,1,[num2str(PropMinEarly),' %']);
text(1,2,[num2str(PropMinLate),' %']);
set(gca,'ytick',[])
yticks([1 2]);
yticklabels({'Early trials','Late trials'});
title('Estimated reward distance (local minimum)');

saveas(gcf,[dirName(1:end-4),'plots\EarlyLateLocalMin.png'])
close;

%dividing following opto and probe
afterProbeEstimate = estimate;
afterOptoEstimate = estimate;

for i = 1:length(estimate)
    if any(probeAfterProbeTrials == WorkingProbeTrials(i)) 
        afterProbeEstimate(i) = estimate(i);
        afterOptoEstimate(i) = 1000;
    else
        afterProbeEstimate(i) = 1000;
        afterOptoEstimate(i) = estimate(i);
    end
end

afterProbeEstimate(afterProbeEstimate == 1000) = [];
afterOptoEstimate(afterOptoEstimate == 1000) = [];

AllEstimatesTrialType = [afterProbeEstimate,afterOptoEstimate];
Groups = [ones(1,length(afterProbeEstimate)),2*(ones(1,length(afterOptoEstimate)))];

PropMinProbe = sum(~isnan(afterProbeEstimate))*100/length(probeAfterProbeTrials);
PropMinOpto = sum(~isnan(afterOptoEstimate))*100/length(probeAfterOptoTrials);

figure,
boxplot(AllEstimatesTrialType(~isnan(AllEstimatesTrialType)),Groups(~isnan(AllEstimatesTrialType)),'orientation','horizontal')
hold on
jitter = [ones(1,length(afterProbeEstimate)).*(1+(rand(1,length(afterProbeEstimate))-0.5)/10),2*ones(1,length(afterOptoEstimate)).*(1+(rand(1,length(afterOptoEstimate))-0.5)/10)];
scatter(AllEstimatesTrialType(~isnan(AllEstimatesTrialType)),jitter(~isnan(AllEstimatesTrialType)),'k','filled')
xlim([1 13]); ylim([0 3]);
xlabel('y dimension');
xline(7,'--r','lineWidth',2)
text(1,1,[num2str(PropMinProbe),' %']);
text(1,2,[num2str(PropMinOpto),' %']);
set(gca,'ytick',[])
yticks([1 2]);
yticklabels({'After probe trials','After opto trials'});
title('Estimated reward distance (local minimum)');

saveas(gcf,[dirName(1:end-4),'plots\PrecedingTrialTypeLocalMin.png'])
close;

save([dirName(1:end-4),'dataFromAnalysis\distanceEstimates'],'AllEstimates','AllEstimatesTrialType','estimate','PropMinProbe','PropMinOpto','PropMinEarly','PropMinLate','PropMin');

%% Adding a similar analysis for maxima of angular velocity

for i = 1:length(WorkingProbeTrials)
    [localmax{i},Prom{i}] = islocalmax(angularStepMeans(WorkingProbeTrials(i),:)); %find local minima  
    if sum(localmax{i}) > 0
        [M{i},I{i}] = max(Prom{1,i}); %save the index of the biggest valley
    else
        I{i} = [];
    end
%      figure,
%      plot(ydimension,angularStepMeans(probeTrials(i),:))
%      hold on
%      plot(ydimension(I{i}),angularStepMeans(probeTrials(i),I{i}),'ro')
%      xlabel('y dimension'); ylabel('Ang speed (deg/s)');
    
    if sum(localmax{i}) > 0
        angEstimate(i) = ydimension(I{i});
    else
        angEstimate(i) = nan;
    end
end

PropMax = sum(~isnan(angEstimate))*100/length(WorkingProbeTrials);

figure,
boxplot(angEstimate,'orientation','horizontal')
hold on
scatter(angEstimate,ones(1,length(angEstimate)).*(1+(rand(1,length(angEstimate))-0.5)/10),'k','filled')
xlim([1 13]); ylim([0 2]);
xlabel('y dimension');
xline(7,'--r','lineWidth',2)
set(gca,'ytick',[])
text(1,1,[num2str(PropMax),' %']);
title('Estimated reward distance (angular speed maximum)');

saveas(gcf,[dirName(1:end-4),'plots\AllAngMax.png'])
close;

%Dividing early and late trials
EarlyAngEstimate = angEstimate(1:floor(length(angEstimate)/2));
LateAngEstimate = angEstimate(floor(length(angEstimate)/2)+1:end);

AllAngEstimates = [EarlyAngEstimate,LateAngEstimate];
AngGroups = [ones(1,length(EarlyAngEstimate)),2*(ones(1,length(LateAngEstimate)))];

PropMaxEarly = sum(~isnan(EarlyAngEstimate))*100/length(earlyProbeTrials);
PropMaxLate = sum(~isnan(LateAngEstimate))*100/length(lateProbeTrials);

figure,
boxplot(AllAngEstimates(~isnan(AllAngEstimates)),Groups(~isnan(AllAngEstimates)),'orientation','horizontal')
hold on
jitter = [ones(1,length(EarlyAngEstimate)).*(1+(rand(1,length(EarlyAngEstimate))-0.5)/10),2*ones(1,length(LateAngEstimate)).*(1+(rand(1,length(LateAngEstimate))-0.5)/10)];
scatter(AllAngEstimates(~isnan(AllAngEstimates)),jitter(~isnan(AllAngEstimates)),'k','filled')
xlim([1 13]); ylim([0 3]);
xlabel('y dimension');
xline(7,'--r','lineWidth',2)
text(1,1,[num2str(PropMaxEarly),' %']);
text(1,2,[num2str(PropMaxLate),' %']);
set(gca,'ytick',[])
yticks([1 2]);
yticklabels({'Early trials','Late trials'});
title('Estimated reward distance (angular speed max)');

saveas(gcf,[dirName(1:end-4),'plots\EarlyLateAngMax.png'])
close;

%dividing following opto and probe
afterProbeAngEstimate = angEstimate;
afterOptoAngEstimate = angEstimate;

for i = 1:length(angEstimate)
    if any(probeAfterProbeTrials == probeTrials(i)) 
        afterProbeAngEstimate(i) = angEstimate(i);
        afterOptoAngEstimate(i) = 1000;
    else
        afterProbeAngEstimate(i) = 1000;
        afterOptoAngEstimate(i) = angEstimate(i);
    end
end

afterProbeAngEstimate(afterProbeAngEstimate == 1000) = [];
afterOptoAngEstimate(afterOptoAngEstimate == 1000) = [];

AllAngEstimatesTrialType = [afterProbeAngEstimate,afterOptoAngEstimate];
AngGroups = [ones(1,length(afterProbeAngEstimate)),2*(ones(1,length(afterOptoAngEstimate)))];

PropMaxProbe = sum(~isnan(afterProbeAngEstimate))*100/length(probeAfterProbeTrials);
PropMaxOpto = sum(~isnan(afterOptoAngEstimate))*100/length(probeAfterOptoTrials);

figure,
boxplot(AllAngEstimatesTrialType(~isnan(AllAngEstimatesTrialType)),AngGroups(~isnan(AllAngEstimatesTrialType)),'orientation','horizontal')
hold on
jitter = [ones(1,length(afterProbeAngEstimate)).*(1+(rand(1,length(afterProbeAngEstimate))-0.5)/10),2*ones(1,length(afterOptoAngEstimate)).*(1+(rand(1,length(afterOptoAngEstimate))-0.5)/10)];
scatter(AllAngEstimatesTrialType(~isnan(AllAngEstimatesTrialType)),jitter(~isnan(AllAngEstimatesTrialType)),'k','filled')
xlim([1 13]); ylim([0 3]);
xlabel('y dimension');
xline(7,'--r','lineWidth',2)
text(1,1,[num2str(PropMaxProbe),' %']);
text(1,2,[num2str(PropMaxOpto),' %']);
set(gca,'ytick',[])
yticks([1 2]);
yticklabels({'After probe trials','After opto trials'});
title('Estimated reward distance (angular speed maximum)');

saveas(gcf,[dirName(1:end-4),'plots\PrecedingTrialTypeLocalMax.png'])
close;

save([dirName(1:end-4),'dataFromAnalysis\angDistanceEstimates'],'AllAngEstimates','AllAngEstimatesTrialType','angEstimate','PropMaxProbe','PropMaxOpto','PropMaxEarly','PropMaxLate','PropMax');


end