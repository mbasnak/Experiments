%Code to analyze Jenny's data with trials of darkness

%% Fly 1 - 60D05, 6f

%Darkness
DARKdata = load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\JennyPanelsOffExp\data\20190329_60D05_6f\dark trial\analysis_sid_8_tid_0.mat');

%get timepoint when the panels turn off and on
panelsOFF = find(DARKdata.data.fr_y_ds == 1);
StartDark = panelsOFF(1);
EndDark = panelsOFF(end);

figure('Position',[100 100 1000 800]),
subplot(3,1,1)
imagesc(DARKdata.data.dff_matrix)
colormap(gray)
hold on
xline(StartDark,'lineWidth',3,'color','r')
xline(EndDark,'lineWidth',3,'color','r')
ylabel('PB glomerulus');
set(gca,'xticklabel',{[]});
title('Fly 1 (60D05-6f), darkness');

subplot(3,1,2)
%stimulus position
stim = rad2deg(DARKdata.data.panel_angle);
changeStim = [0,diff(stim)];
stimToPlot = smooth(stim);
stimToPlot(changeStim>40) = NaN;
plot(DARKdata.data.time, stimToPlot,'lineWidth',1.5)
hold on
%Phase
phase = rad2deg(DARKdata.data.phase);
changePhase = [0,diff(phase)];
phaseToPlot = smooth(phase);
phaseToPlot(changePhase>40) = NaN;
plot(DARKdata.data.time,phaseToPlot,'color',[0.3 0.8 0.4],'lineWidth',1.5)
xline(DARKdata.data.time(StartDark),'lineWidth',3,'color','r')
xline(DARKdata.data.time(EndDark),'lineWidth',3,'color','r')
ylim([-180 180]);
xlim([0 DARKdata.data.time(end)]);
set(gca,'xticklabel',{[]})
ylabel('Deg');
legend('Fly heading', 'EPG phase');

subplot(3,1,3)
%Offset
offset = rad2deg(DARKdata.data.offset);
changeOffset = [0,diff(offset)];
offsetToPlot = smooth(offset);
offsetToPlot(changeOffset > 40) = NaN;
plot(DARKdata.data.time,offsetToPlot,'k','lineWidth',1.5)
hold on
xline(DARKdata.data.time(StartDark),'lineWidth',3,'color','r')
xline(DARKdata.data.time(EndDark),'lineWidth',3,'color','r')
ylim([-180 180]);
xlim([0 DARKdata.data.time(end)]);
xlabel('Time (s)');
ylabel('Offset (deg)');

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\fly1.png')

%% Fly 2 - 60D05, 7f

% darkness
DARKdata = load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\JennyPanelsOffExp\data\20190411_60D05_7f_2\dark trial\analysis_sid_2_tid_0.mat');

%get timepoint when the panels turn off and on
panelsOFF = find(DARKdata.data.fr_y_ds == 1);
StartDark = panelsOFF(1);
EndDark = panelsOFF(end);

figure('Position',[100 100 1000 800]),
subplot(3,1,1)
imagesc(DARKdata.data.dff_matrix)
colormap(gray)
hold on
xline(StartDark,'lineWidth',3,'color','r');
xline(EndDark,'lineWidth',3,'color','r');
ylabel('PB glomerulus');
set(gca,'xticklabel',{[]});
title('Fly 2 (60D05-7f), darkness');

subplot(3,1,2)
%stimulus position
stim = rad2deg(DARKdata.data.panel_angle);
changeStim = [0,diff(stim)];
stimToPlot = smooth(stim);
stimToPlot(changeStim>40) = NaN;
plot(DARKdata.data.time, stimToPlot,'lineWidth',1.5)
hold on
%Phase
phase = rad2deg(DARKdata.data.phase);
changePhase = [0,diff(phase)];
phaseToPlot = smooth(phase);
phaseToPlot(changePhase>40) = NaN;
plot(DARKdata.data.time,phaseToPlot,'color',[0.3 0.8 0.4],'lineWidth',1.5)
xline(DARKdata.data.time(StartDark),'lineWidth',3,'color','r')
xline(DARKdata.data.time(EndDark),'lineWidth',3,'color','r')
ylim([-180 180]);
xlim([0 DARKdata.data.time(end)]);
set(gca,'xticklabel',{[]})
ylabel('Deg');
legend('Fly heading', 'EPG phase');

subplot(3,1,3)
%Offset
offset = rad2deg(DARKdata.data.offset);
changeOffset = [0,diff(offset)];
offsetToPlot = smooth(offset);
offsetToPlot(changeOffset > 40) = NaN;
plot(DARKdata.data.time,offsetToPlot,'k','lineWidth',1.5)
hold on
xline(DARKdata.data.time(StartDark),'lineWidth',3,'color','r')
xline(DARKdata.data.time(EndDark),'lineWidth',3,'color','r')
ylim([-180 180]);
xlim([0 DARKdata.data.time(end)]);
xlabel('Time (s)');
ylabel('Offset (deg)');

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\fly2.png')

%% Fly 3 - 60D05, 7f

% darkness
DARKdata = load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\JennyPanelsOffExp\data\20190415_60D05_7f\dark trial\analysis_sid_4_tid_0.mat');

%get timepoint when the panels turn off and on
panelsOFF = find(DARKdata.data.fr_y_ds == 1);
StartDark = panelsOFF(1);
EndDark = panelsOFF(end);

figure('Position',[100 100 1000 800]),
subplot(3,1,1)
imagesc(DARKdata.data.dff_matrix)
colormap(gray)
hold on
xline(StartDark,'lineWidth',3,'color','r');
xline(EndDark,'lineWidth',3,'color','r');
ylabel('PB glomerulus');
set(gca,'xticklabel',{[]});
title('Fly 3 (60D05-7f), darkness');

subplot(3,1,2)
%stimulus position
stim = rad2deg(DARKdata.data.panel_angle);
changeStim = [0,diff(stim)];
stimToPlot = smooth(stim);
stimToPlot(changeStim>40) = NaN;
plot(DARKdata.data.time, stimToPlot,'lineWidth',1.5)
hold on
%Phase
phase = rad2deg(DARKdata.data.phase);
changePhase = [0,diff(phase)];
phaseToPlot = smooth(phase);
phaseToPlot(changePhase>40) = NaN;
plot(DARKdata.data.time,phaseToPlot,'color',[0.3 0.8 0.4],'lineWidth',1.5)
xline(DARKdata.data.time(StartDark),'lineWidth',3,'color','r')
xline(DARKdata.data.time(EndDark),'lineWidth',3,'color','r')
ylim([-180 180]);
xlim([0 DARKdata.data.time(end)]);
set(gca,'xticklabel',{[]})
ylabel('Deg');
legend('Fly heading', 'EPG phase');

subplot(3,1,3)
%Offset
offset = rad2deg(DARKdata.data.offset);
changeOffset = [0,diff(offset)];
offsetToPlot = smooth(offset);
offsetToPlot(changeOffset > 40) = NaN;
plot(DARKdata.data.time,offsetToPlot,'k','lineWidth',1.5)
hold on
xline(DARKdata.data.time(StartDark),'lineWidth',3,'color','r')
xline(DARKdata.data.time(EndDark),'lineWidth',3,'color','r')
ylim([-180 180]);
xlim([0 DARKdata.data.time(end)]);
xlabel('Time (s)');
ylabel('Offset (deg)');

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\fly3.png')

%% Group bump analyses

% files = dir('Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\JennyData\*\dark trial\*');
% 
% Data = struct;
% Data.StartDarkness = [];
% Data.EndDarkness = [];
% Data.offset = [];
% Data.dff_matrix = [];
% Data.phase_value = [];
% Data.time = [];
% Data.panel_angle = [];
% Data.ang_vel = [];
% Data.for_vel = [];
% 
% %load darkness data from each fly in a loop
% for i = 1:3
%     data = load([files(i).folder,'\',files(i).name]);
% %save start and end of darkness
%     panelsOFF = find(data.data.fr_y_ds == 1);
%     Data(i).StartDarkness = panelsOFF(1);
%     Data(i).EndDarkness = panelsOFF(end);
% %save offset
%     Data(i).offset = data.data.offset;
% %save time
%     Data(i).time = data.data.time;
% %save dff_matrix
%     Data(i).dff_matrix = data.data.dff_matrix;
% %save phase
%     Data(i).phase_value = data.data.phase;
% %save panel angle
%     Data(i).panel_angle = data.data.panel_angle;
% %save angular velocity
%     Data(i).ang_vel = data.data.vel_yaw_ds;
% %save forwards velocity
%     Data(i).for_vel = data.data.vel_for_ds;
% %save ft power
%     Data(i).ftpower = data.data.ftpower;
% %save PVA
%     Data(i).pva = data.data.dff_pva_rad;
% end
% 
% save('Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\allDarkData.mat','Data');


%load group data
load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\JennyPanelsOffExp\data\allDarkData');

%% 1) Offset distribution in each case

figure('Position',[100 100 1000 800]),
subplot(3,3,1)
polarhistogram(Data(1).offset(1:Data(1).StartDarkness-1),15,'FaceColor',[0.2, 0.4, 1])
title('First certain state');
subplot(3,3,2)
polarhistogram(Data(1).offset(Data(1).StartDarkness:Data(1).EndDarkness),15,'FaceColor',[0.2, 0.4, 1])
title('Uncertain state');
subplot(3,3,3)
polarhistogram(Data(1).offset(Data(1).EndDarkness+1:end),15,'FaceColor',[0.2, 0.4, 1])
title('Second certain state');

subplot(3,3,4)
polarhistogram(Data(2).offset(1:Data(2).StartDarkness-1),15, 'FaceColor', [0.6 0.1 0.5])
subplot(3,3,5)
polarhistogram(Data(2).offset(Data(2).StartDarkness:Data(2).EndDarkness),15, 'FaceColor', [0.6 0.1 0.5])
subplot(3,3,6)
polarhistogram(Data(2).offset(Data(2).EndDarkness+1:end),15, 'FaceColor', [0.6 0.1 0.5])

subplot(3,3,7)
polarhistogram(Data(3).offset(1:Data(3).StartDarkness-1),15, 'FaceColor', [0.8 0.3 0.3])
subplot(3,3,8)
polarhistogram(Data(3).offset(Data(3).StartDarkness:Data(3).EndDarkness),15, 'FaceColor', [0.8 0.3 0.3])
subplot(3,3,9)
polarhistogram(Data(3).offset(Data(3).EndDarkness+1:end),15, 'FaceColor', [0.8 0.3 0.3])

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\OffsetDistribution.png')

%% variation of offset distribution per state
offset_variation = {};

for i = 1:3
    offset_variation{i}(1) = circ_std(Data(i).offset(1:Data(i).StartDarkness-1),[],[],2);
    offset_variation{i}(2) = circ_std(Data(i).offset(Data(i).StartDarkness:Data(i).EndDarkness),[],[],2);
    offset_variation{i}(3) = circ_std(Data(i).offset(Data(i).EndDarkness+1:end),[],[],2);    
end

figure('Position',[300 300 800 600]),
plotColors(1,:) = [0.2, 0.4, 1];
plotColors(2,:) = [0.6 0.1 0.5];
plotColors(3,:) = [0.8 0.3 0.3];
for i = 1:3
    plot(offset_variation{1,i},'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    hold on
end
xlim([0 4]);
ylabel('Circular std of offset');
set(gca,'xticklabel',{[]});
xticks([1 2 3]);
xticklabels({'First certain state','Uncertain state','Second certain state'})
legend('Fly 1', 'Fly 2', 'Fly 3');
title('Offset variation');

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\OffsetVariation.png')


%% 2) Bump magnitude in each case
    %I will define the bump magnitude as the difference between the max and
    %min fluorescence at each time point, although I'm not sure this is
    %correct
BumpMagnitude = {};

for i = 1:3
   BumpMagnitude{i,1} = max(Data(i).dff_matrix(:,1:Data(i).StartDarkness-1))-min(Data(1).dff_matrix(:,1:Data(i).StartDarkness-1));
   BumpMagnitude{i,2} = max(Data(i).dff_matrix(:,Data(i).StartDarkness:Data(i).EndDarkness))-min(Data(1).dff_matrix(:,Data(i).StartDarkness:Data(i).EndDarkness));
   BumpMagnitude{i,3} = max(Data(i).dff_matrix(:,Data(i).EndDarkness+1:end))-min(Data(1).dff_matrix(:,Data(i).EndDarkness+1:end));   
end

figure('Position',[400 200 1000 800]),
for i = 1:3
    subplot(3,1,i)
    plot(Data(i).time(1:Data(i).StartDarkness-1),BumpMagnitude{i,1},'k.')
    hold on
    plot(Data(i).time(Data(i).StartDarkness:Data(i).EndDarkness),BumpMagnitude{i,2},'r.')
    plot(Data(i).time(Data(i).EndDarkness+1:end),BumpMagnitude{i,3},'k.')
    title(['Fly ',num2str(i)]);
    ylabel('Bump magnitude');
end

%% Plot median and iqr bump magnitude for each fly

medianBumpMagnitude = cellfun(@median,BumpMagnitude);
iqrBumpMagnitude = cellfun(@iqr,BumpMagnitude);

figure('Position',[300 300 800 600]),
for i = 1:3
    errorbar(medianBumpMagnitude(i,:),iqrBumpMagnitude(i,:),'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    hold on
end
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'First certain state','Uncertain state','Second certain state'})
legend('Fly 1', 'Fly 2', 'Fly 3');
title('Median bump magnitude across states');
ylabel('Bump magnitude')

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\MedianBumpMagnitude.png')

%% Scatter of bump magnitude vs fly angular speed

figure('Position',[200, 200, 1000, 800]),
for i = 1:3
    subplot(1,3,i)
    scatter(abs(Data(i).ang_vel(1:Data(i).StartDarkness-1)),BumpMagnitude{i,1},[],'k')
    hold on
    scatter(abs(Data(i).ang_vel(Data(i).StartDarkness:Data(i).EndDarkness)),BumpMagnitude{i,2},[],'r')
    scatter(abs(Data(i).ang_vel(Data(i).EndDarkness+1:end)),BumpMagnitude{i,3},[],'k')
    ylabel('Bump magnitude');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('Certain state', 'Uncertain state');
end

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\BumpMagnitudevsAngSpeed.png')

%% Binning the velocity

velBins = [0:0.5:4]; %I creat a vector with my bins
velBins = [velBins,max(abs(Data(1).ang_vel(Data(1).EndDarkness+1:end)))+1];

%getting binned medians by state
for i = 1:3
    for j = 1:length(velBins)-1
        medianBin{i}(j) = median(BumpMagnitude{i,1}((abs(Data(i).ang_vel(1:Data(i).StartDarkness-1)))>velBins(j) & (abs(Data(i).ang_vel(1:Data(i).StartDarkness-1)))<velBins(j+1)));
        medianBin2{i}(j) = median(BumpMagnitude{i,3}((abs(Data(i).ang_vel(Data(i).EndDarkness+1:end)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).EndDarkness+1:end)))<velBins(j+1)));
        medianCertainBin{i}(j) = (medianBin{i}(j)+medianBin2{i}(j))/2;
        medianUncertainBin{i}(j) = median(BumpMagnitude{i,2}((abs(Data(i).ang_vel(Data(i).StartDarkness:Data(i).EndDarkness)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).StartDarkness:Data(i).EndDarkness)))<velBins(j+1)));        
    end   
end

velAxes = velBins - 0.5;
velAxes = velAxes(2:end);
velAxes(end) = velAxes(end-1)+0.5;
figure('Position',[200, 200, 1000, 800]),
for i = 1:3
    subplot(1,3,i)
    plot(velAxes, medianCertainBin{1,i},'-ko')
    hold on
    plot(velAxes, medianUncertainBin{1,i},'-ro')
    ylabel('Bump magnitude');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('Certain state', 'Uncertain state');
    ylim([1 6]);
end

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\BinnedBumpMagnitudevsAngSpeed.png')

%%  Bump magnitude in time
BM = {};
BM{1} = [BumpMagnitude{1,1},BumpMagnitude{1,2},BumpMagnitude{1,3}];
BM{2} = [BumpMagnitude{2,1},BumpMagnitude{2,2},BumpMagnitude{2,3}];
BM{3} = [BumpMagnitude{3,1},BumpMagnitude{3,2},BumpMagnitude{3,3}];

figure('Position',[100 100 1600 1000]),
for i = 1:3
    subplot(3,1,i)
    plot(Data(i).time(1:Data(i).StartDarkness-1),BM{i}(1:Data(i).StartDarkness-1),'k')
    hold on
    plot(Data(i).time(Data(i).StartDarkness:Data(i).EndDarkness),BM{i}(Data(i).StartDarkness:Data(i).EndDarkness),'r')
    plot(Data(i).time(Data(i).EndDarkness+1:end),BM{i}(Data(i).EndDarkness+1:end),'k') 
    xlim([0 Data(i).time(end)]);
    ylabel('DeltaF/F');
end
xlabel('Time (s)');
suptitle('Change in bump magnitude');

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\BumpMagnitudeInTime.png')


% Close-up to the change
figure('Position',[100 100 800 1000]),
for i = 1:3
    subplot(3,1,i)
    plot(Data(i).time(1:Data(i).StartDarkness),BM{i}(1:Data(i).StartDarkness),'k')
    hold on
    plot(Data(i).time(Data(i).StartDarkness:Data(i).EndDarkness),BM{i}(Data(i).StartDarkness:Data(i).EndDarkness),'r')
    plot(Data(i).time(Data(i).EndDarkness:end),BM{i}(Data(i).EndDarkness:end),'k') 
    xlim([Data(i).time(Data(i).StartDarkness-100), Data(i).time(Data(i).StartDarkness+150)]);
    ylabel('DeltaF/F');
end
xlabel('Time (s)');
suptitle('Change in bump magnitude');

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\BumpMagnitudeInTimeCloseUp.png')

%%  Bump width at half max in each case

    %to find the width at half max, example:
%I'm going to linearly interpolate to generate more datapoints and ensure
%the half max is met
exampleData = interp1([1:16],Data(1).dff_matrix(:,1),linspace(1,16,1000));
figure, plot(exampleData)
set(gca,'xticklabel',{[]});
xlabel('PB glomerulus');
ylabel('DFF');
[M,I] = max(exampleData);
hold on
plot(I,M,'ro')
%find closest values to half max
[M1,FirstHalf] = min(abs(exampleData(1:I)-(M/2)));
[M2,SecondHalf] = min(abs(exampleData(I+1:end)-(M/2)));
line([FirstHalf SecondHalf+I],[exampleData(FirstHalf) exampleData(FirstHalf)],'color','r','LineWidth',2)
exampleWidth = (SecondHalf + I - FirstHalf)*16/1000;

    %Doing this for each time point across flies for the different states
interpData = {};
for i = 1:3    
    interpData{i} = interp1([1:16],Data(i).dff_matrix(:,:),linspace(1,16,1000));
    [Max{i},Ind{i}] = max(interpData{i});
    [Max1{i},FirstHalfMax{i}] = min(abs(interpData{i}(1:Ind{i},:)-(Max{i}/2)));
    [Max2{i},SecondHalfMax{i}] = min(abs(interpData{i}(Ind{i}+1:end,:)-(Max{i}/2)));
    HalfWidth{i} = (SecondHalfMax{i} + Ind{i} - FirstHalfMax{i})*16/1000;
end
   
%divide the data by state
HalfMaxWidth = {};
for i = 1:3
   HalfMaxWidth{i,1} = HalfWidth{1,i}(1:Data(i).StartDarkness-1);
   HalfMaxWidth{i,2} = HalfWidth{1,i}(Data(i).StartDarkness:Data(i).EndDarkness);
   HalfMaxWidth{i,3} = HalfWidth{1,i}(Data(i).EndDarkness+1:end);   
end

medianHalfMaxWidth = cellfun(@median,HalfMaxWidth);
iqrHalfMaxWidth = cellfun(@iqr,HalfMaxWidth);

figure('Position',[300 300 800 600]),
for i = 1:3
    errorbar(medianHalfMaxWidth(i,:),iqrHalfMaxWidth(i,:),'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    hold on
end
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'First certain state','Uncertain state','Second certain state'})
legend('Fly 1', 'Fly 2', 'Fly 3');
title('Median bump width at half max across states');
ylabel('Bump width at half max')

%saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\WidthHalfMax.png')

%%  Bump angular speed distribution (and median?) in each case - I should
%compare with the fly's angular speed distribution to see that differences
%are not due to that.

%median bump angular speed
for i = 1:3
    unwrappedPhase{i} = unwrap(Data(i).phase_value);
    diffUnwrapped{i} = diff(unwrappedPhase{i});   
    smoothedDiff{i} = abs(smoothdata(diffUnwrapped{i}));
end

bumpAngSpeed = {};
for i = 1:3
   bumpAngSpeed{i,1} = smoothedDiff{1,i}(1:Data(i).StartDarkness-1);
   bumpAngSpeed{i,2} = smoothedDiff{1,i}(Data(i).StartDarkness:Data(i).EndDarkness);
   bumpAngSpeed{i,3} = smoothedDiff{1,i}(Data(i).EndDarkness+1:end);        
end

medianBumpAngSpeed = cellfun(@median,bumpAngSpeed);
iqrBumpAngSpeed = cellfun(@iqr,bumpAngSpeed);

figure('Position',[300 300 800 600]),
for i = 1:3
    errorbar(medianBumpAngSpeed(i,:),iqrBumpAngSpeed(i,:),'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    hold on
end
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'First certain state','Uncertain state','Second certain state'})
legend('Fly 1', 'Fly 2', 'Fly 3');
title('Median bump angular speed');
ylabel('Angular speed')

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\medianBumpAngSpeed.png')

%% add scatter of bump magnitude vs fly angular speed

figure('Position',[200, 200, 1000, 800]),
for i = 1:3
    subplot(1,3,i)
    scatter(abs(Data(i).ang_vel(1:Data(i).StartDarkness-1)),bumpAngSpeed{i,1},[],'k')
    hold on
    scatter(abs(Data(i).ang_vel(Data(i).StartDarkness:Data(i).EndDarkness)),bumpAngSpeed{i,2},[],'r')
    scatter(abs(Data(i).ang_vel(Data(i).EndDarkness+1:end-1)),bumpAngSpeed{i,3},[],'k')
    ylabel('Bump angular speed');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('Certain state', 'Uncertain state');
end

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\BumpSpeedvsFlySpeed.png')


%Binning the velocity
%getting binned medians by state
for i = 1:3
    for j = 1:length(velBins)-1
        MedianBin{i}(j) = median(bumpAngSpeed{i,1}((abs(Data(i).ang_vel(1:Data(i).StartDarkness-1)))>velBins(j) & (abs(Data(i).ang_vel(1:Data(i).StartDarkness-1)))<velBins(j+1)));
        MedianBin2{i}(j) = median(bumpAngSpeed{i,3}((abs(Data(i).ang_vel(Data(i).EndDarkness+1:end-1)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).EndDarkness+1:end-1)))<velBins(j+1)));
        MedianCertainBin{i}(j) = (MedianBin{i}(j)+ MedianBin2{i}(j))/2;
        MedianUncertainBin{i}(j) = median(bumpAngSpeed{i,2}((abs(Data(i).ang_vel(Data(i).StartDarkness:Data(i).EndDarkness)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).StartDarkness:Data(i).EndDarkness)))<velBins(j+1)));        
    end   
end

figure('Position',[200, 200, 1000, 800]),
for i = 1:3
    subplot(1,3,i)
    plot(velAxes, MedianCertainBin{1,i},'-ko')
    hold on
    plot(velAxes, MedianUncertainBin{1,i},'-ro')
    ylabel('Bump angular speed');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('Certain state', 'Uncertain state');
    ylim([0 0.3]);
end

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\BinnedBumpSpeedvsFlySpeed.png')

%% Dividing in walking and standing bouts

allWalkingVelocities = [Data(1).for_vel;Data(2).for_vel;Data(3).for_vel];

%Determine the threshold to say if the fly is walking or not.
figure,
histogram(allWalkingVelocities)

edges = [-5:0.2:10];
[counts] = histcounts(allWalkingVelocities,edges);
probabilities = counts./sum(counts);
degs = linspace(-5,10,length(counts));
figure,
plot(degs,probabilities,'k')
set(0, 'DefaulttextInterpreter', 'none')
[M,I] = max(probabilities);
[M2,I2] = min(abs(probabilities(I+1:end)-(M/2)));
thresholdVel = degs(I+I2);

%% Get and plot bump magnitude in walking and standing bouts for every state

AllVel = {};
for i = 1:3  
    AllVel{i} = Data(i).for_vel;
end
allVel = {};
for i = 1:3
   allVel{i,1} = AllVel{1,i}(1:Data(i).StartDarkness-1);
   allVel{i,2} = AllVel{1,i}(Data(i).StartDarkness:Data(i).EndDarkness);
   allVel{i,3} = AllVel{1,i}(Data(i).EndDarkness+1:end);        
end

%find data points above the vel threshold
walkingBouts = cellfun(@(x) x>thresholdVel, allVel, 'UniformOutput', 0);

%compare bump magnitude in walking and standing bouts for every state
walkingBM = {};
for i = 1:3
    for j = 1:3
        walkingBM{i,j} = [median(BumpMagnitude{i,j}(walkingBouts{i,j}));median(BumpMagnitude{i,j}(~walkingBouts{i,j}))];
    end
end

figure('Position',[200 200 1000 800]),
subplot(3,3,1)
plot(walkingBM{1,1},'LineWidth',2,'color',[0.2,0.4,1])
ylim([2 5.5]);
xlim([0 3]);
set(gca,'xticklabel',{[]});
title('First certain state');

subplot(3,3,2)
plot(walkingBM{1,2},'LineWidth',2,'color',[0.2,0.4,1])
ylim([2 5.5]);
xlim([0 3]);
set(gca,'xticklabel',{[]});
title('Uncertain state');

subplot(3,3,3)
plot(walkingBM{1,3},'LineWidth',2,'color',[0.2,0.4,1])
ylim([2 5.5]);
xlim([0 3]);
set(gca,'xticklabel',{[]});
title('Second certain state');

for i = 1:3
    subplot(3,3,i+3)
    plot(walkingBM{2,i},'LineWidth',2,'color',[0.6,0.1,0.5])
    ylim([2 5.5]);
    xlim([0 3]);
    set(gca,'xticklabel',{[]});
end
for i = 1:3
    subplot(3,3,i+6)
    plot(walkingBM{3,i},'LineWidth',2,'color',[0.8,0.3,0.3])
    ylim([2 5.5]);
    xlim([0 3]);
    xticks([1 2 ]);
    xticklabels({'Walking','Standing'})
end

saveas(gcf,'Z:\Wilson Lab\Mel\Lab_meetings\meetingWithJanAndAnna\plots\WalkingVsStanding.png')