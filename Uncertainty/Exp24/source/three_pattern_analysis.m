
clear all; close all;
%% Group bump analyses

% files = dir('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\*\*.mat');
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
% Data.flyPosRad = [];
% 
% %load data from each fly in a loop
% for i = 1:length(files)
%     data = load([files(i).folder,'\',files(i).name]);
% %save identity of stimulus for each frame
%     Data(i).singleBar = find(data.data.fr_y_ds<4);
%     Data(i).ambiguous = find(data.data.fr_y_ds>6 & data.data.fr_y_ds<8);
%     Data(i).darkness = find(data.data.fr_y_ds>9);
% %save offset
%     Data(i).offset = data.data.offset3;
% %save time
%     Data(i).time = data.data.time;
% %save dff_matrix
%     Data(i).dff_matrix = data.data.dff_matrix;
% %save phase
%     Data(i).phase_value = data.data.phase;
% %save panel angle I THINK THIS ONE IS CURRENTLY INCORRECT
%     Data(i).panel_angle = data.data.panel_angle;
% %save angular velocity
%     Data(i).ang_vel = data.data.vel_yaw_ds;
% %save forwards velocity
%     Data(i).for_vel = data.data.vel_for_ds;
% %save fly angular position
%     Data(i).flyPosRad = data.data.flyPosRad;
% 
% end
% %remove extra data from fly 4
% Data(4).offset = Data(4).offset(1,1:5510);
% Data(4).dff_matrix = Data(4).dff_matrix(:,1:5510);
% Data(4).phase_value = Data(4).phase_value(:,1:5510);
% Data(4).time = Data(4).time(1:5510,:);
% Data(4).panel_angle = Data(4).panel_angle(:,1:5510);
% Data(4).ang_vel = Data(4).ang_vel(:,1:5510);
% Data(4).for_vel = Data(4).for_vel(:,1:5510);
% Data(4).flyPosRad = Data(4).flyPosRad(1:5510,:);
% Data(4).singleBar = Data(4).singleBar(Data(4).singleBar<5510);
% Data(4).ambiguous = Data(4).ambiguous(Data(4).ambiguous<5510);
% Data(4).darkness = Data(4).darkness(Data(4).darkness<5510);
% 
% save('Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\allData.mat','Data');


%load group data
load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\ThreePatternExp\data\allData.mat');


%% variation of offset distribution per state
offset_variation = {};

for i = 1:length(Data)
    [~,offset_variation{i}(1)] = circ_std(Data(i).offset(Data(i).singleBar),[],[],2);
    [~,offset_variation{i}(2)] = circ_std(Data(i).offset(Data(i).ambiguous),[],[],2);
    [~,offset_variation{i}(3)] = circ_std(Data(i).offset(Data(i).darkness),[],[],2);    
end

figure('Position',[300 300 800 600]),
for i = 1:length(Data)
    plot(offset_variation{1,i},'-o')
    hold on
end
xlim([0 4]);
ylabel('Circular std of offset');
set(gca,'xticklabel',{[]});
xticks([1 2 3]);
xticklabels({'Vertical bar','Horizontal bar','Panels off'})
legend('Fly 1', 'Fly 2', 'Fly 3','Fly4','Fly5');
title('Offset variation');

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\OffsetVariation.png')


%% 2) Bump magnitude in each case
    %I will define the bump magnitude as the difference between the max and
    %min fluorescence at each time point, although I'm not sure this is
    %correct
BumpMagnitude = {};

for i = 1:length(Data)
   BumpMagnitude{i,1} = max(Data(i).dff_matrix(:,Data(i).singleBar))-min(Data(i).dff_matrix(:,Data(i).singleBar));
   BumpMagnitude{i,2} = max(Data(i).dff_matrix(:,Data(i).ambiguous))-min(Data(i).dff_matrix(:,Data(i).ambiguous));
   BumpMagnitude{i,3} = max(Data(i).dff_matrix(:,Data(i).darkness))-min(Data(i).dff_matrix(:,Data(i).darkness));   
end

figure('Position',[400 200 1000 800]),
for i = 1:length(Data)
    subplot(length(Data),1,i)
    plot(Data(i).time(Data(i).singleBar),BumpMagnitude{i,1},'b.')
    hold on
    plot(Data(i).time(Data(i).ambiguous),BumpMagnitude{i,2},'r.')
    plot(Data(i).time(Data(i).darkness),BumpMagnitude{i,3},'k.')
    title(['Fly ',num2str(i)]);
    ylabel('Bump magnitude');
    xlim([0 580]);
end

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\BumpMagInTime.png')


%% Plot median and iqr bump magnitude for each fly

medianBumpMagnitude = cellfun(@median,BumpMagnitude);
iqrBumpMagnitude = cellfun(@iqr,BumpMagnitude);

figure('Position',[300 300 800 600]),
for i = 1:length(Data)
    plot(medianBumpMagnitude(i,:),'-o')
    hold on
end
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'vertical bar','horizontal bar','panels off'})
legend('Fly 1', 'Fly 2', 'Fly 3', 'Fly 4', 'Fly 5');
title('Median bump magnitude across states');
ylabel('Bump magnitude')

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\MedianBumpMag.png')


figure('Position',[300 300 800 600]),
for i = 1:length(Data)
    errorbar(medianBumpMagnitude(i,:),iqrBumpMagnitude(i,:),'-o')
    hold on
end
xlim([0 4]);
xticks([1 2 3]);
xticklabels({'vertical bar','horizontal bar','panels off'})
legend('Fly 1', 'Fly 2', 'Fly 3', 'Fly 4', 'Fly 5');
title('Median bump magnitude across states');
ylabel('Bump magnitude')


%% Scatter of bump magnitude vs fly angular speed

figure('Position',[200, 200, 1000, 800]),
for i = 1:length(Data)
    subplot(1,length(Data),i)
    scatter(abs(Data(i).ang_vel(Data(i).singleBar)),BumpMagnitude{i,1},[],'b')
    hold on
    scatter(abs(Data(i).ang_vel(Data(i).ambiguous)),BumpMagnitude{i,2},[],'r')
    scatter(abs(Data(i).ang_vel(Data(i).darkness)),BumpMagnitude{i,3},[],'k')
    ylabel('Bump magnitude');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('vertical bar','horizontal bar','panels off');
end

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\BumpMagVsFlySpeed.png')


%% Binning the velocity

velBins = [0:0.5:4]; %I creat a vector with my bins
velBins = [velBins,5000];

%getting binned medians by state
for i = 1:length(Data)
    for j = 1:length(velBins)-1
        medianBin{i}(j) = median(BumpMagnitude{i,1}((abs(Data(i).ang_vel(Data(i).singleBar)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).singleBar)))<velBins(j+1)));
        medianBin2{i}(j) = median(BumpMagnitude{i,2}((abs(Data(i).ang_vel(Data(i).ambiguous)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).ambiguous)))<velBins(j+1)));
        medianBin3{i}(j) = median(BumpMagnitude{i,3}((abs(Data(i).ang_vel(Data(i).darkness)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).darkness)))<velBins(j+1)));        
    end   
end

velAxes = velBins - 0.5;
velAxes = velAxes(2:end);
velAxes(end) = velAxes(end-1)+0.5;
figure('Position',[200, 200, 1000, 800]),
for i = 1:length(Data)
    subplot(1,length(Data),i)
    plot(velAxes, medianBin{1,i},'-bo')
    hold on
    plot(velAxes, medianBin2{1,i},'-ro')
    plot(velAxes, medianBin3{1,i},'-ko')
    ylabel('Bump magnitude');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('vertical bar', 'horizontal bar', 'panels off');
    %ylim([1 4]);
end

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\BumpMagVsBinnedFlySpeed.png')


%%  Bump angular speed distribution (and median?) in each case - I should
%compare with the fly's angular speed distribution to see that differences
%are not due to that.

%median bump angular speed
for i = 1:length(Data)
    unwrappedPhase{i} = unwrap(Data(i).phase_value);
    diffUnwrapped{i} = diff(unwrappedPhase{i});   
    smoothedDiff{i} = [abs(smoothdata(diffUnwrapped{i})),abs(smoothdata(diffUnwrapped{i}(end)))];
end

bumpAngSpeed = {};
for i = 1:length(Data)
   bumpAngSpeed{i,1} = smoothedDiff{1,i}(Data(i).singleBar);
   bumpAngSpeed{i,2} = smoothedDiff{1,i}(Data(i).ambiguous);
   bumpAngSpeed{i,3} = smoothedDiff{1,i}(Data(i).darkness);        
end

medianBumpAngSpeed = cellfun(@median,bumpAngSpeed);
iqrBumpAngSpeed = cellfun(@iqr,bumpAngSpeed);

figure('Position',[300 300 800 600]),
for i = 1:length(Data)
    %errorbar(medianBumpAngSpeed(i,:),iqrBumpAngSpeed(i,:),'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    plot(medianBumpAngSpeed(i,:),'-o')    
    hold on
end
xlim([0 4]);
xticks([1 2 3 4 5]);
xticklabels({'vertical bar','horizontal bar','panels off'})
legend('Fly 1', 'Fly 2', 'Fly 3', 'Fly4', 'Fly5');
title('Median bump angular speed');
ylabel('Angular speed')

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\BumpAngSpeed.png')

%% add scatter of bump magnitude vs fly angular speed

figure('Position',[200, 200, 1000, 800]),
for i = 1:length(Data)
    subplot(1,length(Data),i)
    scatter(abs(Data(i).ang_vel(Data(i).singleBar)),bumpAngSpeed{i,1},[],'b')
    hold on
    scatter(abs(Data(i).ang_vel(Data(i).ambiguous)),bumpAngSpeed{i,2},[],'r')
    scatter(abs(Data(i).ang_vel(Data(i).darkness)),bumpAngSpeed{i,3},[],'k')
    ylabel('Bump angular speed');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('vertical bar','horizontal bar','panels off');
end


%Binning the velocity
%getting binned medians by state
for i = 1:length(Data)
    for j = 1:length(velBins)-1
        MedianBin{i}(j) = median(bumpAngSpeed{i,1}((abs(Data(i).ang_vel(Data(i).singleBar)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).singleBar)))<velBins(j+1)));
        MedianBin2{i}(j) = median(bumpAngSpeed{i,2}((abs(Data(i).ang_vel(Data(i).ambiguous)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).ambiguous)))<velBins(j+1)));
        MedianBin3{i}(j) = median(bumpAngSpeed{i,3}((abs(Data(i).ang_vel(Data(i).darkness)))>velBins(j) & (abs(Data(i).ang_vel(Data(i).darkness)))<velBins(j+1)));        
    end   
end

figure('Position',[200, 200, 1000, 800]),
for i = 1:length(Data)
    subplot(1,length(Data),i)
    plot(velAxes, MedianBin{1,i},'-bo')
    hold on
    plot(velAxes, MedianBin2{1,i},'-ro')
    plot(velAxes, MedianBin3{1,i},'-ko')
    ylabel('Bump angular speed');
    xlabel('Fly angular speed');
    title(['Fly ',num2str(i)]);
    legend('vertical bar','horizontal bar','panels off');
    %ylim([0 0.3]);
end

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\BumpAngSpeedVsFlySpeed.png')



%% Dividing in walking and standing bouts

allWalkingVelocities = [Data(1).for_vel;Data(2).for_vel;Data(3).for_vel;Data(5).for_vel];

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
for i = 1:length(Data)  
    AllVel{i} = Data(i).for_vel;
end
allVel = {};
for i = 1:length(Data)
   allVel{i,1} = AllVel{1,i}(Data(i).singleBar);
   allVel{i,2} = AllVel{1,i}(Data(i).ambiguous);
   allVel{i,3} = AllVel{1,i}(Data(i).darkness);        
end

%find data points above the vel threshold
walkingBouts = cellfun(@(x) x>thresholdVel, allVel, 'UniformOutput', 0);

%compare bump magnitude in walking and standing bouts for every state
walkingBM = {};
for i = 1:length(Data)
    for j = 1:3
        walkingBM{i,j} = [median(BumpMagnitude{i,j}(walkingBouts{i,j}));median(BumpMagnitude{i,j}(~walkingBouts{i,j}))];
    end
end

%% 

bumpMagWalking1 =[3.2721,4.1522,1.6678,1.5872,1.6633];
bumpMagStanding1 =[2.3701,3.7680,1.4547,1.2562,1.2477];
bumpMagWalking2 =[3.0842,2.9703,1.3966,1.4529,1.4061];
bumpMagStanding2 =[2.5312,2.8572,1.1800,1.0873,1.1596];
bumpMagWalking3 =[3.5682,3.4546,1.6257,1.4509,1.4045];
bumpMagStanding3 =[2.2452,3.3437,1.2427,0.9581,1.1511];
medianBMState(1,1) = median(bumpMagStanding1);
medianBMState(1,2) = median(bumpMagWalking1);
medianBMState(2,1) = median(bumpMagStanding2);
medianBMState(2,2) = median(bumpMagWalking2);
medianBMState(3,1) = median(bumpMagStanding3);
medianBMState(3,2) = median(bumpMagWalking3);

plotColors(1,:) = [0, 0, 1];
plotColors(2,:) = [1 0 0];
plotColors(3,:) = [0 0 0];

figure,
for i = 1:3
    plot(medianBMState(i,:),'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    hold on
end
ylabel('Bump magnitude');
xticks([1 2]);
xlim([0 3]);
xticklabels({'Standing','Walking'})
legend('vertical bar','horizontal bar', 'panels off')

saveas(gcf,'Z:\Wilson Lab\Mel\FlyOnTheBall\data\Experiment17\three_patterns\useful trials\plots\WalkingVsStanding.png')


%% Bump speed vs walking state

%compare bump magnitude in walking and standing bouts for every state
walkingBSpeed = {};
for i = 1:length(Data)
    for j = 1:3
        walkingBSpeed{i,j} = [median(bumpAngSpeed{i,j}(walkingBouts{i,j}));median(bumpAngSpeed{i,j}(~walkingBouts{i,j}))];
    end
end

bumpASWalking1 =[0.0908 0.0671 0.1068 0.0664 0.0578];
bumpASStanding1 =[0.0851 0.0545 0.0833 0.0667 0.0678];
bumpASWalking2 =[0.0484 0.1155 0.1427 0.0995 0.0518];
bumpASStanding2 =[0.00530 0.1276 0.1704 0.1267 0.0677];
bumpASWalking3 =[0.0670 0.0978 0.1301 0.0985 0.0531];
bumpASStanding3 =[0.0884 0.1058 0.1650 0.1209 0.0857];
medianBSState(1,1) = median(bumpASStanding1);
medianBSState(1,2) = median(bumpASWalking1);
medianBSState(2,1) = median(bumpASStanding2);
medianBSState(2,2) = median(bumpASWalking2);
medianBSState(3,1) = median(bumpASStanding3);
medianBSState(3,2) = median(bumpASWalking3);


figure,
for i = 1:3
    plot(medianBSState(i,:),'-o','color',plotColors(i,:),'MarkerFaceColor',plotColors(i,:))
    hold on
end
ylabel('Bump angular speed');
xticks([1 2]);
xlim([0 3]);
xticklabels({'Standing','Walking'})
legend('vertical bar','horizontal bar', 'panels off')

%% von Mises fit

%We will next perform a von Mises fit for find the bump magnitude and width
%at half max

clear all; close all

leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];

load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\ThreePatternExp\data\allData.mat');

for fly = 1:length(Data)
    combined_full_dff{fly} = (Data(fly).dff_matrix(leftPB,:) + Data(fly).dff_matrix(rightPB,:))/2;
    half_width{fly} = zeros(1,length(combined_full_dff{fly}));
    bump_mag{fly} = zeros(1,length(combined_full_dff{fly}));

    for timepoint = 1:length(combined_full_dff{fly})    
        extendedData = interp1([1:8],combined_full_dff{fly}(:,timepoint),linspace(1,8,1000));
        angles = linspace(0,2*pi,length(extendedData));
        [vonMises, rescaledVonMises] = fitTuningCurveToVonMises(extendedData, angles);
        bump_mag{fly}(timepoint) = max(rescaledVonMises);
        fitData = circ_vmpdf(linspace(-pi,pi,100),0,max(rescaledVonMises));
        halfMax = (min(fitData) + max(fitData)) / 2;
        index1 = find(fitData >= halfMax, 1, 'first');
        index2 = find(fitData >= halfMax, 1, 'last');
        half_width{fly}(timepoint) = index2-index1 + 1;
    
    end

end

%Get and plot the median values
BumpMag = {};
HalfWidth = {};

for fly = 1:length(bump_mag)
    BumpMag{fly,1} = bump_mag{fly}(Data(fly).singleBar);
    BumpMag{fly,2} = bump_mag{fly}(Data(fly).ambiguous);
    BumpMag{fly,3} = bump_mag{fly}(Data(fly).darkness); 
    
    HalfWidth{fly,1} = half_width{fly}(Data(fly).singleBar);
    HalfWidth{fly,2} = half_width{fly}(Data(fly).ambiguous);
    HalfWidth{fly,3} = half_width{fly}(Data(fly).darkness); 
    
end

medianBumpMag = cellfun(@median,BumpMag);
medianHalfWidth = cellfun(@median,HalfWidth);

figure('Position',[100 300 1400 600]),
subplot(1,2,1)
plot(medianBumpMag','-o')
xlim([0 4]);
xticks([1 2 3]); 
xticklabels({'vertical bar' ,'horizontal bar', 'panels off'});
ylabel('Median bump magnitude')
legend({'fly1','fly2','fly3','fly4','fly5'});

subplot(1,2,2)
plot(medianHalfWidth','-o')
ylabel('Median bump half width');
xlim([0 4]); 
xticks([1 2 3]);
xticklabels({'vertical bar' ,'horizontal bar', 'panels off'});
legend({'fly1','fly2','fly3','fly4','fly5'});


bump_mag = bump_mag';
half_width = half_width';
bump_mag = [bump_mag{1,1},bump_mag{2,1},bump_mag{3,1},bump_mag{4,1},bump_mag{5,1}];
half_width = [half_width{1,1},half_width{2,1},half_width{3,1},half_width{4,1},half_width{5,1}];

save('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\ThreePatternExp\data\fitData.mat','bump_mag','half_width')