%This code analyzes the first block of Experiment 25 (i.e., the block in
%which the fly gets 5 min darkness - 15 min blue bar in closed loop with
%bar jumps - 5 min darkness


%% Load data

clear all; close all;

[file,path] = uigetfile('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
load([path,file])

%determine the sid we're working with, to save the plots to specific folder
%later
sid = file(10:14);
%% Make directory to save plots

%Move to the analysis folder
cd(path)
%List the contents
contents = dir();
%if there isn't a 'plots' folder already, create one
if (contains([contents.name],'plots') == 0)
   mkdir(path,'plots'); 
end
%List the contents of the 'plots' folder
cd([path,'plots\'])
plotContents = dir();
%if there isn't a folder for this sid already, create one
if (contains([plotContents.name],sid) == 0)
   mkdir([path,'plots\'],sid); 
end

%% Identify changes in stim (darkness to stim)

%Plot the signal from the y dimension of the panels
figure,
subplot 121,
plot(data.fr_y_ds)
title('Panels y signal');
xlabel('Time (downsampled frames)');
ylabel('Voltage');
xlim([0 length(data.fr_y_ds)]);
subplot 122,
plot(abs(diff(data.fr_y_ds)))
title('Frames where the stimulus changes');
ylim([-1 6]);
xlim([0 length(data.fr_y_ds)]);
xlabel('Time (downsampled frames)');


%Obtain the frames where the stim changes using the derivative of the
%signal from the y dimension of the panels
changeStim = find(abs(diff(data.fr_y_ds))>1);
singleBar = [changeStim(1):changeStim(end)];
panelsOff = [1:changeStim(1),changeStim(end):length(data.fr_y_ds)];


%% Identify jumps (using the ypanels signal)

%This section needs to be run for all the flies except that of 20200828

%The jumps are the odd changes in y dim voltage
jumps = changeStim(2:end-1);

%Recover the y dimensions using the voltage
ydimensions = 24;
maxYVoltage = 10;
yDim = data.fr_y_ds*ydimensions/maxYVoltage;

figure,
plot(yDim)
xlim([0 length(yDim)]);
title('y dimension');
ylabel('Y dim'); xlabel('Time');

Jumps = [90,90,90,-90,-90,90,-90,-90];

% Fixing the data relative to the bar jumps
%heading = rad2deg(data.flyPosRad');
heading = rad2deg(data.heading');
correctedHeading = zeros(1,length(heading));
correctedHeading(1:jumps(1)) = heading(1:jumps(1));
for section = 1:length(Jumps)-1
    correctedHeading(jumps(section)+1:jumps(section+1)) = heading(jumps(section)+1:jumps(section+1)) + sum(Jumps(1:section)); 
end
correctedHeading(jumps(end):end) = heading(jumps(end):end) + sum(Jumps(1:end));
correctedHeading = wrapTo180(correctedHeading);

%Plot
figure,
subplot(2,1,1)
plot(data.time,heading)
for jump = 1:length(jumps)
    line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180], 'color', 'r','LineWidth',2)
end
ylim([-180 180]);

subplot(2,1,2)
plot(data.time,correctedHeading)
for jump = 1:length(jumps)
    line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180], 'color', 'r','LineWidth',2)
end
ylim([-180 180]);


%% Plot PB heatmap with overlayed EPG phase in very first part

%heading = rad2deg(data.flyPosRad');
phase = rad2deg(data.phase);
smoothPhase = smooth(phase);

%Overlay the EPG phase to the heatmap to see if the EPG phase we're getting
%makes sense

% Remap phase values to 1-16 to fit the plot
phaseInPBcoordinates = wrapTo360(phase)*16/360;
phase1 = wrapTo16(phaseInPBcoordinates-8);
phase2 = wrapTo16(phaseInPBcoordinates);
% Determine timepoints when the phase wraps around to remove those from the
%plot and avoid jumps in the lines
changePhase1 = abs([0,diff(phase1)]);
changePhase2 = abs([0,diff(phase2)]);
phase1(changePhase1>5==1) = NaN;
phase2(changePhase2>5==1) = NaN;

% Plot
figure('Position',[300 200 1800 800]),
subplot(2,1,1)
imagesc(data.dff_matrix)
colormap(gray)
hold on
%I have to do 16- phase 1 and phase 2 because glomerulus 1 is at the top of
%the heatmap
plot(16-phase1,'r')
plot(16-phase2,'r')
xlim([0 2000]);
ylabel('PB glomerulus');
[~, hobj, ~, ~] = legend('EPG phase');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title('Using the Fourier transform');


%Let's try with the PVA
pva = rad2deg(data.dff_pva');
PVAInPBcoordinates = wrapTo360(pva)*16/360;
PVA1 = wrapTo16(PVAInPBcoordinates-8);
PVA2 = wrapTo16(PVAInPBcoordinates);
% Determine timepoints when the phase wraps around to remove those from the
%plot and avoid jumps in the lines
changePVA1 = abs([0,diff(PVA1)]);
changePVA2 = abs([0,diff(PVA2)]);
PVA1(changePVA1>5==1) = NaN;
PVA2(changePVA2>5==1) = NaN;

% Plot
subplot(2,1,2)
imagesc(data.dff_matrix)
title('Using the PVA');
colormap(gray)
hold on
plot(PVA1,'b')
plot(PVA2,'b')
xlabel('Time (frames)');
ylabel('PB glomerulus');
xlim([0 2000])
[~, hobj, ~, ~] = legend('EPG phase');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

% Save figure
saveas(gcf,[path,'plots\',sid,'\EPGphaseCloseUp.png']);

%% Plot PB heatmap with overlayed EPG phase in full trial

% Plot
figure('Position',[300 200 1800 800]),
subplot(2,1,1)
imagesc(data.dff_matrix)
colormap(gray)
hold on
plot(-phase1,'r')
plot(-phase2,'r')
ylabel('PB glomerulus');
[~, hobj, ~, ~] = legend('EPG phase');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);
title('Using the Fourier transform');

subplot(2,1,2)
imagesc(data.dff_matrix)
title('Using the PVA');
colormap(gray)
hold on
plot(PVA1,'b')
plot(PVA2,'b')
xlabel('Time (frames)');
ylabel('PB glomerulus');
[~, hobj, ~, ~] = legend('EPG phase');
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',1.5);

% Save figure
saveas(gcf,[path,'plots\',sid,'\EPGphase.png']);


%% Full PB heatmap with lines denoting change in stimulus, plus EPG phase and fly's position

offset = wrapTo180(rad2deg(circ_dist(deg2rad(heading),deg2rad(phase))));

%Set the values for which the heading and phase jump to NaNs such that they
%won't appear in the code
changePhase = abs([0;diff(smooth(phase))]);
phaseToPlot = smooth(phase);
phaseToPlot(changePhase>40==1) = NaN;

changeHeading = abs([0;diff(smooth(heading))]);
headingToPlot = smooth(heading);
headingToPlot(changeHeading>40==1) = NaN;

%Plot
% Plot the heatmap of EPG activity
figure('Position',[300 200 1800 800]),
subplot(4,3,[1 3])
imagesc(data.dff_matrix)
colormap(gray)
hold on
%add the changes in stim
line([changeStim(1) changeStim(1)], [0 17], 'LineWidth', 2, 'color', 'r');
line([changeStim(end) changeStim(end)], [0 17], 'LineWidth', 2, 'color', 'r');
ylabel('PB glomerulus');
set(gca,'XTickLabel',[]);
legend('Change in stimulus');
title('EPG activity throughout the experiment');

% Plot the heading and the EPG phase
subplot(4,3,[4 6])
plot(data.time,headingToPlot,'color',[0.6 0.2 0.4],'LineWidth',1.5)
hold on
plot(data.time,phaseToPlot,'color',[0.2 0.6 0.4],'LineWidth',1.5)
line([data.time(changeStim(1)) data.time(changeStim(1))], [-180 180], 'LineWidth', 2, 'color', 'r');
line([data.time(changeStim(end)) data.time(changeStim(end))], [-180 180], 'LineWidth', 2, 'color', 'r');
for jump = 1:length(jumps)
    line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180], 'LineWidth', 2, 'color', 'k');
end
legend('Fly heading', 'EPG phase');
ylim([-180, 180]);
xlim([0,data.time(end)]);
ylabel('PB glomerulus');
set(gca,'XTickLabel',[]);

% Plot the offset
changeOffset = abs([0,diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
subplot(4,3,[7 9])
plot(data.time(1:changeStim(1)),offsetToPlot(1:changeStim(1)),'LineWidth',1.5,'color','k')
hold on
plot(data.time(changeStim(1):changeStim(end)),offsetToPlot(changeStim(1):changeStim(end)),'LineWidth',1.5,'color','b')
plot(data.time(changeStim(end):end),offsetToPlot(changeStim(end):end),'LineWidth',1.5,'color','k')
line([data.time(changeStim(1)) data.time(changeStim(1))], [-180 180], 'LineWidth', 2, 'color', 'r');
line([data.time(changeStim(end)) data.time(changeStim(end))], [-180 180], 'LineWidth', 2, 'color', 'r');
for jump = 1:length(jumps)
    line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180], 'LineWidth', 2, 'color', 'k');
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
legend('Offset');

% Plot the polar histograms of offset for the three bouts
subplot(4,3,10)
polarhistogram(deg2rad(offset(1:changeStim(1))),20,'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
subplot(4,3,11)
polarhistogram(deg2rad(offset(changeStim(1):changeStim(2))),20,'FaceColor','b')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
subplot(4,3,12)
polarhistogram(deg2rad(offset(changeStim(2):end)),20,'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');

% Save figure
saveas(gcf,[path,'plots\',sid,'\HeatmapAndOffset.png']);

%% Get times of jumps

figure('Position',[200 200 1600 500]),
plot(data.time, headingToPlot)
hold on
for i = [86:90:data.time(end)]
    line([i i],[-180 180],'LineWidth',2,'color','k')    
end
ylim([-180 180]); ylabel('Heading (deg)');
title('Heading data with bar jumps');
legend('Heading','bar jumps');

%Store the bar jump times (in sec and frames)
barJumpSec = [86:90:data.time(end)];
barJumpFrames = [];
for i = 1:length(barJumpSec)
    [M,I] = min(abs(data.time-barJumpSec(i)));
    barJumpFrames(i) = I;
end

%Keep only the jumps that happen when the panels are on
barJumpSec = barJumpSec(4:13);
barJumpFrames = barJumpFrames(4:13);

%% Get angular velocity using heading data

figure('Position',[100 100 1600 800]),
subplot(4,1,1)
plot(data.time,heading,'color',[0.3 0.5 0.1])
ylim([-180 180]);
set(gca,'XTickLabel',[]);
ylabel('Deg');
title('Heading');

subplot(4,1,2)
plot(data.time,smooth(heading),'color',[0.8 0.2 0.1])
ylim([-180 180]);
set(gca,'XTickLabel',[]);
ylabel('Deg');
title('Smoothed heading');

subplot(4,1,3)
plot(data.time,[0;diff(smooth(heading))],'color',[0.1 0.5 0.6])
ylim([-100 100]);
set(gca,'XTickLabel',[]);
ylabel('Deg/s');
title('Angular velocity');

subplot(4,1,4)
angVel = [0;smooth(diff(smooth(heading)))];
plot(data.time,angVel,'color',[0.7 0.2 0.6])
ylim([-100 100]);
xlabel('Time (sec)');
ylabel('Deg/s');
title('Smoothed angular velocity');

% Save figure
saveas(gcf,[path,'plots\',sid,'\AngularVel.png']);

%% Plot angular velocity and heatmap around jump times

jumpOrder = [-90, 90, 90, -90, -90, 90, -90, 90, 90, 90];

for jump = 1:length(barJumpSec)
    figure,
    subplot(2,1,1)
    imagesc(data.dff_matrix(:,barJumpFrames(jump)-50:barJumpFrames(jump)+50)) 
    colormap(gray)
    hold on
    line([50 50] , [0 17], 'LineWidth', 2, 'color', 'r');
    xlim([0 100]);
    set(gca,'XTickLabel',[]);
    ylabel('PB glomerulus'); 
    title(['Jump of ',num2str(jumpOrder(jump)),' deg']);
    
    subplot(2,1,2)
    plot(data.time((barJumpFrames(jump)-50:barJumpFrames(jump)+50)),angVel(barJumpFrames(jump)-50:barJumpFrames(jump)+50)) 
    minVel = min(angVel(barJumpFrames(jump)-50:barJumpFrames(jump)+50));
    maxVel = max(angVel(barJumpFrames(jump)-50:barJumpFrames(jump)+50));
    hold on
    line([data.time(barJumpFrames(jump)) data.time(barJumpFrames(jump))] , [minVel-1 maxVel+1], 'LineWidth', 2, 'color', 'r');
    line([data.time(barJumpFrames(jump)-50) data.time(barJumpFrames(jump)+50)], [0 0], 'color', 'k')
    xlim([data.time(barJumpFrames(jump)-50) data.time(barJumpFrames(jump)+50)]);
    ylim([minVel-1 maxVel+1]);
    ylabel('Angular velocity (deg/s)'); xlabel('Time (sec)');
    
    % Save figure
    saveas(gcf,[path,'plots\',sid,'\AroundJump',num2str(jump),'.png']);
end

%% Mean angular vecocity around jumps of similar magnitude

%Store together the angular velocity around all the jumps of similar
%magnitude
angVel90Jumps = [];
angVelNeg90Jumps = [];

for jump = 1:length(barJumpSec)
   if (jumpOrder(jump) == 90) 
       angVel90Jumps = [angVel90Jumps,angVel(barJumpFrames(jump)-50:barJumpFrames(jump)+50)];
   else
       angVelNeg90Jumps = [angVelNeg90Jumps,angVel(barJumpFrames(jump)-50:barJumpFrames(jump)+50)];       
   end   
end

%Calculate the mean angular velocity around each jump type
mean90Jumps = mean(angVel90Jumps,2);
meanNeg90Jumps = mean(angVelNeg90Jumps,2);

%Plot
figure,
plot(mean90Jumps)
hold on
plot(meanNeg90Jumps)
line([50 50], [-20 20], 'color', 'r','LineWidth',2)
ylabel('Deg/s');
xlim([0 100]);
title('Mean angular velocity around jumps');
legend('90 deg','-90 deg');

%Save figure
saveas(gcf,[path,'plots\',sid,'\MeanArounJumpAngVel.png']);

%% Plot fly's heading

figure('Position',[100 100 1400 600]),

subplot(1,3,1)
polarhistogram(deg2rad(headingToPlot(1:changeStim(1))),20,'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    
subplot(1,3,2)
polarhistogram(deg2rad(headingToPlot(changeStim(1):changeStim(2))),20,'FaceColor','b')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    
subplot(1,3,3)
polarhistogram(deg2rad(headingToPlot(changeStim(2):end)),20,'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
suptitle('Fly heading distribution during the different blocks');

%Save figure
saveas(gcf,[path,'plots\',sid,'\HeadingDistribution.png']);

%% Compute bump mag and with at half max throughout

%Combine PB data in the prober two halves
leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];
combined_full_dff = (data.dff_matrix(leftPB,:) + data.dff_matrix(rightPB,:))/2;

%Bump mag as min-max
BumpMagnitude = max(combined_full_dff) - min(combined_full_dff);

figure('Position',[200 200 1600 400]),
subplot(1,4,[1 3]),
plot(data.time,BumpMagnitude, 'k')
hold on
line([data.time(changeStim(1)) data.time(changeStim(1))], [0 5], 'LineWidth', 2, 'color', 'r');
line([data.time(changeStim(end)) data.time(changeStim(end))], [0 5], 'LineWidth', 2, 'color', 'r');
xlabel('Time (sec)'); ylabel('Bump magnitude (max-min)');

medianBM = [median(BumpMagnitude(1:changeStim(1))),median(BumpMagnitude(changeStim(1):changeStim(end))),median(BumpMagnitude(changeStim(end):end))];

subplot(1,4,4)
plot(medianBM,'-ko','MarkerFaceColor','k')
xlim([0 4]);
title('Median bump magnitude');
ylabel('Bump mag (max-min)');

%Save figure
saveas(gcf,[path,'plots\',sid,'\BMMinMax.png']);

%% Using Von Mises fit

half_width = zeros(1,length(combined_full_dff));
bump_mag = zeros(1,length(combined_full_dff));

%Apply von Mises fit and get bump mag and width at half max from there
for timepoint = 1:length(combined_full_dff)
    extendedData{timepoint} = interp1([1:8],combined_full_dff(:,timepoint),linspace(1,8,1000));
    angles = linspace(0,2*pi,length(extendedData{timepoint}));
    [vonMises{timepoint}, rescaledVonMises{timepoint}] = fitTuningCurveToVonMises(extendedData{timepoint}, angles);
    %get bump mag as max - min from fit
    bump_mag(timepoint) = max(rescaledVonMises{timepoint})-min(rescaledVonMises{timepoint});
    fitData = circ_vmpdf(linspace(-pi,pi,1000),0,max(rescaledVonMises{timepoint}));
    halfMax = (min(fitData) + max(fitData)) / 2;
    index1 = find(fitData >= halfMax, 1, 'first');
    index2 = find(fitData >= halfMax, 1, 'last');
    half_width(timepoint) = index2-index1 + 1;  
end


% Plot
figure('Position',[200 200 1600 400]),
subplot(2,4,[1 3]),
plot(data.time,bump_mag, 'k')
hold on
line([data.time(changeStim(1)) data.time(changeStim(1))], [0 5], 'LineWidth', 2, 'color', 'r');
line([data.time(changeStim(end)) data.time(changeStim(end))], [0 5], 'LineWidth', 2, 'color', 'r');
xlabel('Time (sec)'); ylabel('Bump magnitude (max-min)');

subplot(2,4,4)
median_bump_mag = [median(bump_mag(1:changeStim(1))),median(bump_mag(changeStim(1):changeStim(end))),median(bump_mag(changeStim(end):end))];
plot(median_bump_mag,'-ko','MarkerFaceColor','k')
xlim([0 4]);
title('Median bump magnitude');
ylabel('Bump mag (max-min)');

subplot(2,4,[5 7]),
plot(data.time,half_width, 'k')
hold on
line([data.time(changeStim(1)) data.time(changeStim(1))], [0 500], 'LineWidth', 2, 'color', 'r');
line([data.time(changeStim(end)) data.time(changeStim(end))], [0 500], 'LineWidth', 2, 'color', 'r');
xlabel('Time (sec)'); ylabel('Bump magnitude (max-min)');

subplot(2,4,8)
median_half_width = [median(half_width(1:changeStim(1))),median(half_width(changeStim(1):changeStim(end))),median(half_width(changeStim(end):end))];
plot(median_half_width,'-ko','MarkerFaceColor','k')
xlim([0 4]);
title('Median half width');
ylabel('Bump width at half max (max-min)');

%Save figure
saveas(gcf,[path,'plots\',sid,'\VonMisesBump.png']);

%% Flies velocity distribution to see how 'well' they were walking

figure('Position',[200 200 1600 400]),
subplot(2,3,[1 2])
plot(smoothed.xVel,'color',[0.7 0.2 0.6])
ylabel('Forward vel (mm/s)');
ylim([0 8]); xlim([0 length(smoothed.xVel)]);
subplot(2,3,3)
histogram(smoothed.xVel,'FaceColor',[0.7 0.2 0.6])
title('Foward vel distribution');

subplot(2,3,[4 5])
plot(smoothed.angularVel)
xlabel('Time'); ylabel('Angular vel (deg/s)');
ylim([-100 100]); xlim([0 length(smoothed.xVel)]);
subplot(2,3,6)
histogram(smoothed.angularVel)
title('Angular vel distribution');

%Save figure
saveas(gcf,[path,'plots\',sid,'\VelocityDistributions.png']);


%% Bump mag and forward velocity binning the velocity

%Let's always define a fixed number of vel bins
nbins = 10;
maxBin = min(ceil(max(data.vel_for_ds)),15); %I'm adding this line because very high forward velocities are probably more artifact than reality
binWidth = maxBin/nbins;
velBins = [0:binWidth:maxBin]; %I creat a vector with my bins
velBins = [velBins,maxBin+1];

%getting binned medians 
for bin = 1:length(velBins)-1
    medianBin(bin) = median(bump_mag((abs(data.vel_for_ds))>velBins(bin) & (abs(data.vel_for_ds))<velBins(bin+1)));
end


velAxes = velBins - binWidth;
velAxes = velAxes(2:end);
velAxes(end) = velAxes(end-1)+binWidth;

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time,bump_mag, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot forward velocity in time) 
subplot(2,4,[5 7])
plot(data.time,data.vel_for_ds,'k')
xlabel('Time (sec)');
ylabel('Forward velocity (mm/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(velAxes,medianBin,'-ko')
ylabel('Median bump magnitude'); xlabel('Forward velocity (mm/s)');
ylim([0 (max(medianBin)+0.5)]);

%Save figure
saveas(gcf,[path,'plots\',sid,'\BumpVsFwdVel.png']);

%% Relationship with the angular speed

%Resample my smoothed angular velocity to have the same number of
%datapoints as the bump mag
resampledAngVel = resample(smoothed.angularVel,length(bump_mag),length(smoothed.angularVel));


%Plot the relationship between bump and angular speed with my smoothed
%velocity
nbins = 10;
angBinWidth = ceil(max(abs(resampledAngVel)))/nbins;

angVelBins2 = [0:angBinWidth:ceil(max(abs(resampledAngVel)))]; %I create a vector with my bins
angVelBins2 = [angVelBins2,max(abs(resampledAngVel))+10];

%getting binned medians 
for bin = 1:length(angVelBins2)-1
    medianBinAng2(bin) = median(bump_mag((abs(resampledAngVel))>angVelBins2(bin) & (abs(resampledAngVel))<angVelBins2(bin+1)));
end

angVelAxes2 = angVelBins2 - angBinWidth;
angVelAxes2 = angVelAxes2(2:end);
angVelAxes2(end) = angVelAxes2(end-1)+angBinWidth;

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time(1:length(bump_mag)),bump_mag, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot angular speed in time) 
subplot(2,4,[5 7])
plot(data.time(1:length(bump_mag)),abs(resampledAngVel(1:length(bump_mag))),'k')
xlabel('Time (sec)');
ylabel('Angular speed (deg/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(angVelAxes2,medianBinAng2,'-ko')
ylabel('Median bump magnitude'); xlabel('Angular speed (deg/s)');
ylim([0 (max(medianBinAng2)+0.5)]);

%Save figure
saveas(gcf,[path,'plots\',sid,'\BumpVsAngSpeedSmoothed.png']);
