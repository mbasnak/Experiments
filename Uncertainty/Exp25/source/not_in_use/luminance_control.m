%Luminance control

%This code analyses data coming from the luminance control: a horizontal
%set of equidistant dots that change contrast every 100 sec and match the
%changes in luminance from the 'bar contrast' block
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
changeStim = find(abs(diff(data.fr_y_ds))>0.2);

%% Get contrasts

%Recover the y dimensions using the voltage
ydimensions = 6;
maxYVoltage = 10;
yDim = data.fr_y_ds*ydimensions/maxYVoltage;

figure,
plot(yDim)
xlim([0 length(yDim)]);
title('y dimension');
ylabel('Y dim'); xlabel('Time');

ContrastOrder = [6,5,4,1,2,3];

%% Heatmap

phase = rad2deg(data.phase(1:changeStim(6)));
smoothPhase = smooth(phase);

%Overlay the EPG phase to the heatmap to see if the EPG phase we're getting
%makes sense

% Remap phase values to 1-16 to fit the plot
phaseInPBcoordinates = wrapTo360(phase)*16/360;
% Add and subtract 4 glomeruli to get 2 phases to plot on the PB (I'm not
% sure if this is what I should be going)
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
imagesc(data.dff_matrix(:,1:changeStim(6)))
colormap(gray)
hold on
plot(phase1,'r')
plot(phase2,'r')
xlim([0 2000])
xlabel('Time (frames)');
ylabel('PB glomerulus');
legend('EPG phase');
title('Using the Fourier transform');


%Let's try with the PVA
pva = rad2deg(data.dff_pva_rad(1:changeStim(6)));
PVAInPBcoordinates = wrapTo360(pva)*16/360;
% Add and subtract 4 glomeruli to get 2 phases to plot on the PB (I'm not
% sure if this is what I should be going)
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
imagesc(data.dff_matrix(:,1:changeStim(6)))
title('Using the PVA');
colormap(gray)
hold on
plot(PVA1,'b')
plot(PVA2,'b')
xlabel('Time (frames)');
ylabel('PB glomerulus');
xlim([0 2000])
legend('EPG phase');

%Save figure
saveas(gcf,[path,'plots\',sid,'\EPGphaseCloseUp.png']);

%% Full PB heatmap with lines denoting change in stimulus, plus EPG phase and fly's position

heading = rad2deg(data.flyPosRad(1:changeStim(6))');
offset = smooth(wrapTo180(heading-phase));

%Set the values for which the heading and phase jump to NaNs such that they
%won't appear in the code
changePhase = abs([0;diff(smooth(phase))]);
phaseToPlot = smooth(-phase);
phaseToPlot(changePhase>40==1) = NaN;

changeHeading = abs([0;diff(smooth(heading))]);
headingToPlot = smooth(heading);
headingToPlot(changeHeading>40==1) = NaN;

%Plot
% Plot the heatmap of EPG activity
figure('Position',[300 200 1800 800]),
subplot(4,6,[1 6])
imagesc(data.dff_matrix(:,1:changeStim(6)))
colormap(gray)
hold on
%add the changes in stim
for change = 1:5
   line([changeStim(change) changeStim(change)], [0 17], 'LineWidth', 2, 'color', 'r'); 
end
ylabel('PB glomerulus');
set(gca,'XTickLabel',[]);
legend('Change in for_velulus');
title('EPG activity throughout the experiment');

% Plot the heading and the EPG phase
subplot(4,6,[7 12])
plot(data.time(1:changeStim(6)),headingToPlot,'color',[0.6 0.2 0.4],'LineWidth',1.5)
hold on
plot(data.time(1:changeStim(6)),phaseToPlot,'color',[0.2 0.6 0.4],'LineWidth',1.5)
%add the changes in for_vel
for change = 1:5
   line([data.time(changeStim(change)) data.time(changeStim(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
legend('Fly heading', 'EPG phase');
ylim([-180, 180]);
xlim([0,data.time(changeStim(6))]);
ylabel('PB glomerulus');
set(gca,'XTickLabel',[]);

% Plot the offset
changeOffset = abs([0;diff(offset)]);
offsetToPlot = offset;
offsetToPlot(changeOffset>100==1) = NaN;
subplot(4,6,[13 18])
plot(data.time(1:changeStim(6)),offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
for change = 1:5
   line([data.time(changeStim(change)) data.time(changeStim(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
ylim([-180 180]);
xlim([0,data.time(changeStim(6))]);
ylabel('Deg'); xlabel('Time (sec)');
legend('Offset');

%Polar histograms of offset
subplot(4,6,19)
polarhistogram(deg2rad(offset(1:changeStim(1))))
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
for Stim = 1:5
   subplot(4,6,19+Stim)
   polarhistogram(deg2rad(offset(changeStim(Stim):changeStim(Stim+1))))
   set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
end

%Save figure
saveas(gcf,[path,'plots\',sid,'\HeatmapAndOffset.png']);


%% Plot fly's heading

figure('Position',[100 100 1400 600]),

subplot(1,6,1)
polarhistogram(deg2rad(headingToPlot(1:changeStim(1))),20,'FaceColor','k')
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');

for contrast = 1:5
    subplot(1,6,contrast)
    polarhistogram(deg2rad(headingToPlot(changeStim(contrast):changeStim(contrast+1))),20,'FaceColor','b')
    set(gca,'ThetaZeroLocation','top',...
            'ThetaDir','counterclockwise');
end
suptitle('Fly heading distribution during the different blocks');

%Save figure
saveas(gcf,[path,'plots\',sid,'\HeadingDistribution.png']);


%% Compute bump mag and with at half max throughout

%Combine PB data in the prober two halves
leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];
combined_full_dff = (data.dff_matrix(leftPB,1:changeStim(6)) + data.dff_matrix(rightPB,1:changeStim(6)))/2;

%Bump mag as min-max
BumpMagnitude = max(combined_full_dff) - min(combined_full_dff);

figure('Position',[200 200 1600 400]),
subplot(1,4,[1 3]),
plot(data.time(1:changeStim(6)),BumpMagnitude, 'k')
hold on
medianBM = median(BumpMagnitude(1:changeStim(1)));
for contrast = 1:5
    line([data.time(changeStim(contrast)) data.time(changeStim(contrast))], [0 5], 'LineWidth', 2, 'color', 'r');
    medianBM = [medianBM, median(BumpMagnitude(changeStim(contrast):changeStim(contrast+1)))];
end
xlabel('Time (sec)'); ylabel('Bump magnitude (max-min)');

subplot(1,4,4)
plot(medianBM,'-ko','MarkerFaceColor','k')
%xlim([0 4]);
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
plot(data.time(1:changeStim(6)),bump_mag, 'k')
hold on
median_bump_mag = median(bump_mag(1:changeStim(1)));
median_half_width = median(half_width(1:changeStim(1)));
for contrast = 1:5
    line([data.time(changeStim(contrast)) data.time(changeStim(contrast))], [0 5], 'LineWidth', 2, 'color', 'r');
    median_bump_mag = [median_bump_mag, median(bump_mag(changeStim(contrast):changeStim(contrast+1)))];
    median_half_width = [median_half_width, median(half_width(changeStim(contrast):changeStim(contrast+1)))];
end
xlabel('Time (sec)'); ylabel('Bump magnitude');

subplot(2,4,4)
plot(median_bump_mag,'-ko','MarkerFaceColor','k')
title('Median bump magnitude');
ylabel('Bump mag');

subplot(2,4,[5 7]),
plot(data.time(1:changeStim(6)),half_width, 'k')
hold on
for contrast = 1:5
    line([data.time(changeStim(contrast)) data.time(changeStim(contrast))], [0 500], 'LineWidth', 2, 'color', 'r');
end
xlabel('Time (sec)'); ylabel('Bump width at half max');

subplot(2,4,8)
plot(median_half_width,'-ko','MarkerFaceColor','k')
title('Median half width');
ylabel('Bump width at half max');

%Save figure
saveas(gcf,[path,'plots\',sid,'\VonMisesBump.png']);

%% Bump magnitude and width at half max sorted by contrast

figure,
subplot(1,2,1)
plot(median_bump_mag(ContrastOrder),'-ko','MarkerFaceColor','k')
title('Median bump magnitude');
ylabel('Bump mag');
xlabel('Contrast level');
xlim([1 6]);

subplot(1,2,2)
plot(median_half_width(ContrastOrder),'-ko','MarkerFaceColor','k')
title('Median half width');
ylabel('Bump width at half max');
xlabel('Contrast level');
xlim([1 6]);

%Save figure
saveas(gcf,[path,'plots\',sid,'\VonMisesBumpOrdered.png']);


%% Bump mag and forward velocity

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time(1:length(bump_mag)),bump_mag, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot forward velocity in time) 
subplot(2,4,[5 7])
plot(data.time(1:length(bump_mag)),data.vel_for_ds(1:length(bump_mag)),'k')
xlabel('Time (sec)');
ylabel('Forward velocity (mm/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
scatter(data.vel_for_ds(1:length(bump_mag)),bump_mag(1:length(bump_mag)))
ylabel('Bump magnitude'); xlabel('Forward velocity');

%The relationship looks weird this way: it seems like we should bin it

%% Bump mag and forward velocity binning the velocity

%Let's always define a fixed number of vel bins
nbins = 10;
binWidth = ceil(max(data.vel_for_ds))/nbins;
velBins = [0:binWidth:ceil(max(data.vel_for_ds))]; %I creat a vector with my bins
velBins = [velBins,max(data.vel_for_ds)+1];

%getting binned medians 
for bin = 1:length(velBins)-1
    medianBin(bin) = median(bump_mag((abs(data.vel_for_ds(1:length(bump_mag))))>velBins(bin) & (abs(data.vel_for_ds(1:length(bump_mag))))<velBins(bin+1)));
end


velAxes = velBins - binWidth;
velAxes = velAxes(2:end);
velAxes(end) = velAxes(end-1)+binWidth;

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time(1:length(bump_mag)),bump_mag, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot forward velocity in time) 
subplot(2,4,[5 7])
plot(data.time(1:length(bump_mag)),data.vel_for_ds(1:length(bump_mag)),'k')
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
