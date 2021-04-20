%Bar contrast block


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
changeContrast = changeStim(2:2:end); %every other jump is a change in contrast

%% Identify jumps (using the ypanels signal)

%The jumps are the odd changes in y dim voltage
jumps = changeStim(1:2:end);

%Recover the y dimensions using the voltage
ydimensions = 24;
maxYVoltage = 10;
yDim = data.fr_y_ds*ydimensions/maxYVoltage;

figure,
plot(yDim)
xlim([0 length(yDim)]);
title('y dimension');
ylabel('Y dim'); xlabel('Time');

%Use the pos function number to define the direction of the jumps
%load the run object
run_obj_dir = dir([path(1:end-9),'ball\runobj']);
for fileName = 1:length(run_obj_dir)
    if (contains(run_obj_dir(fileName).name,file(10:14)) == 1)
        run_obj_file_name = run_obj_dir(fileName).name;
        run_obj_folder_name = run_obj_dir(fileName).folder;
    end
end
load([run_obj_folder_name,'\',run_obj_file_name])
pos_function = run_obj.function_number;

%Assign the bar jumps and the contrast changes based on the pos function
%The vector 'ContrastOrder' reflects the order in which the intensities
%appear such that ContrastOrder(1) shows in which position of the
%experiment the bar of intensity 1 appears.
%The vector 'Intensities' reflects the intensity that appears in each
%position, such that Intensities(1) returns the intensity level of the bar
%that went first in the experiment
if pos_function == 190   
    Jumps = [90, -90, 90, 90, -90, -90];
    ContrastOrder = [2,6,5,3,1,4]; 
    Intensities = [5,1,4,6,3,2];
elseif pos_function == 185  %This is for the fly from 8/28 with the previous protocol
    Jumps = [-90, 90, -90, 90, 90, 180];
    ContrastOrder = [6,4,3,1,5,2];
    Intensities = [4,6,3,2,5,1];
elseif pos_function == 191   
    Jumps = [90, -90, -90, -90, 90, 90]; 
    ContrastOrder = [5,3,6,2,1,4];
    Intensities = [5,4,2,6,1,3];
elseif pos_function == 192 
    Jumps = [-90, -90, -90, 90, 90, 90];
    ContrastOrder = [1,4,6,3,2,5];
    Intensities = [1,5,4,2,6,3];
else
    Jumps = [90, -90, 90, 90, -90, -90];
    ContrastOrder = [2,5,1,4,3,6];
    Intensities = [3,1,5,4,2,6];
end

% Fixing the data relative to the bar jumps
heading = rad2deg(data.flyPosRad');
correctedHeading = zeros(1,length(heading));
correctedHeading(1:jumps(1)) = heading(1:jumps(1));
for section = 1:length(Jumps)-1
    correctedHeading(jumps(section)+1:jumps(section+1)) = heading(jumps(section)+1:jumps(section+1)) + sum(Jumps(1:section)); 
end
correctedHeading(jumps(end):end) = heading(jumps(end):end) + sum(Jumps(1:end));
correctedHeading = wrapTo180(correctedHeading);

%Plot
figure('Position',[100 100 1600 600]),
subplot(2,1,1)
plot(data.time,heading)
for jump = 1:length(jumps)
    line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180], 'color', 'r','LineWidth',2)
end
ylim([-180 180]);
ylabel('Deg'); title('Heading');
legend('Fly heading', 'Bar jumps');

subplot(2,1,2)
plot(data.time,correctedHeading)
for jump = 1:length(jumps)
    line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180], 'color', 'r','LineWidth',2)
end
ylim([-180 180]); xlabel('Time (sec)');
ylabel('Deg'); title('Corrected heading');
legend('Fly heading', 'Bar jumps');


%% Show a close-up of the changes in heading around the jump times

for jump = 1:length(jumps)
   figure,
   subplot(2,1,1)
   plot(data.time(jumps(jump)-50:jumps(jump)+50),heading(jumps(jump)-50:jumps(jump)+50))
   line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180],'color','r','LineWidth',2)
   ylim([-180 180]); xlim([data.time(jumps(jump)-50) data.time(jumps(jump)+50)])
   xlabel('Time (sec)'); ylabel('Deg');
   title('Heading');
   
   subplot(2,1,2)
   plot(data.time(jumps(jump)-50:jumps(jump)+50),correctedHeading(jumps(jump)-50:jumps(jump)+50))
   line([data.time(jumps(jump)) data.time(jumps(jump))], [-180 180],'color','r','LineWidth',2)
   ylim([-180 180]); xlim([data.time(jumps(jump)-50) data.time(jumps(jump)+50)])
   xlabel('Time (sec)'); ylabel('Deg');
   title('Corrected heading');
   suptitle(['Jump ',num2str(Jumps(jump)),' deg']);
end

%% Heatmap with bar jumps and stim changes

phase = rad2deg(data.phase);
smoothPhase = smooth(phase);

%Overlay the EPG phase to the heatmap to see if the EPG phase we're getting
%makes sense

% Remap phase values to 1-16 to fit the plot
phaseInPBcoordinates = wrapTo360(phase)*16/360;
% Add and subtract 4 glomeruli to get 2 phases to plot on the PB (I'm not
% sure if this is what I should be going)
phase1 = wrapTo16(phaseInPBcoordinates-5);
phase2 = wrapTo16(phaseInPBcoordinates+3);
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
plot(phase1,'r')
plot(phase2,'r')
xlim([0 2000])
xlabel('Time (frames)');
ylabel('PB glomerulus');
legend('EPG phase');
title('Using the Fourier transform');


%Let's try with the PVA
pva = rad2deg(data.dff_pva_rad);
PVAInPBcoordinates = wrapTo360(pva)*16/360;
% Add and subtract 4 glomeruli to get 2 phases to plot on the PB (I'm not
% sure if this is what I should be going)
PVA1 = wrapTo16(PVAInPBcoordinates-6);
PVA2 = wrapTo16(PVAInPBcoordinates+2);
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
legend('EPG phase');

%Save figure
saveas(gcf,[path,'plots\',sid,'\EPGphaseCloseUp.png']);

%% Full heatmap

% Plot
figure('Position',[300 200 1800 800]),
subplot(2,1,1)
imagesc(data.dff_matrix)
colormap(gray)
hold on
plot(phase1,'r')
plot(phase2,'r')
xlabel('Time (frames)');
ylabel('PB glomerulus');
legend('EPG phase');
title('Using the Fourier transform');

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
legend('EPG phase');

%Save figure
saveas(gcf,[path,'plots\',sid,'\EPGphase.png']);


%% Offset

%Set the values for which the heading and phase jump to NaNs such that they
%won't appear in the code
if ~contains(path,'20200904_60D05_7f');
    offset = wrapTo180(-correctedHeading-phase);
    changeHeading = abs([0;diff(smooth(correctedHeading))]);
    headingToPlot = smooth(correctedHeading);
else
    offset = wrapTo180(correctedHeading-phase);
    changeHeading = abs([0;diff(smooth(-correctedHeading))]);
    headingToPlot = smooth(-correctedHeading);
end
headingToPlot(changeHeading>60==1) = NaN;

changePhase = abs([0;diff(smooth(phase))]);
phaseToPlot = smooth(-phase);
phaseToPlot(changePhase>60==1) = NaN;



%Plot
% Plot the heatmap of EPG activity
figure('Position',[300 200 1800 800]),
subplot(4,6,[1 6])
imagesc(data.dff_matrix)
colormap(gray)
hold on
%add the changes in stim
for change = 1:length(changeContrast)
   line([changeContrast(change) changeContrast(change)], [0 17], 'LineWidth', 2, 'color', 'r'); 
end
ylabel('PB glomerulus');
set(gca,'XTickLabel',[]);
legend('Change in stimulus');
title('EPG activity throughout the experiment');

% Plot the heading and the EPG phase
subplot(4,6,[7 12])
plot(data.time,headingToPlot,'color',[0.6 0.2 0.4],'LineWidth',1.5)
hold on
plot(data.time,phaseToPlot,'color',[0.2 0.6 0.4],'LineWidth',1.5)
%add the changes in stim
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
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
subplot(4,6,[13 18])
plot(data.time,offsetToPlot,'LineWidth',1.5,'color','k')
%add the changes in stim
jumpTimes = data.time(jumps);
for change = 1:length(changeContrast)
   line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r'); 
end
%add the jumps
for jump = 1:length(jumps)
   line([jumpTimes(jump) jumpTimes(jump)], [-180 180], 'color', 'b','LineStyle','--', 'LineWidth', 2); 
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
legend('Offset');

%Polar histograms of offset
%Color histograms acording to the intensity level of the bar, using the
%vector Intensities
color_gradient = {[0,0,0],[0.5 0 0], [1.0 0.0156 0],[1.0 0.2094 0],[1.0 0.7 0],[1.0 1 0]};
subplot(4,6,19)
polarhistogram(deg2rad(offset(1:changeContrast(1))),'FaceColor',color_gradient{Intensities(1)})
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title({'Offset distribution',['Contrast ',num2str(Intensities(1))]});

for contrast = 1:length(changeContrast)-1
   subplot(4,6,19+contrast)
   polarhistogram(deg2rad(offset(changeContrast(contrast):changeContrast(contrast+1))),'FaceColor',color_gradient{Intensities(contrast+1)})
   set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
   title({'Offset distribution',['Contrast ',num2str(Intensities(contrast+1))]});
end

%Save figure
saveas(gcf,[path,'plots\',sid,'\HeatmapAndOffset.png']);


%% Plot the EPG phase around the jumps 

for jump = 1:length(jumps)
   figure,
   %Plot the PB heatmap around the jumps
   subplot(3,1,1)
   imagesc(data.dff_matrix(:,jumps(jump)-100:jumps(jump)+100))
   colormap(gray)
   hold on
   line([100 100],[0 17],'color','r','LineWidth',1.5);
   set(gca,'XTickLabel',[]);
   
   %Plot the EPG phase
   subplot(3,1,2)
   plot(data.time(jumps(jump)-100:jumps(jump)+100),-phase(jumps(jump)-100:jumps(jump)+100),'k','LineWidth',1.5) 
   hold on
   line([data.time(jumps(jump)) data.time(jumps(jump))],[-180 180],'color','r','LineWidth',1.5);
   ylim([-180 180]);
   xlim([data.time(jumps(jump)-100) data.time(jumps(jump)+100)]);
   legend('EPG phase', 'Bar jump');
   ylabel('Deg');
   title('EPG phase');
   
   %Plot the fly's angular velocity
   subplot(3,1,3)
   %Get the flie's angular velocity
   %angVel = [0;smooth(diff(smooth(heading)))];  %I'm calculating this with
   %the uncorrected heading because I think if I use the corrected one I will get spurious jumps
   resampledAngVel = resample(smoothed.angularVel,length(data.time),length(smoothed.angularVel));
   plot(data.time(jumps(jump)-100:jumps(jump)+100),resampledAngVel(jumps(jump)-100:jumps(jump)+100),'k','LineWidth',1.5) 
   hold on
   line([data.time(jumps(jump)) data.time(jumps(jump))],[-40 40],'color','r','LineWidth',1.5);
   ylim([-40 40]);
   xlim([data.time(jumps(jump)-100) data.time(jumps(jump)+100)]);
   legend('Fly angular velocity', 'Bar jump');
   ylabel('Angular velocity (deg/s)'); xlabel('Time (sec)');
   title('Fly angular velocity');
   
   suptitle(['Jump ',num2str(Jumps(jump)),' deg'])
   
   %Save
   saveas(gcf,[path,'plots\',sid,'\PhaseAndAngVelAJ_',num2str(jump),'.png']);

end


%% Compare offset before and after the jump

%I will use the methods in Seelig 2015 (Figure 3b, see methods): they computed the offset in the 100
%frames previous to the jump as the mean_circ_dist between PVA and
%heading. They compared that to the offset in the 100 frames previous to a
%second bar jump, which occurred 50 sec later! I will then take 100 frames
%preceding the jumps, 100 frames 50 sec after the jump, and compare both

%Initialize cell array
offset_change = {};
figure('Position',[100, 100, 1400, 800]),
for jump = 1:length(jumps)
    %Pre-jump
    offset_change{jump}(1) = mean(offset(jumps(jump)-105:jumps(jump)-5));
    pre_jump_offset(jump) = offset_change{jump}(1);
    
    %Post-jump
    offset_change{jump}(2) = mean(offset(jumps(jump)+505:jumps(jump)+555));
    post_jump_offset(jump) = offset_change{jump}(2);
    
    %Plot
    subplot(2,2,1)
    plot(offset_change{jump},'-o','color',color_gradient{Intensities(jump)},'MarkerFaceColor',color_gradient{Intensities(jump)})
    hold on
    ylim([-180,180]); xlim([0,3]);
    xticks([1 2]);
    xticklabels({'Offset before the jump', 'Offset after the jump'});
    ylabel('Degrees');
    
end

%Create table to fit regression
offsetChange = array2table([pre_jump_offset',post_jump_offset'], 'VariableNames', {'pre_jump_offset', 'post_jump_offset'});
%Fit a model to all contrasts excluding darkness 
offsetChangeModel = fitlm(offsetChange(Intensities~=1,:));

%Add boxplots condensing changes
subplot(2,2,2)
boxplot(offsetChange.post_jump_offset-offsetChange.pre_jump_offset)
%get ordered pre and post jump data
ordered_pre_jump_offset = offsetChange.pre_jump_offset(ContrastOrder);
ordered_post_jump_offset = offsetChange.post_jump_offset(ContrastOrder);
title('Change in offset');
ylim([-180 180]);
yline(0,'--r');
set(gca,'xticklabel',{[]});

%Plot changes in offset vs contrast of stim
subplot(2,2,3)
orderedOffsetChange = offsetChange.post_jump_offset-offsetChange.pre_jump_offset;
orderedOffsetChange = orderedOffsetChange(ContrastOrder);
plot(orderedOffsetChange,'o')
ylim([-180 180]);
yline(0,'--r');
xlim([0 7]);
xticks([1:6]);
xlabel('Contrast level');
ylabel('Offset change around bar jump (deg)');

%Plot offset before vs after jump, along with regression fit to it.
subplot(2,2,4)
scatter(pre_jump_offset(Intensities~=1),post_jump_offset(Intensities~=1))
xlim([-180 180]); ylim([-180 180]);
xlabel('Offset before cue jump'); ylabel('Offset after cue jump');
%add regression line
hold on
mock_offsets = -180:10:180;
regression_line = offsetChangeModel.Coefficients.Estimate(1) + offsetChangeModel.Coefficients.Estimate(2)*mock_offsets;
plot(mock_offsets,regression_line,'r')
%Add R and p-value as text
text(120,120,['R =',num2str(offsetChangeModel.Rsquared.Ordinary)]);
text(120,140,['p-val =',num2str(offsetChangeModel.Coefficients.pValue(2))]);

%Save
saveas(gcf,[path,'plots\',sid,'\OffsetChangesAJ.png']);
   
%% Set the contrast limits

%First find the start of each contrast
contrastStart = [];
if Intensities(1) == 1
    contrastStart(1) = 1;
else
    contrastStart(1) = changeContrast(find(Intensities==1)-1);
end
for contrast = 1:length(changeContrast)-1
    if find(Intensities==contrast+1) == 1
        contrastStart(contrast+1) = 1;
    else
        contrastStart(contrast+1) = changeContrast(find(Intensities==contrast+1)-1);
    end
end

%Then find the end
contrastEnd = [];
for contrast = 1:length(contrastStart)
    if contrastStart(contrast)~= 1
        contrastEnd(contrast) = changeContrast(find(changeContrast == contrastStart(contrast)) + 1);
    else
        contrastEnd(contrast) = changeContrast(1);
    end
end

%Then store both values
contrastLimits = {};
for contrast = 1:6
    contrastLimits{contrast} = [contrastStart(contrast),contrastEnd(contrast)];
end

%%  Offset distribution

%use circular sd of offset to determine how stable it was and see of that
%relates to the contrast

offset_var = [];
for contrast = 1:length(changeContrast)
    [s(contrast) offset_var(contrast)] = circ_std(deg2rad(offset(contrastLimits{contrast}(1):contrastLimits{contrast}(2))),[],[],2);
end

figure('Position',[200,200,1200,400]),
for contrast = 1:6
    subplot(4,6,contrast)
    polarhistogram(deg2rad(offset(contrastLimits{contrast}(1):contrastLimits{contrast}(2))),'FaceColor',color_gradient{contrast})
    set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    Ax = gca; % current axes
    Ax.RTickLabel = []; 
    Ax.ThetaTick = [0,90,180,270];
    Ax.ThetaTickLabel = [0,90,180,270];
end
subplot(4,6,[7 24])
plot(offset_var,'-ko','MarkerFaceColor','k')
xlabel('Contrast level');
ylabel('Offset variation (circular std)');
xlim([0.5 6.5]);
xticks(1:6);

%Save figure
saveas(gcf,[path,'plots\',sid,'\OffsetVariation.png']);

%% Flies velocity distribution to see how 'well' they were walking

figure('Position',[200 200 1600 400]),
subplot(2,3,[1 2])
plot(smoothed.xVel,'color',[0.7 0.2 0.6])
ylabel('Forward vel (mm/s)');
ylim([0 max(smoothed.xVel)+1]); xlim([0 length(smoothed.xVel)]);
subplot(2,3,3)
histogram(smoothed.xVel,'FaceColor',[0.7 0.2 0.6])
title('Foward vel distribution');

subplot(2,3,[4 5])
plot(smoothed.angularVel)
xlabel('Time'); ylabel('Angular vel (deg/s)');
ylim([-max(abs(smoothed.angularVel))-10 max(abs(smoothed.angularVel))+10]); xlim([0 length(smoothed.xVel)]);
subplot(2,3,6)
histogram(smoothed.angularVel)
title('Angular vel distribution');

%Save figure
saveas(gcf,[path,'plots\',sid,'\VelocityDistributions.png']);


%% Median angular speed and fwd vel per condition

%Resample velocities to match the data length
resampledAngSpeed = resample(abs(smoothed.angularVel),length(phase),length(smoothed.angularVel));
resampledFwdVel = resample(smoothed.xVel,length(phase),length(smoothed.xVel));

%Obtain the mean velocities per condition
for contrast = 1:6
    sortedMeanFwdVel(contrast) = nanmean(resampledFwdVel(contrastLimits{contrast}(1):contrastLimits{contrast}(2)));
    sortedMeanAngSpeed(contrast) = nanmean(resampledAngSpeed(contrastLimits{contrast}(1):contrastLimits{contrast}(2)));
end


%Plot
figure('Position',[100, 200, 500, 800]),
subplot(2,1,1)
plot(sortedMeanFwdVel,'-o','color',[0.7 0.2 0.6])
xlim([1 6]);
ylabel('Median forward velocity (mm/s)');
subplot(2,1,2)
plot(sortedMeanAngSpeed,'-o')
xlim([1 6]);
ylabel('Median angular speed (deg/s)');
xlabel('Contrast level');

%Save figure
saveas(gcf,[path,'plots\',sid,'\MeanVelPerConditions.png']);

%This gives me weird results
%maybe I should be looking at instantaneous velocities but rather at
%sliding window analyses?

%% Plot fly's heading

figure('Position',[100 100 1400 300]),

subplot(1,6,1)
polarhistogram(deg2rad(headingToPlot(1:changeContrast(1))),'FaceColor',color_gradient{Intensities(1)})
set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
title(['Contrast ',num2str(Intensities(1))]);

for contrast = 1:length(changeContrast)-1
    subplot(1,6,contrast+1)
    polarhistogram(deg2rad(headingToPlot(changeContrast(contrast):changeContrast(contrast+1))),'FaceColor',color_gradient{Intensities(contrast+1)})
    set(gca,'ThetaZeroLocation','top',...
            'ThetaDir','counterclockwise');
    title(['Contrast ',num2str(Intensities(contrast+1))]);
end
suptitle('Fly heading distribution during the different blocks');

%Save figure
saveas(gcf,[path,'plots\',sid,'\HeadingDistribution.png']);


%% Ordered flies heading

figure('Position',[100 100 1400 300]),
for contrast = 1:6
    subplot(1,6,contrast)
    polarhistogram(deg2rad(correctedHeading(contrastLimits{contrast}(1):contrastLimits{contrast}(2))),'FaceColor',color_gradient{contrast})
    set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    title(['Contrast ',num2str(contrast)]); 
end
suptitle('Ordered fly heading distribution during the different blocks');

%Save figure
saveas(gcf,[path,'plots\',sid,'\OrderedHeadingDistribution.png']);

%% Heading distribution

%use circular sd of heading to determine how stable it was and see of that
%relates to the contrast

heading_var = [];
for contrast = 1:length(changeContrast)
    [s(contrast) heading_var(contrast)] = circ_std(deg2rad(correctedHeading(contrastLimits{contrast}(1):contrastLimits{contrast}(2))),[],[],2);
end

figure('Position',[200,200,1000,400]),
for contrast = 1:6
    subplot(4,6,contrast)
    polarhistogram(deg2rad(correctedHeading(contrastLimits{contrast}(1):contrastLimits{contrast}(2))),'FaceColor',color_gradient{contrast})
    set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    Ax = gca; % current axes
    Ax.RTickLabel = []; 
    Ax.ThetaTick = [0,90,180,270];
    Ax.ThetaTickLabel = [0,90,180,270];
end
subplot(4,6,[7 24])
plot(heading_var,'-ko','MarkerFaceColor','k')
xlabel('Contrast level');
ylabel('Heading variation (circular std)');
xlim([0.5,6.5]);
xticks(1:6);

%Save figure
saveas(gcf,[path,'plots\',sid,'\HeadingVariation.png']);

%% Compute bump mag and with at half max throughout

%Combine PB data in the prober two halves
leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];
combined_full_dff = (data.dff_matrix(leftPB,:) + data.dff_matrix(rightPB,:))/2;

%Bump mag as min-max
BumpMagnitude = max(combined_full_dff) - min(combined_full_dff);

figure('Position',[200 200 1600 400]),
subplot(2,4,[1 3]),
imagesc(data.dff_matrix)
colormap(gray);
hold on
for contrast = 1:length(changeContrast)-1
    line([changeContrast(contrast) changeContrast(contrast)], [0 17], 'LineWidth', 2, 'color', 'r');
end
ylabel('PB glomerulus');
set(gca,'xticklabel',{[]});

subplot(2,4,[5 7]),
plot(data.time,BumpMagnitude, 'k')
hold on
medianBM = median(BumpMagnitude(1:changeContrast(1)));
for contrast = 1:length(changeContrast)-1
    line([data.time(changeContrast(contrast)) data.time(changeContrast(contrast))], [0 5], 'LineWidth', 2, 'color', 'r');
    medianBM = [medianBM, median(BumpMagnitude(changeContrast(contrast):changeContrast(contrast+1)))];
end
xlabel('Time (sec)'); ylabel('Bump magnitude (max-min)');

subplot(2,4,[4,8])
plot(medianBM,'-ko','MarkerFaceColor','k')
xlabel('Contrast level');
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
subplot(3,4,[1 3]),
imagesc(data.dff_matrix)
colormap(gray);
hold on
for contrast = 1:length(changeContrast)-1
    line([changeContrast(contrast) changeContrast(contrast)], [0 17], 'LineWidth', 2, 'color', 'r');
end
ylabel('PB glomerulus');
set(gca,'xticklabel',{[]});

subplot(3,4,[5 7]),
plot(data.time,bump_mag, 'k')
hold on
median_bump_mag = median(bump_mag(1:changeContrast(1)));
median_half_width = median(half_width(1:changeContrast(1)));
for contrast = 1:length(changeContrast)-1
    line([data.time(changeContrast(contrast)) data.time(changeContrast(contrast))], [0 5], 'LineWidth', 2, 'color', 'r');
    median_bump_mag = [median_bump_mag, median(bump_mag(changeContrast(contrast):changeContrast(contrast+1)))];
    median_half_width = [median_half_width, median(half_width(changeContrast(contrast):changeContrast(contrast+1)))];
end
xlabel('Time (sec)'); ylabel('Bump magnitude');
ylim([0 5]);

subplot(3,4,8)
plot(median_bump_mag,'-ko','MarkerFaceColor','k')
title('Median bump magnitude');
ylabel('Bump mag');
set(gca,'xticklabel',{[]});

subplot(3,4,[9 11]),
plot(data.time,half_width, 'k')
hold on
for contrast = 1:length(changeContrast)-1
    line([data.time(changeContrast(contrast)) data.time(changeContrast(contrast))], [0 500], 'LineWidth', 2, 'color', 'r');
end
xlabel('Time (sec)'); ylabel('Bump width at half max');
ylim([0 500]);

subplot(3,4,12)
plot(median_half_width,'-ko','MarkerFaceColor','k')
set(gca,'xticklabel',{[]});
title('Median half width');
ylabel('Bump width at half max');

suptitle('Bump properties with Von Mises fit');

%Save figure
saveas(gcf,[path,'plots\',sid,'\VonMisesBump.png']);


%% Using mean z-score (like Jenny)

mean_zscore = mean(data.z_matrix);

%Compare the 3 different bump magnitudes obtained
figure('Position',[300,300,1000,600]),
subplot(3,1,1)
plot(data.time,BumpMagnitude)
title('Max-min');
set(gca,'xticklabel',{[]});

subplot(3,1,2)
plot(data.time,bump_mag)
title('Von Mises fit');
ylabel('Bump magnitude');
set(gca,'xticklabel',{[]});

subplot(3,1,3)
plot(data.time,mean_zscore)
title('Mean z-score');
xlabel('Time (sec)');

suptitle('Bump magnitude comparison');

%Save figure
saveas(gcf,[path,'plots\',sid,'\BumpComparison.png']);

%% Bump magnitude and width at half max sorted by contrast

figure('Position',[100 100 600 600]),
subplot(1,2,1)
ordered_median_bump_mag = median_bump_mag(ContrastOrder);
plot(ordered_median_bump_mag,'-ko','MarkerFaceColor','k')
title('Median bump magnitude');
ylabel('Bump mag');
xlabel('Contrast level');
xlim([1 6]);
xticks(1:6);

subplot(1,2,2)
ordered_median_half_width = median_half_width(ContrastOrder);
plot(ordered_median_half_width,'-ko','MarkerFaceColor','k')
title('Median half width');
ylabel('Bump width at half max');
xlabel('Contrast level');
xlim([1 6]);
xticks(1:6);

suptitle('Ordered bump parameters with von Mises fit');

%Save figure
saveas(gcf,[path,'plots\',sid,'\VonMisesBumpOrdered.png']);

%% Bump magnitude and width at half max sorted by contrast with z-scored data

median_zscore_bump_mag = [];

for contrast = 1:6
    median_zscore_bump_mag = [median_zscore_bump_mag, median(mean_zscore(contrastLimits{contrast}(1):contrastLimits{contrast}(2)))];
end

figure
plot(median_zscore_bump_mag,'-ko','MarkerFaceColor','k')
title('Median z-scored bump magnitude');
ylabel('Bump mag');
xlabel('Contrast level');
xlim([1 6]);
xticks(1:6);

%Save figure
saveas(gcf,[path,'plots\',sid,'\ZscoredBumpOrdered.png']);


%% Changes in bump magnitude around the jump

%I have noticed that sometimes after the bar jump, the bump magnitude
%decreases.
%I want to see if this depends on the bar contrast
%I'm going to compare the mean bump magnitude in the 30 frames (~3 sec)
%preceding the jumps with the bump magnitude in the 30 frames following the
%jump

for jump = 1:length(jumps)
   figure,
   %Plot the PB heatmap around the jumps
   subplot(3,1,1)
   imagesc(data.dff_matrix(:,jumps(jump)-30:jumps(jump)+30))
   colormap(gray)
   hold on
   line([31 31],[0 17],'color','r','LineWidth',1.5);
   set(gca,'XTickLabel',[]);
   
   %Plot bump magnitude during that time
   subplot(3,1,2)
   plot(data.time(jumps(jump)-30:jumps(jump)+30),bump_mag(jumps(jump)-30:jumps(jump)+30)) 
   line([data.time(jumps(jump)) data.time(jumps(jump))],[0 4],'color','r','LineWidth',1.5);
   ylim([0 4]); 
   xlim([data.time(jumps(jump)-30) , data.time(jumps(jump)+30)]);
   
   %Plot change in mean bump magnitude
   subplot(3,1,3)
   pre_jump_bump_mag(jump) = mean(bump_mag(jumps(jump)-33:jumps(jump)-3));
   post_jump_bump_mag(jump) = mean(bump_mag(jumps(jump)+3:jumps(jump)+33));
   change_in_bump_mag{jump} = [pre_jump_bump_mag(jump),post_jump_bump_mag(jump)];
   plot(change_in_bump_mag{jump},'-o')
   hold on
   xlim([0 3]);
   xticks(1:2);
   xticklabels({'Pre jump','Post jump'});
   ylabel('Bump magnitude (von Mises fit)');
   
   suptitle(['Jump ',num2str(Jumps(jump)),' deg'])
   
   %Save
   saveas(gcf,[path,'plots\',sid,'\BumpMagAJ_',num2str(jump),'.png']);

end

figure,
changes_in_bump_mag = [post_jump_bump_mag-pre_jump_bump_mag];
changes_in_bump_mag = changes_in_bump_mag(ContrastOrder);
plot(1:6,changes_in_bump_mag,'o')
xlim([0,7]);
xticks(1:6);
xlabel('Contrast level');
ylabel('Changes in bump magnitude');

%Save
saveas(gcf,[path,'plots\',sid,'\ChangeInBM_',num2str(jump),'.png']);

%% Bump mag and forward velocity

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
scatter(data.vel_for_ds,bump_mag)
ylabel('Bump magnitude'); xlabel('Forward velocity');

%Save
saveas(gcf,[path,'plots\',sid,'\BumpAndFwdVel_',num2str(jump),'.png']);

%% Bump mag and forward velocity binning the velocity

%Let's always define a fixed number of vel bins
nbins = 10;
maxBin = min(ceil(max(data.vel_for_ds)),15); %I'm adding this line because very high forward velocities are probably more artifact than reality
binWidth = maxBin/nbins;
velBins = [0:binWidth:maxBin]; %I creat3 a vector with my bins
%velBins = [velBins,maxBin+1];

%getting binned medians 
for bin = 1:length(velBins)-1
    meanBin(bin) = mean(bump_mag((abs(data.vel_for_ds))>velBins(bin) & (abs(data.vel_for_ds))<velBins(bin+1)));
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
plot(velAxes,meanBin,'-ko')
ylabel('Mean bump magnitude'); xlabel('Forward velocity (mm/s)');
ylim([0 (max(meanBin)+0.5)]);

%Save figure
saveas(gcf,[path,'plots\',sid,'\BumpVsFwdVel.png']);

%% Using the z-scored data

%getting binned medians 
for bin = 1:length(velBins)-1
    meanBin2(bin) = mean(mean_zscore((abs(data.vel_for_ds))>velBins(bin) & (abs(data.vel_for_ds))<velBins(bin+1)));
end

velAxes = velBins - binWidth;
velAxes = velAxes(2:end);
velAxes(end) = velAxes(end-1)+binWidth;

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time,mean_zscore, 'k')
ylabel('Bump magnitude');
%ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot forward velocity in time) 
subplot(2,4,[5 7])
plot(data.time,data.vel_for_ds,'k')
xlabel('Time (sec)');
ylabel('Forward velocity (mm/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(velAxes,meanBin2,'-ko')
ylabel('Mean bump magnitude as z-score'); xlabel('Forward velocity (mm/s)');

%Save figure
saveas(gcf,[path,'plots\',sid,'\ZscoredBumpVsFwdVel.png']);


%% Trying with my smoothed angular speed, since this relationship looks weird

%Resample my smoothed angular velocity to have the same number of
%datapoints as the bump mag
resampledAngVel = resample(smoothed.angularVel,length(bump_mag),length(smoothed.angularVel));


%Plot the relationship between bump and angular speed with my smoothed
%velocity
angVelBins2 = [0:5:ceil(max(abs(resampledAngVel)))]; %I create a vector with my bins
%angVelBins2 = [angVelBins2,max(abs(resampledAngVel))+10];

%getting binned medians 
for bin = 1:length(angVelBins2)-1
    meanBinAng2(bin) = mean(bump_mag((abs(resampledAngVel))>angVelBins2(bin) & (abs(resampledAngVel))<angVelBins2(bin+1)));
end

angVelAxes2 = angVelBins2 - 5;
angVelAxes2 = angVelAxes2(2:end);
angVelAxes2(end) = angVelAxes2(end-1)+5;

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time,bump_mag, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot angular speed in time) 
subplot(2,4,[5 7])
plot(data.time,abs(resampledAngVel),'k')
xlabel('Time (sec)');
ylabel('Angular speed (deg/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(angVelAxes2,meanBinAng2,'-ko')
ylabel('Mean bump magnitude'); xlabel('Angular speed (deg/s)');
ylim([0 (max(meanBinAng2)+0.5)]);

%Save figure
saveas(gcf,[path,'plots\',sid,'\BumpVsAngSpeedSmoothed.png']);

%% Using z-scored data

%getting binned medians 
for bin = 1:length(angVelBins2)-1
    meanBinAng3(bin) = mean(mean_zscore((abs(resampledAngVel))>angVelBins2(bin) & (abs(resampledAngVel))<angVelBins2(bin+1)));
end

figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time,mean_zscore, 'k')
ylabel('Bump magnitude');
%ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot angular speed in time) 
subplot(2,4,[5 7])
plot(data.time,abs(resampledAngVel),'k')
xlabel('Time (sec)');
ylabel('Angular speed (deg/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(angVelAxes2,meanBinAng3,'-ko')
ylabel('Mean bump magnitude'); xlabel('Angular speed (deg/s)');

%Save figure
saveas(gcf,[path,'plots\',sid,'\ZscoredBumpVsAngSpeed.png']);
 

%% Make binned version of this and/or 3D plot

%getting binned medians 
for angSpeedBin = 1:length(angVelBins2)-1
    for fwdVelBin = 1:length(velBins)-1
        meanBinBin(angSpeedBin,fwdVelBin) = nanmean(bump_mag((abs(resampledAngVel))>angVelBins2(angSpeedBin) & (abs(resampledAngVel))<angVelBins2(angSpeedBin+1) & (abs(data.vel_for_ds'))>velBins(fwdVelBin) & (abs(data.vel_for_ds'))<velBins(fwdVelBin+1)));
    end
end

figure('Position',[100,100,1600,600]),
subplot(1,2,1)
double_binned_bump_mag = flip(meanBinBin);
imagesc(double_binned_bump_mag);
c = colorbar;
ylabel(c, 'Bump magnitude','FontSize',12)
ylabel('Angular speed (deg/s)');
xlabel('Forward velocity (mm/s)');
xticks(1:2:size(meanBinBin,2))
yticks(1:2:size(meanBinBin,1))
set(gca,'xticklabel',{velAxes(1:2:end)});
set(gca,'yticklabel',{round(angVelAxes2(end:-2:1))});
title('von Mises fit bump');


%Repeat for z-scored bump

%getting binned medians 
for angSpeedBin = 1:length(angVelBins2)-1
    for fwdVelBin = 1:length(velBins)-1
        meanBinBin2(angSpeedBin,fwdVelBin) = nanmean(mean_zscore((abs(resampledAngVel))>angVelBins2(angSpeedBin) & (abs(resampledAngVel))<angVelBins2(angSpeedBin+1) & (abs(data.vel_for_ds'))>velBins(fwdVelBin) & (abs(data.vel_for_ds'))<velBins(fwdVelBin+1)));
    end
end

subplot(1,2,2)
double_binned_zscore_bump = flip(meanBinBin2);
imagesc(double_binned_zscore_bump);
c = colorbar;
ylabel(c, 'Bump magnitude (z-score)','FontSize',12)
ylabel('Angular speed (deg/s)');
xlabel('Forward velocity (mm/s)');
xticks(1:2:size(meanBinBin2,2))
yticks(1:2:size(meanBinBin2,1))
set(gca,'xticklabel',{velAxes(1:2:end)});
set(gca,'yticklabel',{round(angVelAxes2(end:-2:1))});
title('mean z-score bump');

%Save figure
saveas(gcf,[path,'plots\',sid,'\BumpVsFwdVelAndAngSpeed.png']);


%save data
fileName = [path,'bar_contrast_data_',sid,'.mat'];
save(fileName,'orderedOffsetChange','ordered_pre_jump_offset','ordered_post_jump_offset','offset_var','heading_var','ordered_median_bump_mag','ordered_median_half_width','median_zscore_bump_mag','changes_in_bump_mag','velAxes','angVelAxes2','double_binned_bump_mag','double_binned_zscore_bump');


%% Linear model to fit the bump magnitude with several variables

%Create contrast level vector populated with contrast level values for each
%timepoint
for timepoint = 1:length(data.time)
    if timepoint < changeContrast(1)
        contrastLevel(timepoint) = Intensities(1);
    elseif (timepoint < changeContrast(2)) & (timepoint >= changeContrast(1))
        contrastLevel(timepoint) = Intensities(2);
    elseif (timepoint < changeContrast(3)) & (timepoint >= changeContrast(2))
        contrastLevel(timepoint) = Intensities(3);
    elseif (timepoint < changeContrast(4)) & (timepoint >= changeContrast(3))
        contrastLevel(timepoint) = Intensities(4);
    elseif (timepoint < changeContrast(5)) & (timepoint >= changeContrast(4))
        contrastLevel(timepoint) = Intensities(5);
    else
        contrastLevel(timepoint) = Intensities(6);       
    end
end

%Create a vector with the actual contrast levels (0,1,2,4,8,15)
actualContrast = contrastLevel;
actualContrast(contrastLevel == 1) = 0;
actualContrast(contrastLevel == 2) = 1;
actualContrast(contrastLevel == 3) = 2;
actualContrast(contrastLevel == 5) = 8;
actualContrast(contrastLevel == 6) = 15;

%figure, plot(contrastLevel)

%Create table with all the variables for the model
ModelVariables = table(data.vel_for_ds,resampledAngSpeed',contrastLevel',data.time,bump_mag', 'VariableNames', {'FwdVel','AngSpeed','ContrastLevel','Time','BumpMagnitude'});

%Fit a regression
mdl = fitlm(ModelVariables)
figure('Position',[100,100,1000,800]),
subplot(3,2,1)
plotDiagnostics(mdl)
subplot(3,2,2)
plotDiagnostics(mdl,'cookd')
subplot(3,2,3)
plotResiduals(mdl,'probability')
subplot(3,2,4)
plotResiduals(mdl,'lagged')
%This shows a very large correlation!
subplot(3,2,5)
plotResiduals(mdl,'fitted')


%Save figure
saveas(gcf,[path,'plots\',sid,'\ModelAnalysis.png']);

%% Bump angular speed

% %unwrap phase value
% unwrapped_phase = unwrap(phase);
% 
% %obtain the derivative
% bump_vel = gradient(unwrapped_phase).*(length(phase)/data.time(end));
% 
% %smooth
% smooth_bump_vel = smooth(bump_vel);
% 
% %get speed
% bump_speed = abs(smooth_bump_vel);
% 
% figure,
% subplot(4,1,1)
% plot(unwrap(correctedHeading))
% hold on
% plot(unwrapped_phase)
% legend('Unwrapped heading', 'Unwrapped phase');
% 
% subplot(4,1,2)
% plot(gradient(unwrap(correctedHeading)).*(length(phase)/data.time(end)))
% hold on
% plot(bump_vel)
% legend('Fly angular velocity', 'Bump angular velocity');
% 
% subplot(4,1,3)
% plot(smooth(gradient(unwrap(correctedHeading)).*(length(phase)/data.time(end))))
% hold on
% plot(smooth_bump_vel)
% legend('Smooth fly angular velocity', 'Smooth bump angular velocity');
% 
% subplot(4,1,4)
% plot(abs(smooth(gradient(unwrap(correctedHeading)).*(length(phase)/data.time(end)))))
% hold on
% plot(bump_speed)
% legend('Fly angular speed', 'Bump angular speed');
% 
% 
% figure,
% subplot(1,2,1)
% scatter(resampledAngSpeed,bump_speed)
% %We will make a second plot remove certain points where the ftpower is very
% %low and theerfore this analysis might not make sense
%     %set threshold as lower 25%
% ft_thresh = prctile(data.ftpower,25);
% %Plot only those data points above that
% subplot(1,2,2)
% scatter(resampledAngSpeed(data.ftpower>ft_thresh),bump_speed(data.ftpower>ft_thresh))
% 
% 
% %Plotting this for the different contrasts
% 
% %Save fly and bump speed per contrast
% contrastFlySpeed = {};
% contrastBumpSpeed = {};
% for contrast = 1:6
%     contrastFlySpeed{contrast} = resampledAngSpeed(contrastLimits{1,contrast}(1):contrastLimits{1,contrast}(2)); 
%     contrastBumpSpeed{contrast} = bump_speed(contrastLimits{1,contrast}(1):contrastLimits{1,contrast}(2)); 
% end
% 
% figure('Position',[100, 100, 1200, 1000]),
% for contrast = 1:6
%     subplot(1,6,contrast)
%     scatter(contrastFlySpeed{contrast},contrastBumpSpeed{contrast},[],color_gradient{contrast}) 
%     title(['Contrast level ',num2str(contrast)]);
%     ylabel('Bump angular speed (deg/s)');
%     xlabel('Fly angular speed');
% end
% 
% %remove low ftpowers
% 
% 
% %bin the data
% figure,
% for contrast = 1:6    
%     angSpeedBins{contrast} = [0:5:ceil(max(contrastFlySpeed{contrast}))]; %I create a vector with my bins
%     angSpeedBins{contrast} = [angSpeedBins{contrast},max(contrastFlySpeed{contrast})+10];
%     
%     %getting binned medians
%     for bin = 1:length(angSpeedBins{contrast})-1
%         medianAngSpeedBin{contrast}(bin) = median(contrastBumpSpeed{contrast}(contrastFlySpeed{contrast}>angSpeedBins{contrast}(bin) & contrastFlySpeed{contrast}<angSpeedBins{contrast}(bin+1)));
%     end
%     
%     angSpeedAxes{contrast} = angSpeedBins{contrast} - 5;
%     angSpeedAxes{contrast} = angSpeedAxes{contrast}(2:end);
%     angSpeedAxes{contrast}(end) = angSpeedAxes{contrast}(end-1)+5;
%     
%     %Plot
%     subplot(1,6,contrast),
%     plot(angSpeedAxes{contrast},medianAngSpeedBin{contrast},'-ko')
%     %ylabel('Median bump magnitude'); xlabel('Angular speed (deg/s)');
%     ylim([0 (max(medianAngSpeedBin{contrast})+0.5)]);
% end
%The bump speed values should be in the same order as fly speed values!
%I should follow the same procedure to compute th bump speed that I follow
%for fly speed.
