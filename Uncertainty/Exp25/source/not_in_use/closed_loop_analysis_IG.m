%analysis for the closed-loop bouts of the change in contrast experiment
%for the inverted-gain flies


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
changeContrast = find(abs(diff(data.fr_y_ds))>1);

%% Identify contrast order

%Check the position function
pos_function = data.run_obj.function_number;

%Define the order of intensities according to the function used
if pos_function == 195  
    Intensities = [1,2,1,3,2,3];
elseif pos_function == 196  
    Intensities = [2,1,3,1,2,3];
else
    Intensities = [3,1,2,1,2,3];
end


%% Offset


%Plot
% Plot the heatmap of EPG activity
figure('Position',[100 200 1800 800]),
subplot(4,6,[1 6])
%I'm flipping the dff matrix for it to make sense along with the fly's
%heading
imagesc(flip(data.dff_matrix))
colormap(gray)
hold on
%add the changes in stim
for change = 1:length(changeContrast)
    line([changeContrast(change) changeContrast(change)], [0 17], 'LineWidth', 2, 'color', 'r');
end
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
ylabel('PB glomerulus');
title('EPG activity in the PB');
set(gca,'XTickLabel',[]);
legend('Change in stimulus');

% Plot the heading and the EPG phase
subplot(4,6,[7 12])
%Get heading to plot
%The gain was inverted here, so I'm negating the heading to obtain the
%correct one
heading = wrapTo180(-data.heading_deg);
change_heading = abs([0;diff(smooth(heading))]);
heading_to_plot = smooth(heading);
heading_to_plot(change_heading>40==1) = NaN;
plot(heading_to_plot,'LineWidth',1.5)
plot(data.time,heading_to_plot,'color',[0.6 0.2 0.4],'LineWidth',1.5)
hold on
%Get EPG phase to plot
%I'm now going to negate the phase, since I'm plotting heading instead of
%bar position, to the bump moves in the other direction
phase = wrapTo180(rad2deg(-data.dff_pva));
change_phase = abs([0;diff(smooth(phase))]);
phase_to_plot = smooth(phase);
phase_to_plot(change_phase>40==1) = NaN;
plot(data.time,phase_to_plot,'color',[0.2 0.6 0.4],'LineWidth',1.5)
%add the changes in stim
for change = 1:length(changeContrast)
    line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r');
end
legend('Fly heading', 'EPG phase');
ylim([-180, 180]);
xlim([0,data.time(end)]);
ylabel('Deg');
set(gca,'XTickLabel',[]);

% Plot the offset
subplot(4,6,[13 18])
%Get offset to plot
offset = wrapTo180(rad2deg(circ_dist(deg2rad(heading),deg2rad(phase))));
change_offset = abs([0;diff(smooth(offset))]);
offset_to_plot = smooth(offset);
offset_to_plot(change_offset>30==1) = NaN;
plot(data.time,offset_to_plot,'LineWidth',1.5,'color','k')
%add the changes in stim
for change = 1:length(changeContrast)
    line([data.time(changeContrast(change)) data.time(changeContrast(change))], [-180 180], 'LineWidth', 2, 'color', 'r');
end
ylim([-180 180]);
ylabel('Deg'); xlabel('Time (sec)');
legend('Offset');
%improve offset!

% Polar histograms of offset
% Color histograms acording to the intensity level of the bar, using the
% vector Intensities
color_gradient = {[0,0,0],[0 0 0.6],[ 0.5 0.8 0.9]};
subplot(4,6,19)
polarhistogram(deg2rad(offset(1:changeContrast(1))),15,'FaceColor',color_gradient{Intensities(1)})
set(gca,'ThetaZeroLocation','top',...
    'ThetaDir','counterclockwise');
if Intensities(1) == 1
    title({'Offset distribution','darkness'});
elseif Intensities(1) == 2
    title({'Offset distribution','low contrast'});
else
    title({'Offset distribution','high contrast'});
end

for contrast = 1:length(changeContrast)-1
    subplot(4,6,19+contrast)
    polarhistogram(deg2rad(offset(changeContrast(contrast):changeContrast(contrast+1))),15,'FaceColor',color_gradient{Intensities(contrast+1)})
    set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    if Intensities(contrast+1) == 1
        title({'Offset distribution','darkness'});
    elseif Intensities(contrast+1) == 2
        title({'Offset distribution','low contrast'});
    else
        title({'Offset distribution','high contrast'});
    end
end

%add last contrast for pos function 197, which appears to have one less
%change in stim
if pos_function == 197
    subplot(4,6,24)
    polarhistogram(deg2rad(offset(changeContrast(contrast+1):end)),15,'FaceColor',color_gradient{Intensities(6)})
    set(gca,'ThetaZeroLocation','top',...
        'ThetaDir','counterclockwise');
    title({'Offset distribution','high contrast'});
end


%save figure
saveas(gcf,[path,'plots\closedLoopHeatmapAndOffset.png']);

%% Set block limits

blockLimits{1} = [1,changeContrast(1)-1];

if length(changeContrast) == 5
    for block = 2:5
        blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
    end
    blockLimits{6} = [changeContrast(5),length(heading)];
elseif length(changeContrast) == 6
    for block = 2:6
        blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
    end
else
    for block = 2:4
        blockLimits{block} = [changeContrast(block-1),changeContrast(block)-1];
    end
     blockLimits{5} = [changeContrast(4),length(heading)];
end

%% Offset distribution

figure('Position',[200 200 1000 800]),
subplot(2,1,1)
for block = 1:length(blockLimits)
   [~, offset_var(block)] = circ_std(deg2rad(offset(blockLimits{block}(1):blockLimits{block}(2))),[],[],1); 
   plot(block,offset_var(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
   hold on
end
if (contains(path, '20201020_60D05_7f') & ~contains(path,'fly2'))
    xlim([0 6]);
    xticks(1:5);
else
    xlim([0 7]);
    xticks(1:6);
end
title('Offset variation per block');
ylabel({'Circular standard deviation','of the offset'});
%Add custom legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
legend(h, 'Darkness','Low contrast','High Contrast');

subplot(2,1,2)
%Create table with contrast level and offset variation
for block = 1:length(blockLimits)
   if Intensities(block) == 1
       contrasts{block} = 'Darkness';
   elseif Intensities(block) == 2
       contrasts{block} = 'Low contrast';
   else
       contrasts{block} = 'High contrast';
   end       
end
summary_data = table(contrasts',offset_var','VariableNames',{'contrast','offset_var'});
%Get mean offset var by contrast
mean_data = varfun(@mean,summary_data,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast'});
%Plot
newOrder = [1,3,2];
for contrastLevel = 1:3
   plot(contrastLevel,mean_data.mean_offset_var(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
   hold on
end
xlim([0 4]);
ylim([min(mean_data.mean_offset_var)-0.5 max(mean_data.mean_offset_var)+0.5]); 
xticks(1:3);
xticklabels({'Darkness', 'Low contrast', 'High contrast'});
title('Mean offset variation per contrast');
ylabel({'Circular standard deviation','of the offset'});

%save figure
saveas(gcf,[path,'plots\closedLoopOffsetVariation.png']);

%% Get offset value from last bout of high contrast

%Find limits of last high contrast bout
high_contrast_bouts = find(Intensities==3);
last_high_contrast = high_contrast_bouts(end);

%This will be the 'visual' offset, so we need to compute it as the circular
%distance between the stimulus and the phase
visual_offset = circ_dist(data.dff_pva,deg2rad(wrapTo180(data.panel_angle)));
mean_reference_offset = rad2deg(circ_mean(visual_offset(blockLimits{high_contrast_bouts(1)}(1):blockLimits{high_contrast_bouts(1)}(2)),[],1));


%this reference offset represents the difference phase-visualstimulus
figure,
subplot(3,1,1)
plot(rad2deg(data.dff_pva))
title('Bump position');
subplot(3,1,2)
plot(wrapTo180(data.panel_angle))
title('Bar position');
subplot(3,1,3)
plot(wrapTo180(rad2deg(data.dff_pva)+mean_reference_offset))
title('Bump estimate of bar position');



%% Calculate and plot bump magnitude in time

%compute bump magnitude as max-min
for block = 1:length(blockLimits)
   bump_mag{block} = max(data.mean_dff_EB(:,blockLimits{block}(1):blockLimits{block}(2)))-min(data.mean_dff_EB(:,blockLimits{block}(1):blockLimits{block}(2))); 
end

figure('Position',[200 200 1600 600]),
%plot EPG activity
subplot(2,1,1)
imagesc(flip(data.dff_matrix))
colormap(gray)
hold on
%add the changes in stim
for change = 1:length(changeContrast)
   line([changeContrast(change) changeContrast(change)], [0 17], 'LineWidth', 2, 'color', 'r'); 
end
yticks(1:2:16);
yticklabels({'8L','6L','4L','2L','1R','3R','5R','7R'});
ylabel('PB glomerulus');
title('EPG activity in the PB');
set(gca,'XTickLabel',[]);
legend('Change in stimulus');

%plot bump magnitude
subplot(2,1,2)
for block = 1:length(blockLimits)
    if contains(contrasts(block),'Dark')
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{1})
    elseif contains(contrasts(block),'Low')
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{2})
    else
        plot(data.time(blockLimits{block}(1):blockLimits{block}(2)),bump_mag{block},'color',color_gradient{3})
    end
    hold on
end
title('Bump magnitude');
ylabel('Bump magnitude (max-min)');
xlabel('Time (sec)');
xlim([0 data.time(end)]);

%save figure
saveas(gcf,[path,'plots\closedLoopBMinTime.png']);

%% Compute and plot mean bump magnitude per block

figure('Position',[200 200 1000 800]),
subplot(2,1,1)
for block = 1:length(blockLimits)
   mean_bump_mag(block) = mean(bump_mag{block}); 
   plot(block,mean_bump_mag(block),'ko','MarkerFaceColor',color_gradient{Intensities(block)},'MarkerSize',8)
   hold on
end
xlim([0 7]);
xticks(1:6);

title('Mean bump magnitude per block');
ylabel('Bump magnitude (max-min)');
%Add custom legend
h = zeros(3, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{1});
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{2});
h(3) = plot(NaN,NaN,'ko','MarkerFaceColor',color_gradient{3});
legend(h, 'Darkness','Low contrast','High Contrast');

subplot(2,1,2)
%add bump mag data to table
if size(summary_data,2) == 2
    summary_data = addvars(summary_data, mean_bump_mag','NewVariableNames','mean_bump_mag');
end
%Get mean bump mag by contrast
mean_bump_data = varfun(@mean,summary_data,'InputVariables','mean_bump_mag',...
       'GroupingVariables',{'contrast'});
%Plot
for contrastLevel = 1:3
   plot(contrastLevel,mean_bump_data.mean_mean_bump_mag(newOrder(contrastLevel)),'ko','MarkerFaceColor',color_gradient{contrastLevel},'MarkerSize',8)
   hold on
end
xlim([0 4]);
ylim([min(mean_bump_data.mean_mean_bump_mag)-0.5 max(mean_bump_data.mean_mean_bump_mag)+0.5]); 
xticks(1:3);
xticklabels({'Darkness', 'Low contrast', 'High contrast'});
title('Mean bump magnitude per contrast');
ylabel('Bump magnitude (max-min)');

%save figure
saveas(gcf,[path,'plots\closedLoopMeanBM.png']);

%% Plot fly's velocity in all 3 axes  with all vel in deg/s

figure('Position',[100, 100, 1600, 1000]),
subplot(4,5,[1 4])
plot(data.time,data.vel_for_deg_ds)
title('Forward velocity');
ylabel('Forward velocity (deg/s)');

subplot(4,5,5)
histogram(data.vel_for_deg_ds)
title('Velocity distributions');
xlabel('Forward velocity (deg/s)');
ylabel('Counts');

subplot(4,5,[6 9])
plot(data.time,data.vel_side_deg_ds,'color',[0.8 0.2 0.6])
title('Side velocity');
ylabel('Side velocity (deg/s)');

subplot(4,5,10)
histogram(data.vel_side_deg_ds,'FaceColor',[0.8 0.2 0.6])
xlabel('Side velocity (deg/s)');
ylabel('Counts');

subplot(4,5,[11 14])
plot(data.time,data.vel_yaw_ds,'color',[0.4 0.2 0.8])
title('Angular velocity');
ylabel('Angular velocity (deg/s)');

subplot(4,5,15)
histogram(data.vel_yaw_ds,'FaceColor',[0.4 0.2 0.8])
xlabel('Angular velocity (deg/s)');
ylabel('Counts');

subplot(4,5,[16 19])
plot(data.time,data.total_mvt_ds,'k')
title('Total movement');
ylabel('Total movement (deg/s)');
xlabel('Time (sec)');

subplot(4,5,20)
histogram(data.total_mvt_ds,'FaceColor','k')
xlabel('Total movement (deg/s)');
ylabel('Counts');

%save figure
saveas(gcf,[path,'plots\closedLoopVelDistributions.png']);

%% Relationship between bump magnitude and velocity

allBumpMag = [];
for block = 1:length(blockLimits)
    allBumpMag = [allBumpMag,bump_mag{block}];
end

figure('Position',[200 200 1600 800]),
subplot(1,4,1)
scatter(data.vel_for_deg_ds(1:length(allBumpMag)),allBumpMag)
xlim([0 max(data.vel_for_deg_ds)]);
xlabel('Forward velocity (deg/s)');
ylabel('Bump magnitude (max-min)');

subplot(1,4,2)
scatter(abs(data.vel_side_deg_ds(1:length(allBumpMag))),allBumpMag,[],[0.8 0.2 0.6])
xlim([0 max(abs(data.vel_side_deg_ds))]);
xlabel('Side speed (deg/s)');

subplot(1,4,3)
scatter(abs(data.vel_yaw_ds(1:length(allBumpMag))),allBumpMag,[],[0.4 0.2 0.8])
xlim([0 max(abs(data.vel_yaw_ds))]);
xlabel('Angular speed (deg/s)');

subplot(1,4,4)
scatter(data.total_mvt_ds(1:length(allBumpMag)),allBumpMag,[],'k')
xlim([0 max(data.total_mvt_ds)]);
xlabel('Total movement (deg/s)');

%save figure
saveas(gcf,[path,'plots\closedLoopBMvsVel.png']);

%% Binning it

%Let's always define a fixed number of vel bins
nbins = 20;
maxBin = max(data.total_mvt_ds);
binWidth = maxBin/nbins;
mvtBins = [0:binWidth:maxBin]; 

%getting binned means 
for bin = 1:length(mvtBins)-1
    meanBin(bin) = mean(allBumpMag((data.total_mvt_ds(1:length(allBumpMag)) > mvtBins(bin)) & (data.total_mvt_ds(1:length(allBumpMag)) < mvtBins(bin+1))));
end

%create axes for plot
mvtAxes = mvtBins - binWidth;
mvtAxes = mvtAxes(2:end);
mvtAxes(end) = mvtAxes(end-1)+binWidth;

%Plot
figure('Position',[200 200 1400 600]),
%Plot bump magnitude in time
subplot(2,4,[1 3])
plot(data.time(1:length(allBumpMag)),allBumpMag, 'k')
ylabel('Bump magnitude');
ylim([0 5]);
set(gca,'xticklabel',{[]});

%Plot total movement 
subplot(2,4,[5 7])
plot(data.time(1:length(allBumpMag)),data.total_mvt_ds(1:length(allBumpMag)),'k')
xlabel('Time (sec)');
ylabel('Total movement (deg/s)');

%Plot relationship between both parameters
subplot(2,4,[4,8]);
plot(mvtAxes,meanBin,'-ko')
ylabel('Mean bump magnitude'); xlabel('Total movement (deg/s)');
ylim([0 (max(meanBin)+0.5)]);

%save figure
saveas(gcf,[path,'plots\closedLoopBMvsVelBinned.png']);

%% Model bump magnitude as a function of contrastLevel and total movement

%create vector with the contrast level for each timepoint
all_contrast_levels = [];
for block = 1:length(blockLimits);
    contrast_level{block} = repelem(Intensities(block),blockLimits{block}(2)+1-blockLimits{block}(1));
    all_contrast_levels = [all_contrast_levels,contrast_level{block}];
end

%crate table with the model's variables
modelTable = table(all_contrast_levels',data.total_mvt_ds(1:length(allBumpMag))',data.time(1:length(allBumpMag)),allBumpMag','VariableNames',{'ContrastLevel','TotalMovement','Time','BumpMagnitude'});

%fit linear model using contrast level as a categorical variable
mdl_BM = fitlm(modelTable,'BumpMagnitude~ContrastLevel+TotalMovement+Time','CategoricalVars',1);
%This results in very small weigths being assigned to total movement, but
%that probably is because being in deg/s, total movement has large units.

%fit linear model using contrast level as a numerical variable
mdl_BM2 = fitlm(modelTable,'BumpMagnitude~ContrastLevel+TotalMovement+Time');

%% Save useful data

save([path,'\summary_data.mat'],'summary_data','mean_reference_offset','modelTable')

close all; clc;