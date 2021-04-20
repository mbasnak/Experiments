
%load data
load('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\JennyPanelsOffExp\data\20190329_60D05_6f\dark trial\analysis_sid_8_tid_0.mat');

%% Visualize DF/F data

%average both hemispheres
dff_matrix = data.dff_matrix;

%Plot to visualize
figure('Position', [100 500 1200 300]),
imagesc(dff_matrix)
colormap(gray)
ylabel('PB glomerulus');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


%% Let's focus on the first part of the experiment, where the bump is very clear

panelOFF = find(data.fr_y_ds ==1);
StartDark = panelOFF(1);
cropped_dff_matrix = data.dff_matrix(:,1:StartDark);

figure('Position', [100 500 1200 300]),
imagesc(cropped_dff_matrix)
colormap(gray)
ylabel('PB glomerulus');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);

%% Looking at the bump in one time frame

exampleData = interp1([1:16],cropped_dff_matrix(:,4),linspace(1,16,1000));

figure,
plot(exampleData)
xlabel('PB glomerulus');
ylabel('DF/F');
set(gca,'XTickLabel',[]);


%% Let's combine the two halves of the PB averaging DF/F as Tots suggested, and plot that same timepoint

leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];

combined_dff = (cropped_dff_matrix(leftPB,:) + cropped_dff_matrix(rightPB,:))/2;
exampleData2 = interp1([1:8],combined_dff(:,4),linspace(1,8,1000));

figure,
plot(exampleData2)
xlabel('PB glomerulus');
ylabel('DF/F');
set(gca,'XTickLabel',[]);


%% Let's see what the heatmap looks like

figure('Position', [100 500 1200 300]),
imagesc(combined_dff)
colormap(gray)
ylabel('Averaged PB glomerulus');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);


%% Do von Mises fit with one example data point

[mu,kappa,Value,p,Rsquare] = watson_test(exampleData2);


%% Find the width at half max for this datapoint and plot

fitData = circ_vmpdf(linspace(-pi,pi,100),0,kappa);

% Find the half max value.
halfMax = (min(fitData) + max(fitData)) / 2;
% Find where the data first drops below half the max.
index1 = find(fitData >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(fitData >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes.

figure,
plot(fitData)
hold on
line([index1 index2],[fitData(index1) fitData(index1)]);
ylabel('Probability');
title('Example data for one frame');


%% Use all thew data from the fly to see the fit

combined_full_dff = (dff_matrix(leftPB,:) + dff_matrix(rightPB,:))/2;

half_width = zeros(1,length(combined_full_dff));
bump_mag = zeros(1,length(combined_full_dff));

for i = 1:length(combined_full_dff)
    
    extendedData = interp1([1:8],combined_full_dff(:,i),linspace(1,8,1000));
    [mu(i),kappa(i),Value(i),p2{i},Rsquare(i)] = watson_test(extendedData);
    bump_mag(i) = kappa(i);
    fitData = circ_vmpdf(linspace(-pi,pi,100),0,kappa(i));
    halfMax = (min(fitData) + max(fitData)) / 2;
    index1 = find(fitData >= halfMax, 1, 'first');
    index2 = find(fitData >= halfMax, 1, 'last');
    half_width(i) = index2-index1 + 1;
    
end

%% Plot the results

EndDark = panelOFF(end);

figure('Position',[100 300 1000 800]),
subplot(2,2,1)
plot(bump_mag)
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
title('Bump magnitude (kappa from von Mises fit)');

subplot(2,2,2)
plot(half_width)
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
title('Bump width at half max');

subplot(2,2,3)
plot(smoothdata(bump_mag))
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
xlabel('Time');

subplot(2,2,4)
plot(smoothdata(half_width))
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
xlabel('Time');

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\vonMisesFit\allDataFitMatlab.png');

%% When we look at these results, they show the opposite trend from what we
%would expect. However, if we plot the error of the von Mises fit, we get
%the following

figure,
plot(Rsquare,'color',[0.5 0 0])
xlabel('Time');
title('Error of fit');
xline(StartDark,'lineWidth',3,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\vonMisesFit\fitError.png');

%% Let's now look at the data from the three flies in Jenny's experiments

clear all; close all

leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];

load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp23\data\allDarkData.mat');

for fly = 1:length(Data)
    combined_full_dff{fly} = (Data(fly).dff_matrix(leftPB,:) + Data(fly).dff_matrix(rightPB,:))/2;
    half_width{fly} = zeros(1,length(combined_full_dff{fly}));
    bump_mag{fly} = zeros(1,length(combined_full_dff{fly}));

    for timepoint = 1:length(combined_full_dff{fly})
    
        extendedData = interp1([1:8],combined_full_dff{fly}(:,timepoint),linspace(1,8,1000));
        [mu{fly}(timepoint),kappa{fly}(timepoint),Value{fly}(timepoint),p{fly,timepoint},Rsquare{fly}(timepoint)] = watson_test(extendedData);
        bump_mag{fly}(timepoint) = kappa{fly}(timepoint);
        fitData = circ_vmpdf(linspace(-pi,pi,100),0,kappa{fly}(timepoint));
        halfMax = (min(fitData) + max(fitData)) / 2;
        index1 = find(fitData >= halfMax, 1, 'first');
        index2 = find(fitData >= halfMax, 1, 'last');
        half_width{fly}(timepoint) = index2-index1 + 1;
    
    end

end

%% Plot everything

StartDark = Data(1).StartDarkness;
EndDark = Data(1).EndDarkness;

figure('Position',[100 300 1400 800]),

subplot(3,3,1)
plot(bump_mag{1})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
title('Bump magnitude (kappa from von Mises fit)');

subplot(3,3,2)
plot(half_width{1})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
title('Bump width at half max');

subplot(3,3,3)
plot(Rsquare{1})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
title('Error of the fit');

subplot(3,3,4)
plot(bump_mag{2})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);

subplot(3,3,5)
plot(half_width{2})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);

subplot(3,3,6)
plot(Rsquare{2})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);

subplot(3,3,7)
plot(bump_mag{3})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
xlabel('Time');

subplot(3,3,8)
plot(half_width{3})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
xlabel('Time');

subplot(3,3,9)
plot(Rsquare{3})
hold on
xline(StartDark,'lineWidth',2,'color',[0 0.5 0]);xline(EndDark,'lineWidth',2,'color',[0 0.5 0]);
xlabel('Time');

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\vonMisesFit\allFliesFitMatlab.png');

%% Taking the median values for each state

BumpMag = {};
HalfWidth = {};

for fly = 1:length(bump_mag)
    BumpMag{fly,1} = bump_mag{fly}(1:StartDark);
    BumpMag{fly,2} = bump_mag{fly}(StartDark + 1:EndDark);
    BumpMag{fly,3} = bump_mag{fly}(EndDark:end); 
    
    HalfWidth{fly,1} = half_width{fly}(1:StartDark);
    HalfWidth{fly,2} = half_width{fly}(StartDark + 1:EndDark);
    HalfWidth{fly,3} = half_width{fly}(EndDark:end); 
    
end

medianBumpMag = cellfun(@median,BumpMag);
medianHalfWidth = cellfun(@median,HalfWidth);

figure('Position',[100 300 1400 600]),
subplot(1,2,1)
plot(medianBumpMag,'-o')
title('Median bump magnitude across states');
xlim([0 4]); ylim([0 4]);
xticks([1 2 3]); 
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

subplot(1,2,2)
plot(medianHalfWidth,'-o')
title('Median bump half width across states');
xlim([0 4]);  ylim([15 40]);
xticks([1 2 3]);
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\vonMisesFit\allFliesMedianFit.png');

%% One thing to try could be to fit a von Mises distribution separately to each half of the PB and then average the results.


%First half
for fly = 1:length(Data)
    first_dff{fly} = Data(fly).dff_matrix(1:8,:);
    first_half_width{fly} = zeros(1,length(first_dff{fly}));
    first_bump_mag{fly} = zeros(1,length(first_dff{fly}));
    second_dff{fly} = Data(fly).dff_matrix(9:16,:);
    second_half_width{fly} = zeros(1,length(second_dff{fly}));
    second_bump_mag{fly} = zeros(1,length(second_dff{fly}));

    for timepoint = 1:length(first_dff{fly})
    
        extendedData = interp1([1:8],first_dff{fly}(:,timepoint),linspace(1,8,1000));
        [mu{fly}(timepoint),kappa{fly}(timepoint),Value{fly}(timepoint),p{fly,timepoint},Rsquare{fly}(timepoint)] = watson_test(extendedData);
        first_bump_mag{fly}(timepoint) = kappa{fly}(timepoint);
        fitData = circ_vmpdf(linspace(-pi,pi,100),0,kappa{fly}(timepoint));
        halfMax = (min(fitData) + max(fitData)) / 2;
        index1 = find(fitData >= halfMax, 1, 'first');
        index2 = find(fitData >= halfMax, 1, 'last');
        first_half_width{fly}(timepoint) = index2-index1 + 1;
        
        extendedData = interp1([1:8],second_dff{fly}(:,timepoint),linspace(1,8,1000));
        [mu{fly}(timepoint),kappa{fly}(timepoint),Value{fly}(timepoint),p{fly,timepoint},Rsquare{fly}(timepoint)] = watson_test(extendedData);
        second_bump_mag{fly}(timepoint) = kappa{fly}(timepoint);
        fitData = circ_vmpdf(linspace(-pi,pi,100),0,kappa{fly}(timepoint));
        halfMax = (min(fitData) + max(fitData)) / 2;
        index1 = find(fitData >= halfMax, 1, 'first');
        index2 = find(fitData >= halfMax, 1, 'last');
        second_half_width{fly}(timepoint) = index2-index1 + 1;
        
        total_half_width{fly}(timepoint) = (first_half_width{fly}(timepoint) + second_half_width{fly}(timepoint))/2;
        total_bump_mag{fly}(timepoint) = (first_bump_mag{fly}(timepoint) + second_bump_mag{fly}(timepoint))/2;
        
    end

end

%% Plot this

TotalBumpMag = {};
TotalHalfWidth = {};

for fly = 1:length(bump_mag)
    TotalBumpMag{fly,1} = total_bump_mag{fly}(1:StartDark);
    TotalBumpMag{fly,2} = total_bump_mag{fly}(StartDark + 1:EndDark);
    TotalBumpMag{fly,3} = total_bump_mag{fly}(EndDark:end); 
    
    TotalHalfWidth{fly,1} = total_half_width{fly}(1:StartDark);
    TotalHalfWidth{fly,2} = total_half_width{fly}(StartDark + 1:EndDark);
    TotalHalfWidth{fly,3} = total_half_width{fly}(EndDark:end); 
    
end

medianTotalBumpMag = cellfun(@median,TotalBumpMag);
medianTotalHalfWidth = cellfun(@median,TotalHalfWidth);

figure('Position',[100 300 1400 600]),
subplot(1,2,1)
plot(medianTotalBumpMag,'-o')
title('Median bump magnitude across states');
xlim([0 4]); ylim([0 4]);
xticks([1 2 3]); 
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

subplot(1,2,2)
plot(medianTotalHalfWidth,'-o')
title('Median bump half width across states');
xlim([0 4]);  ylim([15 40]);
xticks([1 2 3]);
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\vonMisesFit\allFliesTotalMedianFit.png');


%% Compute von Mises with matlab circ_stats package

for fly = 1:length(Data)
    combined_full_dff{fly} = (Data(fly).dff_matrix(1:8,:) + Data(fly).dff_matrix(9:16,:))/2;
    half_width{fly} = zeros(1,length(combined_full_dff{fly}));
    bump_mag{fly} = zeros(1,length(combined_full_dff{fly}));

    for timepoint = 1:length(combined_full_dff{fly})
    
        extendedData = interp1([1:8],combined_full_dff{fly}(:,timepoint),linspace(1,8,1000));
        [mu{fly}(timepoint),kappa{fly}(timepoint)] = circ_vmpar(extendedData);
        bump_mag{fly}(timepoint) = kappa{fly}(timepoint);
        fitData = circ_vmpdf(linspace(-pi,pi,100),0,kappa{fly}(timepoint));
        halfMax = (min(fitData) + max(fitData)) / 2;
        index1 = find(fitData >= halfMax, 1, 'first');
        index2 = find(fitData >= halfMax, 1, 'last');
        half_width{fly}(timepoint) = index2-index1 + 1;
    
    end

end

BumpMag = {};
HalfWidth = {};

for fly = 1:length(bump_mag)
    BumpMag{fly,1} = bump_mag{fly}(1:StartDark);
    BumpMag{fly,2} = bump_mag{fly}(StartDark + 1:EndDark);
    BumpMag{fly,3} = bump_mag{fly}(EndDark:end); 
    
    HalfWidth{fly,1} = half_width{fly}(1:StartDark);
    HalfWidth{fly,2} = half_width{fly}(StartDark + 1:EndDark);
    HalfWidth{fly,3} = half_width{fly}(EndDark:end); 
    
end

medianBumpMag = cellfun(@median,BumpMag);
medianHalfWidth = cellfun(@median,HalfWidth);

figure('Position',[100 300 1400 600]),
subplot(1,2,1)
plot(medianBumpMag,'-o')
title('Median bump magnitude across states');
xlim([0 4]); ylim([0 4]);
xticks([1 2 3]); 
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

subplot(1,2,2)
plot(medianHalfWidth,'-o')
title('Median bump half width across states');
xlim([0 4]);  ylim([15 40]);
xticks([1 2 3]);
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

%I get the same results using this function...

%% Using Yvette's function


clear all; close all

leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];

load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp23\data\allDarkData.mat');

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
StartDark = Data(1).StartDarkness;
EndDark = Data(1).EndDarkness;

BumpMag = {};
HalfWidth = {};

for fly = 1:length(bump_mag)
    BumpMag{fly,1} = bump_mag{fly}(1:StartDark);
    BumpMag{fly,2} = bump_mag{fly}(StartDark + 1:EndDark);
    BumpMag{fly,3} = bump_mag{fly}(EndDark:end); 
    
    HalfWidth{fly,1} = half_width{fly}(1:StartDark);
    HalfWidth{fly,2} = half_width{fly}(StartDark + 1:EndDark);
    HalfWidth{fly,3} = half_width{fly}(EndDark:end); 
    
end

medianBumpMag = cellfun(@median,BumpMag);
medianHalfWidth = cellfun(@median,HalfWidth);

figure('Position',[100 300 1400 600]),
subplot(1,2,1)
plot(medianBumpMag,'-o')
xlim([0 4]);
xticks([1 2 3]); 
xticklabels({'certain' ,'uncertain', 'certain'});
ylabel('Median bump magnitude')
legend({'fly1','fly2','fly3'});

subplot(1,2,2)
plot(medianHalfWidth','-o')
ylabel('Median bump half width');
xlim([0 4]); ylim([15 30]);
xticks([1 2 3]);
xticklabels({'certain' ,'uncertain', 'certain'});
legend({'fly1','fly2','fly3'});

bump_mag = bump_mag';
half_width = half_width';
bump_mag = [bump_mag{1,1},bump_mag{2,1},bump_mag{3,1}];
half_width = [half_width{1,1},half_width{2,1},half_width{3,1}];

save('C:\Users\Melanie\Dropbox (HMS)\UncertaintyProject\JennyPanelsOffExp\data\fitData.mat','bump_mag','half_width')
saveas(gcf,'Z:\Wilson Lab\Mel\Experiments\Uncertainty\vonMisesFit\allFliesYvetteFit.png');
