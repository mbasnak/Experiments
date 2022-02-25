%Code for the group analysis of the closed-loop bouts using the continuous
%method

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
%Get folder names
folderContents = dir(path);

%Load the summary data of the folder that correspond to experimental flies
for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\continuous_summary_data.mat']);
   end
end

%% Clean and combine data

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%Combine the tables into one
allSummaryData = array2table(zeros(0,5),'VariableNames', {'contrast','offset_var','mean_bump_mag','mean_bump_width','heading_var'});
flyNumber= [];
for fly = 1:length(data)
    flyNumber = [flyNumber,repelem(fly,length(data(fly).summary_data.contrast))];
    allSummaryData = [allSummaryData;data(fly).summary_data]; 
end

%Add the fly ID as a variable
allSummaryData = addvars(allSummaryData,nominal(flyNumber'));
allSummaryData.Properties.VariableNames{'Var6'} = 'Fly';

%% Compare offset variability statistically

%Reorder data to bring the darkness observation to the first row, such that the model
%will use 'darkness' as the default
reordered_data = allSummaryData;
reordered_data(height(reordered_data)+1,:) = reordered_data(1,:);
reordered_data(1,:) = reordered_data(2,:);
reordered_data(2,:) = reordered_data(height(reordered_data),:);
reordered_data(height(reordered_data),:) = [];

%Run model
mdl_offset = fitlme(reordered_data,'offset_var~contrast+(1|Fly)')
disp(anova(mdl_offset))

%Add post-hoc comparison between low and high contrast
for y = 1:numel(mdl_offset.CoefficientNames)
    fprintf('Contrast position %i: %s\n', y, char(mdl_offset.CoefficientNames{y}));
end
bonfCorrection = 1;
[pVal F df1 df2] = coefTest(mdl_offset, [0  1 -1]);   pVal = min(1,pVal*bonfCorrection);           fprintf('Planned comparison low conrast vs high contrast,  F(%i, %i)=%0.2f,\t p=%0.7f (contrast 0  1 -1)\n', df1, df2, F, pVal);

%% Get and plot mean offset variation per fly and in total

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per fly
mean_offset_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast','Fly'});
%Plot
allFlies = [];
for fly = 1:length(data)
    fly_data{fly} = [mean_offset_data_per_fly.mean_offset_var(fly),mean_offset_data_per_fly.mean_offset_var(fly+2*length(data)),mean_offset_data_per_fly.mean_offset_var(fly+length(data))];
    allFlies = [allFlies;fly_data{fly}];
end
plot(1:3,allFlies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset variation by contrast
mean_offset_data = varfun(@mean,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast'});
%Get error in offset variation
std_offset_data = varfun(@std,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast'});
hold on
%Add mean and error
errorbar(1:3,mean(allFlies),std(allFlies)/sqrt(length(allFlies)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)

%Add significance 
% for contrast = 2
%    if (mdl_offset.Coefficients.pValue(contrast)<0.05 & mdl_offset.Coefficients.pValue(contrast)>=0.01)
%       text(contrast-0.58,2,'*','fontsize',26);
%    elseif (mdl_offset.Coefficients.pValue(contrast)<0.01 & mdl_offset.Coefficients.pValue(contrast)>=0.001)
%        text(contrast-0.58,2,'**','fontsize',26);
%    elseif mdl_offset.Coefficients.pValue(contrast)<0.001
%        text(contrast-0.58,2,'***','fontsize',26);
%    else
%       text(contrast-0.58,2,'ns');
%    end
% end
% for contrast = 3
%    if (mdl_offset.Coefficients.pValue(contrast)<0.05 & mdl_offset.Coefficients.pValue(contrast)>=0.01)
%       text(contrast-1.1,2.2,'*','fontsize',26);
%    elseif (mdl_offset.Coefficients.pValue(contrast)<0.01 & mdl_offset.Coefficients.pValue(contrast)>=0.001)
%        text(contrast-1.1,2.2,'**','fontsize',26);
%    elseif mdl_offset.Coefficients.pValue(contrast)<0.001
%        text(contrast-1.1,2.2,'***','fontsize',26);
%    else
%       text(contrast-1.1,2.2,'ns');
%    end
% end
% for contrast = 2
%    if (pVal<0.05 & pVal>=0.01)
%       text(contrast+0.43,1.78,'*','fontsize',26);
%    elseif (pVal<0.01 & pVal>=0.001)
%        text(contrast+0.43,1.78,'**','fontsize',26);
%    elseif pVal<0.001
%        text(contrast+0.43,1.78,'***','fontsize',26);
%    else
%       text(contrast+0.43,1.78,'ns');
%    end
% end
% 
% %add line between asterisks
% %contrasts 1 and 2
% line([1 1.4],[2.02 2.02],'color','k');
% line([1.66 2],[2.02 2.02],'color','k');
% %contrasts 2 and 3
% line([2 2.4],[1.8 1.8],'color','k');
% line([2.66 3],[1.8 1.8],'color','k');
% %contrasts 1 and 3
% line([1 1.88],[2.22 2.22],'color','k');
% line([2.15 3],[2.22 2.22],'color','k');

xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Offset variability (rad)','FontSize',12);
ylim([min(mean_offset_data_per_fly.mean_offset_var)-0.3 max(mean_offset_data_per_fly.mean_offset_var)+0.3]);

%save figure
saveas(gcf,[path,'\continuousGroupPlots\mean_closed_loop_offset_var.png']);
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\CueBrightness-Experiment\mean_closed_loop_offset_var.png');

%% Run model for heading variation

%Run model
mdl_heading = fitlme(reordered_data,'heading_var~contrast+(1|Fly)')
disp(anova(mdl_heading))

%Add post-hoc comparison between low and high contrast
for y = 1:numel(mdl_heading.CoefficientNames)
    fprintf('Contrast position %i: %s\n', y, char(mdl_heading.CoefficientNames{y}));
end
bonfCorrection = 1;
[pVal F df1 df2] = coefTest(mdl_heading, [0  1 -1]);pVal = min(1,pVal*bonfCorrection); fprintf('Planned comparison low conrast vs high contrast,  F(%i, %i)=%0.2f,\t p=%0.7f (contrast 0  1 -1)\n', df1, df2, F, pVal);


%% Get and plot mean heading variability per fly and in total

%Get mean heading var by contrast per fly
mean_heading_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','heading_var',...
       'GroupingVariables',{'contrast','Fly'});
   
figure('Position',[200 200 1000 800]),
%Plot
allflies = [];
for fly = 1:length(data)
    fly_data_heading{fly} = [mean_heading_data_per_fly.mean_heading_var(fly),mean_heading_data_per_fly.mean_heading_var(fly+2*length(data)),mean_heading_data_per_fly.mean_heading_var(fly+length(data))];
    allflies = [allflies;fly_data_heading{fly}];
end
plot(1:3,allflies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean heading var by contrast in total
mean_heading_data = varfun(@mean,allSummaryData,'InputVariables','heading_var',...
       'GroupingVariables',{'contrast'});
hold on
errorbar(1:3,mean(allflies),std(allflies)/sqrt(length(allFlies)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)

% %add significance 
% for contrast = 2
%    if (mdl_heading.Coefficients.pValue(contrast)<0.05 & mdl_heading.Coefficients.pValue(contrast)>=0.01)
%       text(contrast-0.58,2,'*','fontsize',26);
%    elseif (mdl_heading.Coefficients.pValue(contrast)<0.01 & mdl_heading.Coefficients.pValue(contrast)>=0.001)
%        text(contrast-0.58,2,'**','fontsize',26);
%    elseif mdl_heading.Coefficients.pValue(contrast)<0.001
%        text(contrast-0.58,2,'***','fontsize',26);
%    else
%       text(contrast-0.58,2,'ns');
%    end
% end
% for contrast = 3
%    if (mdl_heading.Coefficients.pValue(contrast)<0.05 & mdl_heading.Coefficients.pValue(contrast)>=0.01)
%       text(contrast-1.1,2.2,'*','fontsize',26);
%    elseif (mdl_heading.Coefficients.pValue(contrast)<0.01 & mdl_heading.Coefficients.pValue(contrast)>=0.001)
%        text(contrast-1.1,2.2,'**','fontsize',26);
%    elseif mdl_heading.Coefficients.pValue(contrast)<0.001
%        text(contrast-1.1,2.2,'***','fontsize',26);
%    else
%       text(contrast-1.1,2.2,'ns');
%    end
% end
% for contrast = 2
%    if (pVal<0.05 & pVal>=0.01)
%       text(contrast+0.43,1.78,'*','fontsize',26);
%    elseif (pVal<0.01 & pVal>=0.001)
%        text(contrast+0.43,1.78,'**','fontsize',26);
%    elseif pVal<0.001
%        text(contrast+0.43,1.78,'***','fontsize',26);
%    else
%       text(contrast+0.48,1.81,'ns','fontsize',16);
%    end
% end
% %add line between asterisks
% %contrasts 1 and 2
% line([1 1.4],[2.02 2.02],'color','k');
% line([1.66 2],[2.02 2.02],'color','k');
% %contrasts 2 and 3
% line([2 2.4],[1.8 1.8],'color','k');
% line([2.66 3],[1.8 1.8],'color','k');
% %contrasts 1 and 3
% line([1 1.88],[2.22 2.22],'color','k');
% line([2.15 3],[2.22 2.22],'color','k');

xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Heading variability (rad)','FontSize',12);
ylim([min(mean_heading_data_per_fly.mean_heading_var)-0.3 max(mean_heading_data_per_fly.mean_heading_var)+0.3]);

%save figure
saveas(gcf,[path,'\continuousGroupPlots\mean_closed_loop_heading_var.png']);

%% Combine bump parameters tables into one

%Add zscored bump magnitude and bump width to the tables for each
%individual fly
for fly = 1:length(data)   
    zscoredBM = zscore(data(fly).modelTable.BumpMagnitude);
    zscoredBW = zscore(data(fly).modelTable.BumpWidth);
    if size(data(fly).modelTable,2) < 15
        data(fly).modelTable = addvars(data(fly).modelTable,zscoredBM,'NewVariableNames','zscoredBM');
        data(fly).modelTable = addvars(data(fly).modelTable,zscoredBW,'NewVariableNames','zscoredBW');
    end
end


%Combine the tables into one
allModelData = array2table(zeros(0,15),'VariableNames', {'ContrastLevel','ForVelocity','ZscoredForVel','SideSpeed','ZscoredSideSpeed','YawSpeed','ZscoredYawSpeed','TotalMovement','ZscoredTotalMovement','Time','BumpMagnitude','AdjRSquare','BumpWidth','zscoredBM','zscoredBW'});
Fly = [];
for fly = 1:length(data)
    Fly = [Fly,repelem(fly,length(data(fly).modelTable.ContrastLevel))];
    allModelData = [allModelData;data(fly).modelTable]; 
end
%Add fly id as a categorical variable
if size(data(fly).modelTable,2) < 16
    allModelData = addvars(allModelData,nominal(Fly'));
    allModelData.Properties.VariableNames{'Var16'} = 'Fly';
end

%Change contrast level as a categorical variable
allModelData.ContrastLevel = nominal(allModelData.ContrastLevel);

%% Compute model for bump magnitude
adj_rs = allModelData.AdjRSquare;

%Fit different models
mdl_BM{1} = fitlme(allModelData(adj_rs >= 0.5,:),'BumpMagnitude~ContrastLevel+TotalMovement+(1|Fly)');
mdl_BM{2} = fitlme(allModelData(adj_rs >= 0.5,:),'BumpMagnitude~ContrastLevel+ForVelocity+SideSpeed+YawSpeed+(1|Fly)');
mdl_BM{3} = fitlme(allModelData(adj_rs >= 0.5,:),'BumpMagnitude~ContrastLevel+ForVelocity+SideSpeed+YawSpeed+Time+(1|Fly)');
%zscored data
mdl_BM{4} = fitlme(allModelData(adj_rs >= 0.5,:),'BumpMagnitude~ContrastLevel+ZscoredTotalMovement+(1|Fly)');
mdl_BM{5} = fitlme(allModelData(adj_rs >= 0.5,:),'BumpMagnitude~ContrastLevel+ZscoredForVel+ZscoredSideSpeed+ZscoredYawSpeed+(1|Fly)');
mdl_BM{6} = fitlme(allModelData(adj_rs >= 0.5,:),'BumpMagnitude~ContrastLevel+ZscoredForVel+ZscoredSideSpeed+ZscoredYawSpeed+Time+(1|Fly)');

%Plot model comparison
figure,
%Plot R squared
subplot(1,2,1)
for model = 1:length(mdl_BM)
   Rsquared(model) = mdl_BM{1,model}.Rsquared.Adjusted; 
end
plot(Rsquared,'-ko')
xlabel('Model #');ylim([0 1]);
xlim([1 length(mdl_BM)]);
ylabel('Proportion of variance explained');
title('R squared');

%Plot fit statistics
subplot(1,2,2)
for model = 1:length(mdl_BM)
   AIC(model) = mdl_BM{1,model}.ModelCriterion.AIC; 
end
plot(AIC,'-ko')
xlabel('Model #');
xlim([1 length(mdl_BM)]);
title('AIC');
ylabel('AIC');

suptitle('Model comparison for bump magnitude');

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\model_comparison_bump_mag.png']);


for y = 1:numel(mdl_BM{1}.CoefficientNames)
    fprintf('Contrast position %i: %s\n', y, char(mdl_BM{1}.CoefficientNames{y}));
end
bonfCorrection = 1;
[pVal F df1 df2] = coefTest(mdl_BM{1}, [0  1 -1 0]);   pVal = min(1,pVal*bonfCorrection);           fprintf('Planned comparison low conrast vs high contrast,  F(%i, %i)=%0.2f,\t p=%0.7f (contrast 0  1 -1)\n', df1, df2, F, pVal);


%% Get and plot mean bump magnitude

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per fly
mean_bump_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','mean_bump_mag',...
       'GroupingVariables',{'contrast','Fly'});
%Plot
AllFlies = [];
for fly = 1:length(data)
    bump_data{fly} = [mean_bump_data_per_fly.mean_mean_bump_mag(fly),mean_bump_data_per_fly.mean_mean_bump_mag(fly+2*length(data)),mean_bump_data_per_fly.mean_mean_bump_mag(fly+length(data))];
    AllFlies = [AllFlies;bump_data{fly}];
end
plot(1:3,AllFlies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
hold on
errorbar(1:3,mean(AllFlies),std(AllFlies)/sqrt(length(allFlies)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)

xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel({'Mean bump magnitude';'(DF/F)'},'FontSize',12);
ylim([min(mean_bump_data_per_fly.mean_mean_bump_mag)-0.3 max(mean_bump_data_per_fly.mean_mean_bump_mag)+0.3]);

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\mean_closed_loop_bump_mag.png']);

%% Repeat excluding timepoints with the fly standing still, and poor goodness of fit


%% Plot mean bump magnitude vs movement parameters, parsed by contrast

all_contrast_levels = str2num(char(allModelData.ContrastLevel));
allBumpMag = allModelData.BumpMagnitude;
nbins = 20;

figure('Position',[200 200 1400 600]),

%Forward velocity
subplot(1,4,1)

%Define bin limits
ZForVel = allModelData.ZscoredForVel;
maxBinFV = prctile(ZForVel,97.5); %upper limit
minBinFV = prctile(ZForVel,2.5);
binWidthFV = (maxBinFV-minBinFV)/nbins;
forVelBins = [minBinFV:binWidthFV:maxBinFV];

%Create axes for plot, centering them in the middle of the bin
forVelAxes = forVelBins-binWidthFV/2;
forVelAxes = forVelAxes(2:end);

color_gradient = {[0,0,0],[0 0 0.6],[ 0.5 0.8 0.9]};
%Get binned means
for contrast = 1:3
    for bin = 1:length(forVelBins)-1
        doubleBin(bin,contrast) = nanmean(allBumpMag((ZForVel(1:length(allBumpMag)) > forVelBins(bin)) & (ZForVel(1:length(allBumpMag)) < forVelBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(forVelAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump magnitude'); xlabel('Zscored forward velocity (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinFV-binWidthFV/2 maxBinFV+binWidthFV/2]);
legend('Darkness','Low contrast','High contrast');


%Side speed
subplot(1,4,2)

%Define bin limits
ZSideSpeed = allModelData.ZscoredSideSpeed;
maxBinSS = prctile(ZSideSpeed,97.5); %upper limit
minBinSS = prctile(ZSideSpeed,2.5);
binWidthSS = (maxBinSS-minBinSS)/nbins;
sideSpeedBins = [minBinSS:binWidthSS:maxBinSS];

%Create axes for plot, centering them in the middle of the bin
sideSpeedAxes = sideSpeedBins-binWidthSS/2;
sideSpeedAxes = sideSpeedAxes(2:end);

%Get binned means
for contrast = 1:3
    for bin = 1:length(sideSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allBumpMag((ZSideSpeed(1:length(allBumpMag)) > sideSpeedBins(bin)) & (ZSideSpeed(1:length(allBumpMag)) < sideSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(sideSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump magnitude'); xlabel('Zscored side speed (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinSS-binWidthSS/2 maxBinSS+binWidthSS/2]);
legend('Darkness','Low contrast','High contrast');


%Yaw speed
subplot(1,4,3)

%Define bin limits
ZYawSpeed = allModelData.ZscoredYawSpeed;
maxBinYS = prctile(ZYawSpeed,97.5); %upper limit
minBinYS = prctile(ZYawSpeed,2.5);
binWidthYS = (maxBinYS-minBinYS)/nbins;
yawSpeedBins = [minBinYS:binWidthYS:maxBinYS];

%Create axes for plot, centering them in the middle of the bin
yawSpeedAxes = yawSpeedBins-binWidthYS/2;
yawSpeedAxes = yawSpeedAxes(2:end);

%Get binned means
for contrast = 1:3
    for bin = 1:length(yawSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allBumpMag((ZYawSpeed(1:length(allBumpMag)) > yawSpeedBins(bin)) & (ZYawSpeed(1:length(allBumpMag)) < yawSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(yawSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump magnitude'); xlabel('Zscored yaw speed (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinYS-binWidthYS/2 maxBinYS+binWidthYS/2]);
legend('Darkness','Low contrast','High contrast');

%Total movement
subplot(1,4,4)

%Define bin limits
ZTotalMvt = allModelData.ZscoredTotalMovement;
maxBinTM = prctile(ZTotalMvt,97.5); %upper limit
minBinTM = prctile(ZTotalMvt,2.5);
binWidthTM = (maxBinTM-minBinTM)/nbins;
totalMvtBins = [minBinTM:binWidthTM:maxBinTM];

%Create axes for plot, centering them in the middle of the bin
totalMvtAxes = totalMvtBins-binWidthTM/2;
totalMvtAxes = totalMvtAxes(2:end);

%Get binned means
for contrast = 1:3
    for bin = 1:length(totalMvtBins)-1
        doubleBin(bin,contrast) = nanmean(allBumpMag((ZTotalMvt(1:length(allBumpMag)) > totalMvtBins(bin)) & (ZTotalMvt(1:length(allBumpMag)) < totalMvtBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(totalMvtAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump magnitude'); xlabel('Zscored total movement (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinTM-binWidthTM/2 maxBinTM+binWidthTM/2]);
legend('Darkness','Low contrast','High contrast');

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\bump_mag_vs_mvt_and_contrast.png']);

%% Repeat using zscored BM data

allzBumpMag = allModelData.zscoredBM;

figure('Position',[200 200 1400 600]),

%Forward velocity
subplot(1,4,1)

%Get binned means
for contrast = 1:3
    for bin = 1:length(forVelBins)-1
        doubleBin(bin,contrast) = nanmean(allzBumpMag((ZForVel(1:length(allzBumpMag)) > forVelBins(bin)) & (ZForVel(1:length(allzBumpMag)) < forVelBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(forVelAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump magnitude'); xlabel('Zscored forward velocity (deg/s)');
xlim([minBinFV-binWidthFV/2 maxBinFV+binWidthFV/2]);
legend('Darkness','Low contrast','High contrast');


%Side speed
subplot(1,4,2)

%Get binned means
for contrast = 1:3
    for bin = 1:length(sideSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allzBumpMag((ZSideSpeed(1:length(allzBumpMag)) > sideSpeedBins(bin)) & (ZSideSpeed(1:length(allzBumpMag)) < sideSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(sideSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump magnitude'); xlabel('Zscored side speed (deg/s)');
xlim([minBinSS-binWidthSS/2 maxBinSS+binWidthSS/2]);
legend('Darkness','Low contrast','High contrast');


%Yaw speed
subplot(1,4,3)

%Get binned means
for contrast = 1:3
    for bin = 1:length(yawSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allzBumpMag((ZYawSpeed(1:length(allzBumpMag)) > yawSpeedBins(bin)) & (ZYawSpeed(1:length(allzBumpMag)) < yawSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(yawSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump magnitude'); xlabel('Zscored yaw speed (deg/s)');
xlim([minBinYS-binWidthYS/2 maxBinYS+binWidthYS/2]);
legend('Darkness','Low contrast','High contrast');

%Total movement
subplot(1,4,4)

%Get binned means
for contrast = 1:3
    for bin = 1:length(totalMvtBins)-1
        doubleBin(bin,contrast) = nanmean(allzBumpMag((ZTotalMvt(1:length(allzBumpMag)) > totalMvtBins(bin)) & (ZTotalMvt(1:length(allzBumpMag)) < totalMvtBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(totalMvtAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump magnitude'); xlabel('Zscored total movement (deg/s)');
xlim([minBinTM-binWidthTM/2 maxBinTM+binWidthTM/2]);
legend('Darkness','Low contrast','High contrast');

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\z_bump_mag_vs_mvt_and_contrast.png']);

%% Compute model for bump width at half max

%Fit different models
mdl_new_HW{1} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+(1|Fly)');
mdl_new_HW{2} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+TotalMovement+(1|Fly)');
mdl_new_HW{3} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+ForVelocity+SideSpeed+YawSpeed+(1|Fly)');
mdl_new_HW{4} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+ForVelocity+SideSpeed+YawSpeed+Time+(1|Fly)');
%Zscored data
mdl_new_HW{5} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+(1|Fly)');
mdl_new_HW{6} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+ZscoredTotalMovement+(1|Fly)');
mdl_new_HW{7} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+ZscoredForVel+ZscoredSideSpeed+ZscoredYawSpeed+(1|Fly)');
mdl_new_HW{8} = fitlme(allModelData(adj_rs >= 0.5,:),'NewBumpWidth~ContrastLevel+ZscoredForVel+ZscoredSideSpeed+ZscoredYawSpeed+Time+(1|Fly)');

%Plot model comparison
figure,
%Plot R squared
subplot(1,2,1)
for model = 1:length(mdl_new_HW)
   Rsquared(model) = mdl_new_HW{1,model}.Rsquared.Adjusted; 
end
plot(Rsquared,'-ko')
xlabel('Model #');ylim([0 1]);
xlim([1 length(mdl_new_HW)]);
ylabel('Proportion of variance explained');
title('R squared');

%Plot fit statistics
subplot(1,2,2)
for model = 1:length(mdl_new_HW)
   AIC(model) = mdl_new_HW{1,model}.ModelCriterion.AIC; 
end
plot(AIC,'-ko')
xlabel('Model #');
xlim([1 length(mdl_new_HW)]);
title('AIC');
ylabel('AIC');

suptitle('Model comparison for bump width');

%save figure
saveas(gcf,[path,'\continuousGroupPlots\modelComparisonNewBW.png']);


%% Get and plot mean bump width at half max

figure('Position',[200 200 1000 800]),
mean_new_hw_per_fly = varfun(@mean,allSummaryData,'InputVariables','mean_new_half_width',...
       'GroupingVariables',{'contrast','Fly'});
%Plot
AllFliesHW = [];
for fly = 1:length(data)
    new_hw_bump_data{fly} = [mean_new_hw_per_fly.mean_mean_new_half_width(fly),mean_new_hw_per_fly.mean_mean_new_half_width(fly+2*length(data)),mean_new_hw_per_fly.mean_mean_new_half_width(fly+length(data))];
    AllFliesHW = [AllFliesHW;new_hw_bump_data{fly}];
end
plot(1:3,AllFliesHW','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
hold on
errorbar(1:3,mean(AllFliesHW),std(AllFliesHW)/sqrt(length(AllFliesHW)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Mean bump width at half max','FontSize',12);
ylim([0 max(mean_new_hw_per_fly.mean_mean_new_half_width)+0.8]);

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\mean_closed_loop_new_half_width.png']);

%% Plot mean bump width vs movement parameters, parsed by contrast

allNewBumpWidth = allModelData.NewBumpWidth;

figure('Position',[200 200 1400 600]),

%Forward velocity
subplot(1,4,1)

%Get binned means
for contrast = 1:3
    for bin = 1:length(forVelBins)-1
        doubleBin(bin,contrast) = nanmean(allNewBumpWidth((ZForVel(1:length(allNewBumpWidth)) > forVelBins(bin)) & (ZForVel(1:length(allNewBumpWidth)) < forVelBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(forVelAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump width'); xlabel('Zscored forward velocity (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinFV-binWidthFV/2 maxBinFV+binWidthFV/2]);
legend('Darkness','Low contrast','High contrast');


%Side speed
subplot(1,4,2)

%Get binned means
for contrast = 1:3
    for bin = 1:length(sideSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allNewBumpWidth((ZSideSpeed(1:length(allNewBumpWidth)) > sideSpeedBins(bin)) & (ZSideSpeed(1:length(allNewBumpWidth)) < sideSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(sideSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump width'); xlabel('Zscored side speed (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinSS-binWidthSS/2 maxBinSS+binWidthSS/2]);
legend('Darkness','Low contrast','High contrast');


%Yaw speed
subplot(1,4,3)

%Get binned means
for contrast = 1:3
    for bin = 1:length(yawSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allNewBumpWidth((ZYawSpeed(1:length(allNewBumpWidth)) > yawSpeedBins(bin)) & (ZYawSpeed(1:length(allNewBumpWidth)) < yawSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(yawSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump width'); xlabel('Zscored yaw speed (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinYS-binWidthYS/2 maxBinYS+binWidthYS/2]);
legend('Darkness','Low contrast','High contrast');

%Total movement
subplot(1,4,4)

%Get binned means
for contrast = 1:3
    for bin = 1:length(totalMvtBins)-1
        doubleBin(bin,contrast) = nanmean(allNewBumpWidth((ZTotalMvt(1:length(allNewBumpWidth)) > totalMvtBins(bin)) & (ZTotalMvt(1:length(allNewBumpWidth)) < totalMvtBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(totalMvtAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean bump width'); xlabel('Zscored total movement (deg/s)');
ylim([0 (max(max(doubleBin))+0.5)]);
xlim([minBinTM-binWidthTM/2 maxBinTM+binWidthTM/2]);
legend('Darkness','Low contrast','High contrast');

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\NewBumpWidth_vs_mvt_and_contrast.png']);

%% Repeat using zscored data

allzNewBumpWidth = allModelData.zscoredNewBW;

figure('Position',[200 200 1400 600]),

%Forward velocity
subplot(1,4,1)

%Get binned means
for contrast = 1:3
    for bin = 1:length(forVelBins)-1
        doubleBin(bin,contrast) = nanmean(allzNewBumpWidth((ZForVel(1:length(allzNewBumpWidth)) > forVelBins(bin)) & (ZForVel(1:length(allzNewBumpWidth)) < forVelBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(forVelAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump width'); xlabel('Zscored forward velocity (deg/s)');
xlim([minBinFV-binWidthFV/2 maxBinFV+binWidthFV/2]);
legend('Darkness','Low contrast','High contrast');


%Side speed
subplot(1,4,2)

%Get binned means
for contrast = 1:3
    for bin = 1:length(sideSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allzNewBumpWidth((ZSideSpeed(1:length(allzNewBumpWidth)) > sideSpeedBins(bin)) & (ZSideSpeed(1:length(allzNewBumpWidth)) < sideSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(sideSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump width'); xlabel('Zscored side speed (deg/s)');
xlim([minBinSS-binWidthSS/2 maxBinSS+binWidthSS/2]);
legend('Darkness','Low contrast','High contrast');


%Yaw speed
subplot(1,4,3)

%Get binned means
for contrast = 1:3
    for bin = 1:length(yawSpeedBins)-1
        doubleBin(bin,contrast) = nanmean(allzNewBumpWidth((ZYawSpeed(1:length(allzNewBumpWidth)) > yawSpeedBins(bin)) & (ZYawSpeed(1:length(allzNewBumpWidth)) < yawSpeedBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(yawSpeedAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump width'); xlabel('Zscored yaw speed (deg/s)');
xlim([minBinYS-binWidthYS/2 maxBinYS+binWidthYS/2]);
legend('Darkness','Low contrast','High contrast');

%Total movement
subplot(1,4,4)

%Get binned means
for contrast = 1:3
    for bin = 1:length(totalMvtBins)-1
        doubleBin(bin,contrast) = nanmean(allzNewBumpWidth((ZTotalMvt(1:length(allzNewBumpWidth)) > totalMvtBins(bin)) & (ZTotalMvt(1:length(allzNewBumpWidth)) < totalMvtBins(bin+1)) & (adj_rs >= 0.5) & (all_contrast_levels == contrast)));
    end
    plot(totalMvtAxes,doubleBin(:,contrast),'-o','color',color_gradient{contrast})
    hold on
end
ylabel('Mean zscored bump width'); xlabel('Zscored total movement (deg/s)');
xlim([minBinTM-binWidthTM/2 maxBinTM+binWidthTM/2]);
legend('Darkness','Low contrast','High contrast');

%Save figure
saveas(gcf,[path,'\continuousGroupPlots\zNewBumpWidth_vs_mvt_and_contrast.png']);

%% Get and plot total movement per contrast

figure('Position',[100 100 1400 800]),

subplot(1,2,1)
mean_total_mvt_per_fly = varfun(@mean,allModelData,'InputVariables','TotalMovement',...
       'GroupingVariables',{'ContrastLevel','Fly'});
AllFliesTM = [];
for fly = 1:length(data)
    total_mvt_data{fly} = [mean_total_mvt_per_fly.mean_TotalMovement(fly),mean_total_mvt_per_fly.mean_TotalMovement(fly+2*length(data)),mean_total_mvt_per_fly.mean_TotalMovement(fly+length(data))];
    AllFliesTM = [AllFliesTM;total_mvt_data{fly}];
end
plot(1:3,AllFliesTM','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
hold on
errorbar(1:3,mean(AllFliesTM),std(AllFliesTM)/sqrt(length(AllFliesTM)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Total movement (deg/s)','FontSize',12);
ylim([0 max(mean_total_mvt_per_fly.mean_TotalMovement)+10]);

subplot(1,2,2)
%zscored movement
mean_ztotal_mvt_per_fly = varfun(@mean,allModelData,'InputVariables','ZscoredTotalMovement',...
       'GroupingVariables',{'ContrastLevel','Fly'});
AllFlieszTM = [];
for fly = 1:length(data)
    ztotal_mvt_data{fly} = [mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement(fly),mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement(fly+2*length(data)),mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement(fly+length(data))];
    AllFlieszTM = [AllFlieszTM;ztotal_mvt_data{fly}];
end
plot(1:3,AllFlieszTM','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
hold on
errorbar(1:3,mean(AllFlieszTM),std(AllFlieszTM)/sqrt(length(AllFlieszTM)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Zscored total movement (deg/s)','FontSize',12);
ylim([min(mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement)-1 max(mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement)+1]);

%save figure
saveas(gcf,[path,'\continuousGroupPlots\mean_closed_loop_total_mvt.png']);

%%
clear all; close all; clc