%Code for the group analysis of the closed-loop bouts

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');
%Get folder names
folderContents = dir(path);

%Load the summary data of the folder that correspond to experimental flies
for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\summary_data.mat']);
   end
end

%% Clean and combine data

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%Combine the tables into one
allSummaryData = array2table(zeros(0,5),'VariableNames', {'contrast','offset_var','mean_bump_mag','mean_half_width','heading_var'});
flyNumber= [];
for fly = 1:length(data)
    flyNumber = [flyNumber,repelem(fly,length(data(fly).summary_data.contrast))];
    allSummaryData = [allSummaryData;data(fly).summary_data]; 
end

%Add the fly ID as a variable
allSummaryData = addvars(allSummaryData,nominal(flyNumber'));
allSummaryData.Properties.VariableNames{'Var6'} = 'Fly';

%% Compare offset variation statistically

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
ylabel('Circular standard deviation of offset','FontSize',12);
ylim([min(mean_offset_data_per_fly.mean_offset_var)-0.3 max(mean_offset_data_per_fly.mean_offset_var)+0.3]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_closed_loop_offset_var.png']);

%% Run model for heading variation

%Run model
mdl_heading = fitlme(reordered_data,'heading_var~contrast+(1|Fly)')
disp(anova(mdl_heading))

%Add post-hoc comparison between low and high contrast
for y = 1:numel(mdl_heading.CoefficientNames)
    fprintf('Contrast position %i: %s\n', y, char(mdl_heading.CoefficientNames{y}));
end
bonfCorrection = 1;
[pVal F df1 df2] = coefTest(mdl_heading, [0  1 -1]);   pVal = min(1,pVal*bonfCorrection);           fprintf('Planned comparison low conrast vs high contrast,  F(%i, %i)=%0.2f,\t p=%0.7f (contrast 0  1 -1)\n', df1, df2, F, pVal);


%% Get and plot mean heading variation per fly and in total

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
ylabel('Circular standard deviation of heading','FontSize',12);
ylim([min(mean_heading_data_per_fly.mean_heading_var)-0.3 max(mean_heading_data_per_fly.mean_heading_var)+0.3]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_closed_loop_heading_var.png']);

%% Compute model for bump magnitude

%Combine the tables into one
allModelData = array2table(zeros(0,5),'VariableNames', {'ContrastLevel','TotalMovement','ZscoredTotalMovement','Time','BumpMagnitude'});
Fly = [];
for fly = 1:length(data)
    Fly = [Fly,repelem(fly,length(data(fly).modelTable.ContrastLevel))];
    allModelData = [allModelData;data(fly).modelTable]; 
end
%Add fly id as a categorical variable
allModelData = addvars(allModelData,nominal(Fly'));
allModelData.Properties.VariableNames{'Var6'} = 'Fly';

%Change contrast level as a categorical variable
allModelData.ContrastLevel = nominal(allModelData.ContrastLevel);

%Fit different models
mdl_BM = fitlme(allModelData,'BumpMagnitude~ContrastLevel+TotalMovement+(1|Fly)')
mdl_BM2 = fitlme(allModelData,'BumpMagnitude~ContrastLevel+ZscoredTotalMovement+(1|Fly)')


%Add post-hoc comparison between low and high contrast
for y = 1:numel(mdl_BM2.CoefficientNames)
    fprintf('Contrast position %i: %s\n', y, char(mdl_BM2.CoefficientNames{y}));
end
bonfCorrection = 1;
[pVal F df1 df2] = coefTest(mdl_BM2, [0  1 -1 0]);   pVal = min(1,pVal*bonfCorrection);           fprintf('Planned comparison low conrast vs high contrast,  F(%i, %i)=%0.2f,\t p=%0.7f (contrast 0  1 -1)\n', df1, df2, F, pVal);


%% Get and plot mean bump magnitude

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per fly
mean_bump_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','mean_bump_mag',...
       'GroupingVariables',{'contrast','Fly'});
mean_bump_data_per_fly = addvars(mean_bump_data_per_fly,contrastLevel');
mean_bump_data_per_fly.Properties.VariableNames{'Var5'} = 'contrastLevel';
%Plot
%plot(mean_bump_data_per_fly.contrastLevel,mean_bump_data_per_fly.mean_mean_bump_mag,'o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
%change the code above to be able to get lines per fly
AllFlies = [];
for fly = 1:length(data)
    bump_data{fly} = [mean_bump_data_per_fly.mean_mean_bump_mag(fly),mean_bump_data_per_fly.mean_mean_bump_mag(fly+2*length(data)),mean_bump_data_per_fly.mean_mean_bump_mag(fly+length(data))];
    AllFlies = [AllFlies;bump_data{fly}];
end
plot(1:3,AllFlies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_bump_data = varfun(@mean,allSummaryData,'InputVariables','mean_bump_mag',...
       'GroupingVariables',{'contrast'});
hold on
mean_bump_data = addvars(mean_bump_data,contrastLevel2');
mean_bump_data.Properties.VariableNames{'Var4'} = 'contrastLevel';
%plot(mean_bump_data.contrastLevel,mean_bump_data.mean_mean_bump_mag,'ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
errorbar(1:3,mean(AllFlies),std(AllFlies)/sqrt(length(allFlies)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)

%add significance 
for contrast = 2
   if (mdl_BM6.Coefficients.pValue(contrast)<0.05 & mdl_BM6.Coefficients.pValue(contrast)>=0.01)
      text(contrast-0.58,1,'*','fontsize',26);
   elseif (mdl_BM6.Coefficients.pValue(contrast)<0.01 & mdl_BM6.Coefficients.pValue(contrast)>=0.001)
       text(contrast-0.58,1,'**','fontsize',26);
   elseif mdl_BM6.Coefficients.pValue(contrast)<0.001
       text(contrast-0.58,1,'***','fontsize',26);
   else
      text(contrast-0.58,1,'ns');
   end
end
for contrast = 3
   if (mdl_BM6.Coefficients.pValue(contrast)<0.05 & mdl_BM6.Coefficients.pValue(contrast)>=0.01)
      text(contrast-1.1,1.2,'*','fontsize',26);
   elseif (mdl_BM6.Coefficients.pValue(contrast)<0.01 & mdl_BM6.Coefficients.pValue(contrast)>=0.001)
       text(contrast-1.1,1.2,'**','fontsize',26);
   elseif mdl_BM6.Coefficients.pValue(contrast)<0.001
       text(contrast-1.1,1.2,'***','fontsize',26);
   else
      text(contrast-1.1,1.2,'ns');
   end
end
for contrast = 2
   if (pVal<0.05 & pVal>=0.01)
      text(contrast+0.43,1.1,'*','fontsize',26);
   elseif (pVal<0.01 & pVal>=0.001)
       text(contrast+0.43,1.1,'**','fontsize',26);
   elseif pVal<0.001
       text(contrast+0.43,1.1,'***','fontsize',26);
   elseif (pVal>0.05 & pVal<0.1)
       text(contrast+0.39,1.11,['p = ',num2str(round(pVal,2))],'fontsize',12);
   else
      text(contrast+0.48,1.1,'ns','fontsize',16);
   end
end

%add line between asterisks
%contrasts 1 and 2
line([1 1.4],[1.01 1.01],'color','k');
line([1.66 2],[1.01 1.01],'color','k');
%contrasts 2 and 3
line([2 2.35],[1.1 1.1],'color','k');
line([2.7 3],[1.1 1.1],'color','k');
%contrasts 1 and 3
line([1 1.88],[1.21 1.21],'color','k');
line([2.15 3],[1.21 1.21],'color','k');


xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel({'Mean bump magnitude';'(amplitude of Fourier component'},'FontSize',12);
ylim([min(mean_bump_data_per_fly.mean_mean_bump_mag)-0.3 max(mean_bump_data_per_fly.mean_mean_bump_mag)+0.3]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_closed_loop_bump_mag.png']);



%% Compute model for bump width at half max

%combine the tables into one
allModelDataHW = array2table(zeros(0,5),'VariableNames', {'ContrastLevel','TotalMovement','ZscoredTotalMovement','Time','BumpWidth'});
Fly = [];
for fly = 1:length(data)
    Fly = [Fly,repelem(fly,length(data(fly).modelTable_HW.ContrastLevel))];
    allModelDataHW = [allModelDataHW;data(fly).modelTable_HW]; 
end
allModelDataHW = addvars(allModelDataHW,Fly');
allModelDataHW.Properties.VariableNames{'Var6'} = 'Fly';

%fit mixed linear model using fly number as a random variable
mdl_HW = fitlme(allModelDataHW,'BumpWidth~ContrastLevel+TotalMovement+Time+(1|Fly)');

%force the intercept to be 0
mdl_HW2 = fitlme(allModelDataHW,'BumpWidth~-1+ContrastLevel+TotalMovement+Time+(1|Fly)');

%with contrast level as a categorical variable
allModelDataHWCat = allModelDataHW;
allModelDataHWCat.ContrastLevel = nominal(allModelDataHWCat.ContrastLevel);

mdl_HW3 = fitlme(allModelDataHWCat,'BumpWidth~ContrastLevel+TotalMovement+Time+(1|Fly)');

%using zscored time data
%1) per fly
mdl_HW4 = fitlme(allModelDataHWCat,'BumpWidth~ContrastLevel+ZscoredTotalMovement+Time+(1|Fly)');

%2) across flies
zscored_mvt_across = zscore(allModelDataHW.TotalMovement);
allModelDataHWCat = addvars(allModelDataHWCat,zscored_mvt_across,'NewVariableNames','ZscoredMovementAcrossFlies');

mdl_HW5 = fitlme(allModelDataHWCat,'BumpWidth~ContrastLevel+ZscoredMovementAcrossFlies+(1|Fly)')
%add post-hoc comparison between low and high contrast
for y = 1:numel(mdl_HW5.CoefficientNames)
    fprintf('Contrast position %i: %s\n', y, char(mdl_HW5.CoefficientNames{y}));
end
bonfCorrection = 1;
[pVal F df1 df2] = coefTest(mdl_HW5, [0  1 -1 0]);   pVal = min(1,pVal*bonfCorrection);           fprintf('Planned comparison low conrast vs high contrast,  F(%i, %i)=%0.2f,\t p=%0.7f (contrast 0  1 -1)\n', df1, df2, F, pVal);


%% Get and plot mean bump width at half max

figure('Position',[200 200 1000 800]),
mean_hw_per_fly = varfun(@mean,allSummaryData,'InputVariables','mean_half_width',...
       'GroupingVariables',{'contrast','Fly'});
mean_hw_per_fly = addvars(mean_hw_per_fly,contrastLevel');
mean_hw_per_fly.Properties.VariableNames{'Var5'} = 'contrastLevel';
%Plot
AllFliesHW = [];
for fly = 1:length(data)
    hw_bump_data{fly} = [mean_hw_per_fly.mean_mean_half_width(fly),mean_hw_per_fly.mean_mean_half_width(fly+2*length(data)),mean_hw_per_fly.mean_mean_half_width(fly+length(data))];
    AllFliesHW = [AllFliesHW;hw_bump_data{fly}];
end
plot(1:3,AllFliesHW','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_hw = varfun(@mean,allSummaryData,'InputVariables','mean_half_width',...
       'GroupingVariables',{'contrast'});
hold on
mean_hw = addvars(mean_hw,contrastLevel2');
mean_hw.Properties.VariableNames{'Var4'} = 'contrastLevel';
errorbar(1:3,mean(AllFliesHW),std(AllFliesHW)/sqrt(length(AllFliesHW)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)

%add significance 
for contrast = 2
   if (mdl_HW5.Coefficients.pValue(contrast)<0.05 & mdl_HW5.Coefficients.pValue(contrast)>=0.01)
      text(contrast-0.58,3.5,'*','fontsize',26);
   elseif (mdl_HW5.Coefficients.pValue(contrast)<0.01 & mdl_HW5.Coefficients.pValue(contrast)>=0.001)
       text(contrast-0.58,3.5,'**','fontsize',26);
   elseif mdl_HW5.Coefficients.pValue(contrast)<0.001
       text(contrast-0.58,3.5,'***','fontsize',26);
   else
      text(contrast-0.58,3.5,'ns');
   end
end
for contrast = 3
   if (mdl_HW5.Coefficients.pValue(contrast)<0.05 & mdl_HW5.Coefficients.pValue(contrast)>=0.01)
      text(contrast-1.1,3.9,'*','fontsize',26);
   elseif (mdl_HW5.Coefficients.pValue(contrast)<0.01 & mdl_HW5.Coefficients.pValue(contrast)>=0.001)
       text(contrast-1.1,3.9,'**','fontsize',26);
   elseif mdl_HW5.Coefficients.pValue(contrast)<0.001
       text(contrast-1.1,3.9,'***','fontsize',26);
   else
      text(contrast-1.1,3.9,'ns');
   end
end
for contrast = 2
   if (pVal<0.05 & pVal>=0.01)
      text(contrast+0.43,3.4,'*','fontsize',26);
   elseif (pVal<0.01 & pVal>=0.001)
       text(contrast+0.43,3.4,'**','fontsize',26);
   elseif pVal<0.001
       text(contrast+0.43,3.4,'***','fontsize',26);
   elseif (pVal>0.05 & pVal<0.1)
       text(contrast+0.39,3.4,['p = ',num2str(round(pVal,2))],'fontsize',12);
   else
      text(contrast+0.48,3.4,'ns','fontsize',16);
   end
end

%add line between asterisks
%contrasts 1 and 2
line([1 1.4],[3.52 3.52],'color','k');
line([1.66 2],[3.52 3.52],'color','k');
%contrasts 2 and 3
line([2 2.4],[3.42 3.42],'color','k');
line([2.66 3],[3.42 3.42],'color','k');
%contrasts 1 and 3
line([1 1.88],[3.92 3.92],'color','k');
line([2.15 3],[3.92 3.92],'color','k');

xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Mean bump width at half max','FontSize',12);
ylim([0 max(mean_hw_per_fly.mean_mean_half_width)+0.8]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_closed_loop_half_width.png']);

%% Get and plot total movement per contrast


figure('Position',[100 100 1400 800]),

subplot(1,2,1)
mean_total_mvt_per_fly = varfun(@mean,allModelDataCat,'InputVariables','TotalMovement',...
       'GroupingVariables',{'ContrastLevel','Fly'});
AllFliesTM = [];
for fly = 1:length(data)
    total_mvt_data{fly} = [mean_total_mvt_per_fly.mean_TotalMovement(fly),mean_total_mvt_per_fly.mean_TotalMovement(fly+2*length(data)),mean_total_mvt_per_fly.mean_TotalMovement(fly+length(data))];
    AllFliesTM = [AllFliesTM;total_mvt_data{fly}];
end
plot(1:3,AllFliesTM','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_tm = varfun(@mean,allModelDataCat,'InputVariables','TotalMovement',...
       'GroupingVariables',{'ContrastLevel'});
hold on
mean_tm = addvars(mean_tm,contrastLevel2');
mean_tm.Properties.VariableNames{'Var4'} = 'contrastLevel';
errorbar(1:3,mean(AllFliesTM),std(AllFliesTM)/sqrt(length(AllFliesTM)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('total movement (deg/s)','FontSize',12);
ylim([0 max(mean_total_mvt_per_fly.mean_TotalMovement)+10]);

subplot(1,2,2)
%zscored movement
mean_ztotal_mvt_per_fly = varfun(@mean,allModelDataCat,'InputVariables','ZscoredTotalMovement',...
       'GroupingVariables',{'ContrastLevel','Fly'});
AllFlieszTM = [];
for fly = 1:length(data)
    ztotal_mvt_data{fly} = [mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement(fly),mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement(fly+2*length(data)),mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement(fly+length(data))];
    AllFlieszTM = [AllFlieszTM;ztotal_mvt_data{fly}];
end
plot(1:3,AllFlieszTM','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_ztm = varfun(@mean,allModelDataCat,'InputVariables','ZscoredTotalMovement',...
       'GroupingVariables',{'ContrastLevel'});
hold on
mean_ztm = addvars(mean_ztm,contrastLevel2');
mean_ztm.Properties.VariableNames{'Var4'} = 'contrastLevel';
errorbar(1:3,mean(AllFlieszTM),std(AllFlieszTM)/sqrt(length(AllFlieszTM)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'Darkness','Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('zscored total movement (deg/s)','FontSize',12);
ylim([min(mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement)-1 max(mean_ztotal_mvt_per_fly.mean_ZscoredTotalMovement)+1]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_closed_loop_total_mvt.png']);