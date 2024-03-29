%Code for the group analysis of the open-loop bouts

clear all; close all;

%% Load data

path = uigetdir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental');

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\open_loop_data.mat']);
   end
end

%% Clean and combine data

%remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));

%combine the tables into one
allSummaryData = array2table(zeros(0,6),'VariableNames', {'offset_var','bump_mag','half_width','contrast_level','stim_vel','mvt_offset_var'});
flyNumber= [];
for fly = 1:length(data)
    flyNumber = [flyNumber,repelem(fly,length(data(fly).summarydata.contrast_level))];
    allSummaryData = [allSummaryData;data(fly).summarydata]; 
end

allSummaryData = addvars(allSummaryData,flyNumber');
allSummaryData.Properties.VariableNames{'Var7'} = 'Fly';


%% Compute model for offset variation

%fit mixed linear model using fly number as a random variable
mdl_offset = fitlme(allSummaryData,'offset_var~contrast_level+(1|Fly)')

%% Get and plot mean offset variation per fly and in total

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per fly
mean_offset_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast_level','Fly'});
%plot(mean_offset_data_per_fly.contrast_level,mean_offset_data_per_fly.mean_offset_var,'o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
%change the code above to be able to get lines per fly
allFlies = [];
for fly = 1:length(data)
    fly_data{fly} = [mean_offset_data_per_fly.mean_offset_var(fly),mean_offset_data_per_fly.mean_offset_var(fly+length(data))];
    allFlies = [allFlies;fly_data{fly}];
end
plot(56:57,allFlies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_offset_data = varfun(@mean,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast_level'});
hold on
plot(56:57,mean(allFlies),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
if (mdl_offset.Coefficients.pValue(2)<0.05 & mdl_offset.Coefficients.pValue(2)>=0.01)
    text(56.4,1.8,'*','fontsize',26);
elseif (mdl_offset.Coefficients.pValue(2)<0.01 & mdl_offset.Coefficients.pValue(2)>=0.001)
    text(56.4,1.8,'**','fontsize',26);
elseif mdl_offset.Coefficients.pValue(2)<0.001
    text(56.4,1.8,'***','fontsize',26);
else
    text(56.4,1.8,'ns');
end

xlim([55 58]);
xticks(56:57);
xticklabels({'Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Circular standard deviation of offset','FontSize',12);
ylim([min(mean_offset_data_per_fly.mean_offset_var)-0.3 max(mean_offset_data_per_fly.mean_offset_var)+0.3]);


%save figure
saveas(gcf,[path,'\globalPlots\mean_open_loop_stim_offset_var.png']);


%% Offset variation per speed per contrast

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per stim speed
mean_offset_data_per_speed = varfun(@mean,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast_level','stim_vel'});
plot(mean_offset_data_per_speed.stim_vel(mean_offset_data_per_speed.contrast_level==56),mean_offset_data_per_speed.mean_offset_var(mean_offset_data_per_speed.contrast_level==56),'-o','lineWidth',2,'color',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',8)
%add lines per fly
mean_offset_data_per_speed_and_fly = varfun(@mean,allSummaryData,'InputVariables','offset_var',...
       'GroupingVariables',{'contrast_level','stim_vel','Fly'});
allFlies_data{1} = [];
allFlies_data{2} = [];
for fly = 1:length(data)
    fly_offset_data{fly} = [mean_offset_data_per_speed_and_fly.mean_offset_var(mean_offset_data_per_speed_and_fly.contrast_level==56 & mean_offset_data_per_speed_and_fly.Fly==fly),mean_offset_data_per_speed_and_fly.mean_offset_var(mean_offset_data_per_speed_and_fly.contrast_level==57 & mean_offset_data_per_speed_and_fly.Fly==fly)];
    allFlies_data{1} = [allFlies_data{1},fly_offset_data{fly}(:,1)]; 
    allFlies_data{2} = [allFlies_data{1},fly_offset_data{fly}(:,2)]; 
end

hold on
plot(mean_offset_data_per_speed.stim_vel(mean_offset_data_per_speed.contrast_level==57),mean_offset_data_per_speed.mean_offset_var(mean_offset_data_per_speed.contrast_level==57),'-o','lineWidth',2,'color',[ 0.5 0.8 0.9],'MarkerFaceColor',[ 0.5 0.8 0.9],'MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'20 deg/s','30 deg/s','60 deg/s'});
%add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
a = get(gca,'XTickLabel');  
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Circular standard deviation of offset','FontSize',12);
ylim([min(mean_offset_data_per_speed.mean_offset_var)-0.3 max(mean_offset_data_per_speed.mean_offset_var)+0.3]);

%save figure
saveas(gcf,[path,'\globalPlots\open_loop_offset_var_per_speed.png']);

%% Compute model for bump magnitude

%combine the tables into one
allModelData = array2table(zeros(0,5),'VariableNames', {'ContrastLevel','Time','TotalMovement','ZscoredMvt','BumpMagnitude'});
Fly = [];
for fly = 1:length(data)
    Fly = [Fly,repelem(fly,length(data(fly).bump_mag_data.ContrastLevel))];
    allModelData = [allModelData;data(fly).bump_mag_data]; 
end
allModelData = addvars(allModelData,Fly');
allModelData.Properties.VariableNames{'Var6'} = 'Fly';

%fit mixed linear model using fly number as a random variable
mdl_BM = fitlme(allModelData,'BumpMagnitude~ContrastLevel+TotalMovement+Time+(1|Fly)');

%force the intercept to be 0
mdl_BM2 = fitlme(allModelData,'BumpMagnitude~-1+ContrastLevel+TotalMovement+Time+(1|Fly)');

%removing the time as variable
mdl_BM3 = fitlme(allModelData,'BumpMagnitude~-1+ContrastLevel+TotalMovement+(1|Fly)');

%using zscored mvt data
%1) per fly
mdl_BM4 = fitlme(allModelData,'BumpMagnitude~-1+ContrastLevel+ZscoredMvt+(1|Fly)');

%2) across flies
zscored_mvt_across = zscore(allModelData.TotalMovement);
allModelData = addvars(allModelData,zscored_mvt_across,'NewVariableNames','ZscoredMovementAcrossFlies');
mdl_BM5 = fitlme(allModelData,'BumpMagnitude~ContrastLevel+ZscoredMovementAcrossFlies+(1|Fly)')

%% Get and plot mean bump magnitude

figure('Position',[200 200 1000 800]),
%Get mean bump mag by contrast per fly
mean_bump_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','bump_mag',...
       'GroupingVariables',{'contrast_level','Fly'});
%Plot
%plot(mean_bump_data_per_fly.contrast_level,mean_bump_data_per_fly.mean_bump_mag,'o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
%change the code above to be able to get lines per fly
AllFlies = [];
for fly = 1:length(data)
    bump_data{fly} = [mean_bump_data_per_fly.mean_bump_mag(fly),mean_bump_data_per_fly.mean_bump_mag(fly+length(data))];
    AllFlies = [AllFlies;bump_data{fly}];
end
plot(56:57,AllFlies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_bump_data = varfun(@mean,allSummaryData,'InputVariables','bump_mag',...
    'GroupingVariables',{'contrast_level'});
hold on
errorbar(56:57,mean(AllFlies),std(AllFlies)/sqrt(length(AllFlies)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
if (mdl_BM5.Coefficients.pValue(2)<0.05 & mdl_BM5.Coefficients.pValue(2)>=0.01)
    text(56.5,1,'*','fontsize',26);
elseif (mdl_BM5.Coefficients.pValue(2)<0.01 & mdl_BM5.Coefficients.pValue(2)>=0.001)
    text(56.5,1,'**','fontsize',26);
elseif mdl_BM5.Coefficients.pValue(2)<0.001
    text(56.5,1,'***','fontsize',26);
else
    text(56.5,1,'ns');
end
xlim([55 58]);
xticks(56:57);
xticklabels({'Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel({'Mean bump magnitude';'(amplitude of Fourier component)'},'FontSize',12);
ylim([min(mean_bump_data_per_fly.mean_bump_mag)-0.3 max(mean_bump_data_per_fly.mean_bump_mag)+0.3]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_open_loop_bump_mag.png']);


%% Bump magnitude per speed per contrast

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per stim speed
mean_bump_mag_data_per_speed = varfun(@mean,allSummaryData,'InputVariables','bump_mag',...
       'GroupingVariables',{'contrast_level','stim_vel'});
plot(mean_bump_mag_data_per_speed.stim_vel(mean_bump_mag_data_per_speed.contrast_level==56),mean_bump_mag_data_per_speed.mean_bump_mag(mean_bump_mag_data_per_speed.contrast_level==56),'-o','lineWidth',2,'color',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',8)
%add lines per fly
mean_bump_mag_data_per_speed_and_fly = varfun(@mean,allSummaryData,'InputVariables','bump_mag',...
       'GroupingVariables',{'contrast_level','stim_vel','Fly'});
hold on
plot(mean_bump_mag_data_per_speed.stim_vel(mean_bump_mag_data_per_speed.contrast_level==57),mean_bump_mag_data_per_speed.mean_bump_mag(mean_bump_mag_data_per_speed.contrast_level==57),'-o','lineWidth',2,'color',[ 0.5 0.8 0.9],'MarkerFaceColor',[ 0.5 0.8 0.9],'MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'20 deg/s','30 deg/s','60 deg/s'});
%add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
a = get(gca,'XTickLabel');  
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel({'Bump magnitude','(amplitude of Fourier component)'},'FontSize',12);
ylim([min(mean_bump_mag_data_per_speed.mean_bump_mag)-0.3 max(mean_bump_mag_data_per_speed.mean_bump_mag)+0.3]);
saveas(gcf,[path,'\globalPlots\bump_mag_per_speed.png']);


%% Compute model for half width

%combine the tables into one
allModelDataHW = array2table(zeros(0,5),'VariableNames', {'ContrastLevel','Time','TotalMovement','ZscoredMvt','HalfWidth'});
Fly = [];
for fly = 1:length(data)
    Fly = [Fly,repelem(fly,length(data(fly).half_width_data.ContrastLevel))];
    allModelDataHW = [allModelDataHW;data(fly).half_width_data]; 
end
allModelDataHW = addvars(allModelDataHW,Fly');
allModelDataHW.Properties.VariableNames{'Var6'} = 'Fly';

%fit mixed linear model using fly number as a random variable
mdl_HW = fitlme(allModelDataHW,'HalfWidth~ContrastLevel+TotalMovement+Time+(1|Fly)');

%force the intercept to be 0
mdl_HW2 = fitlme(allModelDataHW,'HalfWidth~-1+ContrastLevel+TotalMovement+Time+(1|Fly)');

%removing the time as variable
mdl_HW3 = fitlme(allModelDataHW,'HalfWidth~-1+ContrastLevel+TotalMovement+(1|Fly)');

%using zscored mvt data
%1) per fly
mdl_HW4 = fitlme(allModelDataHW,'HalfWidth~-1+ContrastLevel+ZscoredMvt+Time+(1|Fly)');

%2) across flies
zscored_mvt_across = zscore(allModelDataHW.TotalMovement);
allModelDataHW = addvars(allModelDataHW,zscored_mvt_across,'NewVariableNames','ZscoredMovementAcrossFlies');
mdl_HW5 = fitlme(allModelDataHW,'HalfWidth~ContrastLevel+ZscoredMovementAcrossFlies+(1|Fly)')

%% Get and plot mean half width

figure('Position',[200 200 1000 800]),
%Get mean bump mag by contrast per fly
mean_half_width_data_per_fly = varfun(@mean,allSummaryData,'InputVariables','half_width',...
       'GroupingVariables',{'contrast_level','Fly'});
%Plot
AllFlies = [];
for fly = 1:length(data)
    half_width_data{fly} = [mean_half_width_data_per_fly.mean_half_width(fly),mean_half_width_data_per_fly.mean_half_width(fly+length(data))];
    AllFlies = [AllFlies;half_width_data{fly}];
end
plot(56:57,AllFlies','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])

%Get mean offset var by contrast in total
mean_half_width_data = varfun(@mean,allSummaryData,'InputVariables','half_width',...
       'GroupingVariables',{'contrast_level'});
hold on
errorbar(56:57,mean(AllFlies),std(AllFlies)/sqrt(length(AllFlies)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
if (mdl_HW5.Coefficients.pValue(2)<0.05 & mdl_HW5.Coefficients.pValue(2)>=0.01)
    text(56.5,3.1,'*','fontsize',26);
elseif (mdl_HW5.Coefficients.pValue(2)<0.01 & mdl_HW5.Coefficients.pValue(2)>=0.001)
    text(56.5,3.1,'**','fontsize',26);
elseif mdl_HW5.Coefficients.pValue(2)<0.001
    text(56.5,3.1,'***','fontsize',26);
else
    text(56.5,3.1,'ns');
end
xlim([55 58]);
xticks(56:57);
xticklabels({'Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Mean half width','FontSize',12);
ylim([min(mean_half_width_data_per_fly.mean_half_width)-0.3 max(mean_half_width_data_per_fly.mean_half_width)+0.3]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_open_loop_half_width.png']);


%% Half width per speed per contrast

figure('Position',[200 200 1000 800]),
%Get mean offset var by contrast per stim speed
mean_half_width_data_per_speed = varfun(@mean,allSummaryData,'InputVariables','half_width',...
       'GroupingVariables',{'contrast_level','stim_vel'});
plot(mean_half_width_data_per_speed.stim_vel(mean_half_width_data_per_speed.contrast_level==56),mean_half_width_data_per_speed.mean_half_width(mean_half_width_data_per_speed.contrast_level==56),'-o','lineWidth',2,'color',[0 0 0.6],'MarkerFaceColor',[0 0 0.6],'MarkerSize',8)
%add lines per fly
mean_half_width_data_per_speed_and_fly = varfun(@mean,allSummaryData,'InputVariables','half_width',...
       'GroupingVariables',{'contrast_level','stim_vel','Fly'});
hold on
plot(mean_half_width_data_per_speed.stim_vel(mean_half_width_data_per_speed.contrast_level==57),mean_half_width_data_per_speed.mean_half_width(mean_half_width_data_per_speed.contrast_level==57),'-o','lineWidth',2,'color',[ 0.5 0.8 0.9],'MarkerFaceColor',[ 0.5 0.8 0.9],'MarkerSize',8)
xlim([0 4]);
xticks(1:3);
xticklabels({'20 deg/s','30 deg/s','60 deg/s'});
%add custom legend
h = zeros(2, 1);
h(1) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0 0 0.6]);
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor',[ 0.5 0.8 0.9]);
legend(h, 'Low contrast','High Contrast');
a = get(gca,'XTickLabel');  
[ax,h2]=suplabel('Stimulus angular velocity','x');
set(h2,'FontSize',12,'FontWeight','bold')
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Bump half width','FontSize',12);
ylim([min(mean_half_width_data_per_speed.mean_half_width)-0.3 max(mean_half_width_data_per_speed.mean_half_width)+0.3]);
saveas(gcf,[path,'\globalPlots\half_width_per_speed.png']);

%% Add total movement comparison

figure('Position',[100 100 800 800]),

%Get mean total mvt by contrast per fly
mean_total_mvt_data_per_fly = varfun(@mean,allModelData,'InputVariables','TotalMovement',...
       'GroupingVariables',{'ContrastLevel','Fly'});
%Plot
AllFliesTM = [];
for fly = 1:length(data)
    total_mvt_data{fly} = [mean_total_mvt_data_per_fly.mean_TotalMovement(fly),mean_total_mvt_data_per_fly.mean_TotalMovement(fly+length(data))];
    AllFliesTM = [AllFliesTM;total_mvt_data{fly}];
end
plot(56:57,AllFliesTM','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
hold on
errorbar(56:57,mean(AllFliesTM),std(AllFliesTM)/sqrt(length(AllFliesTM)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
xlim([55 58]);
xticks(56:57);
xticklabels({'Low contrast','High contrast'});
a = get(gca,'XTickLabel');  
set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
ylabel('Mean total movement (deg/s)','FontSize',12);
ylim([0 max(mean_total_mvt_data_per_fly.mean_TotalMovement)+10]);

% subplot(1,2,2)
% %Get mean total mvt by contrast per fly
% mean_ztotal_mvt_data_per_fly = varfun(@mean,allModelData,'InputVariables','ZscoredMvt',...
%        'GroupingVariables',{'ContrastLevel','Fly'});
% %Plot
% AllFlieszTM = [];
% for fly = 1:length(data)
%     ztotal_mvt_data{fly} = [mean_ztotal_mvt_data_per_fly.mean_ZscoredMvt(fly),mean_ztotal_mvt_data_per_fly.mean_ZscoredMvt(fly+length(data))];
%     AllFlieszTM = [AllFlieszTM;ztotal_mvt_data{fly}];
% end
% plot(56:57,AllFlieszTM','-o','color',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6])
% hold on
% errorbar(56:57,mean(AllFlieszTM),std(AllFlieszTM)/sqrt(length(AllFlieszTM)),'-ko','LineWidth',2,'MarkerFaceColor','k','MarkerSize',8)
% xlim([55 58]);
% xticks(56:57);
% xticklabels({'Low contrast','High contrast'});
% a = get(gca,'XTickLabel');  
% set(gca,'XTickLabel',a,'fontsize',12,'FontWeight','bold')
% ylabel('Mean zscored total movement','FontSize',12);
% %ylim([min(mean_ztotal_mvt_data_per_fly.mean_ZscoredMvt)-0.1 max(mean_ztotal_mvt_data_per_fly.mean_ZscoredMvt)+0.1]);

%save figure
saveas(gcf,[path,'\globalPlots\mean_open_loop_total_mvt.png']);