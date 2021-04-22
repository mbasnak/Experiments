%Code to plot the coefficients of the open-loop model for bump velocity

clear all; close all;

%% Load data

%Get directory you're interested in
folderNames = dir('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts');

for content = 1:length(folderNames)
    if contains(folderNames(content).name,'60D05')
        flyData{content} = [folderNames(content).folder,'\',folderNames(content).name];
    end
end

%remove empty cells
data = flyData(~cellfun(@isempty,flyData));

%% Load data for 'good' flies for the model

%load data on flies to include
load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts\flies_for_model');

for fly = 1:length(data)
    if any(flies_for_model==fly)
        modelData{fly} = load([data{fly},'\analysis\open_loop_data.mat']);
    end
end

%remove empty cells
modelData = modelData(~cellfun(@isempty,modelData));


%% Combine data from weights

beta_visual = [];
beta_internal = [];

for fly = 1:length(modelData)
   %store the low contrast weigths in the first column
   beta_visual(fly,1) = modelData{1,fly}.bump_vel_model_LC.Coefficients.Estimate(2);
   beta_visual_error(fly,1) = modelData{1,fly}.bump_vel_model_LC.Coefficients.SE(2);
   beta_internal(fly,1) = modelData{1,fly}.bump_vel_model_LC.Coefficients.Estimate(1);
   beta_internal_error(fly,1) = modelData{1,fly}.bump_vel_model_LC.Coefficients.SE(1);
   
   %store the high contrast weigths in the second column
   beta_visual(fly,2) = modelData{1,fly}.bump_vel_model_HC.Coefficients.Estimate(2);
   beta_visual_error(fly,2) = modelData{1,fly}.bump_vel_model_HC.Coefficients.SE(2);
   beta_internal(fly,2) = modelData{1,fly}.bump_vel_model_HC.Coefficients.Estimate(1); 
   beta_internal_error(fly,2) = modelData{1,fly}.bump_vel_model_HC.Coefficients.SE(1);
   
   %store the goodness of fit metrics for both models to compare
   RMSE(fly,1) = modelData{1,fly}.bump_vel_model_LC.RMSE;
   RMSE(fly,2) = modelData{1,fly}.bump_vel_model_HC.RMSE;
end


%% Run model to see if the differences are statistically significant

contrast = [repelem(1,length(beta_visual)),repelem(2,length(beta_visual))];
Fly = [1:11,1:11];

all_beta_visual = [beta_visual(:,1);beta_visual(:,2)];
visual_data = [all_beta_visual,contrast',Fly'];
model_visual_data = array2table(visual_data, 'VariableNames',{'beta_visual','contrast','fly'});

mdl_visual = fitlme(model_visual_data,'beta_visual~contrast+(1|fly)')


all_beta_internal = [beta_internal(:,1);beta_internal(:,2)];
internal_data = [all_beta_internal,contrast',Fly'];
model_internal_data = array2table(internal_data, 'VariableNames',{'beta_internal','contrast','fly'});

mdl_internal = fitlme(model_internal_data,'beta_internal~contrast+(1|fly)')

%% Add mean plot

mean_beta_internal = mean(beta_internal);
mean_beta_visual = mean(beta_visual);

sd_beta_internal = std(beta_internal);
sd_beta_visual = std(beta_visual);

figure('Position',[200 200 1400 800]),
subplot(1,2,1)
for fly = 1:length(modelData)
    plot([1:2],[beta_visual(fly),beta_visual(fly+length(modelData))],'color',[0.5 0.5 0.5])
    hold on
end
errorbar(mean_beta_visual,sd_beta_visual/sqrt(length(modelData)),'k','lineWidth',3)
if mdl_visual.Coefficients.pValue(2)<0.005
    text(1.4,0.4,'***','fontsize',26);
else
    text(1.4,0.4,'ns');
end
xlim([0 3]);
ylabel('Beta visual','fontSize',14);
ylim([0 0.6]);
xticks(1:2);
xticklabels({'Low contrast','High contrast'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
for fly = 1:length(modelData)
    plot([1:2],[beta_internal(fly),beta_internal(fly+length(modelData))],'color',[0.5 0.5 0.5])
    hold on
end
errorbar(mean_beta_internal,sd_beta_internal/sqrt(length(modelData)),'k','LineWidth',3)
if mdl_internal.Coefficients.pValue(2)<0.005
    text(1.4,0.4,'***','fontsize',26);
else
    text(1.4,0.4,'ns');
end
xlim([0 3]);
ylabel('Beta internal','fontSize',14);
ylim([0 0.6]);
xticks(1:2);
xticklabels({'Low contrast','High contrast'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
%save figure
path = folderNames(1).folder;
saveas(gcf,[path,'\globalPlots\mean_bump_vel_model_OL.png']);


%% Compare goodness of fit of models

figure,
plot(RMSE','color',[0.5 0.5 0.5]);
%add mean
hold on
errorbar(1:2,mean(RMSE),std(RMSE)/sqrt(length(RMSE)),'k','linewidth',2)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

xlim([0 3]);
xticks(1:2);
xticklabels({'Low contrast','High contrast'});
ylim([0 100]);
ylabel('RMSE (fit error)')

saveas(gcf,[path,'\globalPlots\fit_error.png']);

%% Run the model with all the data
%for the previous plots, I had actually run the model on each individual
%fly and recovered the beta parameters from each.
%I will now instead pool all the data to run a unique model

all_LC_data = array2table(zeros(0,5));
all_LC_data.Properties.VariableNames = {'FlyAngVel','VisualCueDrive','BumpAngVel','contrast','fly_ID'};
all_LC_data.contrast = nominal(all_LC_data.contrast);
all_HC_data = array2table(zeros(0,5));
all_HC_data.Properties.VariableNames = {'FlyAngVel','VisualCueDrive','BumpAngVel','contrast','fly_ID'};
all_HC_data.contrast = nominal(all_LC_data.contrast);
all_data = array2table(zeros(0,5));
all_data.Properties.VariableNames = {'FlyAngVel','VisualCueDrive','BumpAngVel','contrast','fly_ID'};
all_data.contrast = nominal(all_data.contrast);

for fly = 1:length(modelData)
   %get the low and high contrast data, and add 1 column for the contrast
   %level, and 1 column for the fly ID
   LC_data = [modelData{1,fly}.bump_vel_data_LC];
   fly_ID = repelem(fly,size(modelData{1,fly}.bump_vel_data_LC,1),1);
   LC_data = addvars(LC_data,fly_ID);
   %combine all the low contrast data
   all_LC_data = [all_LC_data;LC_data];
   
   %repeat all for high contrast
   HC_data = [modelData{1,fly}.bump_vel_data_HC];
   fly_ID = repelem(fly,size(modelData{1,fly}.bump_vel_data_HC,1),1);
   HC_data = addvars(HC_data,fly_ID);
   %combine all the high contrast data
   all_HC_data = [all_HC_data;HC_data];
      
end

%combine all the data
all_data = [all_LC_data;all_HC_data];
%convert contrast level and fly id to categorical
all_data.contrast = nominal(all_data.contrast);
all_data.fly_ID = nominal(all_data.fly_ID);

%fit model
bump_vel_model = fitlme(all_data,'BumpAngVel ~ -1 + VisualCueDrive:contrast + FlyAngVel:contrast + VisualCueDrive + FlyAngVel + (1|fly_ID)')
