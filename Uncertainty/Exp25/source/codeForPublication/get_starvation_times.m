%Code to see starvation time for each fly in this experiment:

clear all; close all;

%% Load data

load('Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp25\data\Experimental\two_ND_filters_3_contrasts\starvation_times.mat')

%% Plot

figure,
boxplot(starvation_times,'color','k')
hold on
yline(0);
scatter(repelem(1,length(starvation_times),1),starvation_times,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Starvation time');

max(starvation_times)
min(starvation_times)