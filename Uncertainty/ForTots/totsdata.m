%Analysis for Tot's data

clear all; close all;


%% Fly 1
% 1 day old fly, cold collected, starved overnight, 60D05_7f

% Before I ran this experiment, I ran two 10 min menotaxis session with a
% bright bar. She walked very little in the 1st session, so I didn't analyze it. In the second she walked a lot better.
% The corresponding heatmaps and offset plots for this session can be found in
% the 'fly1' folder and are called 'heatmap_menotaxis_trial' and
% 'offset_dist_menotaxis_trial' respectively.

%I then ran 2 sessions of the 'three pattern experiment' that I was running.
%The plots for one of these sessions are also in the folder.

%Finally, I lifted the fly off the ball and ran this experiment, in which
%the animal got a vertical stripe at a fixed location and then the panels
%turned off, before turning on again, and off one final time (the switch here happens every 120 sec). Next is the
%code to analyze this last part

load('Z:\Wilson lab\Tots\ForTots\fly1\20200211_analysis_sid_4_tid_0.mat')


%get the times of stim change
stimChange = [1;diff(data.fr_y_ds)];
stimChangePoints = find(abs(stimChange)>1);

figure('Position',[100,300, 1600, 800]),
subplot(2,1,1)
imagesc(data.dff_matrix)
hold on
for i = 1:length(stimChangePoints)
    xline(stimChangePoints(i),'lineWidth',2,'color','r')
end
ylabel('PB glomerulus');
colormap(bone)
set(gca,'xticklabel',{[]})

%averaging activity between the two PB halves
leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];

combinedDFF = (data.dff_matrix(leftPB,:) + data.dff_matrix(rightPB,:))/2;

subplot(2,1,2)
imagesc(combinedDFF)
hold on
for i = 1:length(stimChangePoints)
    xline(stimChangePoints(i),'lineWidth',2,'color','r')
end
ylabel('PB glomerulus');
colormap(bone)
set(gca,'xticklabel',{[]})
xlabel('Time');

saveas(gcf,('Z:\Wilson Lab\Mel\Experiments\Uncertainty\ForTots\heatmapsFly1.png'));

%% I eliminated fly 2 because from my records she wasn't moving at all

%% Fly 3

clear all;
% 1 day old fly, starved overnight, 60D05_7f

% Before I ran this experiment, I ran two 10 min menotaxis session with a bright bar.
% The corresponding heatmaps and offset plots for one of these sessions can be found in
% the 'fly1' folder and are called 'heatmap_menotaxis_trial' and
% 'offset_dist_menotaxis_trial' respectively.

%I then ran 2 sessions of the 'three pattern experiment' that I was running.
%The plots for one of these sessions are also in the folder.

%Finally, I lifted the fly off the ball and ran this experiment, in which
%I alternate giving the animal 1 bar at a fixed location with panels off (the switch here happens every 120 sec).
%Next is the code to analyze this last part

load('Z:\Wilson Lab\Tots\ForTots\fly3\20200220_analysis_sid_4_tid_0.mat')


%get the times of stim change
stimChange = [1;diff(data.fr_y_ds)];
stimChangePoints = find(abs(stimChange)>1);

figure('Position',[100,300, 1600, 800]),
subplot(2,1,1)
imagesc(data.dff_matrix)
hold on
for i = 1:length(stimChangePoints)
    xline(stimChangePoints(i),'lineWidth',2,'color','r')
end
ylabel('PB glomerulus');
colormap(bone)
set(gca,'xticklabel',{[]})

%averaging activity between the two PB halves
leftPB = [1,2,3,4,5,6,7,8];
rightPB = [10,11,12,13,14,15,16,9];

combinedDFF = (data.dff_matrix(leftPB,:) + data.dff_matrix(rightPB,:))/2;

subplot(2,1,2)
imagesc(combinedDFF)
hold on
for i = 1:length(stimChangePoints)
    xline(stimChangePoints(i),'lineWidth',2,'color','r')
end
ylabel('PB glomerulus');
colormap(bone)
set(gca,'xticklabel',{[]})
xlabel('Time');

saveas(gcf,('Z:\Wilson Lab\Mel\Experiments\Uncertainty\ForTots\heatmapsFly3.png'));
