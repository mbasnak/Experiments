%Group analysis code for the cue combination experiment

clear all; close all;

%% Load data

path = 'Z:\Wilson Lab\Mel\Experiments\Uncertainty\Exp38\data\third_version';

folderContents = dir(path);

for content = 1:length(folderContents)
   if contains(folderContents(content).name,'60D05')
       data(content) = load([folderContents(content).folder,'\',folderContents(content).name,'\analysis\data.mat']);
   end
end

%Remove empty rows
data = data(all(~cellfun(@isempty,struct2cell(data))));


%% Combine all data

%Initialize variables
pref_index = [];
pref_index_2 = [];
heading_pref_index = [];
stick_index = [];
PI = [];
heading_PI = [];
SI = [];

initial_bar_bm = [];
initial_wind_bm = [];
initial_bar_bw = [];
initial_wind_bw = [];
initial_bar_offset = [];
initial_wind_offset = [];

short_bm_bar_jump = [];
short_bw_bar_jump = [];
short_bm_wind_jump = [];
short_bw_wind_jump = [];

long_bm_bar_jump = [];
long_bw_bar_jump = [];
long_bm_wind_jump = [];
long_bw_wind_jump = [];

short_bm_pref_cue = [];
short_bw_pref_cue = [];
short_bm_non_pref_cue = [];
short_bw_non_pref_cue = [];

bar_offset_variability_abj_75 = [];
bar_offset_variability_awj_75 = [];
wind_offset_variability_abj_75 = [];
wind_offset_variability_awj_75 = [];
bar_offset_variability_abj_100 = [];
bar_offset_variability_awj_100 = [];
wind_offset_variability_abj_100 = [];
wind_offset_variability_awj_100 = [];

configuration = [];

for fly = 1:length(data)
    
    pref_index = [pref_index;data(fly).pref_index];
    pref_index_2 = [pref_index_2;data(fly).pref_index_2];    
    heading_pref_index = [heading_pref_index;data(fly).heading_pref_index];
    stick_index = [stick_index;data(fly).stick_index];
    PI = [PI,data(fly).PI];
    heading_PI = [heading_PI,data(fly).PI_heading];
    SI = [SI,data(fly).SI];
    
    initial_bar_bm{fly} = data(fly).initial_bar_bm;
    initial_wind_bm{fly} = data(fly).initial_wind_bm;
    initial_bar_bw{fly} = data(fly).initial_bar_bw;
    initial_wind_bw{fly} = data(fly).initial_wind_bw;
    initial_bar_offset{fly} = data(fly).initial_bar_offset;
    initial_wind_offset{fly} = data(fly).initial_wind_offset;
    
    short_bm_bar_jump = [short_bm_bar_jump;data(fly).short_bm_bar_jump];
    short_bw_bar_jump = [short_bw_bar_jump;data(fly).short_bw_bar_jump];
    short_bm_wind_jump = [short_bm_wind_jump;data(fly).short_bm_wind_jump];
    short_bw_wind_jump = [short_bw_wind_jump;data(fly).short_bw_wind_jump];
    
    long_bm_bar_jump = [long_bm_bar_jump;data(fly).long_bm_bar_jump];
    long_bw_bar_jump = [long_bw_bar_jump;data(fly).long_bw_bar_jump];
    long_bm_wind_jump = [long_bm_wind_jump;data(fly).long_bm_wind_jump];
    long_bw_wind_jump = [long_bw_wind_jump;data(fly).long_bw_wind_jump];
    
    short_bm_pref_cue = [short_bm_pref_cue;data(fly).short_bm_pref_cue];
    short_bw_pref_cue = [short_bw_pref_cue;data(fly).short_bw_pref_cue];
    short_bm_non_pref_cue = [short_bm_non_pref_cue;data(fly).short_bm_non_pref_cue];
    short_bw_non_pref_cue = [short_bw_non_pref_cue;data(fly).short_bw_non_pref_cue];    
    
    bar_offset_variability_abj_75 = [bar_offset_variability_abj_75;data(fly).bar_offset_variability_abj_75];
    bar_offset_variability_awj_75 = [bar_offset_variability_awj_75;data(fly).bar_offset_variability_awj_75];
    wind_offset_variability_abj_75 = [wind_offset_variability_abj_75;data(fly).wind_offset_variability_abj_75];
    wind_offset_variability_awj_75 = [wind_offset_variability_awj_75;data(fly).wind_offset_variability_awj_75];
    bar_offset_variability_abj_100 = [bar_offset_variability_abj_100;data(fly).bar_offset_variability_abj_100];
    bar_offset_variability_awj_100 = [bar_offset_variability_awj_100;data(fly).bar_offset_variability_awj_100];
    wind_offset_variability_abj_100 = [wind_offset_variability_abj_100;data(fly).wind_offset_variability_abj_100];
    wind_offset_variability_awj_100 = [wind_offset_variability_awj_100;data(fly).wind_offset_variability_awj_100];
    
    configuration = [configuration;data(fly).configuration];
end


%% Short timescale raster plots

% %1) As they are
% %Plot bump magnitude
% figure('Position',[100 100 600 600]),
% subplot(2,1,1)
% imagesc(short_bm_bar_jump)
% hold on
% xline(19,'r','linewidth',2)
% colormap(flipud(gray))
% title('Bar jumps');
% 
% subplot(2,1,2)
% imagesc(short_bm_wind_jump)
% hold on
% xline(19,'r','linewidth',2)
% colormap(flipud(gray))
% title('Wind jumps');
% 
% suptitle('Bump magnitude');
% 
% 
% %Plot bump width
% figure('Position',[100 100 600 600]),
% subplot(2,1,1)
% imagesc(short_bw_bar_jump)
% hold on
% xline(19,'r','linewidth',2)
% colormap(flipud(gray))
% title('Bar jumps');
% 
% subplot(2,1,2)
% imagesc(short_bw_wind_jump)
% hold on
% xline(19,'r','linewidth',2)
% colormap(flipud(gray))
% title('Wind jumps');
% 
% suptitle('Bump width');


%2) z-scored
%zscore
zscored_short_bm_bar_jump = zscore(short_bm_bar_jump,[],2);
zscored_short_bm_wind_jump = zscore(short_bm_wind_jump,[],2);
zscored_short_bw_bar_jump = zscore(short_bw_bar_jump,[],2);
zscored_short_bw_wind_jump = zscore(short_bw_wind_jump,[],2);

%Plot bump magnitude
figure('Position',[100 100 600 600]),
subplot(2,1,1)
imagesc(zscored_short_bm_bar_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_short_bm_wind_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump magnitude');

saveas(gcf,[path,'\groupPlots\bump_mag_rasterplots.png'])


%Plot bump width
figure('Position',[100 100 600 600]),
subplot(2,1,1)
imagesc(zscored_short_bw_bar_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_short_bw_wind_jump)
hold on
xline(19,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump width');

saveas(gcf,[path,'\groupPlots\bump_width_rasterplots.png'])

%% Long raster plots

%zscore
zscored_long_bm_bar_jump = zscore(long_bm_bar_jump,[],2);
zscored_long_bm_wind_jump = zscore(long_bm_wind_jump,[],2);
zscored_long_bw_bar_jump = zscore(long_bw_bar_jump,[],2);
zscored_long_bw_wind_jump = zscore(long_bw_wind_jump,[],2);

%Plot bump magnitude
figure('Position',[100 100 600 600]),
subplot(2,1,1)
imagesc(zscored_long_bm_bar_jump)
hold on
xline(1103,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_long_bm_wind_jump)
hold on
xline(1103,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump magnitude');

%saveas(gcf,[path,'\groupPlots\bump_mag_long_rasterplots.png'])


%Plot bump width
figure('Position',[100 100 600 600]),
subplot(2,1,1)
imagesc(zscored_long_bw_bar_jump)
hold on
xline(1103,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_long_bw_wind_jump)
hold on
xline(1103,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump width');

%saveas(gcf,[path,'\groupPlots\bump_width_long_rasterplots.png'])

%% Mean bump parameters around the jumps, all cues

figure('Position',[100 100 1400 1000]),

%bump magnitude around cue jumps
subplot(1,2,1)
mean_bm_around_cue_jump = nanmean([short_bm_bar_jump;short_bm_wind_jump]);
std_bm_around_cue_jump = std([short_bm_bar_jump;short_bm_wind_jump]);
ste_bm_around_cue_jump = std_bm_around_cue_jump/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(short_bm_bar_jump,2)],mean_bm_around_cue_jump,ste_bm_around_cue_jump,'cmap',gray)
hold on
line([size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_bar_jump,2)]);
xticks([0 size(short_bm_bar_jump,2)/4 size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)*(3/4) size(short_bm_bar_jump,2)]);
xticklabels({'-2','-1','0','1','2'});
ylabel('Bump magnitude');
xlabel('Time (s)');
title('Cue jumps');

%bump width around cue jumps
subplot(1,2,2)
mean_bw_around_cue_jump = nanmean([short_bw_bar_jump;short_bw_wind_jump]);
std_bw_around_cue_jump = std([short_bw_bar_jump;short_bw_wind_jump]);
ste_bw_around_cue_jump = std_bw_around_cue_jump/sqrt(size(short_bw_bar_jump,1)*2);
boundedline([1:size(short_bw_bar_jump,2)],mean_bw_around_cue_jump,ste_bw_around_cue_jump,'cmap',gray)
hold on
line([size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_bar_jump,2)]);
xticks([0 size(short_bm_bar_jump,2)/4 size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)*(3/4) size(short_bm_bar_jump,2)]);
xticklabels({'-2','-1','0','1','2'});
ylabel('Bump width');
xlabel('Time (s)');
title('Cue jumps');

saveas(gcf,[path,'\groupPlots\mean_bump_pars_aj_all_cues.png'])

%% Mean bump parameters around the jumps, by cue type

figure('Position',[100 100 1400 1000]),

%bump magnitude around bar jumps
subplot(2,2,1)
mean_bm_around_bar_jump = nanmean(short_bm_bar_jump);
std_bm_around_bar_jump = std(short_bm_bar_jump);
ste_bm_around_bar_jump = std_bm_around_bar_jump/sqrt(size(short_bm_bar_jump,1));
boundedline([1:size(short_bm_bar_jump,2)],mean_bm_around_bar_jump,ste_bm_around_bar_jump,'cmap',gray)
hold on
line([size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_bar_jump,2)]);
xticks([0 size(short_bm_bar_jump,2)/4 size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)*(3/4) size(short_bm_bar_jump,2)]);
xticklabels({'-2','-1','0','1','2'});
%ylim([1 2]);
ylabel('Bump magnitude');
title('Bar jumps');

%bump magnitude around wind jumps
subplot(2,2,2)
mean_bm_around_wind_jump = nanmean(short_bm_wind_jump);
std_bm_around_wind_jump = std(short_bm_wind_jump);
ste_bm_around_wind_jump = std_bm_around_wind_jump/sqrt(size(short_bm_wind_jump,1));
boundedline([1:size(short_bm_wind_jump,2)],mean_bm_around_wind_jump,ste_bm_around_wind_jump,'cmap',gray)
hold on
line([size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_bar_jump,2)]);
xticks([0 size(short_bm_bar_jump,2)/4 size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)*(3/4) size(short_bm_bar_jump,2)]);
xticklabels({'-2','-1','0','1','2'});
%ylim([1 2]);
title('Wind jumps');

%bump width around bar jumps
subplot(2,2,3)
mean_bw_around_bar_jump = nanmean(short_bw_bar_jump);
std_bw_around_bar_jump = std(short_bw_bar_jump);
ste_bw_around_bar_jump = std_bw_around_bar_jump/sqrt(size(short_bw_bar_jump,1));
boundedline([1:size(short_bw_bar_jump,2)],mean_bw_around_bar_jump,ste_bw_around_bar_jump,'cmap',gray)
hold on
line([size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_bar_jump,2)]);
xticks([0 size(short_bm_bar_jump,2)/4 size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)*(3/4) size(short_bm_bar_jump,2)]);
xticklabels({'-2','-1','0','1','2'});
ylabel('Bump width');
%ylim([2 3]);
xlabel('Time (s)');

%bump magnitude around wind jumps
subplot(2,2,4)
mean_bw_around_wind_jump = nanmean(short_bw_wind_jump);
std_bw_around_wind_jump = std(short_bw_wind_jump);
ste_bw_around_wind_jump = std_bw_around_wind_jump/sqrt(size(short_bw_wind_jump,1));
boundedline([1:size(short_bw_wind_jump,2)],mean_bw_around_wind_jump,ste_bw_around_wind_jump,'cmap',gray)
hold on
line([size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_bar_jump,2)]);
xticks([0 size(short_bm_bar_jump,2)/4 size(short_bm_bar_jump,2)/2 size(short_bm_bar_jump,2)*(3/4) size(short_bm_bar_jump,2)]);
xticklabels({'-2','-1','0','1','2'});
xlabel('Time (s)');
%ylim([2 3]);

saveas(gcf,[path,'\groupPlots\mean_bump_pars_aj_by_cue_type.png'])

%% By preferred cue type

figure('Position',[100 100 700 1000]),

%bump magnitude around preferred cue
subplot(2,2,1)
mean_bm_around_pref_cue = nanmean(short_bm_pref_cue);
std_bm_around_pref_cue = std(short_bm_pref_cue);
ste_bm_around_pref_cue = std_bm_around_pref_cue/sqrt(size(short_bm_pref_cue,1));
boundedline([1:size(short_bm_pref_cue,2)],mean_bm_around_pref_cue,ste_bm_around_pref_cue,'cmap',gray)
hold on
line([size(short_bm_pref_cue,2)/2 size(short_bm_pref_cue,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_pref_cue,2)]);
xticks([0 size(short_bm_pref_cue,2)/4 size(short_bm_pref_cue,2)/2 size(short_bm_pref_cue,2)*(3/4) size(short_bm_pref_cue,2)]);
xticklabels({'-2','-1','0','1','2'});
%ylim([1 2]);
ylabel('Bump magnitude');
title('Preferred cue');

%bump magnitude around non preferred cue
subplot(2,2,2)
mean_bm_around_non_pref_cue = nanmean(short_bm_non_pref_cue);
std_bm_around_non_pref_cue = std(short_bm_non_pref_cue);
ste_bm_around_non_pref_cue = std_bm_around_non_pref_cue/sqrt(size(short_bm_non_pref_cue,1));
boundedline([1:size(short_bm_non_pref_cue,2)],mean_bm_around_non_pref_cue,ste_bm_around_non_pref_cue,'cmap',gray)
hold on
line([size(short_bm_non_pref_cue,2)/2 size(short_bm_non_pref_cue,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bm_non_pref_cue,2)]);
xticks([0 size(short_bm_non_pref_cue,2)/4 size(short_bm_non_pref_cue,2)/2 size(short_bm_non_pref_cue,2)*(3/4) size(short_bm_non_pref_cue,2)]);
xticklabels({'-2','-1','0','1','2'});
%ylim([1 2]);
title('Non-preferred cue');

%bump width around preferred cue
subplot(2,2,3)
mean_bw_around_pref_cue = nanmean(short_bw_pref_cue);
std_bw_around_pref_cue = std(short_bw_pref_cue);
ste_bw_around_pref_cue = std_bw_around_pref_cue/sqrt(size(short_bw_pref_cue,1));
boundedline([1:size(short_bw_pref_cue,2)],mean_bw_around_pref_cue,ste_bw_around_pref_cue,'cmap',gray)
hold on
line([size(short_bw_pref_cue,2)/2 size(short_bw_pref_cue,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bw_pref_cue,2)]);
xticks([0 size(short_bw_pref_cue,2)/4 size(short_bw_pref_cue,2)/2 size(short_bw_pref_cue,2)*(3/4) size(short_bw_pref_cue,2)]);
xticklabels({'-2','-1','0','1','2'});
ylabel('Bump width');
%ylim([2 3]);
xlabel('Time (s)');

%bump magnitude around non preferred cue
subplot(2,2,4)
mean_bw_around_non_pref_cue = nanmean(short_bw_non_pref_cue);
std_bw_around_non_pref_cue = std(short_bw_non_pref_cue);
ste_bw_around_non_pref_cue = std_bw_around_non_pref_cue/sqrt(size(short_bw_non_pref_cue,1));
boundedline([1:size(short_bw_non_pref_cue,2)],mean_bw_around_non_pref_cue,ste_bw_around_non_pref_cue,'cmap',gray)
hold on
line([size(short_bw_non_pref_cue,2)/2 size(short_bw_non_pref_cue,2)/2],[0 3],'color','r','linewidth',3);
xlim([0 size(short_bw_non_pref_cue,2)]);
xticks([0 size(short_bw_non_pref_cue,2)/4 size(short_bw_non_pref_cue,2)/2 size(short_bw_non_pref_cue,2)*(3/4) size(short_bw_non_pref_cue,2)]);
xticklabels({'-2','-1','0','1','2'});
xlabel('Time (s)');
%ylim([2 3]);

saveas(gcf,[path,'\groupPlots\mean_bump_pars_aj_by_pref_cue.png'])
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\CueJump-Experiment\mean_bump_pars_aj_by_pref_cue.png')

%% Longer mean bump parameters around the jumps

figure('Position',[100 100 1400 1000]),

%bump magnitude around cue jumps
subplot(1,2,1)
long_mean_bm_around_cue_jump = nanmean([long_bm_bar_jump;long_bm_wind_jump]);
long_std_bm_around_cue_jump = std([long_bm_bar_jump;long_bm_wind_jump]);
long_ste_bm_around_cue_jump = long_std_bm_around_cue_jump/sqrt(size(long_bm_bar_jump,1)*2);
boundedline([1:size(long_bm_bar_jump,2)],long_mean_bm_around_cue_jump,long_ste_bm_around_cue_jump,'cmap',gray)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
ylabel('Bump magnitude');
xlabel('Time (s)');
title('Cue jumps');

%bump width around cue jumps
subplot(1,2,2)
long_mean_bw_around_cue_jump = nanmean([long_bw_bar_jump;long_bw_wind_jump]);
long_std_bw_around_cue_jump = std([long_bw_bar_jump;long_bw_wind_jump]);
long_ste_bw_around_cue_jump = long_std_bw_around_cue_jump/sqrt(size(long_bw_bar_jump,1)*2);
boundedline([1:size(long_bw_bar_jump,2)],long_mean_bw_around_cue_jump,long_ste_bw_around_cue_jump,'cmap',gray)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
ylabel('Bump width');
xlabel('Time (s)');
title('Cue jumps');

saveas(gcf,[path,'\groupPlots\long_mean_bump_pars_aj_all_cues.png'])

%% Mean bump parameters around the jumps, by cue type

figure('Position',[100 100 1400 1000]),

%bump magnitude around bar jumps
subplot(2,2,1)
long_mean_bm_around_bar_jump = nanmean(long_bm_bar_jump);
long_std_bm_around_bar_jump = std(long_bm_bar_jump);
long_ste_bm_around_bar_jump = long_std_bm_around_bar_jump/sqrt(size(long_bm_bar_jump,1));
boundedline([1:size(long_bm_bar_jump,2)],long_mean_bm_around_bar_jump,long_ste_bm_around_bar_jump,'cmap',gray)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-120','-60','0','60','120'});
ylabel('Bump magnitude');
title('Bar jumps');

%bump magnitude around wind jumps
subplot(2,2,2)
long_mean_bm_around_wind_jump = nanmean(long_bm_wind_jump);
long_std_bm_around_wind_jump = std(long_bm_wind_jump);
long_ste_bm_around_wind_jump = long_std_bm_around_wind_jump/sqrt(size(long_bm_wind_jump,1));
boundedline([1:size(long_bm_wind_jump,2)],long_mean_bm_around_wind_jump,long_ste_bm_around_wind_jump,'cmap',gray)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
title('Wind jumps');

%bump width around bar jumps
subplot(2,2,3)
long_mean_bw_around_bar_jump = nanmean(long_bw_bar_jump);
long_std_bw_around_bar_jump = std(long_bw_bar_jump);
long_ste_bw_around_bar_jump = long_std_bw_around_bar_jump/sqrt(size(long_bw_bar_jump,1));
boundedline([1:size(long_bw_bar_jump,2)],long_mean_bw_around_bar_jump,long_ste_bw_around_bar_jump,'cmap',gray)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
ylabel('Bump width');
xlabel('Time (s)');

%bump magnitude around wind jumps
subplot(2,2,4)
long_mean_bw_around_wind_jump = nanmean(long_bw_wind_jump);
long_std_bw_around_wind_jump = std(long_bw_wind_jump);
long_ste_bw_around_wind_jump = long_std_bw_around_wind_jump/sqrt(size(long_bw_wind_jump,1));
boundedline([1:size(long_bw_wind_jump,2)],long_mean_bw_around_wind_jump,long_ste_bw_around_wind_jump,'cmap',gray)
hold on
line([length(long_bm_bar_jump)/2 length(long_bm_bar_jump)/2],[0 3],'color','r','linewidth',3);
xlim([0 length(long_bm_bar_jump)]);
xticks([0 length(long_bm_bar_jump)/4 length(long_bm_bar_jump)/2 length(long_bm_bar_jump)*(3/4) length(long_bm_bar_jump)]);
xticklabels({'-2','-1','0','1','2'});
xlabel('Time (s)');

saveas(gcf,[path,'\groupPlots\long_mean_bump_pars_aj_by_cue_type.png'])

%% Mean bump parameters by preferred cue


%% Mean offset variability around the jumps, all cues

%window = 75 frames
figure('Position',[100 100 1400 1000]),
mean_offset_var_around_cue_jump_75 = nanmean([wind_offset_variability_abj_75;bar_offset_variability_awj_75]);
std_offset_var_around_cue_jump_75 = std([wind_offset_variability_abj_75;bar_offset_variability_awj_75]);
ste_offset_var_around_cue_jump_75 = std_offset_var_around_cue_jump_75/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(bar_offset_variability_abj_75,2)],mean_offset_var_around_cue_jump_75,ste_offset_var_around_cue_jump_75,'cmap',gray)
hold on
line([length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)/2],[0 1],'color','r','linewidth',3);
xlim([0 length(bar_offset_variability_abj_75)]);
xticks([0 length(bar_offset_variability_abj_75)/4 length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)*(3/4) length(bar_offset_variability_abj_75)]);
xticklabels({'-120','-60','0','60','120'});
ylim([0 1]);
ylabel('Offset variability (circ std)');
xlabel('Time (s)');
title('Window = 75 frames');
saveas(gcf,[path,'\groupPlots\mean_offset_var_all_cues_75_frames.png'])


%window = 100 frames
figure('Position',[100 100 1400 1000]),
mean_offset_var_around_cue_jump_100 = nanmean([wind_offset_variability_abj_100;bar_offset_variability_awj_100]);
std_offset_var_around_cue_jump_100 = std([wind_offset_variability_abj_100;bar_offset_variability_awj_100]);
ste_offset_var_around_cue_jump_100 = std_offset_var_around_cue_jump_100/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(bar_offset_variability_abj_100,2)],mean_offset_var_around_cue_jump_100,ste_offset_var_around_cue_jump_100,'cmap',gray)
hold on
line([length(bar_offset_variability_abj_100)/2 length(bar_offset_variability_abj_100)/2],[0 1],'color','r','linewidth',3);
xlim([0 length(bar_offset_variability_abj_75)]);
xticks([0 length(bar_offset_variability_abj_75)/4 length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)*(3/4) length(bar_offset_variability_abj_75)]);
xticklabels({'-120','-60','0','60','120'});
ylim([0 1]);
ylabel('Offset variability (circ std)');
xlabel('Time (s)');
title('window = 100 frames');
saveas(gcf,[path,'\groupPlots\mean_offset_var_all_cues_100_frames.png'])

%% Mean offset variability around the jumps, by cue

%window = 75 frames
figure('Position',[100 100 1400 1000]),
subplot(1,2,1)
mean_offset_var_around_bar_jump_75 = nanmean(wind_offset_variability_abj_75);
std_offset_var_around_bar_jump_75 = std(wind_offset_variability_abj_75);
ste_offset_var_around_bar_jump_75 = std_offset_var_around_bar_jump_75/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(bar_offset_variability_abj_75,2)],mean_offset_var_around_bar_jump_75,ste_offset_var_around_bar_jump_75,'cmap',gray)
hold on
line([length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)/2],[0 1],'color','r','linewidth',3);
xlim([0 length(bar_offset_variability_abj_75)]);
xticks([0 length(bar_offset_variability_abj_75)/4 length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)*(3/4) length(bar_offset_variability_abj_75)]);
xticklabels({'-120','-60','0','60','120'});
ylim([0 1]);
ylabel('Offset variability (rad)');
xlabel('Time (s)');
title('Bar jumps');

subplot(1,2,2)
mean_offset_var_around_wind_jump_75 = nanmean(bar_offset_variability_awj_75);
std_offset_var_around_wind_jump_75 = std(bar_offset_variability_awj_75);
ste_offset_var_around_wind_jump_75 = std_offset_var_around_wind_jump_75/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(bar_offset_variability_abj_75,2)],mean_offset_var_around_wind_jump_75,ste_offset_var_around_wind_jump_75,'cmap',gray)
hold on
line([length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)/2],[0 1],'color','r','linewidth',3);
xlim([0 length(bar_offset_variability_abj_75)]);
xticks([0 length(bar_offset_variability_abj_75)/4 length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)*(3/4) length(bar_offset_variability_abj_75)]);
xticklabels({'-120','-60','0','60','120'});
ylim([0 1]);
xlabel('Time (s)');
title('Wind jumps');

saveas(gcf,[path,'\groupPlots\mean_offset_var_aj_by_cue_type_75_frames.png'])
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\CueJump-Experiment\mean_offset_var_aj_by_cue_type_75_frames.svg')


%% Repeat for 100 frames

figure('Position',[100 100 1400 1000]),
subplot(1,2,1)
mean_offset_var_around_bar_jump_100 = nanmean(wind_offset_variability_abj_100);
std_offset_var_around_bar_jump_100 = std(wind_offset_variability_abj_100);
ste_offset_var_around_bar_jump_100 = std_offset_var_around_bar_jump_100/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(bar_offset_variability_abj_100,2)],mean_offset_var_around_bar_jump_100,ste_offset_var_around_bar_jump_100,'cmap',gray)
hold on
line([length(bar_offset_variability_abj_100)/2 length(bar_offset_variability_abj_100)/2],[0 1],'color','r','linewidth',3);
xlim([0 length(bar_offset_variability_abj_75)]);
xticks([0 length(bar_offset_variability_abj_75)/4 length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)*(3/4) length(bar_offset_variability_abj_75)]);
xticklabels({'-120','-60','0','60','120'});
ylim([0 1]);
ylabel('Offset variability (circ std)');
xlabel('Time (s)');
title('Bar jumps');

subplot(1,2,2)
mean_offset_var_around_wind_jump_100 = nanmean(bar_offset_variability_awj_100);
std_offset_var_around_wind_jump_100 = std(bar_offset_variability_awj_100);
ste_offset_var_around_wind_jump_100 = std_offset_var_around_wind_jump_100/sqrt(size(short_bm_bar_jump,1)*2);
boundedline([1:size(bar_offset_variability_abj_100,2)],mean_offset_var_around_wind_jump_100,ste_offset_var_around_wind_jump_100,'cmap',gray)
hold on
line([length(bar_offset_variability_abj_100)/2 length(bar_offset_variability_abj_100)/2],[0 1],'color','r','linewidth',3);
xlim([0 length(bar_offset_variability_abj_75)]);
xticks([0 length(bar_offset_variability_abj_75)/4 length(bar_offset_variability_abj_75)/2 length(bar_offset_variability_abj_75)*(3/4) length(bar_offset_variability_abj_75)]);
xticklabels({'-120','-60','0','60','120'});
ylim([0 1]);
xlabel('Time (s)');
title('Wind jumps');

saveas(gcf,[path,'\groupPlots\mean_offset_var_aj_by_cue_type_100_frames.png'])

%% Preference index (bump)

figure('Position',[100 100 1400 1000]),
boxplot(PI,'color','k')
hold on
yline(0);
scatter(repmat([1:length(data)],8,1),PI,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
xlabel('Fly #');
ylabel('Preference index');
ylim([-1 1]);
    
saveas(gcf,[path,'\groupPlots\PI_across_flies.png'])
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\CueJump-Experiment\PI_across_flies.svg')

%combining all flies
PI_across_flies = reshape(PI,8*length(data),1);
figure,
boxplot(PI_across_flies,'color','k')
hold on
yline(0);
scatter(repelem(1,8*length(data)),PI_across_flies,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Preference index');
ylim([-1 1]);
set(gca,'xticklabel',{[]});

saveas(gcf,[path,'\groupPlots\global_PI.png'])


%% Behavioral preference index (bump)

figure('Position',[100 100 1400 1000]),
boxplot(heading_PI,'color','k')
hold on
yline(0);
scatter(repmat([1:length(data)],8,1),heading_PI,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
xlabel('Fly #');
ylabel('Preference index');
ylim([-1 1]);
    
saveas(gcf,[path,'\groupPlots\heading_PI_across_flies.png'])
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\CueJump-Experiment\heading_PI_across_flies.svg')


%combining all flies
PI_across_flies = reshape(heading_PI,8*length(data),1);
figure,
boxplot(PI_across_flies,'color','k')
hold on
yline(0);
scatter(repelem(1,8*length(data)),PI_across_flies,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Preference index');
ylim([-1 1]);
set(gca,'xticklabel',{[]});

saveas(gcf,[path,'\groupPlots\global_heading_PI.png'])


%% Stickiness index

figure,
boxplot(SI,'color','k')
hold on
yline(0);
scatter(repmat([1:length(data)],8,1),SI,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
xlabel('Fly #');
ylabel('Stickiness index');
ylim([-1 1]);
    
saveas(gcf,[path,'\groupPlots\SI_across_flies.png'])


%combining all flies
SI_across_flies = reshape(SI,8*length(data),1);
figure,
boxplot(SI_across_flies,'color','k')
hold on
yline(0);
scatter(repelem(1,8*length(data)),SI_across_flies,[],[.5 .5 .5],'filled')
set(findobj(gca,'type','line'),'linew',2)
ylabel('Stickiness index');
ylim([-1 1]);
set(gca,'xticklabel',{[]});

saveas(gcf,[path,'\groupPlots\global_SI.png'])

%% Relationship between mean PI and initial cue variability ratio

%1) compute initial offset ratio
for fly = 1:length(data)
    initial_bar_offset_var(fly) = circ_std(deg2rad(initial_bar_offset{fly}));
    initial_wind_offset_var(fly) = circ_std(deg2rad(initial_wind_offset{fly}),[],[],2);
    initial_offset_ratio(fly) = initial_bar_offset_var(fly)/initial_wind_offset_var(fly);
end
%if the ratio < 1, then bar offset var < wind offset var

%2) relate preference index to initial offset variabililty
figure('Position',[100 100 1000 800]),
plot(initial_offset_ratio(configuration==1),pref_index(configuration==1),'ko','MarkerFaceColor','r','MarkerSize',8)
hold on
plot(initial_offset_ratio(configuration==2),pref_index(configuration==2),'ko','MarkerFaceColor','b','MarkerSize',8)
legend('Bar first','Wind first');
xlabel('Initial offset variability ratio (bar offset var / wind offset var)','fontsize',14);
ylabel('Preference index','fontsize',14);
ylim([-1 1]);
xline(1,'HandleVisibility','off'); yline(0,'HandleVisibility','off');
saveas(gcf,[path,'\groupPlots\pref_ind_vs_initial_offset_var.png'])
saveas(gcf,'C:\Users\Melanie\Dropbox (HMS)\Manuscript-Basnak\CueJump-Experiment\pref_ind_vs_initial_offset_var.svg')


%% 3) repeat for behavioral preference index

figure,
plot(initial_offset_ratio(configuration==1),heading_pref_index(configuration==1),'ko','MarkerFaceColor','r')
hold on
plot(initial_offset_ratio(configuration==2),heading_pref_index(configuration==2),'ko','MarkerFaceColor','b')
legend('Bar first','Wind first');xlabel('Initial offset variability ratio (bar offset var / wind offset var)');
ylabel('Behavioral preference index');
ylim([-1 1]);
xline(1); yline(0);
saveas(gcf,[path,'\groupPlots\heading_pref_ind_vs_initial_offset_var.png'])

%% Look at PI evolution to see if there are consistent trends

for fly = 1:length(PI)
    
    figure,
    plot(PI(:,fly),'-ko')
    hold on
    if configuration(fly) == 1
       plot([1:2:8],PI(1:2:8,fly),'ko','MarkerFaceColor','r')
       plot([2:2:8],PI(2:2:8,fly),'ko','MarkerFaceColor','b')       
    else
       plot([1:2:8],PI(1:2:8,fly),'ko','MarkerFaceColor','b')
       plot([2:2:8],PI(2:2:8,fly),'ko','MarkerFaceColor','r')    
    end
    ylim([-1 1]);
    ylabel('Preference index');
    xlabel('Jump #');
    yline(0);
    title(['Fly #',num2str(fly)]);
    
end

%% Using the preference index computed as the ratio of means

figure,
plot(initial_offset_ratio,pref_index_2,'ko')
xlabel('Initial offset variability ratio (bar offset var / wind offset var)');
ylabel('Preference index');
ylim([-1 1]);
xline(1); yline(0);
saveas(gcf,[path,'\groupPlots\pref_ind_2_vs_initial_offset_var.png'])

%% Clear space

clear all; close all;
