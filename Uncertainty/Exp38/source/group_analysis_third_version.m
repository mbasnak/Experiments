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
long_bm_bar_jump = [];
long_bw_bar_jump = [];
long_bm_wind_jump = [];
long_bw_wind_jump = [];
short_bm_bar_jump = [];
short_bw_bar_jump = [];
short_bm_wind_jump = [];
short_bw_wind_jump = [];

for fly = 1:length(data)
    
    pref_index = [pref_index;data(fly).pref_index];
    long_bm_bar_jump = [long_bm_bar_jump;data(fly).long_bm_bar_jump];
    long_bw_bar_jump = [long_bw_bar_jump;data(fly).long_bw_bar_jump];
    long_bm_wind_jump = [long_bm_wind_jump;data(fly).long_bm_wind_jump];
    long_bw_wind_jump = [long_bw_wind_jump;data(fly).long_bw_wind_jump];
    short_bm_bar_jump = [short_bm_bar_jump;data(fly).short_bm_bar_jump];
    short_bw_bar_jump = [short_bw_bar_jump;data(fly).short_bw_bar_jump];
    short_bm_wind_jump = [short_bm_wind_jump;data(fly).short_bm_wind_jump];
    short_bw_wind_jump = [short_bw_wind_jump;data(fly).short_bw_wind_jump];
    
end


%% Short timescale raster plots

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


%% Long timescale raster plots

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
xline(1102,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_long_bm_wind_jump)
hold on
xline(1102,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump magnitude');


%Plot bump width
figure('Position',[100 100 600 600]),
subplot(2,1,1)
imagesc(zscored_long_bw_bar_jump)
hold on
xline(1102,'r','linewidth',2)
colormap(flipud(gray))
title('Bar jumps');

subplot(2,1,2)
imagesc(zscored_long_bw_wind_jump)
hold on
xline(1102,'r','linewidth',2)
colormap(flipud(gray))
title('Wind jumps');

suptitle('Bump width');